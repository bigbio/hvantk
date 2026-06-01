import csv
import gzip
import io
import logging
import os
import os.path as path
import shutil
import struct
import subprocess
import time
import zipfile
import zlib

import requests
from tqdm import tqdm
from hvantk.core.utils.bgzf import BGZF_BLOCK_SIZE, BgzfWriter, make_bgzf_block

logger = logging.getLogger(__name__)

# Backward-compatible aliases for internal helpers used in tests/downstream code.
_BGZF_BLOCK_SIZE = BGZF_BLOCK_SIZE
_make_bgzf_block = make_bgzf_block


def download_file(url: str, out_dir: str, file_name: str):
    """
    Downloads a file from the specified URL to a local directory with a given file name.

    Validates that the file name does not contain path traversal components. Creates the output directory if it does not exist. Raises an exception if the download fails or if the file name is invalid.

    Returns:
        str: The local filesystem path to the downloaded file on success.
    """
    # Validate file_name doesn't contain path traversal
    if path.isabs(file_name) or ".." in file_name:
        raise ValueError(
            f"Invalid file_name: {file_name}. Must be a simple filename without path components."
        )

    os.makedirs(out_dir, exist_ok=True)
    local_path = os.path.join(out_dir, file_name)
    logger.debug(f"Downloading {url} to {local_path}")
    try:
        response = requests.get(url, stream=True, timeout=30)  # Add timeout
        response.raise_for_status()  # Raises exception for 4XX/5XX responses
        logger.debug(f"Download response status code: {response.status_code}")

        total_size = int(response.headers.get("content-length", 0))

        # Optional: Add file size check
        max_size = 1024 * 1024 * 1024  # 1GB
        if total_size > max_size:
            logger.warning(
                f"File is very large ({total_size/1024/1024:.1f} MB), exceeding recommended size of {max_size/1024/1024:.1f} MB"
            )

        with open(local_path, "wb") as f, tqdm(
            desc=file_name,
            total=total_size,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
        ) as bar:
            for chunk in response.iter_content(chunk_size=8192):
                size = f.write(chunk)
                bar.update(size)
        logger.info(f"Download completed: {local_path}")
        return local_path

    except requests.exceptions.RequestException as e:
        # Clean up partial download if it exists
        if os.path.exists(local_path):
            os.remove(local_path)
        logger.exception(f"Failed to download {url}: {str(e)}")
        raise Exception(f"Failed to download {url}: {str(e)}")


def sanitize_tsv(filepath: str, *, delimiter: str = "\t") -> int:
    """Sanitize a delimited text file in-place for Hail compatibility.

    Handles two common problems that ``hl.import_table`` cannot deal with:

    1. **Blank lines** — removed entirely.
    2. **Multiline quoted fields** — newlines inside quoted values are
       replaced with spaces so every logical record occupies exactly one line.

    The file is rewritten atomically (write to temp, then rename).

    Parameters
    ----------
    filepath : str
        Path to the file to sanitize.
    delimiter : str
        Column delimiter (default: tab).

    Returns
    -------
    int
        Number of data rows in the cleaned file (excluding header).
    """
    logger.info("Sanitizing %s (delimiter=%r)", filepath, delimiter)
    tmp_path = filepath + ".sanitized.tmp"
    n_cols = None
    n_data = 0

    with open(filepath, newline="") as fh_in, open(tmp_path, "w", newline="") as fh_out:
        reader = csv.reader(fh_in, delimiter=delimiter, quotechar='"')
        writer = csv.writer(
            fh_out,
            delimiter=delimiter,
            quotechar='"',
            lineterminator="\n",
        )
        for row_index, row in enumerate(reader):
            if not any(cell.strip() for cell in row):
                continue

            if n_cols is None:
                n_cols = len(row)

            # Pad or truncate to match header column count
            if len(row) < n_cols:
                logger.warning(
                    "sanitize_tsv: row %d has %d columns, expected %d; padding",
                    row_index,
                    len(row),
                    n_cols,
                )
                row = row + [""] * (n_cols - len(row))
            elif len(row) > n_cols:
                logger.warning(
                    "sanitize_tsv: row %d has %d columns, expected %d; truncating",
                    row_index,
                    len(row),
                    n_cols,
                )
                row = row[:n_cols]

            # Flatten embedded newlines within cells
            cleaned = [cell.replace("\n", " ").replace("\r", " ") for cell in row]
            writer.writerow(cleaned)
            n_data += 1

    if n_cols is None:
        logger.warning("sanitize_tsv: file %s is empty after cleanup", filepath)
        os.remove(tmp_path)
        return 0

    os.replace(tmp_path, filepath)
    n_data -= 1  # exclude header row
    logger.info("Sanitized %s: %d columns, %d data rows", filepath, n_cols, n_data)
    return n_data


# define a function to test if a given url exists


def url_exists(url: str) -> bool:
    """
    Checks whether the specified URL exists by performing an HTTP HEAD request.

    Returns:
        True if the URL responds with status code 200; otherwise, False.
    """
    try:
        logger.debug(f"Checking if URL exists: {url}")
        response = requests.head(url, allow_redirects=True)
        if response.status_code == 200:
            logger.debug(f"URL {url} exists with status code: {response.status_code}")
            return True
        else:
            logger.warning(
                f"URL {url} does not exist with status code: {response.status_code}"
            )
            return False
    except requests.RequestException as e:
        logger.exception(f"Error checking URL {url}: {e}")
        return False


def compress_files(
    source_dir: str, output_zip: str, remove_originals: bool = False
) -> None:
    """
    Compresses all files in a directory into a ZIP archive, preserving folder structure.

    Args:
        source_dir: Path to the directory whose contents will be compressed.
        output_zip: Path to the output ZIP file.
        remove_originals: If True, deletes the source directory after compression.

    Raises:
        FileNotFoundError: If the source directory does not exist.
        ValueError: If the source path is not a directory.
        PermissionError: If there are file access permission issues.
    """
    if not os.path.exists(source_dir):
        raise FileNotFoundError(f"Source directory '{source_dir}' does not exist")

    if not os.path.isdir(source_dir):
        raise ValueError(f"Source path '{source_dir}' is not a directory")

    # Create parent directory for output_zip if it doesn't exist
    output_dir = os.path.dirname(output_zip)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    with zipfile.ZipFile(output_zip, "w", zipfile.ZIP_DEFLATED) as zipf:
        for root, _, files in os.walk(source_dir):
            for file in files:
                file_path = os.path.join(root, file)
                # Use relative path in the archive to preserve folder structure.
                arcname = os.path.relpath(file_path, source_dir)
                zipf.write(file_path, arcname)
    logger.info(f"Compressed '{source_dir}' into '{output_zip}'")

    if remove_originals:
        shutil.rmtree(source_dir)
        logger.info(f"Removed original directory '{source_dir}' after compression.")


def decompress_files(
    zip_path: str, extract_to: str, remove_originals: bool = False
) -> None:
    """
    Extracts all files from a ZIP archive to a target directory.

    If the ZIP file is corrupted, a BadZipFile exception is raised. Optionally deletes the original ZIP file after extraction if remove_originals is True.

    Args:
        zip_path: Path to the ZIP archive.
        extract_to: Directory where files will be extracted.
        remove_originals: If True, deletes the ZIP file after extraction.

    Raises:
        FileNotFoundError: If the ZIP file does not exist.
        zipfile.BadZipFile: If the ZIP file is invalid or corrupted.
        PermissionError: If there are file access permission issues.
    """
    if not os.path.exists(zip_path):
        raise FileNotFoundError(f"Zip file '{zip_path}' does not exist")

    # Create extract directory if it doesn't exist
    os.makedirs(extract_to, exist_ok=True)

    with zipfile.ZipFile(zip_path, "r") as zipf:
        # Check if the zip file is valid
        if zipf.testzip() is not None:
            raise zipfile.BadZipFile(f"'{zip_path}' is corrupted")

        zipf.extractall(extract_to)
    logger.info(f"Extracted '{zip_path}' into '{extract_to}'")

    if remove_originals:
        os.remove(zip_path)
        logger.info(f"Removed archive file '{zip_path}' after decompression.")


# ---------------------------------------------------------------------------
# BGZF / GZIP detection and conversion
# ---------------------------------------------------------------------------

def is_gzipped(filepath: str) -> bool:
    """Check if a file starts with the gzip magic bytes (``\\x1f\\x8b``).

    Parameters
    ----------
    filepath : str
        Path to the file to check.

    Returns
    -------
    bool
        True if the file begins with the gzip magic number.

    Raises
    ------
    FileNotFoundError
        If *filepath* does not exist.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    with open(filepath, "rb") as f:
        magic = f.read(2)
    return magic == b"\x1f\x8b"


def _is_bgzf_header(header: bytes) -> bool:
    """Return True if *header* (>= 18 bytes) is a valid BGZF block header."""
    return (
        len(header) >= 18
        and header[0:2] == b"\x1f\x8b"  # gzip magic
        and header[2:3] == b"\x08"  # deflate method
        and (header[3] & 0x04) != 0  # FEXTRA flag set
        and header[12:14] == b"BC"  # BGZF subfield marker
    )


def is_bgzf(filepath: str, num_blocks: int = 3) -> bool:
    """Check if a file is block-gzip (BGZF) compressed.

    Validates up to *num_blocks* consecutive BGZF blocks to avoid false
    positives from files whose first block header matches BGZF but whose
    remaining content is standard gzip.

    Parameters
    ----------
    filepath : str
        Path to the file to check.
    num_blocks : int
        Number of consecutive blocks to validate (default 3).

    Returns
    -------
    bool
        True if the file conforms to BGZF format.

    Raises
    ------
    FileNotFoundError
        If *filepath* does not exist.
    """
    if num_blocks < 1:
        raise ValueError(f"num_blocks must be >= 1, got {num_blocks}")
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    file_size = os.path.getsize(filepath)
    checked = 0
    with open(filepath, "rb") as f:
        for _ in range(num_blocks):
            block_start = f.tell()
            header = f.read(18)
            if len(header) < 18:
                # EOF — valid only if we checked at least one block
                return checked > 0
            if not _is_bgzf_header(header):
                return False
            checked += 1
            # BSIZE (little-endian uint16 at offset 16) = total block size - 1
            bsize = int.from_bytes(header[16:18], byteorder="little")
            block_end = block_start + bsize + 1
            if block_end > file_size:
                return False  # truncated/corrupt block
            f.seek(block_end)
    return True


def detect_compression(filepath: str) -> str:
    """Detect the compression type of a file.

    Parameters
    ----------
    filepath : str
        Path to the file to check.

    Returns
    -------
    str
        One of ``'bgzf'``, ``'gzip'``, or ``'none'``.

    Raises
    ------
    FileNotFoundError
        If *filepath* does not exist.
    """
    if is_bgzf(filepath):
        return "bgzf"
    if is_gzipped(filepath):
        return "gzip"
    return "none"


# ---------------------------------------------------------------------------
# Conversion: standard gzip → BGZF
# ---------------------------------------------------------------------------


def _get_conversion_backend() -> str:
    """Detect the best available backend for gz-to-bgz conversion.

    Returns
    -------
    str
        ``'bgzip'`` if the system ``bgzip`` binary is on PATH,
        ``'pysam'`` if the ``pysam`` package is importable, or
        ``'python'`` as the stdlib-only fallback.
    """
    if shutil.which("bgzip"):
        return "bgzip"
    try:
        import pysam  # noqa: F401

        return "pysam"
    except ImportError:
        return "python"


def _convert_with_bgzip(input_path: str, output_path: str, threads: int) -> None:
    """Convert using system ``bgzip`` via subprocess pipe."""
    with open(output_path, "wb") as f_out:
        p1 = subprocess.Popen(
            ["gzip", "-dc", input_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        p2 = subprocess.Popen(
            ["bgzip", "-c", f"-@{threads}"],
            stdin=p1.stdout,
            stdout=f_out,
            stderr=subprocess.PIPE,
        )
        p1.stdout.close()  # allow p1 to receive SIGPIPE if p2 exits
        _, p2_stderr = p2.communicate()
        p1.wait()
        p1_stderr = p1.stderr.read()
        p1.stderr.close()
    errors = []
    if p1.returncode != 0:
        errors.append(
            f"gzip -dc failed (rc={p1.returncode}): {p1_stderr.decode(errors='replace').strip()}"
        )
    if p2.returncode != 0:
        errors.append(
            f"bgzip failed (rc={p2.returncode}): {p2_stderr.decode(errors='replace').strip()}"
        )
    if errors:
        # Remove partial output
        if os.path.exists(output_path):
            os.remove(output_path)
        raise RuntimeError("bgzip conversion failed:\n" + "\n".join(errors))


def _convert_with_pysam(input_path: str, output_path: str) -> None:
    """Convert using pysam's ``BGZFile`` writer.

    Since ``pysam.tabix_compress`` reads the input file as raw bytes (which
    would double-compress a ``.gz`` input), we decompress with ``gzip`` and
    write through ``pysam.BGZFile`` in write mode.
    """
    import pysam

    with gzip.open(input_path, "rb") as fin, pysam.BGZFile(output_path, "wb") as fout:
        while True:
            chunk = fin.read(BGZF_BLOCK_SIZE)
            if not chunk:
                break
            fout.write(chunk)


def _convert_with_python(input_path: str, output_path: str) -> None:
    """Convert using pure Python ``zlib`` BGZF block writing."""
    with gzip.open(input_path, "rb") as fin, open(output_path, "wb") as fout:
        while True:
            chunk = fin.read(BGZF_BLOCK_SIZE)
            if not chunk:
                break
            fout.write(make_bgzf_block(chunk))
        # Write the empty EOF block required by BGZF
        fout.write(make_bgzf_block(b""))


def convert_gz_to_bgz(
    input_path: str,
    output_path: str | None = None,
    threads: int = 4,
) -> str:
    """Convert a standard gzip file to block gzip (BGZF).

    Uses a tiered strategy, automatically selecting the best available
    method:

    1. System ``bgzip`` via subprocess (fastest, multi-threaded).
    2. ``pysam.tabix_compress()`` (fast C implementation).
    3. Pure Python ``zlib`` BGZF block writer (stdlib-only fallback).

    Parameters
    ----------
    input_path : str
        Path to the input ``.gz`` file.
    output_path : str, optional
        Path for the output ``.bgz`` file.  Defaults to *input_path*
        with ``.gz`` replaced by ``.bgz``.
    threads : int, optional
        Number of threads for system ``bgzip`` (default 4).  Ignored by
        other backends.

    Returns
    -------
    str
        Path to the converted ``.bgz`` file.

    Raises
    ------
    FileNotFoundError
        If *input_path* does not exist.
    RuntimeError
        If conversion fails.
    """
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"File not found: {input_path}")

    if output_path is None:
        if input_path.endswith(".gz"):
            output_path = input_path[:-3] + ".bgz"
        else:
            output_path = input_path + ".bgz"

    # Reuse existing output if it's already valid BGZF
    if os.path.exists(output_path) and detect_compression(output_path) == "bgzf":
        logger.info("Reusing existing BGZF file '%s'", output_path)
        return output_path

    backend = _get_conversion_backend()
    t0 = time.time()

    if backend == "bgzip":
        logger.info(
            "Converting '%s' to BGZF using system bgzip with %d threads...",
            input_path,
            threads,
        )
        _convert_with_bgzip(input_path, output_path, threads)
    elif backend == "pysam":
        logger.info(
            "Converting '%s' to BGZF using pysam (system bgzip not found)...",
            input_path,
        )
        _convert_with_pysam(input_path, output_path)
    else:
        logger.info(
            "Converting '%s' to BGZF using pure Python (pysam and bgzip not available)...",
            input_path,
        )
        _convert_with_python(input_path, output_path)

    elapsed = time.time() - t0
    logger.info("Converted to '%s' (%.1fs, backend=%s)", output_path, elapsed, backend)
    return output_path


# ---------------------------------------------------------------------------
# BGZF validation
# ---------------------------------------------------------------------------


def validate_bgzf(filepath: str, full: bool = True) -> tuple[bool, int, str]:
    """Validate all BGZF blocks in *filepath*.

    Parameters
    ----------
    filepath : str
        Path to a file to validate.
    full : bool
        If True, decompress each block and verify CRC32 + ISIZE.
        If False, only check block headers (much faster).

    Returns
    -------
    tuple[bool, int, str]
        (is_valid, block_number, message).  On success *block_number* is the
        total number of validated blocks; on failure it is the 1-based index
        of the failing block.
    """
    filesize = os.path.getsize(filepath)
    block_num = 0

    with open(filepath, "rb") as f:
        while True:
            block_start = f.tell()
            header = f.read(18)

            if len(header) == 0:
                return True, block_num, "OK"

            if len(header) < 18:
                return (
                    False,
                    block_num + 1,
                    f"Truncated header ({len(header)} bytes) at offset {block_start}",
                )

            if header[0:2] != b"\x1f\x8b":
                return (
                    False,
                    block_num + 1,
                    f"Bad gzip magic at offset {block_start} (got {header[0:2].hex()})",
                )
            if header[2] != 0x08:
                return (
                    False,
                    block_num + 1,
                    f"Bad compression method at offset {block_start}",
                )
            if (header[3] & 0x04) == 0:
                return (
                    False,
                    block_num + 1,
                    f"FEXTRA flag not set at offset {block_start}",
                )
            if header[12:14] != b"BC":
                return (
                    False,
                    block_num + 1,
                    f"Missing BC subfield at offset {block_start} (got {header[12:14].hex()})",
                )

            xlen = int.from_bytes(header[10:12], "little")
            if xlen < 6:
                return (
                    False,
                    block_num + 1,
                    (
                        f"Invalid BGZF header at offset {block_start}: "
                        f"xlen too small ({xlen}, expected >= 6)"
                    ),
                )

            bsize = int.from_bytes(header[16:18], "little")
            block_size = bsize + 1

            if block_start + block_size > filesize:
                return (
                    False,
                    block_num + 1,
                    (
                        f"Block extends past EOF at offset {block_start} "
                        f"(block_size={block_size}, remaining={filesize - block_start})"
                    ),
                )

            if full:
                f.seek(block_start)
                block_data = f.read(block_size)
                if len(block_data) != block_size:
                    return (
                        False,
                        block_num + 1,
                        f"Truncated block at offset {block_start}",
                    )

                header_size = 12 + xlen
                if header_size + 8 > block_size:
                    return (
                        False,
                        block_num + 1,
                        (
                            f"Invalid BGZF header at offset {block_start}: "
                            f"header_size ({header_size}) + footer exceeds "
                            f"block_size ({block_size})"
                        ),
                    )
                cdata = block_data[header_size:-8]
                stored_crc = struct.unpack("<I", block_data[-8:-4])[0]
                stored_isize = struct.unpack("<I", block_data[-4:])[0]

                try:
                    decompressed = zlib.decompress(cdata, -15)
                except zlib.error as e:
                    return (
                        False,
                        block_num + 1,
                        f"Decompression error at offset {block_start}: {e}",
                    )

                calc_crc = zlib.crc32(decompressed) & 0xFFFFFFFF
                if calc_crc != stored_crc:
                    return (
                        False,
                        block_num + 1,
                        (
                            f"CRC32 mismatch at offset {block_start} "
                            f"(stored={stored_crc:#010x}, calc={calc_crc:#010x})"
                        ),
                    )

                calc_isize = len(decompressed) & 0xFFFFFFFF
                if calc_isize != stored_isize:
                    return (
                        False,
                        block_num + 1,
                        (
                            f"ISIZE mismatch at offset {block_start} "
                            f"(stored={stored_isize}, calc={calc_isize})"
                        ),
                    )
            else:
                f.seek(block_start + block_size)

            block_num += 1

    return True, block_num, "OK"


# ---------------------------------------------------------------------------
# Smart import wrapper
# ---------------------------------------------------------------------------


def resolve_compression(
    filepath: str,
    force_bgz: bool = True,
    auto_convert: bool = False,
    threads: int = 4,
) -> tuple[str, bool]:
    """Resolve file compression for Hail import.

    Detects whether *filepath* is BGZF, standard gzip, or uncompressed, and
    returns the (possibly converted) path together with the appropriate
    ``force_bgz`` flag for ``hl.import_*``.

    Parameters
    ----------
    filepath : str
        Path to the input file.
    force_bgz : bool
        Desired ``force_bgz`` setting (only honoured when the file is
        actually BGZF).
    auto_convert : bool
        When True, convert all gzip-family files (both standard gzip and
        BGZF) to a clean BGZF file before import.  This ensures Hail
        compatibility even when files pass header-level BGZF checks but
        contain corrupted blocks deeper in the file.
    threads : int
        Thread count passed to :func:`convert_gz_to_bgz` when converting.

    Returns
    -------
    tuple[str, bool]
        ``(filepath, force_bgz)`` ready to pass to ``hl.import_*``.
    """
    compression = detect_compression(filepath)

    if auto_convert and compression in ("gzip", "bgzf"):
        logger.info(
            "File '%s' detected as %s. Converting to BGZF (auto_convert=True)...",
            filepath,
            compression,
        )
        bgz_path = convert_gz_to_bgz(filepath, threads=threads)
        return bgz_path, True

    if compression == "bgzf":
        logger.debug("File '%s' detected as BGZF — no conversion needed.", filepath)
        return filepath, force_bgz

    if compression == "gzip":
        logger.warning(
            "File '%s' is standard gzip, not block gzip (BGZF). "
            "Hail will read this file on a single core, which is significantly slower. "
            "For parallel multi-core reading, convert to BGZF:\n"
            "  hvantk utils convert-bgz %s\n"
            "Or re-run with --auto-convert-bgz to convert automatically.",
            filepath,
            filepath,
        )
        return filepath, False

    # Uncompressed
    return filepath, False

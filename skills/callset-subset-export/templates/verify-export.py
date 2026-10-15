"""Verify an exported subset VCF, write SHA256SUMS and the recipient's README.
Needs no Hail: reads the block-gzipped VCF as gzip. Runs inside the image.

    python verify-export.py <out dir> --prefix <name> --samples <list.txt>
        [--source <callset.vds>] [--commit SHA] [--image-sha256 HEX] [--recipe <text file>]
        [--requested-samples <original list>]

Writes <out>/verify.json, <out>/SHA256SUMS, <out>/README.md; exits 1 if a check fails.
"""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import json
import sys
from pathlib import Path


def read_list(path: Path) -> list[str]:
    return sorted(
        {
            s.strip()
            for s in path.read_text().splitlines()
            if s.strip() and not s.startswith("#")
        }
    )


def inspect_vcf(path: Path) -> dict:
    samples: list[str] = []
    records = 0
    multiallelic = 0
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
            elif not line.startswith("#"):
                records += 1
                if "," in line.split("\t", 5)[4]:
                    multiallelic += 1
    return {
        "samples": sorted(samples),
        "records": records,
        "biallelic": multiallelic == 0,
    }


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--requested-samples")
    ap.add_argument("--source", default="")
    ap.add_argument("--commit", default="")
    ap.add_argument("--image-sha256", default="")
    ap.add_argument("--recipe")
    args = ap.parse_args()
    out = Path(args.out)
    vcf = out / f"{args.prefix}.vcf.bgz"
    if not vcf.is_file():
        print(f"ERROR: {vcf} not found", file=sys.stderr)
        return 1
    wanted = read_list(Path(args.samples))
    info = inspect_vcf(vcf)
    verify = {
        "vcf": str(vcf),
        "n_samples": len(info["samples"]),
        "samples": info["samples"],
        "samples_match": info["samples"] == wanted,
        "records": info["records"],
        "biallelic": info["biallelic"],
        "checked_at": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    (out / "verify.json").write_text(json.dumps(verify, indent=2))

    files = sorted(
        p
        for p in out.iterdir()
        if p.is_file() and p.name not in ("SHA256SUMS", "README.md")
    )
    (out / "SHA256SUMS").write_text(
        "".join(f"{sha256_file(p)}  {p.name}\n" for p in files)
    )

    recipe = Path(args.recipe).read_text() if args.recipe else "(not recorded)"
    requested = (
        read_list(Path(args.requested_samples)) if args.requested_samples else wanted
    )
    readme = f"""# {args.prefix}: sample subset of a joint callset

Produced {verify["checked_at"]} with hvantk.

- Source callset: `{args.source or "(not recorded)"}`
- Code commit: `{args.commit or "(not recorded)"}`; image sha256: `{args.image_sha256 or "(not recorded)"}`
- Samples requested ({len(requested)}): {", ".join(requested)}
- Samples in the callset's spelling ({len(wanted)}): {", ".join(wanted)}
- Records in `{vcf.name}`: {verify["records"]}; biallelic: {verify["biallelic"]}; header samples match the list: {verify["samples_match"]}

## How it was made

The subset was taken from the lossless VDS (alleles no listed sample carries are
dropped), densified, and exported as a biallelic VCF with genotype-quality filtering
and a minimum allele count of 1.

```
{recipe.rstrip()}
```

## Checking the files

`sha256sum -c SHA256SUMS` in this directory must pass. The VCF has no index unless one
was added afterwards; `tabix -p vcf {vcf.name}` or `bcftools index -t {vcf.name}` creates it.
"""
    (out / "README.md").write_text(readme)
    print(
        json.dumps(
            {
                k: verify[k]
                for k in ("n_samples", "samples_match", "records", "biallelic")
            }
        )
    )
    return 0 if (verify["samples_match"] and verify["biallelic"]) else 1


if __name__ == "__main__":
    sys.exit(main())

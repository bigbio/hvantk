"""Helpers shared by INSIDER's datasets (`variants`, `interfaces`).

INSIDER became a multi-dataset provider when `interfaces` was added, and both builders
have to accept the same two input shapes: a plain filesystem path from `reprocess`, and a
``file://`` URI from the snapshot-test adapter (`Path(...).as_uri()`). That normalisation
was written twice, once per builder, which is exactly what the shared/ folder exists to
prevent -- see `hvantk/skills/_conventions/SKILL.md` and the `cptac` provider, which
splits `expression/`, `phospho/` and `shared/` the same way.
"""

_FILE_URI_PREFIX = "file://"


def normalize_hadoop_path(path: str) -> str:
    """Strip a ``file://`` prefix so stdlib/filesystem APIs can use the path.

    Hail accepts URIs; ``open()``, ``os.path.isdir`` and friends do not. Anything that is
    not a local file URI is returned unchanged, so remote URIs (gs://, hdfs://) pass
    through to callers that can handle them rather than being silently mangled.
    """
    return path[len(_FILE_URI_PREFIX):] if path.startswith(_FILE_URI_PREFIX) else path

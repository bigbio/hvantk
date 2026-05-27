"""Errors specific to artifact serialization."""


class ArtifactTypeError(Exception):
    """Raised when an artifact on disk doesn't match the requested type."""


class SchemaIdMismatchError(Exception):
    """Raised when a loaded artifact's schema_id disagrees with its manifest."""

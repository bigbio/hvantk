"""Stub builder for loader tests."""

SENTINEL_OUTPUT = {"called": True}


def build(input_path: str, output_path: str, **kwargs):
    return {"input_path": input_path, "output_path": output_path, **kwargs, **SENTINEL_OUTPUT}

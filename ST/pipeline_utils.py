"""Shared filesystem and Visium validation helpers (no scientific dependencies)."""

from pathlib import Path


def input_path(value, *, directory=False):
    """Resolve inputs relative to the invocation directory, before output setup."""
    path = Path(value).expanduser().resolve()
    valid = path.is_dir() if directory else path.is_file()
    if not valid:
        raise FileNotFoundError(f"Input {'directory' if directory else 'file'} not found: {path}")
    return path


def output_directory(value):
    path = Path(value).expanduser().resolve()
    path.mkdir(parents=True, exist_ok=True)
    return path


def select_library(adata, requested=None):
    """Select a real Visium library instead of assuming its name is 'spatial'."""
    libraries = list(adata.uns.get("spatial", {}))
    if requested is not None:
        if requested not in libraries:
            raise ValueError(f"Unknown library {requested!r}; available: {libraries}")
        return requested
    if len(libraries) != 1:
        raise ValueError(f"Specify --library_id; available libraries: {libraries}")
    return libraries[0]

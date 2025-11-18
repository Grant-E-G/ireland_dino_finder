import h5py
import numpy as np
from typing import Tuple, Dict, Any


def load_wave_file(
    path: str,
    field_name: str = "u",
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Load wave simulation output written by Hdf5WaveWriter.

    Returns
    -------
    u : np.ndarray
        3D array with shape (T, ny, nx)
    meta : dict
        Dictionary with optional keys:
        - 'x', 'y' : coordinate arrays (1D)
        - 'dt', 'dx', 'dy' : floats (attrs)
    """
    meta: Dict[str, Any] = {}

    with h5py.File(path, "r") as f:
        if field_name not in f:
            raise KeyError(f"Dataset '{field_name}' not found in {path}")

        # load main field
        u = f[field_name][...]  # shape: (T, ny, nx)

        # optional coordinate vectors
        if "x" in f:
            meta["x"] = f["x"][...]
        if "y" in f:
            meta["y"] = f["y"][...]

        # attributes
        for name in ("dt", "dx", "dy"):
            if name in f.attrs:
                meta[name] = f.attrs[name]

    return u, meta

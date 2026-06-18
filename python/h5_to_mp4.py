import argparse
from pathlib import Path

import h5py
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np


def render_hdf5_to_mp4(input_path: str, output_path: str, field_name: str = "u", fps: int = 30) -> None:
    with h5py.File(input_path, "r") as h5:
        if field_name not in h5:
            raise KeyError(f"Dataset '{field_name}' not found in {input_path}")
        data = h5[field_name][...]

    if data.ndim != 3:
        raise ValueError(f"Expected field shape (time, ny, nx), got {data.shape}")

    vmax = float(np.max(np.abs(data)))
    if vmax == 0.0:
        vmax = 1.0

    fig, ax = plt.subplots()
    image = ax.imshow(data[0], cmap="seismic", origin="lower", vmin=-vmax, vmax=vmax)
    ax.set_title("pressure")
    fig.colorbar(image, ax=ax)

    def update(frame_index: int):
        image.set_array(data[frame_index])
        ax.set_xlabel(f"frame {frame_index + 1}/{data.shape[0]}")
        return (image,)

    ani = animation.FuncAnimation(fig, update, frames=data.shape[0], blit=True)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    ani.save(output_path, fps=fps)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Render an HDF5 wavefield to MP4.")
    parser.add_argument("input", nargs="?", default="output/wavefield.h5")
    parser.add_argument("output", nargs="?", default="output/wavefield.mp4")
    parser.add_argument("--field", default="u")
    parser.add_argument("--fps", type=int, default=30)
    args = parser.parse_args()

    render_hdf5_to_mp4(args.input, args.output, args.field, args.fps)


if __name__ == "__main__":
    main()

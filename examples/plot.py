"""Plot small open-boundary examples, not a critical-exponent estimate."""
import argparse
import csv
import tomllib
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output already exists")
    with (args.input / "metadata.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    with (args.input / "trace.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.5), layout="constrained")
    for index, L in enumerate(metadata["sizes"]):
        color = ("#008875", "#702e96")[index % 2]
        subset = [r for r in rows if int(float(r["L"])) == L]
        layers = sorted({int(float(r["layer"])) for r in subset})
        if len(layers) != metadata["layers_per_site"] * L:
            raise ValueError("missing layers")
        expected_ids = set(range(1, metadata["trajectories_per_size"] + 1))
        for layer in layers:
            ids = [int(float(r["trajectory"])) for r in subset if int(float(r["layer"])) == layer]
            if len(ids) != len(expected_ids) or set(ids) != expected_ids:
                raise ValueError("missing or duplicate trajectories")
        for ax, key in zip(axes, ("surprisal_nats", "surprisal_per_site_layer")):
            means = [np.mean([float(r[key]) for r in subset if int(float(r["layer"])) == t]) for t in layers]
            if not layers or not np.isfinite(means).all():
                raise ValueError("invalid trace")
            ax.plot(np.array(layers) / L, means, label=f"L = {L}", color=color, linewidth=1.8)
            ax.set_xlabel("Layer / L")
            ax.spines[["top", "right"]].set_visible(False)
            ax.legend(frameon=False)
    axes[0].set_ylabel("Mean record surprisal (nats)")
    axes[1].set_ylabel("Mean surprisal / (L x layer)")
    fig.suptitle(f"Open boundaries: p = {metadata['p']:g}, "
                 f"{metadata['trajectories_per_size']} trajectories per size", fontsize=12)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("xb") as stream:
        fig.savefig(stream, format="png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()

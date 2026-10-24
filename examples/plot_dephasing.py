"""Plot paired Monte Carlo errors against full density-matrix references."""
import argparse
import csv
from pathlib import Path
import tomllib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot(folder, output):
    folder, output = Path(folder), Path(output)
    if output.exists():
        raise FileExistsError(output)
    metadata = tomllib.loads((folder / "metadata.toml").read_text())
    if not metadata["config"]["density_reference"]:
        raise ValueError("This plot requires density_reference=true")
    with (folder / "convergence.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    sizes = sorted({int(row["L"]) for row in rows})
    fig, axes = plt.subplots(1, len(sizes), figsize=(3.8 * len(sizes), 3.6), squeeze=False)
    for ax, L in zip(axes[0], sizes):
        ax.axhline(0, color="0.3", lw=1)
        for kind, color, marker in [("haar", "#7142a0", "o"), ("product", "#128274", "s")]:
            selected = [r for r in rows if int(r["L"]) == L and r["initial_state"] == kind]
            if not selected:
                continue
            values = np.array([[float(r[k]) for k in ("tau", "paired_error", "paired_se")] for r in selected])
            if not np.isfinite(values).all():
                raise ValueError("Non-finite convergence data")
            ax.errorbar(*values[:, :2].T, yerr=values[:, 2], color=color, marker=marker,
                        markersize=3, linewidth=1, capsize=2, label=kind.capitalize())
        ax.set(xscale="log", title=f"L = {L}", xlabel="Hidden histories per record")
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        ax.spines[["top", "right"]].set_visible(False)
    axes[0, 0].set_ylabel(r"$f_{\rm MC} - f_{\rm density}$ (nats / site / layer)")
    axes[0, -1].legend(frameon=False)
    config = metadata["config"]
    fig.suptitle(f"Paired density-matrix benchmark: p = {config['p']}, {config['layers']} layers, "
                 f"{config['records_per_state']} records per initial-state ensemble", fontsize=11)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    plot(args.folder, args.output)

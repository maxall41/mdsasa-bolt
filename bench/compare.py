# Copyright (C) 2025 Maxwell J. Campbell
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr


def load_results(filename: str) -> np.ndarray:
    """Load SASA results from a text file."""
    with Path.open(Path("bench") / Path(filename)) as f:
        results = [float(line.strip()) for line in f]
    return np.array(results)


def main() -> None:
    """Compare SASA results from old and new implementations."""
    # Load results from both implementations
    old_results = load_results("mdakit_results.txt")
    new_results = load_results("rustsasa_results.txt")

    # Calculate Pearson correlation
    correlation, p_value = pearsonr(old_results, new_results)

    # Calculate statistics
    rmse = np.sqrt(np.mean((new_results - old_results) ** 2))

    # Create figure with proper axis styling
    fig, ax = plt.subplots(figsize=(8, 8))

    # Style parameters
    fontsize = 14
    linewidth = 2.0
    colors = ["#4477AA"]  # Color-blind friendly blue

    # Remove top and right spines
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_linewidth(linewidth)
    ax.spines["bottom"].set_linewidth(linewidth)
    ax.xaxis.set_tick_params(width=linewidth, length=8, direction="out")
    ax.yaxis.set_tick_params(width=linewidth, length=8, direction="out")

    # Scatter plot with professional styling
    ax.scatter(
        old_results,
        new_results,
        alpha=0.6,
        s=50,
        color=colors[0],
        edgecolor="none",
    )

    # Add diagonal line (perfect correlation)
    min_val = min(old_results.min(), new_results.min())
    max_val = max(old_results.max(), new_results.max())
    ax.plot([min_val, max_val], [min_val, max_val], color="gray", linestyle="--", linewidth=2)

    # Labels and title
    ax.set_xlabel("mdakit-sasa", fontsize=fontsize + 2)
    ax.set_ylabel("RustSASA", fontsize=fontsize + 2)
    ax.set_title("SASA Results Comparison", fontsize=fontsize + 2, pad=10)

    # Tick styling
    ax.tick_params(axis="both", which="major", labelsize=fontsize)
    ax.ticklabel_format(style="sci", axis="both", scilimits=(0, 0))
    ax.xaxis.get_offset_text().set_fontsize(fontsize - 2)
    ax.yaxis.get_offset_text().set_fontsize(fontsize - 2)

    # Grid
    ax.grid(True, alpha=0.3, linestyle="--")

    # Equal aspect ratio
    ax.set_aspect("equal", adjustable="box")

    # Statistics box
    stats_text = f"n = {len(old_results)}\nr = {correlation:.4f}\np = {p_value:.2e}\nRMSE = {rmse:.2f}"
    ax.text(
        0.05,
        0.95,
        stats_text,
        transform=ax.transAxes,
        va="top",
        fontsize=fontsize - 2,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="gray", alpha=0.8),
    )

    plt.tight_layout()
    plt.savefig("bench/sasa_comparison.png", dpi=300, bbox_inches="tight")

    # Print summary statistics


if __name__ == "__main__":
    main()

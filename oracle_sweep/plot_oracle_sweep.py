#!/usr/bin/env python3
"""Plot saved oracle-sweep results without rerunning the numerical sweep."""

from argparse import ArgumentParser
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_CSV = SCRIPT_DIR / "oracle_sweep_output" / "oracle_sweep_results.csv"

REQUIRED_COLUMNS = {
    "mu",
    "ours_exante_length",
    "lei_exante_length",
    "hpd_stage1_length",
    "lei_stage1_length",
}


def make_figure(df: pd.DataFrame, output_dir: Path) -> tuple[Path, Path]:
    """Draw the four saved population-length curves on a single axis."""
    fig, ax = plt.subplots(
        figsize=(7, 5.2),
        constrained_layout=True,
    )

    membership_color = "#0072B2"
    lei_color = "#D55E00"

    ax.plot(
        df["mu"],
        df["ours_exante_length"],
        color=membership_color,
        linestyle="-",
        linewidth=2.2,
        label="Ex-ante membership level set",
    )
    ax.plot(
        df["mu"],
        df["lei_exante_length"],
        color=lei_color,
        linestyle="-",
        linewidth=2.2,
        label="Ex-ante Lei-Candès nested interval",
    )
    ax.plot(
        df["mu"],
        df["hpd_stage1_length"],
        color=membership_color,
        linestyle="--",
        linewidth=2.2,
        label="First-stage density-level set (HPD)",
    )
    ax.plot(
        df["mu"],
        df["lei_stage1_length"],
        color=lei_color,
        linestyle="--",
        linewidth=2.2,
        label="First-stage Lei-Candès central-quantile interval",
    )

    for mu_ref in (3.0, 6.0):
        ax.axvline(
            mu_ref,
            color="0.65",
            linestyle=":",
            linewidth=1.2,
            zorder=0,
        )

    ax.set_xlim(1.5, 8.5)
    ax.set_ylim(7, 45)
    ax.set_xlabel(
        r"Mixture location parameter $\mu$ ($\sigma=1$)",
        fontsize=12,
        labelpad=8,
    )
    ax.set_ylabel(
        "Population set length",
        fontsize=12,
        labelpad=8,
    )
    ax.tick_params(axis="both", labelsize=13)
    ax.legend(
        loc="upper left",
        fontsize=10.5,
        ncol=1,
        borderpad=0.7,
        labelspacing=0.5,
        handlelength=2.8,
    )
    ax.grid(axis="y", color="0.90", linewidth=0.7)

    output_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = output_dir / "oracle_sweep_figure.pdf"
    png_path = output_dir / "oracle_sweep_figure.png"

    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return pdf_path, png_path


def parse_args():
    parser = ArgumentParser(
        description="Plot an existing oracle_sweep_results.csv file."
    )
    parser.add_argument(
        "csv_path",
        nargs="?",
        type=Path,
        default=DEFAULT_CSV,
        help=(
            "Path to oracle_sweep_results.csv. By default, use "
            "oracle_sweep_output/oracle_sweep_results.csv next to this script."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    csv_path = args.csv_path.expanduser().resolve()

    if not csv_path.is_file():
        raise FileNotFoundError(f"Saved sweep results not found: {csv_path}")

    df = pd.read_csv(csv_path)
    missing_columns = sorted(REQUIRED_COLUMNS.difference(df.columns))
    if missing_columns:
        raise ValueError(
            "Saved CSV is missing required columns: "
            + ", ".join(missing_columns)
        )

    pdf_path, png_path = make_figure(df, csv_path.parent)

    print("Read saved results without rerunning the oracle sweep:")
    print(f"  {csv_path}")
    print("Saved:")
    print(f"  {pdf_path}")
    print(f"  {png_path}")


if __name__ == "__main__":
    main()

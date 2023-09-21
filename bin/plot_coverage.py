#!/usr/bin/env python3

import argparse
import csv


import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import os
import platform
import pandas
import matplotlib.ticker as mtick


def plot_coverage(depths: pd.DataFrame, sample_name: str, threshold: int, rolling_window,) -> matplotlib.figure.Figure:
    """
    """
    percent_coverage_above_threshold = (
        sum(1 if x > threshold else 0 for x in depths.depth)
        / depths.shape[0]
        * 100
    )

    coverage_plot = plt.figure(figsize=(64, 8))
    line_plot = sns.lineplot(data=depths, x="pos", y="depth", linewidth=0.5)
    line_plot.set_xticks(range(0, depths.shape[0], 100000), minor=True)
    threshold_line = plt.axhline(y=threshold, color="red", linestyle="--", linewidth=0.5)
    threshold_line.set_label("Threshold")
    plt.legend(loc="upper right")
    plt.title(
        f"Percent bases with coverage above {threshold}X: {percent_coverage_above_threshold: .1f}% | Rolling window: {rolling_window} nt"
    )
    plt.suptitle(f"Ref: {depths.iloc[0].chrom} | Sample: {sample_name}")
    plt.close()

    return coverage_plot


def main(args):

    all_depths = pd.read_csv(args.input, sep="\t")
    rolling_window_depths = all_depths.assign(
            depth=lambda x: x.depth.rolling(args.rolling_window).mean()
    )
    coverage_plot = plot_coverage(rolling_window_depths, args.sample_id, args.threshold, args.rolling_window)
    coverage_plot.savefig(f"{args.sample_id}_coverage_plot.svg")
    coverage_plot.savefig(f"{args.sample_id}_coverage_plot.png")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot coverage')
    parser.add_argument('-i', '--input', help='Input file', required=True)
    parser.add_argument('-s', '--sample-id', default="Sample", help='Sample ID')
    parser.add_argument('-t', '--threshold', type=int, default=10, help='Threshold for depth of coverage')
    parser.add_argument('-r', '--rolling-window', type=int, default=1000, help='Rolling window for coverage')
    args = parser.parse_args()
    main(args)

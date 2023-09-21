#!/usr/bin/env python3

import argparse
import csv

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

import pandas as pd

import seaborn as sns


def plot_coverage(depths: pd.DataFrame, sample_name: str="Sample", threshold: int=10, rolling_window: int=100, y_limit: int=500, log_scale: bool=False) -> matplotlib.figure.Figure:
    """
    Plot coverage

    :param depths: DataFrame of depths (columns: chrom, pos, depth)
    :param sample_name: Sample name
    :param threshold: Threshold for depth of coverage
    :param rolling_window: Rolling window for coverage
    :param log_scale: Log scale y-axis (default: False)
    :return: matplotlib figure
    :rtype: matplotlib.figure.Figure
    """
    percent_coverage_above_threshold = (
        sum(1 if x > threshold else 0 for x in depths.depth)
        / depths.shape[0]
        * 100
    )

    coverage_plot = plt.figure(figsize=(64, 8))
    line_plot = sns.lineplot(data=depths, x="pos", y="depth", linewidth=0.5)
    if log_scale:
        plt.yscale("log")
        plt.ylim(1, y_limit)
        plt.ylabel("Depth of coverage (log scale)")
    else:
        plt.ylim(0, y_limit)
        plt.ylabel("Depth of coverage")
    plt.xlabel("Position")
    line_plot.set_xticks(range(0, depths.shape[0], 100000), minor=True)
    threshold_line = plt.axhline(y=threshold, color="red", linestyle="--", linewidth=0.5)
    plt.title(
        f"Percent bases with coverage above {threshold}X: {percent_coverage_above_threshold: .1f}% | Rolling window: {rolling_window} nt"
    )
    plt.suptitle(f"Ref: {depths.iloc[0].chrom} | Sample: {sample_name}")
    plt.close()

    return coverage_plot


def main(args):

    all_depths = pd.read_csv(args.input, sep="\t")
    rolling_window_depths = all_depths.assign(
        depth=lambda x: x.depth.rolling(args.rolling_window, center=True).mean()
    )
    coverage_plot = plot_coverage(rolling_window_depths, args.sample_id, args.threshold, args.rolling_window, args.y_limit, args.log_scale)
    coverage_plot.savefig(f"{args.sample_id}_coverage_plot.svg")
    coverage_plot.savefig(f"{args.sample_id}_coverage_plot.png")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot coverage')
    parser.add_argument('-i', '--input', help='Input file', required=True)
    parser.add_argument('-s', '--sample-id', default="Sample", help='Sample ID')
    parser.add_argument('-t', '--threshold', type=int, default=10, help='Threshold for depth of coverage')
    parser.add_argument('-r', '--rolling-window', type=int, default=100, help='Rolling window for coverage')
    parser.add_argument('--log-scale', action='store_true', help='Log scale y-axis')
    parser.add_argument('--y-limit', type=int, default=500, help='Maximum y-axis value in plot (default: 500)')
    args = parser.parse_args()
    main(args)

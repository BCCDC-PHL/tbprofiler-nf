#!/usr/bin/env python3

import argparse
import csv

from pathlib import Path
from typing import Optional

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import pandas as pd
import seaborn as sns

from matplotlib.patches import Rectangle

def parse_bed(bed_path: Path) -> dict:
    """
    Parse BED file

    :param bed: BED file
    :type bed: Path
    :return: Dictionary of genes by name
    :rtype: dict
    """
    genes_by_name = {}
    with open(bed_path, 'r') as f:
        for line in f:
            line_split = line.strip().split('\t')
            chrom = line_split[0]
            start = int(line_split[1])
            end = int(line_split[2])
            name = line_split[3]
            strand = line_split[5]
            genes_by_name[name] = {
                'name': name,
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand,
            }

    return genes_by_name


def plot_coverage(depths: pd.DataFrame, resistance_genes: dict={}, sample_name: str="Sample", threshold: int=10, rolling_window: int=100, y_limit: int=500, log_scale: bool=False) -> matplotlib.figure.Figure:
    """
    Plot coverage

    :param depths: DataFrame of depths (columns: chrom, pos, depth)
    :type depths: pandas.DataFrame
    :param resistance_genes: Dictionary of resistance genes by name
    :type resistance_genes: dict
    :param sample_name: Sample name
    :type sample_name: str
    :param threshold: Threshold for depth of coverage
    :type threshold: int
    :param rolling_window: Rolling window for coverage
    :type rolling_window: int
    :param log_scale: Log scale y-axis (default: False)
    :type log_scale: bool
    :return: matplotlib figure
    :rtype: matplotlib.figure.Figure
    """
    percent_coverage_above_threshold = (
        sum(1 if x > threshold else 0 for x in depths.depth)
        / depths.shape[0]
        * 100
    )

    if resistance_genes:
        coverage_plot, (ax_depth, ax_genes) = plt.subplots(2, gridspec_kw={'height_ratios': [10, 1]}, figsize=(64, 8), sharex=True)
        depth_plot = sns.lineplot(data=depths, x="pos", y="depth", linewidth=0.5, ax=ax_depth)
        threshold_line = ax_depth.axhline(y=threshold, color="red", linestyle="--", linewidth=0.5)
        ax_depth.set_ylabel("Depth of coverage")
        ax_depth.set_title(f"Percent bases with coverage above {threshold}X: {percent_coverage_above_threshold: .1f}% | Rolling window: {rolling_window} nt")
        coverage_plot.suptitle(f"Ref: {depths.iloc[0].chrom} | Sample: {sample_name}")
        if log_scale:
            ax_depth.set_yscale('log')
            ax_depth.set_ylim(1, y_limit)
            ax_depth.set_ylabel("Depth of coverage (log scale)")

        prev_gene_end = 0
        prev_gene_name_vertical_padding = 0
        for gene_name, gene in resistance_genes.items():
            print(gene)
            if gene['strand'] == '+':
                gene_color = 'blue'
                gene_vertical_padding = 3
            else:
                gene_color = 'red'
                gene_vertical_padding = 2
            gene_rectangle = Rectangle((gene['start'], gene_vertical_padding), 1000, 1, color=gene_color, alpha=0.5)
            rx, ry = gene_rectangle.get_xy()
            cx = rx + gene_rectangle.get_width() / 2.0
            cy = ry + gene_rectangle.get_height() / 2.0
            ax_genes.add_patch(gene_rectangle)
            if gene['start'] - prev_gene_end < 10000:
                if gene['strand'] == '+':
                    gene_name_vertical_padding = prev_gene_name_vertical_padding + 1
                else:
                    gene_name_vertical_padding = prev_gene_name_vertical_padding - 1
            else:
                if gene['strand'] == '+':
                    gene_name_vertical_padding = 1
                else:
                    gene_name_vertical_padding = -1
            ax_genes.annotate(gene_name, (cx, cy + gene_name_vertical_padding), color='black', fontsize=6, ha='center', va='center')
            prev_gene_end = gene['start'] + 1000
            prev_gene_name_vertical_padding = gene_name_vertical_padding

        ax_genes.set_ylim(top=6, bottom=0)
        ax_genes.yaxis.set_visible(False)
        ax_genes.xaxis.set_visible(False)
    else:
        coverage_plot = plt.figure(figsize=(64, 8))
        line_plot = sns.lineplot(data=depths, x="pos", y="depth", linewidth=0.5)
        plt.axhline(y=threshold, color="red", linestyle="--", linewidth=0.5)
        plt.ylabel("Depth of coverage")
        if log_scale:
            plt.yscale("log")
            plt.ylim(1, y_limit)
            plt.ylabel("Depth of coverage (log scale)")
        else:
            plt.ylim(0, y_limit)
            plt.ylabel("Depth of coverage")
        plt.title(
            f"Percent bases with coverage above {threshold}X: {percent_coverage_above_threshold: .1f}% | Rolling window: {rolling_window} nt"
        )
        plt.suptitle(f"Ref: {depths.iloc[0].chrom} | Sample: {sample_name}")

    # line_plot.set_xticks(range(0, depths.shape[0], 100000), minor=True)
    
    plt.close()

    return coverage_plot


def main(args):

    all_depths = pd.read_csv(args.input, sep="\t")
    rolling_window_depths = all_depths.assign(
        depth=lambda x: x.depth.rolling(args.rolling_window, center=True).mean()
    )

    resistance_genes = {}
    if args.genes is not None:
        resistance_genes = parse_bed(args.genes)

    coverage_plot = plot_coverage(rolling_window_depths, resistance_genes, args.sample_id, args.threshold, args.rolling_window, args.y_limit, args.log_scale)
    coverage_plot.savefig(f"{args.sample_id}_coverage_plot.svg")
    coverage_plot.savefig(f"{args.sample_id}_coverage_plot.png")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot coverage')
    parser.add_argument('-i', '--input', help='Input file', required=True)
    parser.add_argument('-g', '--genes', type=Path, help='Resistance genes (bed file)')
    parser.add_argument('-s', '--sample-id', default="Sample", help='Sample ID')
    parser.add_argument('-t', '--threshold', type=int, default=10, help='Threshold for depth of coverage')
    parser.add_argument('-r', '--rolling-window', type=int, default=100, help='Rolling window for coverage')
    parser.add_argument('--log-scale', action='store_true', help='Log scale y-axis')
    parser.add_argument('--y-limit', type=int, default=500, help='Maximum y-axis value in plot (default: 500)')
    args = parser.parse_args()
    main(args)

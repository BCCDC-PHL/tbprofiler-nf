#!/usr/bin/env python3

import csv
import argparse

# calculate mean depth coverage for each resistance gene
def calculate_mean_depth(gene_positions, depth):
    total_depth = 0
    total_positions = 0

    for position in range(gene_positions[0], gene_positions[1] + 1):
        if position in depth:
            total_depth += depth[position]
            total_positions += 1
    
    mean_depth = total_depth / total_positions
    return round(mean_depth, 3)

# calculate percent of gene above depth threshold
def calculate_percent_gene_covered(gene_positions, depth, depth_threshold):
    covered_positions = 0

    for position in range(gene_positions[0], gene_positions[1] + 1):
        if position in depth and depth[position] >= depth_threshold:
            covered_positions += 1

    total_positions = gene_positions[1] - gene_positions[0] + 1


    pct_cov = (covered_positions / total_positions) * 100
    return round(pct_cov, 3)

def main(args):
    # Read gene positions from the resistance genes bed file
    gene_data = []
    with open(args.bed, 'r') as bed_file:
        for line in bed_file:
            fields = line.strip().split()
            gene_name = fields[3]
            start_position = int(fields[1]) +1 # add +1 because bed file is 0 indexed and depth file is 1 indexed
            end_position = int(fields[2]) 
            gene_data.append((gene_name, (start_position, end_position)))

    # Read depth data from intermediate mpileup tsv file
    depth_at_position = {}
    with open(args.depth, 'r') as tsv_file:
        for line in tsv_file:
            line = line.strip()
            if not line:
                continue  # Skip empty lines
            
            chrom, pos, ref, depth = line.split('\t')
            
            #skipping header line 
            try:
                position = int(pos)
                depth = float(depth)
            except ValueError:
                continue
            
            depth_at_position[position] = depth

    # Calculate qc metrics and write results to csv
    output_file = args.output
    with open(output_file, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['gene_name', 'gene_position_start', 'gene_position_end', 'mean_depth_coverage', 'percent_of_gene_covered_above_depth_threshold'])
        
        for gene_name, gene_positions in gene_data:
            mean_depth = calculate_mean_depth(gene_positions, depth_at_position)
            percent_covered = calculate_percent_gene_covered(gene_positions, depth_at_position, args.threshold)
            
            writer.writerow([gene_name, gene_positions[0], gene_positions[1], mean_depth, percent_covered])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate resistance gene coverage qc metrics.')
    parser.add_argument('--bed', required=True, help='Path for resistance genes bed file.')
    parser.add_argument('--depth', required=True, help='Path for tsv file with depth for each position (intermediate mpileup output).')
    parser.add_argument('--threshold', type=float, default=10, help='min_depth threshold (default: 10).')
    parser.add_argument('--output',  help='Output CSV file name.')
    args = parser.parse_args()
    main(args)

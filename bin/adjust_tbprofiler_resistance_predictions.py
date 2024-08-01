#!/usr/bin/env python3


import argparse
import csv
import json
import logging
import sys

from pathlib import Path


def parse_tbprofiler_full_report(input_file: Path) -> dict:
    """
    Parse the tbprofiler full report file

    :param input_file: The input tbprofiler full report file (json)
    :type input_file: str
    :return: The parsed tbprofiler full report
    :rtype: dict
    """
    tbprofiler_full_report = {}
    with open(input_file, 'r') as f:
        tbprofiler_full_report = json.load(f)
    return tbprofiler_full_report


def parse_tbprofiler_resistance_predictions(input_file: Path) -> dict:
    """
    Parse the tbprofiler resistance predictions file

    :param input_file: The input tbprofiler resistance predictions file (csv)
    :type input_file: str
    :return: The parsed tbprofiler resistance predictions
    :rtype: dict
    """
    tbprofiler_resistance_prediction_by_drug = {}
    with open(input_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            drug = row['drug'].lower()
            if row['genotypic_resistance'] == 'R':
                row['tbprofiler_resistance_prediction'] = 'Predicted Resistant'
            elif row['genotypic_resistance'] == '':
                row['tbprofiler_resistance_prediction'] = ''
            if drug not in tbprofiler_resistance_prediction_by_drug:
                tbprofiler_resistance_prediction_by_drug[drug] = row

    return tbprofiler_resistance_prediction_by_drug


def parse_tbprofiler_resistance_mutations(input_file: Path) -> dict:
    """
    Parse the tbprofiler resistance mutations file

    :param input_file: The input tbprofiler resistance mutations file (csv)
    :type input_file: str
    :return: The parsed tbprofiler resistance mutations
    :rtype: dict
    """
    tbprofiler_resistance_mutations_by_drug = {}
    float_fields = [
        'estimated_fraction',
    ]
    with open(input_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['drug'] not in tbprofiler_resistance_mutations_by_drug:
                tbprofiler_resistance_mutations_by_drug[row['drug']] = []
            for field in float_fields:
                if row[field] == '':
                    row[field] = None
                else:
                    row[field] = float(row[field])
            if row['mutation'] != '':
                mutation = {
                    'gene': row['gene'],
                    'mutation': row['mutation'],
                    'estimated_fraction': row['estimated_fraction'],
                }
                tbprofiler_resistance_mutations_by_drug[row['drug']].append(mutation)

    return tbprofiler_resistance_mutations_by_drug


def parse_mutation_catalogue(catalogue_path):
    """
    Parse the WHO TB mutation catalogue.

    :param catalogue_path: Path to the WHO TB mutation catalogue
    :type catalogue_path: str
    :return: The parsed mutation catalogue
    :rtype: list
    """
    fieldname_translation = {
        'initial_confidence_grading_all_dataset_dataset_all': 'initial_confidence_grading_dataset_all',
        'initial_confidence_grading_who_dataset_dataset_who': 'initial_confidence_grading_dataset_who',
        'relaxed_thresholds_simulation_bdq_rv0678,_cfz_rv0678,_inh_katg,_dlm_ddn/fbia/fbib/fbic/fgd1/rv2983': 'relaxed_thresholds_simulation',
        'who_guidance_before_cat_ver1_miotto_et_al._pmid_29284687,_who_ngs_guide_2018_source_of_additional_grading_evidence': 'who_guidance_before_cat_ver1_from_miotto_et_al_pmid_29284687',
        'who_guidance_before_cat_ver1_rif_cc_guide_2021_source_of_additional_grading_evidence': 'who_guidance_before_cat_ver1_rif_cc_guide_2021',
    }

    lowercase_fields = [
        'drug',
    ]

    int_fields = [
        'tier',
        'algorithm_pass_dataset_all',
        'present_solo_sr_dataset_all',
        'present_solo_r_dataset_all',
        'present_solo_s_dataset_all',
        'present_r_dataset_all',
        'present_s_dataset_all',
        'absent_r_dataset_all',
        'absent_s_dataset_all',
        'seta_dataset_who',
        'setb_dataset_who',
        'setc_dataset_who',
        'setd1_dataset_who',
        'setd2_dataset_who',
        'v1_dataset_who',
        'literature_dataset_who',
        'algorithm_pass_dataset_who',
        'present_solo_sr_dataset_who',
        'present_solo_r_dataset_who',
        'present_solo_s_dataset_who',
        'present_r_dataset_who',
        'present_s_dataset_who',
        'absent_r_dataset_who',
        'absent_s_dataset_who',
        'present_nopheno_no_pdst',
        'changes_vs_ver1',
    ]

    float_fields = [
        'sens_dataset_all',
        'sens_lb_dataset_all',
        'sens_ub_dataset_all',
        'sens_dataset_who',
        'sens_lb_dataset_who',
        'sens_ub_dataset_who',
        'spec_dataset_all',
        'spec_lb_dataset_all',
        'spec_ub_dataset_all',
        'spec_dataset_who',
        'spec_lb_dataset_who',
        'spec_ub_dataset_who',
        'ppv_dataset_all',
        'ppv_lb_dataset_all',
        'ppv_ub_dataset_all',
        'ppv_dataset_who',
        'ppv_lb_dataset_who',
        'ppv_ub_dataset_who',
        'ppv_solo_dataset_all',
        'ppv_solo_lb_dataset_all',
        'ppv_solo_ub_dataset_all',
        'ppv_solo_dataset_who',
        'ppv_solo_lb_dataset_who',
        'ppv_solo_ub_dataset_who',
        'ppv_conditional_solo_dataset_all',
        'ppv_conditional_solo_lb_dataset_all',
        'ppv_conditional_solo_ub_dataset_all',
        'ppv_conditional_solo_dataset_who',
        'ppv_conditional_solo_lb_dataset_who',
        'ppv_conditional_solo_ub_dataset_who',
        'or_solo_dataset_all',
        'or_solo_exact_lb_dataset_all',
        'or_solo_exact_ub_dataset_all',
        'or_solo_pvalue_dataset_all',
        'or_solo_pval_rank_dataset_all',
        'or_solo_dataset_who',
        'or_solo_exact_lb_dataset_who',
        'or_solo_exact_ub_dataset_who',
        'or_solo_pvalue_dataset_who',
        'or_solo_pval_rank_dataset_who',
        'k_dataset_all',
        'k_dataset_who',
        'or_solo_fe_sig_dataset_all',
        'neutral_masked_dataset_all',
        'or_dataset_all',
        'or_exact_lb_dataset_all',
        'or_exact_ub_dataset_all',
        'or_dataset_who',
        'or_exact_lb_dataset_who',
        'or_exact_ub_dataset_who',
        'sens_solo_dataset_all',
        'sens_solo_lb_dataset_all',
        'sens_solo_ub_dataset_all',
        'sens_solo_dataset_who',
        'sens_solo_lb_dataset_who',
        'sens_solo_ub_dataset_who',
        'spec_solo_dataset_all',
        'spec_solo_lb_dataset_all',
        'spec_solo_ub_dataset_all',
        'spec_solo_dataset_who',
        'spec_solo_lb_dataset_who',
        'spec_solo_ub_dataset_who',
    ]

    percent_fields = [
        'ppv',
        'ppv_lb',
        'ppv_ub',
    ]

    boolean_fields = [
        'or_solo_fe_sig_dataset_all',
        'or_solo_fe_sig_dataset_who',
        'neutral_masked_dataset_all',
        'neutral_masked_dataset_who',
    ]

    yes_no_fields = [
        'listed_in_abridged_tables',
    ]

    silent_mutation_fields = [
        'silent_mutation',
    ]

    grading_fields = [
        'initial_confidence_grading',
        'initial_confidence_grading_dataset_all',
        'initial_confidence_grading_dataset_who',
        'final_confidence_grading',
    ]

    catalogue = []

    with open(catalogue_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            parsed_record = {}
            for k in row.keys():
                cleaned_key = k.lower().replace(' ', '_').replace('(', '').replace(')', '').replace('_|_', '_')
                if cleaned_key in fieldname_translation:
                    cleaned_key = fieldname_translation[cleaned_key]

                if cleaned_key in lowercase_fields:
                    row[k] = row[k].lower()
                elif cleaned_key in int_fields:
                    try:
                        row[k] = int(row[k])
                    except ValueError as e:
                        row[k] = None
                elif cleaned_key in float_fields:
                    try:
                        row[k] = float(row[k])
                    except ValueError as e:
                        row[k] = None
                elif cleaned_key in percent_fields:
                    try:
                        row[k] = float(row[k].rstrip('%'))
                    except ValueError as e:
                        row[k] = None
                elif cleaned_key in boolean_fields:
                    try:
                        row[k] = row[k].lower() == 'true'
                    except ValueError as e:
                        row[k] = None
                elif cleaned_key in yes_no_fields:
                    try:
                        row[k] = row[k].lower() == 'yes'
                    except ValueError as e:
                        row[k] = None
                elif cleaned_key in silent_mutation_fields:
                    try:
                        row[k] = row[k].lower() == 'silent mutation'
                    except ValueError as e:
                        row[k] = None
                elif cleaned_key in grading_fields:
                    try:
                        grading_level = int(row[k].split(')')[0])
                        grading_category = row[k].split(')')[1].strip()
                        if 'Interim' in grading_category:
                            interim_grading = True
                        else:
                            interim_grading = False
                        if 'Assoc w R' in grading_category:
                            associated_with_resistance = True
                        else:
                            associated_with_resistance = False
                    except ValueError as e:
                        grading_number = None
                        grading_category = None
                        interim_grading = None
                    row[k] = {
                        'original_value': row[k],
                        'level': grading_level,
                        'category': grading_category,
                        'interim_grading': interim_grading,
                        'associated_with_resistance': associated_with_resistance,
                    }
                
                parsed_record[cleaned_key] = row[k]

            # The 'genomic_position' field is always "(see \"Genomic_coordinates\" sheet)"
            parsed_record.pop('genomic_position', None)

            catalogue.append(parsed_record)
            
    return catalogue


def index_mutation_catalogue_by_mutation_by_gene_by_drug(mutation_catalogue):
    """
    Index the mutation catalogue by mutation, gene, and drug

    :param mutation_catalogue: The mutation catalogue
    :type mutation_catalogue: list
    :return: The indexed mutation catalogue. Nested by drug, gene, and mutation.
    :rtype: dict
    """
    catalogue_by_mutation_by_gene_by_drug = {}
    for record in mutation_catalogue:
        mutation = record['mutation']
        gene = record['gene']
        drug = record['drug']
        if drug not in catalogue_by_mutation_by_gene_by_drug:
            catalogue_by_mutation_by_gene_by_drug[drug] = {}
        if gene not in catalogue_by_mutation_by_gene_by_drug[drug]:
            catalogue_by_mutation_by_gene_by_drug[drug][gene] = {}
        catalogue_by_mutation_by_gene_by_drug[drug][gene][mutation] = record

    return catalogue_by_mutation_by_gene_by_drug


def index_full_report_resistance_variants_by_mutation_by_gene_by_drug(tbprofiler_full_report_resistance_variants):
    """
    Index the full report resistance variants by mutation, gene, and drug

    :param tbprofiler_full_report_resistance_variants: The full report resistance variants
    :type tbprofiler_full_report_resistance_variants: list
    :return: The indexed full report resistance variants. Nested by drug, gene, and mutation.
    :rtype: dict
    """
    full_report_resistance_variants_by_mutation_by_gene_by_drug = {}
    for variant in tbprofiler_full_report_resistance_variants:
        mutation = variant['change']
        gene = variant['gene_name']
        drugs = []
        for annotation in variant['annotation']:
            if annotation['type'] == 'drug_resistance' or annotation['type'] == 'who_confidence':
                drugs.append(annotation['drug'])
        if len(drugs) == 0:
            continue
        for drug in drugs:
            if drug not in full_report_resistance_variants_by_mutation_by_gene_by_drug:
                full_report_resistance_variants_by_mutation_by_gene_by_drug[drug] = {}
            if gene not in full_report_resistance_variants_by_mutation_by_gene_by_drug[drug]:
                full_report_resistance_variants_by_mutation_by_gene_by_drug[drug][gene] = {}
            full_report_resistance_variants_by_mutation_by_gene_by_drug[drug][gene][mutation] = variant

    return full_report_resistance_variants_by_mutation_by_gene_by_drug


def parse_drug_ppv_thresholds(input_file: Path) -> dict:
    """
    Parse the drug PPV thresholds file. The file should be a CSV with the following columns:
    - drug
    - ppv_threshold

    :param input_file: The input drug PPV thresholds file (csv)
    :type input_file: str
    :return: The parsed drug PPV thresholds
    :rtype: dict
    """
    drug_ppv_thresholds = {}

    with open(input_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            drug = row['drug']
            try:
                ppv_threshold = float(row['ppv_threshold'])
            except ValueError as e:
                continue
            drug_ppv_thresholds[drug] = ppv_threshold

    return drug_ppv_thresholds


def combine_resistance_probabilities(ppvs):
    """
    Calculate the combined resistance probability given a list of PPVs

    :param ppvs: The list of PPVs
    :type ppvs: list
    :return: The combined resistance probability
    :rtype: float
    """

    combined_prob = 0
    for ppv in ppvs:
        if ppv is not None and ppv > 0.0:
            if combined_prob == 0:
                combined_prob = 1
            combined_prob *= (1 - ppv)
    if combined_prob > 0:
        combined_prob = 1 - combined_prob

    return combined_prob


def collect_catalogue_info_for_mutation(tbprofiler_resistance_mutation_record, indexed_catalogue, drug):
    """
    Collect information from the mutation catalogue for a given mutation

    :param mutation: The mutation
    :type mutation: dict
    :param indexed_catalogue: The indexed mutation catalogue
    :type indexed_catalogue: dict
    :param drug: The drug name to look up in the indexed mutation catalogue
    :type drug: str
    :return: The mutation with additional information from the mutation catalogue
    """
    gene = tbprofiler_resistance_mutation_record.get('gene', None)
    if gene is None or gene not in indexed_catalogue[drug]:
        logging.warning(f"Gene {gene} not found in mutation catalogue for drug {drug}")

        return tbprofiler_resistance_mutation_record
    
    mutation = tbprofiler_resistance_mutation_record.get('mutation', None)
    if mutation is None or mutation not in indexed_catalogue[drug][gene]:
        logging.warning(f"Mutation {mutation} not found in mutation catalogue for gene {gene} and drug {drug}")

        return tbprofiler_resistance_mutation_record


    mutation_catalogue_record = indexed_catalogue[drug][gene][mutation]
    catalogue_dataset = "all"
    tbprofiler_resistance_mutation_record['catalogue_ppv_solo_estimate'] = mutation_catalogue_record.get('ppv_dataset_all', None)
    tbprofiler_resistance_mutation_record['catalogue_ppv_solo_lower_bound'] = mutation_catalogue_record.get('ppv_lb_dataset_all', None)
    tbprofiler_resistance_mutation_record['catalogue_ppv_solo_upper_bound'] = mutation_catalogue_record.get('ppv_ub_dataset_all', None)
    tbprofiler_resistance_mutation_record['catalogue_dataset'] = catalogue_dataset
    tbprofiler_resistance_mutation_record['catalogue_final_confidence_grading'] = mutation_catalogue_record.get('final_confidence_grading', {}).get('category', None)

    #print(json.dumps(tbprofiler_resistance_mutation_record, indent=2))
    #exit()

    return tbprofiler_resistance_mutation_record


def collect_full_report_info_for_resistance_mutation(tbprofiler_resistance_mutation_record, indexed_tbprofiler_full_report_resistance_variants, drug):
    """
    Collect information from the full report for a given mutation

    :param tbprofiler_resistance_mutation_record: The mutation
    :type tbprofiler_resistance_mutation_record: dict
    :param indexed_tbprofiler_full_report_resistance_variants: The indexed full report resistance variants
    :type indexed_tbprofiler_full_report_resistance_variants: dict
    :param drug: The drug name to look up in the indexed full report resistance variants
    :type drug: str
    :return: The mutation with additional information from the full report
    :rtype: dict
    """
    gene = tbprofiler_resistance_mutation_record.get('gene', None)
    if gene not in indexed_tbprofiler_full_report_resistance_variants[drug]:
        logging.warning(f"Gene {gene} not found in full report resistance variants for drug {drug}")
        return tbprofiler_resistance_mutation_record
    
    mutation = tbprofiler_resistance_mutation_record.get('mutation', None)
    if mutation not in indexed_tbprofiler_full_report_resistance_variants[drug][gene]:
        logging.warning(f"Mutation {mutation} not found in full report resistance variants for gene {gene} and drug {drug}")
        return tbprofiler_resistance_mutation_record

    full_report_variant = indexed_tbprofiler_full_report_resistance_variants[drug][gene][mutation]
    # print(json.dumps(full_report_variant, indent=2))
    # exit()
    
    for annotation in full_report_variant['annotation']:
        if annotation['type'] == 'drug_resistance' or annotation['type'] == 'who_confidence':
            if annotation['drug'] == drug:
                tbprofiler_resistance_mutation_record['drug'] = drug
                tbprofiler_resistance_mutation_record['ref_genome_accession'] = "NC_000962.3"
                tbprofiler_resistance_mutation_record['position'] = full_report_variant.get('pos', None)
                tbprofiler_resistance_mutation_record['ref_allele'] = full_report_variant.get('ref', None)
                tbprofiler_resistance_mutation_record['alt_allele'] = full_report_variant.get('alt', None)
                tbprofiler_resistance_mutation_record['depth_coverage'] = full_report_variant.get('depth', None)
                tbprofiler_resistance_mutation_record['alt_freq'] = full_report_variant.get('freq', None)
                tbprofiler_resistance_mutation_record['tbprofiler_filter'] = full_report_variant.get('filter', None)
                tbprofiler_resistance_mutation_record['mutation_type'] = full_report_variant.get('type', None)
                tbprofiler_resistance_mutation_record['tbprofiler_confidence'] = annotation.get('confidence', None)
                tbprofiler_resistance_mutation_record['tbprofiler_source'] = annotation.get('source', None)
                tbprofiler_resistance_mutation_record['tbprofiler_comment'] = annotation.get('comment', None)
                tbprofiler_resistance_mutation_record['tbprofiler_variant_type'] = 'resistance'
                break

    return tbprofiler_resistance_mutation_record


def format_full_report_info_for_other_variant(tbprofiler_full_report_other_variant):
    """
    Collect and format relevant information from the full report for a given 'other' variant

    :param tbprofiler_full_report_other_variant: The variant
    :type tbprofiler_full_report_other_variant: dict
    :return: The variant with additional information from the full report, indexed by drug
    :rtype: dict
    """
    formatted_variant_by_drug = {}
    drugs = []
    for annotation in tbprofiler_full_report_other_variant['annotation']:
        if annotation['type'] == 'drug_resistance' or annotation['type'] == 'who_confidence':
            drug = annotation.get('drug', None)
            if drug is not None:
                drugs.append(drug)
                confidence = annotation.get('confidence', None)
                source = annotation.get('source', None)
                comment = annotation.get('comment', None)
                formatted_variant_by_drug[drug] = {
                    'drug': drug,
                    'tbprofiler_confidence': confidence,
                    'tbprofiler_source': source,
                    'tbprofiler_comment': comment,
                    'tbprofiler_variant_type': 'other',
                }

    if len(drugs) == 0:
        return formatted_variant_by_drug

    for drug in drugs:
        formatted_variant = formatted_variant_by_drug[drug]
        formatted_variant['drug'] = drug
        formatted_variant['ref_genome_accession'] = "NC_000962.3"
        formatted_variant['gene'] = tbprofiler_full_report_other_variant.get('gene_name', None)
        formatted_variant['mutation'] = tbprofiler_full_report_other_variant.get('change', None)
        formatted_variant['position'] = tbprofiler_full_report_other_variant.get('pos', None)
        formatted_variant['ref_allele'] = tbprofiler_full_report_other_variant.get('ref', None)
        formatted_variant['alt_allele'] = tbprofiler_full_report_other_variant.get('alt', None)
        formatted_variant['depth_coverage'] = tbprofiler_full_report_other_variant.get('depth', None)
        formatted_variant['alt_freq'] = tbprofiler_full_report_other_variant.get('freq', None)
        formatted_variant['tbprofiler_filter'] = tbprofiler_full_report_other_variant.get('filter', None)
        formatted_variant['mutation_type'] = tbprofiler_full_report_other_variant.get('type', None)

    return formatted_variant_by_drug


def parse_adjusted_resistance_mutations(input_file: Path) -> dict:
    """
    Parse the adjusted resistance mutations file

    :param input_file: The input adjusted resistance mutations file (csv)
    :type input_file: str
    :return: The parsed adjusted resistance mutations
    :rtype: dict
    """
    adjusted_resistance_mutations_by_drug = {}

    int_fields = [
        'position',
    ]

    float_fields = [
        'depth_coverage',
        'alt_freq',
        'catalogue_ppv_solo_estimate',
        'catalogue_ppv_solo_lower_bound',
        'catalogue_ppv_solo_upper_bound',
    ]

    with open(input_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            for field in int_fields:
                if row[field] == '':
                    row[field] = None
                else:
                    try:
                        row[field] = int(row[field])
                    except ValueError as e:
                        row[field] = None
            for field in float_fields:
                if row[field] == '':
                    row[field] = None
                else:
                    try:
                        row[field] = float(row[field])
                    except ValueError as e:
                        row[field] = None

            drug = row['drug']
            if drug not in adjusted_resistance_mutations_by_drug:
                adjusted_resistance_mutations_by_drug[drug] = []
            adjusted_resistance_mutations_by_drug[drug].append(row)

    return adjusted_resistance_mutations_by_drug


def main(args):

    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_format)

    mutation_catalogue = parse_mutation_catalogue(args.input_who_mutation_catalogue)

    indexed_catalogue = index_mutation_catalogue_by_mutation_by_gene_by_drug(mutation_catalogue)

    tbprofiler_resistance_mutations_by_drug = parse_tbprofiler_resistance_mutations(args.input_tbprofiler_resistance_mutations)

    tbprofiler_resistance_prediction_by_drug = parse_tbprofiler_resistance_predictions(args.input_tbprofiler_resistance_predictions)

    tbprofiler_full_report = parse_tbprofiler_full_report(args.input_tbprofiler_full_report)
    tbprofiler_version = tbprofiler_full_report.get('pipeline', {}).get('software_version', None)
    tbdb_commit = tbprofiler_full_report.get('pipeline', {}).get('db_version', {}).get('commit', None)

    tbprofiler_full_report_resistance_variants = tbprofiler_full_report.get('dr_variants', [])

    # Index resistance variants ("dr_variants") from tbprofiler full json report by drug, gene, and mutation
    indexed_tbprofiler_full_report_resistance_variants = index_full_report_resistance_variants_by_mutation_by_gene_by_drug(tbprofiler_full_report_resistance_variants)

    # Index other variants ("other_variants") from tbprofiler full json report by drug, gene, and mutation
    indexed_tbprofiler_full_report_other_variants = {}
    tbprofiler_full_report_other_variants = tbprofiler_full_report.get('other_variants', [])
    for variant in tbprofiler_full_report_other_variants:
        formatted_variant_by_drug = format_full_report_info_for_other_variant(variant)
        for drug, formatted_variant in formatted_variant_by_drug.items():
            if drug not in indexed_tbprofiler_full_report_other_variants:
                indexed_tbprofiler_full_report_other_variants[drug] = {}
            gene = formatted_variant['gene']
            if gene not in indexed_tbprofiler_full_report_other_variants[drug]:
                indexed_tbprofiler_full_report_other_variants[drug][gene] = {}
            mutation = formatted_variant['mutation']
            formatted_variant = collect_catalogue_info_for_mutation(formatted_variant, indexed_catalogue, drug)
            formatted_variant['tbprofiler_version'] = tbprofiler_version
            formatted_variant['tbdb_commit'] = tbdb_commit
            indexed_tbprofiler_full_report_other_variants[drug][gene][mutation] = formatted_variant

    tbprofiler_other_variants_by_drug = {}
    for drug, variants_by_gene in indexed_tbprofiler_full_report_other_variants.items():
        for gene, variants in variants_by_gene.items():
            for mutation, variant in variants.items():
                if drug not in tbprofiler_other_variants_by_drug:
                    tbprofiler_other_variants_by_drug[drug] = []
                tbprofiler_other_variants_by_drug[drug].append(variant)

    drug_ppv_thresholds = {}
    if args.drug_ppv_thresholds:
        drug_ppv_thresholds = parse_drug_ppv_thresholds(args.drug_ppv_thresholds)

    for drug, tbprofiler_resistance_mutations in tbprofiler_resistance_mutations_by_drug.items():
        drug = drug.lower()
        if drug not in indexed_catalogue:
            logging.warning(f"Drug {drug} not found in mutation catalogue")
            continue
        for idx, tbprofiler_resistance_mutation in enumerate(tbprofiler_resistance_mutations):
            resistance_mutation_with_catalogue_info = collect_catalogue_info_for_mutation(tbprofiler_resistance_mutation, indexed_catalogue, drug)
            tbprofiler_resistance_mutations[idx] = resistance_mutation_with_catalogue_info

    for drug, tbprofiler_resistance_mutations in tbprofiler_resistance_mutations_by_drug.items():
        drug = drug.lower()
        if drug not in indexed_tbprofiler_full_report_resistance_variants:
            logging.warning(f"Drug {drug} not found in full report resistance variants")
            continue
        for idx, tbprofiler_resistance_mutation in enumerate(tbprofiler_resistance_mutations):
            resistance_mutation_with_full_report_info = collect_full_report_info_for_resistance_mutation(tbprofiler_resistance_mutation, indexed_tbprofiler_full_report_resistance_variants, drug)
            tbprofiler_resistance_mutations[idx] = resistance_mutation_with_full_report_info

    adjusted_resistance_mutation_output_fieldnames = [
        'sample_id',
        'drug',
        'gene',
        'mutation',
        'mutation_type',
        'ref_genome_accession',
        'position',
        'ref_allele',
        'alt_allele',
        'depth_coverage',
        'alt_freq',
        'tbprofiler_filter',
        'catalogue_dataset',
        'catalogue_ppv_solo_estimate',
        'catalogue_ppv_solo_lower_bound',
        'catalogue_ppv_solo_upper_bound',
        'catalogue_final_confidence_grading',
        'tbprofiler_variant_type',
        'tbprofiler_confidence',
        'tbprofiler_source',
        'tbprofiler_comment',
        'tbprofiler_version',
        'tbdb_commit',
    ]

    rounded_fields = [
        'alt_freq',
        'estimated_fraction',
        'catalogue_ppv_solo_estimate',
        'catalogue_ppv_solo_lower_bound',
        'catalogue_ppv_solo_upper_bound',
    ]
    with open(args.output_adjusted_resistance_mutations, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=adjusted_resistance_mutation_output_fieldnames, extrasaction='ignore')
        writer.writeheader()
        for drug, tbprofiler_resistance_mutations in tbprofiler_resistance_mutations_by_drug.items():
            for mutation in tbprofiler_resistance_mutations:
                for field in rounded_fields:
                    if field in mutation:
                        mutation[field] = round(mutation[field], 6)
                mutation['sample_id'] = args.sample_id
                mutation['tbprofiler_version'] = tbprofiler_version
                mutation['tbdb_commit'] = tbdb_commit
                writer.writerow(mutation)

        for drug, tbprofiler_other_variants in tbprofiler_other_variants_by_drug.items():
            for variant in tbprofiler_other_variants:
                for field in rounded_fields:
                    if field in variant:
                        variant[field] = round(variant[field], 6)
                variant['sample_id'] = args.sample_id
                writer.writerow(variant)

    adjusted_resistance_and_other_mutations_by_drug = parse_adjusted_resistance_mutations(args.output_adjusted_resistance_mutations)

    adjusted_resistance_prediction_output_fieldnames = [
        'sample_id',
        'drug',
        'total_num_putative_drug_related_mutations_detected',
        'total_num_putative_drug_related_mutations_with_ppv_values',
        'sum_all_putative_drug_related_mutation_ppvs',
        'combined_all_putative_drug_related_mutation_ppvs',
        'num_tbprofiler_resistance_mutations_detected',
        'num_tbprofiler_resistance_mutations_with_ppv_values',
        'sum_tbprofiler_resistance_mutation_ppvs',
        'combined_tbprofiler_resistance_mutation_ppvs',
        'ppv_threshold',
        'tbprofiler_resistance_prediction',
        'adjusted_resistance_prediction',
        'tbprofiler_version',
        'tbdb_commit',
    ]

    rounded_fields = [
        'sum_all_putative_drug_related_mutation_ppvs',
        'combined_all_putative_drug_related_mutation_ppvs',
        'sum_tbprofiler_resistance_mutation_ppvs',
        'combined_tbprofiler_resistance_mutation_ppvs',
    ]

    adjusted_resistance_prediction_output = []
    for drug, tbprofiler_resistance_prediction_record in tbprofiler_resistance_prediction_by_drug.items():
        tbprofiler_resistance_prediction = tbprofiler_resistance_prediction_record['tbprofiler_resistance_prediction']
        ppv_threshold = args.global_ppv_threshold
        if drug in drug_ppv_thresholds:
            ppv_threshold = drug_ppv_thresholds[drug]

        resistance_and_other_mutations = adjusted_resistance_and_other_mutations_by_drug.get(drug, [])
        all_ppvs = []
        resistance_ppvs = []
        num_mutations_total = 0
        num_mutations_with_ppv_values = 0
        num_resistance_mutations_total = 0
        num_resistance_mutations_with_ppv_values = 0
        for mutation in resistance_and_other_mutations:
            num_mutations_total += 1
            if mutation['tbprofiler_variant_type'] == 'resistance':
                num_resistance_mutations_total += 1
            if 'catalogue_ppv_solo_estimate' in mutation:
                ppv = mutation['catalogue_ppv_solo_estimate']
                if ppv is not None:
                    num_mutations_with_ppv_values += 1
                    all_ppvs.append(ppv)
                    if mutation['tbprofiler_variant_type'] == 'resistance':
                        num_resistance_mutations_total += 1
                        num_resistance_mutations_with_ppv_values += 1
                        resistance_ppvs.append(ppv)

        num_mutations_with_ppv_values = len(all_ppvs)
        num_resistance_mutations_with_ppv_values = len(resistance_ppvs)

        sum_all_ppv = sum(all_ppvs)
        combined_all_ppv = combine_resistance_probabilities(all_ppvs)
        sum_resistance_ppv = sum(resistance_ppvs)
        combined_resistance_ppv = combine_resistance_probabilities(resistance_ppvs)
        output_row = {
            'sample_id': args.sample_id,
            'drug': drug,
            'total_num_putative_drug_related_mutations_detected': num_mutations_total,
            'total_num_putative_drug_related_mutations_with_ppv_values': num_mutations_with_ppv_values,
            'sum_all_putative_drug_related_mutation_ppvs': sum_all_ppv,
            'combined_all_putative_drug_related_mutation_ppvs': combined_all_ppv,
            'num_tbprofiler_resistance_mutations_detected': num_resistance_mutations_total,
            'num_tbprofiler_resistance_mutations_with_ppv_values': num_resistance_mutations_with_ppv_values,
            'sum_tbprofiler_resistance_mutation_ppvs': sum_resistance_ppv,
            'combined_tbprofiler_resistance_mutation_ppvs': combined_resistance_ppv,
            'ppv_threshold': ppv_threshold,
            'tbprofiler_resistance_prediction': tbprofiler_resistance_prediction,
            'tbprofiler_version': tbprofiler_version,
            'tbdb_commit': tbdb_commit,
        }
        if ppv_threshold is not None:
            if combined_all_ppv >= ppv_threshold:
                output_row['adjusted_resistance_prediction'] = 'Predicted Resistant'
            else:
                output_row['adjusted_resistance_prediction'] = 'Predicted Susceptible'
        else:
            output_row['adjusted_resistance_prediction'] = 'Prediction Undetermined'
            
        adjusted_resistance_prediction_output.append(output_row)


    with open(args.output_adjusted_resistance_predictions, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=adjusted_resistance_prediction_output_fieldnames, extrasaction='ignore')
        writer.writeheader()
        for row in adjusted_resistance_prediction_output:
            for field in rounded_fields:
                if field in row and row[field] is not None:
                    row[field] = round(row[field], 6)
            writer.writerow(row)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Adjust the TBProfiler resistance prediction output')
    parser.add_argument('--sample-id', '-s', type=str, help='The sample ID')
    parser.add_argument('--input-tbprofiler-resistance-mutations', type=Path, help='The input tbprofiler resistance mutations file (csv)')
    parser.add_argument('--input-tbprofiler-resistance-predictions', type=Path, help='The input tbprofiler resistance predictions file (csv)')
    parser.add_argument('--input-tbprofiler-full-report', type=Path, help='The input tbprofiler full report file (json)')
    parser.add_argument('--input-who-mutation-catalogue', type=Path, help='')
    parser.add_argument('--global-ppv-threshold', type=float, default=0.90, help='The global PPV threshold for resistance prediction')
    parser.add_argument('--drug-ppv-thresholds', type=Path, help='The input drug PPV thresholds file (csv)')
    parser.add_argument('--output-adjusted-resistance-predictions', type=Path, help='The output adjusted tbprofiler resistance prediction file (csv)')
    parser.add_argument('--output-adjusted-resistance-mutations', type=Path, help='The output adjusted tbprofiler resistance mutations file (csv)')
    args = parser.parse_args()
    main(args)

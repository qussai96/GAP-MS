#!/usr/bin/env python

import sys
import pandas as pd
import re
from pathlib import Path

output_dir = Path.cwd()/ 'outputs'
output_dir.mkdir(exist_ok=True)


def extract_id(suppl_value, pattern):
    match = re.search(pattern, suppl_value)
    if match:
        return match.group(1)
    return None

def gtf_to_df_with_genes(gtf_file):
    column_names = ['Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Suppl']
    df = pd.read_csv(gtf_file, sep='\t', index_col=False, names=column_names, dtype={'Start': int, 'End': int})
    
    df['Gene'] = df.apply(lambda row: extract_id(row['Suppl'], r'gene_id "([^"]+)"') , axis=1)
    df['Protein'] = df.apply(lambda row: extract_id(row['Suppl'], r'transcript_id "([^"]+)"'), axis=1)
    return df

def get_list(input_list):
    list_df = pd.read_csv(input_list, sep='\t', header=0, index_col=False)
    filter_list = set(list_df['Protein'])
    return filter_list
    

def get_filtered_df(pred_df, filter_list, name):
    # Filter the GTF file to keep only rows where the 'Protein' is in the filter list
    filtered_df = pred_df[pred_df["Type"] == "CDS"]
    filtered_df = filtered_df[filtered_df['Protein'].isin(filter_list)]
    filtered_df.drop(columns=['Gene', 'Protein'], inplace=True)
    filtered_df.to_csv(output_dir / f'filtered_{name}_pred.gtf', sep='\t', index=False, header=False)


def main(argv):
    input_default = argv[0]
    input_list = argv[1]
    name=argv[2]

    pred_df = gtf_to_df_with_genes(input_default)
    filter_list = get_list(input_list)
    
    get_filtered_df(pred_df, filter_list, name)

if __name__ == "__main__":
    main(sys.argv[1:])
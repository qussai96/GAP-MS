import pandas as pd
import re
from Bio.Seq import Seq
from Bio import SeqIO
import sys

def df_with_genes(df):
    df = df[df['Type'] == 'CDS'].copy()  # Ensure you're working on a copy to avoid SettingWithCopyWarning
    gene_ids = []
    for _, row in df.iterrows():
        string = row['Suppl']
        if row['Source'] in [
            'EMBL', 'EVM', 'EuGene', 'Helixer', 'BMAPv1', 'phytozome8_0', 'phytozome9_0',
            'phytozomev10', 'phytozomev11', 'phytozomev12', 'phytozomev13', 'Gmove', 'maker',
            'PCPHybridSetV1', 'Psat_v1a', 'juglans_regia'
        ]:
            extracted_id = extract_id(string, r'Parent=([^;]+)')
        else:
            extracted_id = extract_id(string, r'transcript_id "([^"]+)"')
        gene_ids.append(extracted_id)
    df['Genes'] = gene_ids
    return df

def extract_id(string, pattern):
    match = re.search(pattern, string)
    if match:
        return match.group(1)
    return "unknown_id"  # Return a default value if pattern not found

def add_seqs(input_genome_path, input_gtf_path):
    column_names = ['Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Suppl']
    df = pd.read_csv(input_gtf_path, sep='\t', header=None, names=column_names, dtype={'Start': int, 'End': int})
    df['Seq'] = ""
    df['Frame'] = df['Frame'].replace('.', 0).astype(int)  # Replace '.' and convert to int

    genome_sequences = SeqIO.to_dict(SeqIO.parse(input_genome_path, "fasta"))

    for index, row in df.iterrows():
        seqid = row['Seqid']
        start = row['Start'] - 1
        end = row['End']
        try:
            seq = str(genome_sequences[seqid].seq[start:end])
            df.at[index, 'Seq'] = seq
        except KeyError:
            print(f"Warning: sequence ID '{seqid}' not found in genome.")  # Optional: warning message
            continue
    return df

def get_proteins_from_cds_gtf(input_genome_path, input_gtf_path):
    df = add_seqs(input_genome_path, input_gtf_path)
    df_genes = df_with_genes(df)
    df_groups = df_genes.groupby('Genes')

    output_path = input_gtf_path.with_suffix('.faa')

    with open(output_path, 'w') as out_f:
        for name, sub_df in df_groups:
            df_cds = sub_df[sub_df['Type'] == 'CDS']
            if df_cds.empty:
                continue

            strand = sub_df['Strand'].iloc[0]
            df_cds = df_cds.sort_values('Start')
            cds_regions = [(int(row['Start']), int(row['End'])) for _, row in df_cds.iterrows()]
            cds_regions.sort(key=lambda x: x[0])

            try:
                gene = ''.join(df_cds['Seq'].tolist()).replace(' ', '')
            except TypeError:
                gene = 'TAANNNTTANNNNNNNNNTTANNNTAA'

            first_frame = df_cds['Frame'].iloc[0] if strand == '+' else df_cds['Frame'].iloc[-1]
            frame_offset = int(first_frame) if isinstance(first_frame, int) else 0

            if strand == '-':
                gene = str(Seq(gene).reverse_complement())

            gene = gene[frame_offset:]

            if len(gene) % 3 != 0:
                gene = gene[:len(gene) - (len(gene) % 3)]

            protein = str(Seq(gene).translate(to_stop=False))

            out_f.write(f'>{name}\n{protein}\n')

    check_abnormal_amino_acids(output_path)
    return output_path

def check_abnormal_amino_acids(fasta_file):
    count = 0
    valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWYX")
    abnormal_sequences = {}

    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq).rstrip("*")  # Avoid slicing, just remove the last '*'
        abnormal_amino_acids = set(sequence) - valid_amino_acids

        if abnormal_amino_acids:
            count += 1
            abnormal_sequences[record.id] = abnormal_amino_acids

    if abnormal_sequences:
        print(f"Number of sequences with abnormal amino acids = {count}")
        print("Sequences with abnormal amino acids are:")
        for seq_id, amino_acids in abnormal_sequences.items():
            print(f"ID: {seq_id}, Abnormal Amino Acids: {', '.join(amino_acids)}")
    else:
        print("All sequences are normal.")

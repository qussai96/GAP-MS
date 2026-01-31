import sys
import pandas as pd
import re
from Bio import SeqIO
from .gtf_utils import gtf_to_df_with_genes
from pathlib import Path


def compare_gtf_entries(df1, df2):
    # Drop duplicates from both dataframes and keep only relevant columns
    df1 = df1.drop_duplicates(subset=['Seqid', 'Start', 'End', 'Strand', 'Protein'], keep='first')
    df2 = df2.drop_duplicates(subset=['Seqid', 'Start', 'End', 'Strand', 'Protein'], keep='first')

    def add_protein_rows(df):
        new_rows = []
        for (seqid, strand, protein, gene), group in df.groupby(['Seqid', 'Strand', 'Protein', 'Gene']):
            start = group['Start'].min()
            end = group['End'].max()
            new_rows.append({
                'Seqid': seqid,
                'Source': group.iloc[0].get('Source', 'Gnomon'),
                'Type': 'protein',
                'Start': start,
                'End': end,
                'Prediction_score': '.',
                'Strand': strand,
                'Frame': '0',
                'Suppl': f"ID={protein};Parent={gene}",
                'Gene': gene,
                'Protein': protein
            })
        return pd.concat([df, pd.DataFrame(new_rows)], ignore_index=True)

    df1 = add_protein_rows(df1)
    df2 = add_protein_rows(df2)
    df1 = df1[df1['Type'].isin(['protein'])]
    df2 = df2[df2['Type'].isin(['protein'])]
    # save_gtf_subset(df1, df1['Protein'], Path('/home/students/q.abbas/GAPMS/test/moak_dataset_phaseolus_vulgaris/'),'ref_proteins.gtf')
    # save_gtf_subset(df2, df2['Protein'], Path('/home/students/q.abbas/GAPMS/test/moak_dataset_phaseolus_vulgaris/'),'pred_proteins.gtf')
    df1= df1[['Seqid', 'Start', 'End', 'Gene', 'Protein']]
    df2= df2[['Seqid', 'Start', 'End', 'Gene', 'Protein']]
    
    # Perform an inner merge on seqid and start
    identical_start_transcripts = pd.merge(df1, df2, on=['Seqid', 'Start'], suffixes=('_df1', '_df2'))
    if len(identical_start_transcripts) == 0:
        print("Couldn't merge reference and prediction gtfs, they might have differnet Seqids")
        sys.exit(1)
    # save_gtf_subset(df2, identical_start_transcripts['Protein_df2'], Path('/home/students/q.abbas/GAPMS/test/moak_dataset_phaseolus_vulgaris/'),'identical_start_proteins.gtf')
    
    identical_end_transcripts = df1.merge(df2, on=['Seqid'], suffixes=('_df1', '_df2'))
    identical_end_transcripts = identical_end_transcripts[abs(identical_end_transcripts['End_df1'] - identical_end_transcripts['End_df2']) <= 3]

    # save_gtf_subset(df2, identical_end_transcripts['Protein_df2'], Path('/home/students/q.abbas/GAPMS/test/moak_dataset_phaseolus_vulgaris/'),'identical_end_proteins.gtf')
    
    # Filter rows where the 'end' values are within ±3
    identical = identical_start_transcripts[identical_start_transcripts.apply(lambda row: abs(row['End_df1'] - row['End_df2']) <= 3, axis=1)]
    
    # Extract identical predicted transcripts
    identical_predicted_transcript = identical[['Seqid', 'Start', 'End_df2', 'Protein_df2']].rename(
        columns={'End_df2': 'End', 'Protein_df2': 'Protein'}
    )
    
    # Filter out identical transcripts from df1
    not_identical_ref = df1.merge(
        identical[['Seqid', 'Start', 'End_df1']].rename(columns={'End_df1': 'End'}),
        on=['Seqid', 'Start', 'End'],
        how='left',
        indicator=True
    ).query('_merge == "left_only"').drop(columns=['_merge'])
    
    # Filter out identical transcripts from df2
    not_identical_pred = df2.merge(
        identical[['Seqid', 'Start', 'End_df2']].rename(columns={'End_df2': 'End'}),
        on=['Seqid', 'Start', 'End'],
        how='left',
        indicator=True
    ).query('_merge == "left_only"').drop(columns=['_merge'])

    
    not_identical_ref['Seqid'] = not_identical_ref['Seqid'].astype('category')
    not_identical_pred['Seqid'] = not_identical_pred['Seqid'].astype('category')

    # # Identify overlapped transcript
    # Set the chunk size
    chunk_size = 50000  # You can adjust this based on your data size

    # Get the total number of chunks
    total_chunks = len(not_identical_pred) // chunk_size + 1

    # Initialize an empty DataFrame to store the final result
    result_df = pd.DataFrame()

    # Process each chunk separately
    for chunk_num in range(total_chunks):
        start_idx = chunk_num * chunk_size
        end_idx = (chunk_num + 1) * chunk_size

        # Get the chunk of not_identical_pred
        chunk_not_identical_pred = not_identical_pred.iloc[start_idx:end_idx]

        # Merge with not_identical_ref
        not_identical_transcript = pd.merge(not_identical_ref, chunk_not_identical_pred, on='Seqid', how='inner', suffixes=('_x', '_y'))

        # Filter overlapped transcript based on start and end conditions
        condition = (
            (not_identical_transcript['Start_y'] <= (not_identical_transcript['End_x'] + 3)) &
            (not_identical_transcript['End_y'] >= (not_identical_transcript['Start_x'] - 3))
        )

        overlapped_transcript = not_identical_transcript[condition]

        # Concatenate the result to the final DataFrame
        result_df = pd.concat([result_df, overlapped_transcript], ignore_index=True)
    # Identify overlapped reference transcript and remove duplicates
    
    overlapped_reference_transcript = result_df[['Seqid', 'Start_x', 'End_x', 'Protein_x']].rename(columns={'Start_x': 'Start', 'End_x': 'End', 'Protein_x': 'Protein'}).drop_duplicates()
    # Identify overlapped predicted transcript and remove duplicates
    overlapped_predicted_transcript = result_df[['Seqid', 'Start_y', 'End_y', 'Protein_y']].rename(columns={'Start_y': 'Start', 'End_y': 'End', 'Protein_y': 'Protein'}).drop_duplicates()
    
    # Identify missed transcript in df1
    missed_transcript = not_identical_ref.merge(overlapped_reference_transcript, on=['Seqid', 'Start', 'End', 'Protein'], how='left', indicator=True)
    missed_transcript = missed_transcript.query('_merge == "left_only"').drop(columns=['_merge'])

    # Identify full false positives in df2
    wrong_transcript = not_identical_pred.merge(overlapped_predicted_transcript, on=['Seqid', 'Start', 'End', 'Protein'], how='left', indicator=True)
    wrong_transcript = wrong_transcript.query('_merge == "left_only"').drop(columns=['_merge'])
    identical_start_proteins = set(identical_start_transcripts['Protein_df2'])
    identical_end_proteins = set(identical_end_transcripts['Protein_df2'])
    
    identical_proteins = set(identical_predicted_transcript['Protein'])
    overlapped_predicted_proteins= set(overlapped_predicted_transcript['Protein'])
    new_predicted_proteins= set(wrong_transcript['Protein'])

    # save_gtf_subset(df2, overlapped_predicted_proteins, Path('/home/students/q.abbas/GAPMS/test/moak_dataset_phaseolus_vulgaris/'),'overlapped_predicted_proteins.gtf')
    # save_gtf_subset(df2, new_predicted_proteins, Path('/home/students/q.abbas/GAPMS/test/moak_dataset_phaseolus_vulgaris/'),'new_predicted_proteins.gtf')

    return identical_start_proteins, identical_end_proteins, identical_proteins, overlapped_predicted_proteins, new_predicted_proteins

def get_new_proteins(output_dir, reference_gtf, reference_fasta=None):
    reference_fasta = reference_fasta if reference_fasta else output_dir / "translated_reference_proteins.fa"
    prediction_gtf= output_dir / "supported_proteins.gtf"
    prediction_fasta= output_dir / "supported_proteins.fa"
    all_proteins_scores= output_dir / "all_proteins_scores.tsv"
    
    new_proteins_dir = output_dir / "Novel"

    reference_gtf_df = gtf_to_df_with_genes(reference_gtf)
    reference_gtf_df = reference_gtf_df[reference_gtf_df['Type'].isin(['CDS'])]
    prediction_gtf_df = gtf_to_df_with_genes(prediction_gtf)
    prediction_gtf_df = prediction_gtf_df[prediction_gtf_df['Type'].isin(['CDS'])]
    
    identical_start_proteins, identical_end_proteins, identical_proteins, overlapped_predicted_proteins, new_predicted_proteins  = compare_gtf_entries(reference_gtf_df, prediction_gtf_df)
    
    print (f'Number of reference transcripts = {len(set(reference_gtf_df["Protein"]))}')
    print (f'Number of supported predicted transcripts = {len(set(prediction_gtf_df["Protein"]))}')
    print (f'Number of identical predicted transcripts = {len(set(identical_proteins))}')
    print (f'Number of identical start predicted transcripts = {len(identical_start_proteins)-len(set(identical_proteins))}')
    print (f'Number of identical end predicted transcripts = {len(identical_end_proteins)-len(set(identical_proteins))}')
    print (f'Number of overlapped predicted transcripts = {len(overlapped_predicted_proteins)}')
    print (f'Number of newly predicted transcripts no overlap with any reference = {len(new_predicted_proteins)}')

    overlapped_predicted_proteins_gtf = prediction_gtf_df[prediction_gtf_df['Protein'].isin(overlapped_predicted_proteins)]
    overlapped_predicted_proteins_gtf = overlapped_predicted_proteins_gtf.drop(['Protein', 'Gene'], axis=1)
    overlapped_predicted_proteins_gtf.to_csv(new_proteins_dir / f"overlapped_predicted_proteins.gtf", sep='\t', header=False, index=False)


    new_predicted_proteins_gtf = prediction_gtf_df[prediction_gtf_df['Protein'].isin(new_predicted_proteins)]
    new_predicted_proteins_gtf = new_predicted_proteins_gtf.drop(['Protein', 'Gene'], axis=1)
    new_predicted_proteins_gtf.to_csv(new_proteins_dir / f"new_predicted_proteins.gtf", sep='\t', header=False, index=False)


    prediction_seq_dict = {record.id: str(record.seq) for record in SeqIO.parse(prediction_fasta, "fasta")}
    new_predicted_proteins_dict = {key: prediction_seq_dict[key] for key in new_predicted_proteins if key in prediction_seq_dict}
    overlapped_predicted_proteins_dict = {key: prediction_seq_dict[key] for key in overlapped_predicted_proteins if key in prediction_seq_dict}


    def write_seq_dict_to_file(data, filename):
        with open(filename, "w") as file:
            for key, value in data.items():
                file.write(f">{key}\n{value}\n")
    write_seq_dict_to_file(overlapped_predicted_proteins_dict, new_proteins_dir/"overlapped_proteins.fa")
    write_seq_dict_to_file(new_predicted_proteins_dict, new_proteins_dir/"new_predicted_proteins.fa")


    all_protein_features_df = pd.read_csv(all_proteins_scores, sep='\t', header=0)
    new_predicted_proteins_df = all_protein_features_df[
    all_protein_features_df['Protein'].isin(set(new_predicted_proteins))]
    new_predicted_proteins_df.to_csv(new_proteins_dir/"new_predicted_proteins_scores.tsv", sep='\t', index=False)
        
    # Read sequences from both files
    sequences = []
    for file in [reference_fasta, new_proteins_dir/"new_predicted_proteins.fa"]:
        sequences.extend(SeqIO.parse(file, "fasta"))

    # Write combined sequences to a new file
    with open(new_proteins_dir / "reference+new_proteins.fa", "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")

    print(f"Combined FASTA file saved")
import sys
from bisect import bisect_right
import pandas as pd
from Bio import SeqIO
from pathlib import Path
from .compare_supp_ref import load_peptides_bed
from .gtf_utils import gtf_to_df_with_genes


def _merged_reference_cds(reference_gtf):
    """Build merged CDS intervals per contig, deliberately ignoring strand."""
    reference_df = gtf_to_df_with_genes(reference_gtf)
    reference_df = reference_df[reference_df['Type'].eq('CDS')]
    merged = {}
    for seqid, group in reference_df.groupby('Seqid'):
        intervals = sorted((int(row.Start), int(row.End)) for row in group.itertuples())
        compact = []
        for start, end in intervals:
            if compact and start <= compact[-1][1] + 1:
                compact[-1] = (compact[-1][0], max(compact[-1][1], end))
            else:
                compact.append((start, end))
        merged[seqid] = ([interval[0] for interval in compact], compact)
    return merged


def _overlaps_reference_cds(seqid, start, end, reference_index):
    starts, intervals = reference_index.get(seqid, ([], []))
    if not intervals:
        return False
    index = bisect_right(starts, int(end)) - 1
    return index >= 0 and intervals[index][1] >= int(start)


def parse_gffcompare_tmap(
    output_dir, prediction_gtf, reference_gtf, reference_fasta=None,
    peptides_bed=None, prediction_fasta=None, tmap_file=None,
):
    # 1. Setup Paths
    output_dir = Path(output_dir)
    compare_dir = output_dir / "Compare_to_Reference"
    compare_dir.mkdir(parents=True, exist_ok=True)
    new_proteins_dir = compare_dir / "Novel"
    new_proteins_dir.mkdir(parents=True, exist_ok=True)

    # Set file paths
    # Prefer gffcompare artifacts written into Compare_to_Reference/, but keep a
    # fallback for older layouts.
    if tmap_file is not None:
        tmap_file = Path(tmap_file)
    else:
        preferred_tmap = (
            compare_dir / "gffcmp.peptide_mapped_proteins_cleaned.gtf.tmap"
        )
        tmap_files = sorted(compare_dir.glob("*.tmap")) or sorted(
            output_dir.glob("*.tmap")
        )
        if preferred_tmap.exists():
            tmap_file = preferred_tmap
        elif tmap_files:
            tmap_file = tmap_files[0]
        else:
            tmap_file = preferred_tmap
    
    prediction_fasta = (
        Path(prediction_fasta)
        if prediction_fasta is not None
        else output_dir / "supported_proteins.fa"
    )
    all_proteins_scores = output_dir / "all_proteins_scores.tsv"
    peptides_bed = Path(peptides_bed) if peptides_bed else output_dir / "peptides_mapped.bed"
    
    reference_fasta = Path(reference_fasta) if reference_fasta else output_dir / "translated_reference_proteins.fa"

    # 2. Parse gffcompare TMAP file
    print(f"Parsing TMAP file: {tmap_file}")
    # TMAP is tab-separated. Key columns: 
    # 'class_code' (col 3) and 'qry_id' (col 5, which matches your Protein ID)
    tmap_df = pd.read_csv(tmap_file, sep='\t')
    
    # '=' is an identical intron chain in gffcompare, not proof of an
    # identical protein sequence.
    identical_proteins = set(tmap_df[tmap_df['class_code'] == '=']['qry_id'])

    prediction_seq_dict = {
        record.id: str(record.seq)
        for record in SeqIO.parse(prediction_fasta, "fasta")
    }
    peptide_supported_ids = set()
    if all_proteins_scores.exists():
        support_scores = pd.read_csv(all_proteins_scores, sep='\t')
        if {'Protein', 'mapped_peptides'}.issubset(support_scores.columns):
            peptide_supported_ids = set(
                support_scores.loc[
                    pd.to_numeric(support_scores['mapped_peptides'], errors='coerce').fillna(0) > 0,
                    'Protein',
                ].dropna().astype(str)
            ) & set(prediction_seq_dict)
    if not peptide_supported_ids:
        print("Warning: no peptide-supported proteins were available for novel-gene classification")
    prediction_gtf_df = gtf_to_df_with_genes(
        prediction_gtf,
        valid_protein_ids=prediction_seq_dict,
    )
    prediction_cds = prediction_gtf_df[prediction_gtf_df['Type'].eq('CDS')].copy()
    prediction_cds['Gene'] = prediction_cds['Gene'].where(
        prediction_cds['Gene'].notna(), prediction_cds['Protein']
    )

    # A novel gene is deliberately strict: none of the CDS bases from any of
    # its supported isoforms may overlap any reference CDS, on either strand.
    reference_index = _merged_reference_cds(reference_gtf)
    candidate_novel_genes = []
    supported_genes = set(
        prediction_cds.loc[
            prediction_cds['Protein'].isin(peptide_supported_ids), 'Gene'
        ].dropna()
    )
    for gene_id, gene_cds in prediction_cds.groupby('Gene', dropna=False):
        if gene_id not in supported_genes:
            continue
        has_reference_cds_overlap = any(
            _overlaps_reference_cds(row.Seqid, row.Start, row.End, reference_index)
            for row in gene_cds.itertuples()
        )
        if not has_reference_cds_overlap:
            supported_gene_proteins = sorted(set(
                gene_cds.loc[
                    gene_cds['Protein'].isin(peptide_supported_ids), 'Protein'
                ].dropna().astype(str)
            ))
            candidate_novel_genes.append((gene_id, gene_cds, supported_gene_proteins))

    # A peptide supports a novel gene only when all of its prediction mappings
    # belong to that same gene. Ambiguity among isoforms of the same novel gene
    # is allowed; a mapping to any other gene invalidates that peptide.
    peptide_locations = load_peptides_bed(peptides_bed)
    protein_gene_pairs = (
        prediction_cds[['Protein', 'Gene']]
        .dropna()
        .astype(str)
        .drop_duplicates()
    )
    unambiguous_protein_genes = protein_gene_pairs.groupby('Protein')['Gene'].nunique()
    unambiguous_protein_genes = set(
        unambiguous_protein_genes[unambiguous_protein_genes.eq(1)].index
    )
    protein_to_gene = dict(
        protein_gene_pairs[
            protein_gene_pairs['Protein'].isin(unambiguous_protein_genes)
        ].itertuples(index=False, name=None)
    )

    # Resolve peptide ambiguity once. The previous gene-by-peptide nested scan
    # was unnecessarily quadratic on large annotations.
    unique_peptides_by_gene = {}
    peptide_mapped_proteins = {}
    for peptide, locations in peptide_locations.items():
        mapped_proteins = [location.get('protein') for location in locations]
        mapped_genes = {protein_to_gene.get(protein) for protein in mapped_proteins}
        genomic_footprints = {
            (
                location.get('seqid'), location.get('strand'),
                tuple(location.get('blocks', [])),
            )
            for location in locations
        }
        if (
            locations
            and None not in mapped_genes
            and len(mapped_genes) == 1
            and len(genomic_footprints) == 1
        ):
            unique_peptides_by_gene.setdefault(mapped_genes.pop(), set()).add(peptide)
            peptide_mapped_proteins[peptide] = set(mapped_proteins)

    new_predicted_proteins = set()
    novel_gene_rows = []
    for gene_id, gene_cds, supported_gene_proteins in candidate_novel_genes:
        gene_id = str(gene_id)
        supported_gene_proteins = set(supported_gene_proteins)
        supporting_peptides = {
            peptide
            for peptide in unique_peptides_by_gene.get(gene_id, set())
            if peptide_mapped_proteins[peptide] & supported_gene_proteins
        }

        if not supporting_peptides:
            continue

        new_predicted_proteins.update(supported_gene_proteins)
        novel_gene_rows.append({
            'Gene': gene_id,
            'Supported_proteins': ','.join(sorted(supported_gene_proteins)),
            'Event_supporting_peptides': ','.join(sorted(supporting_peptides)),
            'Seqids': ','.join(sorted(set(gene_cds['Seqid'].dropna().astype(str)))),
            'CDS_start': int(gene_cds['Start'].min()),
            'CDS_end': int(gene_cds['End'].max()),
        })

    pd.DataFrame(
        novel_gene_rows,
        columns=[
            'Gene', 'Supported_proteins', 'Event_supporting_peptides',
            'Seqids', 'CDS_start', 'CDS_end',
        ],
    ).to_csv(new_proteins_dir / 'new_predicted_genes.tsv', sep='\t', index=False)

    # 3. Print Statistics
    # We count FASTA records for accurate file lengths
    reference_fasta_exists = reference_fasta.exists()
    ref_count = len(list(SeqIO.parse(reference_fasta, "fasta"))) if reference_fasta_exists else 0
    supp_count = len(peptide_supported_ids)

    print(f'Number of reference transcripts = {ref_count}')
    print(f'Number of peptide-supported predicted transcripts = {supp_count}')
    print(f'Number of matching-intron-chain predicted transcripts = {len(identical_proteins)}')
    print(f'Number of entirely novel genes = {len(novel_gene_rows)}')
    print(f'Number of supported proteins from entirely novel genes = {len(new_predicted_proteins)}')

    # 4. Generate GTF for New Proteins
    print("Generating GTF for new proteins...")
    # Filter for CDS only
    prediction_gtf_df = prediction_gtf_df[prediction_gtf_df['Type'].isin(['CDS'])]
    
    # Keep only rows belonging to new proteins
    new_predicted_proteins_gtf = prediction_gtf_df[prediction_gtf_df['Protein'].isin(new_predicted_proteins)]
    
    # Remove 'Protein' and 'Gene' columns before saving (as per your snippet)
    cols_to_drop = [c for c in ['Protein', 'Gene'] if c in new_predicted_proteins_gtf.columns]
    new_predicted_proteins_gtf = new_predicted_proteins_gtf.drop(cols_to_drop, axis=1)
    
    new_predicted_proteins_gtf.to_csv(new_proteins_dir / "new_predicted_proteins.gtf", sep='\t', header=False, index=False)

    # 5. Generate FASTA for New Proteins
    print("Generating FASTA for new proteins...")
    # Extract sequences for new proteins only
    new_predicted_proteins_dict = {
        key: prediction_seq_dict[key] 
        for key in new_predicted_proteins 
        if key in prediction_seq_dict
    }

    def write_seq_dict_to_file(data, filename):
        with open(filename, "w") as file:
            for key, value in data.items():
                file.write(f">{key}\n{value}\n")
                
    write_seq_dict_to_file(new_predicted_proteins_dict, new_proteins_dir / "new_predicted_proteins.fa")

    # 6. Extract Scores for New Proteins
    print("Extracting scores for new proteins...")
    if all_proteins_scores.exists():
        all_protein_features_df = pd.read_csv(all_proteins_scores, sep='\t', header=0)
        
        # Ensure the score file has a 'Protein' column to match against
        if 'Protein' in all_protein_features_df.columns:
            new_predicted_proteins_df = all_protein_features_df[
                all_protein_features_df['Protein'].isin(new_predicted_proteins)
            ]
            new_predicted_proteins_df.to_csv(new_proteins_dir / "new_predicted_proteins_scores.tsv", sep='\t', index=False)
        else:
            print("Warning: 'Protein' column missing in scores file. Skipping score extraction.")
    else:
        print(f"Warning: Scores file not found at {all_proteins_scores}")

    # 7. Create Combined FASTA (Reference + New)
    combined_fasta_path = new_proteins_dir / "reference+new_proteins.fa"
    if reference_fasta_exists:
        print("Creating combined FASTA (Reference + New)...")
        sequences = list(SeqIO.parse(reference_fasta, "fasta"))
        sequences.extend(SeqIO.parse(new_proteins_dir / "new_predicted_proteins.fa", "fasta"))
        with open(combined_fasta_path, "w") as output_handle:
            SeqIO.write(sequences, output_handle, "fasta")
        print(f"Combined FASTA file saved at: {combined_fasta_path}")
    else:
        if combined_fasta_path.exists():
            combined_fasta_path.unlink()
        print("Reference protein FASTA not supplied; skipping combined reference+novel FASTA")

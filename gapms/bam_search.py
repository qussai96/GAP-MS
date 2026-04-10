from pathlib import Path
import subprocess

import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO

from .gtf_utils import gtf_to_df_with_genes
from .tools import run_gffcompare


def _run_command(command, tool_name, cwd=None):
    """Run an external BAM-search preparation command with helpful errors."""
    printable = " ".join(str(part) for part in command)
    print(f"Running {tool_name}: {printable}")

    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            cwd=str(cwd) if cwd else None,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"{tool_name} was not found on PATH. Please install it before using --bam."
        ) from exc

    if result.returncode != 0:
        raise RuntimeError(
            f"{tool_name} failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )

    return result


def run_stringtie(bam_path, output_dir):
    """Assemble transcript models from a BAM alignment using StringTie."""
    bam_path = Path(bam_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    transcripts_gtf = output_dir / "transcripts.gtf"
    command = ["stringtie", str(bam_path), "-o", str(transcripts_gtf)]
    _run_command(command, "StringTie")
    return transcripts_gtf


def _refresh_fasta_index(assembly_path):
    """Remove a stale FASTA index and regenerate it with samtools faidx."""
    assembly_path = Path(assembly_path)
    fai_path = Path(f"{assembly_path}.fai")

    if fai_path.exists():
        print(f"Removing stale FASTA index: {fai_path}")
        fai_path.unlink()

    _run_command(["samtools", "faidx", str(assembly_path)], "samtools faidx")


def run_gffread_transcripts(assembly_path, transcripts_gtf, output_dir):
    """Extract transcript FASTA sequences from the assembled GTF."""
    output_dir = Path(output_dir)
    transcripts_fa = output_dir / "transcripts.fa"
    command = [
        "gffread",
        str(transcripts_gtf),
        "-g",
        str(assembly_path),
        "-w",
        str(transcripts_fa),
    ]

    try:
        _run_command(command, "gffread")
    except RuntimeError as exc:
        error_message = str(exc)
        if "improper genomic coordinate" not in error_message:
            raise

        print("Detected possible stale FASTA index; regenerating .fai and retrying gffread...")
        _refresh_fasta_index(assembly_path)
        _run_command(command, "gffread")

    return transcripts_fa


def run_transdecoder(transcripts_fa, output_dir):
    """Predict likely coding ORFs and translated peptide sequences."""
    output_dir = Path(output_dir)
    _run_command(["TransDecoder.LongOrfs", "-t", str(transcripts_fa)], "TransDecoder.LongOrfs", cwd=output_dir)
    _run_command(["TransDecoder.Predict", "-t", str(transcripts_fa)], "TransDecoder.Predict", cwd=output_dir)
    return Path(f"{transcripts_fa}.transdecoder.pep")


def _parse_attributes(attributes_text):
    """Parse GTF/GFF attributes into a dict."""
    parsed = {}
    for raw_part in str(attributes_text).strip().strip(';').split(';'):
        part = raw_part.strip()
        if not part:
            continue
        if '=' in part:
            key, value = part.split('=', 1)
        elif ' ' in part:
            key, value = part.split(' ', 1)
        else:
            continue
        parsed[key.strip()] = value.strip().strip('"')
    return parsed


def _load_transcript_models(transcripts_gtf):
    """Load StringTie transcript models and exon structures keyed by transcript_id."""
    transcript_models = {}

    with open(transcripts_gtf) as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue

            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue

            seqid, source, feature_type, start, end, score, strand, frame, attributes = fields
            feature_type = feature_type.lower()
            if feature_type not in {'transcript', 'mrna', 'exon'}:
                continue

            attrs = _parse_attributes(attributes)
            transcript_id = attrs.get('transcript_id') or attrs.get('ID')
            if not transcript_id:
                continue

            gene_id = attrs.get('gene_id') or attrs.get('gene') or attrs.get('Parent') or transcript_id
            start = int(start)
            end = int(end)

            model = transcript_models.setdefault(
                transcript_id,
                {
                    'seqid': seqid,
                    'source': source,
                    'strand': strand,
                    'score': score,
                    'gene_id': gene_id,
                    'start': start,
                    'end': end,
                    'exons': [],
                },
            )
            model['seqid'] = seqid
            model['strand'] = strand
            model['gene_id'] = gene_id
            if score not in {None, ''}:
                model['score'] = score
            model['start'] = min(model['start'], start)
            model['end'] = max(model['end'], end)

            if feature_type == 'exon':
                model['exons'].append((start, end))

    for model in transcript_models.values():
        if not model['exons']:
            model['exons'] = [(model['start'], model['end'])]
        model['exons'] = sorted(model['exons'], key=lambda exon: exon[0])

    return transcript_models


def _map_transcript_interval_to_genome(feature_start, feature_end, transcript_model):
    """Project a transcript-space interval back onto genomic exon coordinates."""
    feature_start = int(feature_start)
    feature_end = int(feature_end)
    strand = transcript_model['strand']
    ordered_exons = list(transcript_model['exons'])
    if strand == '-':
        ordered_exons = list(reversed(ordered_exons))

    genomic_segments = []
    transcript_cursor = 1

    for exon_start, exon_end in ordered_exons:
        exon_len = exon_end - exon_start + 1
        transcript_exon_start = transcript_cursor
        transcript_exon_end = transcript_cursor + exon_len - 1

        overlap_start = max(feature_start, transcript_exon_start)
        overlap_end = min(feature_end, transcript_exon_end)
        if overlap_start <= overlap_end:
            if strand == '+':
                genomic_start = exon_start + (overlap_start - transcript_exon_start)
                genomic_end = exon_start + (overlap_end - transcript_exon_start)
            else:
                genomic_end = exon_end - (overlap_start - transcript_exon_start)
                genomic_start = exon_end - (overlap_end - transcript_exon_start)

            genomic_segments.append((min(genomic_start, genomic_end), max(genomic_start, genomic_end)))

        transcript_cursor = transcript_exon_end + 1

    return sorted(genomic_segments, key=lambda segment: segment[0])


def convert_transdecoder_to_genome_gtf(transdecoder_gff3, transcripts_gtf, output_gtf):
    """Convert TransDecoder ORF coordinates from transcript space back to genome coordinates."""
    transdecoder_gff3 = Path(transdecoder_gff3)
    transcripts_gtf = Path(transcripts_gtf)
    output_gtf = Path(output_gtf)

    transcript_models = _load_transcript_models(transcripts_gtf)
    if not transcript_models:
        raise RuntimeError(f"No transcript models could be loaded from {transcripts_gtf}")

    orf_records = {}
    with open(transdecoder_gff3) as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue

            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue

            transcript_id, source, feature_type, start, end, score, strand, frame, attributes = fields
            attrs = _parse_attributes(attributes)
            feature_type = feature_type.lower()

            if feature_type in {'mrna', 'transcript'}:
                orf_id = attrs.get('ID') or attrs.get('transcript_id')
                if not orf_id:
                    continue
                entry = orf_records.setdefault(
                    orf_id,
                    {'transcript_id': transcript_id, 'source': source, 'score': score, 'cds_features': []},
                )
                entry['transcript_id'] = transcript_id
                if score not in {None, ''}:
                    entry['score'] = score
            elif feature_type == 'cds':
                orf_id = attrs.get('Parent') or attrs.get('protein_id') or attrs.get('transcript_id')
                if not orf_id:
                    continue
                entry = orf_records.setdefault(
                    orf_id,
                    {'transcript_id': transcript_id, 'source': source, 'score': score, 'cds_features': []},
                )
                entry['transcript_id'] = transcript_id
                entry['cds_features'].append((int(start), int(end), frame if frame != '.' else '0'))
                if score not in {None, ''}:
                    entry['score'] = score

    if not orf_records:
        raise RuntimeError(f"No ORF annotations were found in {transdecoder_gff3}")

    output_gtf.parent.mkdir(parents=True, exist_ok=True)
    converted_orfs = 0
    with open(output_gtf, 'w') as handle:
        for orf_id in sorted(orf_records):
            entry = orf_records[orf_id]
            transcript_model = transcript_models.get(entry['transcript_id'])
            if transcript_model is None or not entry['cds_features']:
                continue

            gene_id = transcript_model['gene_id']
            seqid = transcript_model['seqid']
            strand = transcript_model['strand']
            score = transcript_model.get('score', '.')
            attributes = f'gene_id "{gene_id}"; transcript_id "{orf_id}"; protein_id "{orf_id}";'

            handle.write(
                f"{seqid}\tGAPMS\ttranscript\t{transcript_model['start']}\t{transcript_model['end']}\t{score}\t{strand}\t.\t{attributes}\n"
            )

            wrote_cds = False
            for cds_start, cds_end, frame in sorted(entry['cds_features'], key=lambda feature: feature[0]):
                genomic_segments = _map_transcript_interval_to_genome(cds_start, cds_end, transcript_model)
                for genomic_start, genomic_end in genomic_segments:
                    handle.write(
                        f"{seqid}\tGAPMS\tCDS\t{genomic_start}\t{genomic_end}\t{score}\t{strand}\t{frame}\t{attributes}\n"
                    )
                    wrote_cds = True

            if wrote_cds:
                converted_orfs += 1

    if converted_orfs == 0:
        raise RuntimeError(
            f"No TransDecoder ORFs could be projected back to genome coordinates from {transdecoder_gff3}"
        )

    return output_gtf


def prepare_bam_search_inputs(bam_path, assembly_path, output_dir):
    """Build BAM-derived transcript and peptide-search inputs for the GAP-MS pipeline.

    Workflow:
      1. StringTie -> transcripts.gtf
      2. gffread -> transcripts.fa
      3. TransDecoder.LongOrfs / TransDecoder.Predict -> translated peptide FASTA + ORF annotation

    Returns a dict with the assembled transcript GTF, a genome-coordinate ORF GTF
    derived from the TransDecoder predictions, the transcript FASTA, and the peptide FASTA.
    """
    output_dir = Path(output_dir)
    bam_dir = output_dir if output_dir.name in {"bam_search", "Bam_Search"} else output_dir / "bam_search"
    bam_dir.mkdir(parents=True, exist_ok=True)

    transcripts_gtf = run_stringtie(bam_path, bam_dir)
    transcript_fasta = run_gffread_transcripts(assembly_path, transcripts_gtf, bam_dir)
    protein_fasta = run_transdecoder(transcript_fasta, bam_dir)
    transdecoder_gff3 = Path(f"{transcript_fasta}.transdecoder.gff3")

    if not transdecoder_gff3.exists():
        raise RuntimeError(
            f"TransDecoder annotation file was not created as expected: {transdecoder_gff3}"
        )

    genome_gtf = Path(f"{transcript_fasta}.transdecoder.genome.gtf")
    convert_transdecoder_to_genome_gtf(transdecoder_gff3, transcripts_gtf, genome_gtf)

    return {
        "assembled_gtf": transcripts_gtf,
        "gtf": genome_gtf,
        "transcript_fasta": transcript_fasta,
        "protein_fasta": protein_fasta,
        "orf_gff3": transdecoder_gff3,
    }


def _summarize_loci_from_gtf(gtf_path):
    """Summarize each protein model into a genomic locus row."""
    gtf_path = Path(gtf_path)
    if not gtf_path.exists() or gtf_path.stat().st_size == 0:
        return pd.DataFrame(columns=['protein', 'gene', 'seqid', 'strand', 'start', 'end'])

    gtf_df = gtf_to_df_with_genes(gtf_path)
    cds_df = gtf_df[gtf_df['Type'] == 'CDS'].copy()
    if cds_df.empty:
        return pd.DataFrame(columns=['protein', 'gene', 'seqid', 'strand', 'start', 'end'])

    loci_df = (
        cds_df.groupby('Protein', dropna=False)
        .agg(
            gene=('Gene', lambda values: ','.join(sorted({str(v) for v in values if pd.notna(v) and str(v)}))),
            seqid=('Seqid', 'first'),
            strand=('Strand', 'first'),
            start=('Start', 'min'),
            end=('End', 'max'),
        )
        .reset_index()
        .rename(columns={'Protein': 'protein'})
    )
    loci_df['gene'] = loci_df['gene'].replace('', pd.NA).fillna(loci_df['protein'])
    return loci_df


def _write_candidate_gtf(source_gtf, protein_ids, output_path):
    """Write a GTF subset for the selected candidate proteins."""
    source_gtf = Path(source_gtf)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if not protein_ids or not source_gtf.exists() or source_gtf.stat().st_size == 0:
        output_path.write_text('')
        return output_path

    gtf_df = gtf_to_df_with_genes(source_gtf)
    subset = gtf_df[gtf_df['Protein'].astype(str).isin({str(pid) for pid in protein_ids})].copy()
    subset.drop(columns=['Protein', 'Gene'], errors='ignore').to_csv(
        output_path,
        sep='\t',
        index=False,
        header=False,
    )
    return output_path


def _write_candidate_fasta(source_fasta, protein_ids, output_path):
    """Write a FASTA subset for the selected candidate proteins."""
    source_fasta = Path(source_fasta) if source_fasta else None
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if not protein_ids or source_fasta is None or not source_fasta.exists() or source_fasta.stat().st_size == 0:
        output_path.write_text('')
        return output_path

    wanted = {str(pid) for pid in protein_ids}
    with open(source_fasta) as input_handle, open(output_path, 'w') as output_handle:
        selected = [record for record in SeqIO.parse(input_handle, 'fasta') if record.id.split()[0] in wanted]
        SeqIO.write(selected, output_handle, 'fasta')
    return output_path


def _load_bam_class_codes(bam_tmap_file, accepted_class_codes=None):
    """Load BAM protein IDs carrying the accepted gffcompare overlap classes."""
    accepted_class_codes = list(accepted_class_codes or ['=', 'j', 'c', 'k', 'm', 'n'])
    empty_df = pd.DataFrame(columns=['bam_protein', 'bam_class_codes'])

    bam_tmap_file = Path(bam_tmap_file) if bam_tmap_file else None
    if bam_tmap_file is None or not bam_tmap_file.exists() or bam_tmap_file.stat().st_size == 0:
        return empty_df

    try:
        tmap_df = pd.read_csv(bam_tmap_file, sep='\t')
    except pd.errors.EmptyDataError:
        return empty_df

    if not {'qry_id', 'class_code'}.issubset(tmap_df.columns):
        return empty_df

    tmap_df = tmap_df[['qry_id', 'class_code']].copy()
    tmap_df['qry_id'] = tmap_df['qry_id'].astype(str).str.strip()
    tmap_df['class_code'] = tmap_df['class_code'].fillna('').astype(str).str.strip()
    tmap_df = tmap_df[tmap_df['class_code'].isin(set(accepted_class_codes))]
    if tmap_df.empty:
        return empty_df

    return (
        tmap_df.groupby('qry_id', dropna=False)
        .agg(bam_class_codes=('class_code', lambda values: '{' + ','.join(sorted({str(v) for v in values if pd.notna(v) and str(v)})) + '}'))
        .reset_index()
        .rename(columns={'qry_id': 'bam_protein'})
    )


def _summarize_novel_class_support_from_tmap(
    prediction_novel_gtf,
    bam_annotation_gtf,
    bam_tmap_file,
    accepted_class_codes=None,
):
    """Summarize prediction-branch novel proteins supported by accepted BAM `gffcompare` classes."""
    accepted_class_codes = list(accepted_class_codes or ['=', 'j', 'c', 'k', 'm', 'n'])
    empty_columns = [
        'candidate_id', 'seqid', 'strand', 'prediction_protein', 'prediction_gene',
        'prediction_start', 'prediction_end', 'bam_proteins', 'bam_genes', 'bam_class_codes',
        'bam_match_count',
    ]

    pred_loci = _summarize_loci_from_gtf(prediction_novel_gtf).rename(
        columns={
            'protein': 'prediction_protein',
            'gene': 'prediction_gene',
            'start': 'prediction_start',
            'end': 'prediction_end',
        }
    )
    counts = {
        'accepted_overlap': 0,
        'other_classes': int(len(pred_loci)),
        'total': int(len(pred_loci)),
        'accepted_codes': accepted_class_codes,
        'matched_bam_transcripts': 0,
    }

    if pred_loci.empty:
        return pd.DataFrame(columns=empty_columns), counts

    bam_tmap_file = Path(bam_tmap_file) if bam_tmap_file else None
    if bam_tmap_file is None or not bam_tmap_file.exists() or bam_tmap_file.stat().st_size == 0:
        return pd.DataFrame(columns=empty_columns), counts

    try:
        tmap_df = pd.read_csv(bam_tmap_file, sep='\t')
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=empty_columns), counts

    required_columns = {'qry_id', 'ref_id', 'class_code'}
    if not required_columns.issubset(tmap_df.columns):
        return pd.DataFrame(columns=empty_columns), counts

    pred_ids = set(pred_loci['prediction_protein'].dropna().astype(str))
    accepted_set = set(accepted_class_codes)

    tmap_df = tmap_df[['qry_id', 'ref_id', 'class_code']].copy()
    tmap_df['qry_id'] = tmap_df['qry_id'].astype(str).str.strip()
    tmap_df['ref_id'] = tmap_df['ref_id'].astype(str).str.strip()
    tmap_df['class_code'] = tmap_df['class_code'].fillna('').astype(str).str.strip()

    accepted_df = tmap_df[
        tmap_df['ref_id'].isin(pred_ids) &
        tmap_df['class_code'].isin(accepted_set)
    ].copy()

    counts['accepted_overlap'] = int(accepted_df['ref_id'].nunique())
    counts['other_classes'] = max(counts['total'] - counts['accepted_overlap'], 0)
    counts['matched_bam_transcripts'] = int(accepted_df['qry_id'].nunique())

    if accepted_df.empty:
        return pd.DataFrame(columns=empty_columns), counts

    bam_loci = _summarize_loci_from_gtf(bam_annotation_gtf).rename(
        columns={'protein': 'bam_protein', 'gene': 'bam_gene'}
    ) if bam_annotation_gtf else pd.DataFrame(columns=['bam_protein', 'bam_gene'])

    if not bam_loci.empty:
        accepted_df = accepted_df.merge(
            bam_loci[['bam_protein', 'bam_gene']].drop_duplicates(),
            left_on='qry_id',
            right_on='bam_protein',
            how='left',
        )
    else:
        accepted_df['bam_gene'] = pd.NA

    grouped_matches = (
        accepted_df.groupby('ref_id', dropna=False)
        .agg(
            bam_proteins=('qry_id', lambda values: '{' + ','.join(sorted({str(v) for v in values if pd.notna(v) and str(v)})) + '}'),
            bam_genes=('bam_gene', lambda values: '{' + ','.join(sorted({str(v) for v in values if pd.notna(v) and str(v)})) + '}'),
            bam_class_codes=('class_code', lambda values: '{' + ','.join(sorted({str(v) for v in values if pd.notna(v) and str(v)})) + '}'),
            bam_match_count=('qry_id', 'nunique'),
        )
        .reset_index()
        .rename(columns={'ref_id': 'prediction_protein'})
    )

    candidates_df = pred_loci.merge(grouped_matches, on='prediction_protein', how='inner')
    candidates_df = candidates_df.sort_values(
        ['bam_match_count', 'seqid', 'prediction_start'],
        ascending=[False, True, True],
    ).reset_index(drop=True)
    candidates_df.insert(0, 'candidate_id', [f'HPNGC_{idx:04d}' for idx in range(1, len(candidates_df) + 1)])
    return candidates_df, counts


def report_high_potential_new_gene_candidates(
    prediction_novel_gtf,
    bam_novel_gtf,
    output_dir,
    prediction_protein_fasta=None,
    bam_protein_fasta=None,
    bam_tmap_file=None,
    accepted_class_codes=None,
):
    """Write a TSV of prediction novel proteins with accepted BAM `gffcompare` class-code support."""
    del prediction_protein_fasta, bam_protein_fasta  # kept for backward-compatible calls

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    report_path = output_dir / 'high_potential_new_gene_candidates.tsv'
    legacy_paths = [
        output_dir / 'high_potential_new_gene_candidates_prediction.gtf',
        output_dir / 'high_potential_new_gene_candidates_bam.gtf',
        output_dir / 'high_potential_new_gene_candidates_prediction.fa',
        output_dir / 'high_potential_new_gene_candidates_bam.fa',
    ]
    for legacy_path in legacy_paths:
        if legacy_path.exists():
            legacy_path.unlink()

    candidates_df, counts = _summarize_novel_class_support_from_tmap(
        prediction_novel_gtf=prediction_novel_gtf,
        bam_annotation_gtf=bam_novel_gtf,
        bam_tmap_file=bam_tmap_file,
        accepted_class_codes=accepted_class_codes,
    )
    candidates_df.to_csv(report_path, sep='\t', index=False)

    if candidates_df.empty:
        print(f'No accepted BAM class-code support was found for prediction novel proteins; empty report written to: {report_path}')
    else:
        print(f'High-potential new gene candidate loci written to: {report_path}')
    print(f'Accepted BAM class codes used: {", ".join(counts["accepted_codes"])}')
    print(f'Prediction novel proteins with accepted BAM support: {counts["accepted_overlap"]}/{counts["total"]}')
    return candidates_df


def _summarize_overlap_clusters(prediction_gtf, bam_gtf):
    """Collapse overlapping loci across two branches into prediction-only/shared/BAM-only clusters."""
    empty_loci = pd.DataFrame(columns=['protein', 'gene', 'seqid', 'strand', 'start', 'end'])
    prediction_loci = _summarize_loci_from_gtf(prediction_gtf).copy() if prediction_gtf else empty_loci.copy()
    bam_loci = _summarize_loci_from_gtf(bam_gtf).copy() if bam_gtf else empty_loci.copy()

    frames = []
    if not prediction_loci.empty:
        prediction_loci['branch'] = 'prediction_search'
        frames.append(prediction_loci)
    if not bam_loci.empty:
        bam_loci['branch'] = 'bam_search'
        frames.append(bam_loci)

    counts = {
        'prediction_only': 0,
        'shared': 0,
        'bam_only': 0,
        'prediction_total': int(len(prediction_loci)),
        'bam_total': int(len(bam_loci)),
        'total_clusters': 0,
    }
    if not frames:
        return counts

    loci_df = pd.concat(frames, ignore_index=True).sort_values(
        ['seqid', 'strand', 'start', 'end', 'branch']
    )
    cluster_branches = []

    for (_, _), group in loci_df.groupby(['seqid', 'strand'], dropna=False):
        current_start = None
        current_end = None
        current_branches = set()

        for row in group.itertuples(index=False):
            row_start = int(row.start)
            row_end = int(row.end)
            if current_start is None or row_start > current_end:
                if current_start is not None:
                    cluster_branches.append(set(current_branches))
                current_start = row_start
                current_end = row_end
                current_branches = {row.branch}
            else:
                current_end = max(current_end, row_end)
                current_branches.add(row.branch)

        if current_start is not None:
            cluster_branches.append(set(current_branches))

    counts['total_clusters'] = len(cluster_branches)
    for branches in cluster_branches:
        if branches == {'prediction_search'}:
            counts['prediction_only'] += 1
        elif branches == {'bam_search'}:
            counts['bam_only'] += 1
        else:
            counts['shared'] += 1

    return counts


def _summarize_supported_overlap_from_classes(summary_df, accepted_codes=None):
    """Summarize BAM-vs-prediction overlap using the selected gffcompare class codes."""
    accepted_codes = list(accepted_codes or ['=', 'j', 'c', 'k', 'm', 'n'])
    counts = {
        'accepted_overlap': 0,
        'other_classes': 0,
        'total': 0,
        'accepted_codes': accepted_codes,
    }

    if summary_df is None or getattr(summary_df, 'empty', True):
        return counts
    if not {'class_code', 'count'}.issubset(summary_df.columns):
        return counts

    class_df = summary_df[['class_code', 'count']].copy()
    class_df['class_code'] = class_df['class_code'].astype(str).str.strip()
    class_df['count'] = pd.to_numeric(class_df['count'], errors='coerce').fillna(0).astype(int)

    accepted_set = set(accepted_codes)
    accepted_overlap = int(class_df[class_df['class_code'].isin(accepted_set)]['count'].sum())
    total = int(class_df['count'].sum())

    counts['accepted_overlap'] = accepted_overlap
    counts['other_classes'] = max(total - accepted_overlap, 0)
    counts['total'] = total
    return counts


def _plot_supported_class_pie(ax, counts, title):
    """Draw the supported-overlap pie using accepted gffcompare class codes."""
    ax.set_title(title, fontweight='bold')
    total = counts['total']

    if total <= 0:
        ax.axis('off')
        ax.text(0.5, 0.56, 'No supported transcripts', ha='center', va='center', fontsize=12, transform=ax.transAxes)
        return

    code_label = ', '.join(counts['accepted_codes'])
    labels = [f'Accepted overlap\n({code_label})', 'Other classes']
    values = [counts['accepted_overlap'], counts['other_classes']]
    colors = ['#72B7B2', '#E45756']
    filtered = [(label, value, color) for label, value, color in zip(labels, values, colors) if value > 0]

    ax.pie(
        [item[1] for item in filtered],
        labels=[item[0] for item in filtered],
        colors=[item[2] for item in filtered],
        startangle=90,
        autopct=lambda pct: f"{pct:.1f}%\n({int(round(total * pct / 100.0))})",
        wedgeprops={'edgecolor': 'white'},
        textprops={'fontsize': 10},
    )
    ax.text(
        0.5,
        -0.08,
        f"Total BAM-supported transcripts: {total}",
        ha='center',
        va='top',
        fontsize=9,
        transform=ax.transAxes,
    )


def _plot_novel_class_pie(ax, counts, title):
    """Draw the novel-support pie using accepted BAM `gffcompare` class codes."""
    ax.set_title(title, fontweight='bold')
    total = counts['total']

    if total <= 0:
        ax.axis('off')
        ax.text(0.5, 0.56, 'No prediction novel proteins', ha='center', va='center', fontsize=12, transform=ax.transAxes)
        return

    code_label = ', '.join(counts['accepted_codes'])
    labels = [f'Accepted BAM support\n({code_label})', 'No accepted BAM support']
    values = [counts['accepted_overlap'], counts['other_classes']]
    colors = ['#72B7B2', '#E45756']
    filtered = [(label, value, color) for label, value, color in zip(labels, values, colors) if value > 0]

    ax.pie(
        [item[1] for item in filtered],
        labels=[item[0] for item in filtered],
        colors=[item[2] for item in filtered],
        startangle=90,
        autopct=lambda pct: f"{pct:.1f}%\n({int(round(total * pct / 100.0))})",
        wedgeprops={'edgecolor': 'white'},
        textprops={'fontsize': 10},
    )
    ax.text(
        0.5,
        -0.08,
        f"Prediction novel proteins: {total} | BAM transcripts in accepted classes: {counts['matched_bam_transcripts']}",
        ha='center',
        va='top',
        fontsize=9,
        transform=ax.transAxes,
    )


def plot_bam_gtf_comparison_summary(
    summary_df,
    no_overlap_count,
    output_dir,
    prediction_supported_gtf=None,
    bam_supported_gtf=None,
    prediction_novel_gtf=None,
    bam_novel_gtf=None,
    bam_tmap_file=None,
):
    """Create the requested two-pie summary of supported and novel overlap using `gffcompare` class codes."""
    del no_overlap_count, prediction_supported_gtf, bam_novel_gtf  # retained for API compatibility

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    supported_counts = _summarize_supported_overlap_from_classes(summary_df)
    _, novel_counts = _summarize_novel_class_support_from_tmap(
        prediction_novel_gtf=prediction_novel_gtf,
        bam_annotation_gtf=bam_supported_gtf,
        bam_tmap_file=bam_tmap_file,
        accepted_class_codes=supported_counts['accepted_codes'],
    )

    if summary_df.empty and novel_counts['total'] == 0:
        return None

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.6))
    _plot_supported_class_pie(
        axes[0],
        supported_counts,
        'Supported overlap by gffcompare classes\naccepted set: =, j, c, k, m, n'
    )
    _plot_novel_class_pie(
        axes[1],
        novel_counts,
        'Prediction novel proteins with BAM class-code support'
    )

    fig.suptitle('Comparison summary: bam_search vs prediction_search', fontsize=13, fontweight='bold')
    fig.text(
        0.5,
        0.02,
        'Both pies use accepted `gffcompare` classes: =, j, c, k, m, n  |  Detailed counts remain in bam_vs_input_gtf_summary.tsv',
        ha='center',
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.05, 1, 0.93])
    output_path = output_dir / 'bam_vs_input_gtf_summary.png'
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return output_path


def compare_bam_support_to_input_gtf(
    input_gtf,
    bam_supported_gtf,
    output_dir,
    bam_protein_fasta=None,
    prediction_supported_gtf=None,
    prediction_novel_gtf=None,
    bam_novel_gtf=None,
):
    """Compare BAM-supported transcripts against the user's original input GTF.

    Writes `gffcmp.*` outputs plus a summary of class codes and the subset of BAM-supported
    transcripts with no overlap to any gene in the original input GTF (`class_code == 'u'`).
    """
    input_gtf = Path(input_gtf)
    bam_supported_gtf = Path(bam_supported_gtf)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    run_gffcompare(input_gtf, bam_supported_gtf, output_dir)

    tmap_candidates = sorted(output_dir.glob("*.tmap"))
    if not tmap_candidates:
        raise FileNotFoundError(f"No gffcompare .tmap file was produced in {output_dir}")

    tmap_file = tmap_candidates[0]
    print(f"Loading BAM-vs-input TMAP file: {tmap_file}")
    tmap_df = pd.read_csv(tmap_file, sep='\t')

    code_descriptions = {
        '=': 'exact match to input GTF',
        'c': 'contained within input transcript',
        'j': 'novel isoform with shared splice junctions',
        'e': 'single exon transfrag partially covering a reference intron',
        'i': 'fully contained within a reference intron',
        'k': 'contains a reference transcript (reverse containment)',
        'm': 'retained intron(s); all introns matched or retained',
        'n': 'retained intron(s); not all introns matched/covered',
        'o': 'other same-strand overlap with reference exons',
        'p': 'possible polymerase run-on fragment',
        'r': 'repeat-associated transcript',
        'u': 'no overlap with any input-GTF gene locus',
        'x': 'exonic overlap on opposite strand',
        's': 'intron match on the opposite strand',
        'y': 'contains a reference within its intron(s)',
    }

    summary_df = (
        tmap_df['class_code']
        .fillna('NA')
        .value_counts()
        .rename_axis('class_code')
        .reset_index(name='count')
    )
    total = int(summary_df['count'].sum()) if not summary_df.empty else 0
    summary_df['description'] = summary_df['class_code'].map(code_descriptions).fillna('other/unknown')
    summary_df['fraction'] = (summary_df['count'] / total).round(4) if total else 0.0
    summary_file = output_dir / 'bam_vs_input_gtf_summary.tsv'
    summary_df.to_csv(summary_file, sep='\t', index=False)

    no_overlap_df = tmap_df[tmap_df['class_code'].fillna('').eq('u')].copy()
    no_overlap_file = output_dir / 'bam_supported_no_gene_overlap.tsv'
    no_overlap_df.to_csv(no_overlap_file, sep='\t', index=False)

    no_overlap_ids = set(no_overlap_df['qry_id'].dropna().astype(str))
    subset_gtf_file = output_dir / 'bam_supported_no_gene_overlap.gtf'
    subset_faa_file = output_dir / 'bam_supported_no_gene_overlap.faa'

    bam_gtf_df = gtf_to_df_with_genes(bam_supported_gtf)
    no_overlap_gtf_df = bam_gtf_df[
        bam_gtf_df['Protein'].isin(no_overlap_ids) | bam_gtf_df['Gene'].isin(no_overlap_ids)
    ].copy()
    no_overlap_gtf_df.drop(columns=['Protein', 'Gene'], errors='ignore').to_csv(
        subset_gtf_file,
        sep='\t',
        index=False,
        header=False,
    )

    with open(subset_faa_file, 'w') as fasta_handle:
        if bam_protein_fasta and no_overlap_ids:
            with open(bam_protein_fasta) as protein_handle:
                selected = [
                    record for record in SeqIO.parse(protein_handle, 'fasta')
                    if record.id.split()[0] in no_overlap_ids
                ]
            SeqIO.write(selected, fasta_handle, 'fasta')

    plot_file = plot_bam_gtf_comparison_summary(
        summary_df,
        len(no_overlap_df),
        output_dir,
        prediction_supported_gtf=prediction_supported_gtf or input_gtf,
        bam_supported_gtf=bam_supported_gtf,
        prediction_novel_gtf=prediction_novel_gtf,
        bam_novel_gtf=bam_novel_gtf,
        bam_tmap_file=tmap_file,
    )

    print(f"BAM-vs-input GTF summary written to: {summary_file}")
    if plot_file is not None:
        print(f"BAM-vs-input summary figure written to: {plot_file}")
    print(f"BAM-supported transcripts with no input-GTF overlap: {len(no_overlap_df)}")

    return summary_df, no_overlap_df

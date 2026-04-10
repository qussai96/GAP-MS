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


def plot_bam_gtf_comparison_summary(summary_df, no_overlap_count, output_dir):
    """Create a compact summary figure for BAM-vs-input-GTF comparison results."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if summary_df.empty:
        return None

    plot_df = summary_df.copy()
    plot_df['label'] = plot_df['class_code'].astype(str) + "\n" + plot_df['description'].str.replace(' ', '\n', regex=False)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), gridspec_kw={'width_ratios': [1.6, 1.0]})

    bars = axes[0].bar(plot_df['class_code'], plot_df['count'], color='#4C78A8', edgecolor='black')
    axes[0].set_title('BAM-supported transcripts vs input GTF', fontweight='bold')
    axes[0].set_xlabel('gffcompare class code')
    axes[0].set_ylabel('Transcript count')
    axes[0].grid(axis='y', linestyle='--', alpha=0.35)
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)
    for bar, frac in zip(bars, plot_df['fraction']):
        axes[0].text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{int(bar.get_height())}\n({frac:.1%})",
            ha='center',
            va='bottom',
            fontsize=9,
            fontweight='bold'
        )

    total = int(plot_df['count'].sum())
    overlap_count = total - int(no_overlap_count)
    pie_values = [overlap_count, int(no_overlap_count)]
    pie_labels = ['Overlaps input GTF', 'No gene overlap (u)']
    pie_colors = ['#72B7B2', '#E45756']
    axes[1].pie(
        pie_values,
        labels=pie_labels,
        autopct=lambda pct: f"{pct:.1f}%\n({round(total * pct / 100):d})" if total else '0',
        colors=pie_colors,
        startangle=90,
        wedgeprops={'edgecolor': 'white'}
    )
    axes[1].set_title('Interesting BAM-only support', fontweight='bold')

    fig.suptitle('Comparison summary: bam_search vs prediction_search', fontsize=13, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    output_path = output_dir / 'bam_vs_input_gtf_summary.png'
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return output_path


def compare_bam_support_to_input_gtf(input_gtf, bam_supported_gtf, output_dir, bam_protein_fasta=None):
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
        'e': 'single exon overlap with input transcript',
        'i': 'fully contained in an intron of an input transcript',
        'o': 'generic exonic overlap with input transcript',
        'p': 'possible polymerase run-on fragment',
        'r': 'repeat-associated transcript',
        'u': 'no overlap with any input-GTF gene locus',
        'x': 'exonic overlap on opposite strand',
        's': 'intron overlap on opposite strand',
        'y': 'contains a reference within its intron',
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

    plot_file = plot_bam_gtf_comparison_summary(summary_df, len(no_overlap_df), output_dir)

    print(f"BAM-vs-input GTF summary written to: {summary_file}")
    if plot_file is not None:
        print(f"BAM-vs-input summary figure written to: {plot_file}")
    print(f"BAM-supported transcripts with no input-GTF overlap: {len(no_overlap_df)}")

    return summary_df, no_overlap_df

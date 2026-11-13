"""Strand- and phase-aware TransDecoder projection, retained from GAPMS (MIT)."""
from pathlib import Path

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
            # gffread emits unstranded StringTie transcripts in ascending
            # genomic order, so "." must be projected like a plus-strand
            # transcript here.
            if strand != '-':
                genomic_start = exon_start + (overlap_start - transcript_exon_start)
                genomic_end = exon_start + (overlap_end - transcript_exon_start)
            else:
                genomic_end = exon_end - (overlap_start - transcript_exon_start)
                genomic_start = exon_end - (overlap_end - transcript_exon_start)

            genomic_segments.append((min(genomic_start, genomic_end), max(genomic_start, genomic_end)))

        transcript_cursor = transcript_exon_end + 1

    return sorted(genomic_segments, key=lambda segment: segment[0])

def _genomic_orf_strand(transcript_strand, orf_strand):
    """Compose transcript and TransDecoder ORF orientation on the genome."""
    if orf_strand not in {'+', '-'}:
        return transcript_strand
    if transcript_strand == '-':
        return '+' if orf_strand == '-' else '-'
    # A "." StringTie transcript is extracted by gffread in genomic-plus
    # orientation, so the ORF strand directly determines its genomic strand.
    return orf_strand

def _next_cds_phase(segment_length, phase):
    """Return the GFF phase for the next CDS segment in translation order."""
    coding_length = max(int(segment_length) - int(phase), 0)
    return (3 - (coding_length % 3)) % 3

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
                    {
                        'transcript_id': transcript_id,
                        'source': source,
                        'score': score,
                        'orf_strand': strand,
                        'cds_features': [],
                    },
                )
                entry['transcript_id'] = transcript_id
                if strand in {'+', '-'}:
                    entry['orf_strand'] = strand
                if score not in {None, ''}:
                    entry['score'] = score
            elif feature_type == 'cds':
                orf_id = attrs.get('Parent') or attrs.get('protein_id') or attrs.get('transcript_id')
                if not orf_id:
                    continue
                entry = orf_records.setdefault(
                    orf_id,
                    {
                        'transcript_id': transcript_id,
                        'source': source,
                        'score': score,
                        'orf_strand': strand,
                        'cds_features': [],
                    },
                )
                entry['transcript_id'] = transcript_id
                entry['cds_features'].append((int(start), int(end), frame if frame != '.' else '0'))
                if strand in {'+', '-'}:
                    entry['orf_strand'] = strand
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
            orf_strand = entry.get('orf_strand', '+')
            strand = _genomic_orf_strand(transcript_model['strand'], orf_strand)
            score = transcript_model.get('score', '.')
            attributes = f'gene_id "{gene_id}"; transcript_id "{orf_id}"; protein_id "{orf_id}";'

            cds_features = sorted(
                entry['cds_features'],
                key=lambda feature: feature[0],
                reverse=orf_strand == '-',
            )
            projected_segments = []
            for cds_start, cds_end, _ in cds_features:
                genomic_segments = _map_transcript_interval_to_genome(cds_start, cds_end, transcript_model)
                projected_segments.extend(genomic_segments)

            if not projected_segments:
                continue

            projected_segments.sort(key=lambda segment: segment[0], reverse=strand == '-')
            try:
                phase = int(cds_features[0][2])
            except (TypeError, ValueError):
                phase = 0
            if phase not in {0, 1, 2}:
                phase = 0

            handle.write(
                f"{seqid}\tGAPMS\ttranscript\t{transcript_model['start']}\t{transcript_model['end']}\t{score}\t{strand}\t.\t{attributes}\n"
            )
            for genomic_start, genomic_end in projected_segments:
                handle.write(
                    f"{seqid}\tGAPMS\tCDS\t{genomic_start}\t{genomic_end}\t{score}\t{strand}\t{phase}\t{attributes}\n"
                )
                phase = _next_cds_phase(genomic_end - genomic_start + 1, phase)

            converted_orfs += 1

    if converted_orfs == 0:
        raise RuntimeError(
            f"No TransDecoder ORFs could be projected back to genome coordinates from {transdecoder_gff3}"
        )

    return output_gtf

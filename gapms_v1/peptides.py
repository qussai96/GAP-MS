"""Preserve identification metadata without inventing peptide confidence."""
import csv
import gzip
import json
from pathlib import Path

from gapms._peplocus.io import read_peptides, normalize_peptide, HEADERS
from .common import stable_id, write_table


def load_peptides(args):
    peptides, stats = read_peptides(args.peptides, args.peptide_column)
    kept = [p for p in peptides if len(p) >= args.min_length]
    if not kept:raise ValueError('no peptides meet --min-length')
    metadata = {p:dict(peptide_id=stable_id('P',p), quality='unknown', qvalues=[], records=0) for p in peptides}
    opener = gzip.open if str(args.peptides).endswith('.gz') else open
    with opener(args.peptides, 'rt', encoding='utf-8-sig') as stream:
        first = next((line for line in stream if line.strip() and not line.startswith('#')), '')
    delimiter = '\t' if '\t' in first else ',' if ',' in first else None
    header = [part.strip().lower() for part in first.strip().split(delimiter)] if delimiter else []
    is_table = any(name in HEADERS for name in header)
    rows = []
    if is_table:
        with opener(args.peptides, 'rt', encoding='utf-8-sig') as stream:
            reader = csv.DictReader((line for line in stream if line.strip() and not line.startswith('#')), delimiter=delimiter)
            peptide_key = next((key for key in reader.fieldnames if key.strip().lower() == (args.peptide_column or '').lower()), None)
            peptide_key = peptide_key or next(key for key in reader.fieldnames if key.strip().lower() in HEADERS)
            for index, record in enumerate(reader, 1):
                peptide = normalize_peptide(record[peptide_key]); item=metadata[peptide]; item['records'] += 1
                qkeys = [key for key in record if key.strip().lower().replace('_',' ').replace('-',' ') in {'q value','qvalue','peptide q value','peptide qvalue','global q value'}]
                for key in qkeys:
                    value = record[key].strip()
                    if value:
                        qvalue = float(value)
                        if not 0 <= qvalue <= 1:raise ValueError(f'invalid q-value {value}')
                        item['qvalues'].append(qvalue)
                rows.append(dict(identification_id=stable_id('I',str(args.peptides.resolve()),index),peptide_id=item['peptide_id'],peptide=peptide,original_record=json.dumps(record,sort_keys=True),input_row=index))
    else:
        for index, peptide in enumerate(peptides,1):
            metadata[peptide]['records']=1
            rows.append(dict(identification_id=stable_id('I',str(args.peptides.resolve()),index),peptide_id=metadata[peptide]['peptide_id'],peptide=peptide,original_record='',input_row=index))
    for item in metadata.values():
        if item['qvalues']:
            item['quality'] = 'qvalue_pass' if min(item['qvalues']) <= args.max_peptide_qvalue else 'qvalue_fail'
    kept_set = set(kept)
    directory=args.output/'mapping';directory.mkdir(exist_ok=True)
    write_table(directory/'peptide_identifications.tsv',['identification_id','peptide_id','peptide','input_row','original_record'],rows)
    write_table(directory/'excluded_peptides.tsv',['peptide','length','reason'],(dict(peptide=p,length=len(p),reason='below_min_length') for p in peptides if p not in kept_set))
    normalized=args.output/'mapping/queries.txt';normalized.write_text('\n'.join(kept)+'\n')
    stats.update(search_peptides=len(kept),excluded_peptides=len(peptides)-len(kept))
    return kept, metadata, stats, normalized

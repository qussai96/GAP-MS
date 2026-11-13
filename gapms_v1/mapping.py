"""One disk-backed mapping space shared by all annotation sources and ORFs."""
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import shutil
import sqlite3

from peplocus.sequence import canonical, join_blocks, translate
from peplocus.digestion import protein_context
from .protein_mapper import ProteinMatcher
from .common import stable_id, write_table, read_table
from .models import write_proteins


def locus_key(chrom, strand, blocks):
    merged = []
    for left, right in sorted(blocks):
        if merged and left == merged[-1][1]:
            merged[-1][1] = right
        else:
            merged.append([left, right])
    return (chrom, strand, json.dumps(merged, separators=(',',':')))


class EvidenceStore:
    def __init__(self, path, require_tryptic=False):
        self.require_tryptic = require_tryptic
        self.db = sqlite3.connect(path)
        self.db.execute('PRAGMA journal_mode=WAL')
        self.db.execute('PRAGMA synchronous=NORMAL')
        self.db.executescript('''
          CREATE TABLE hits (peptide TEXT, canonical TEXT, chrom TEXT, strand TEXT, blocks TEXT,
            matched TEXT, evidence TEXT, PRIMARY KEY(peptide,chrom,strand,blocks)) WITHOUT ROWID;
          CREATE TABLE tryptic_hits (canonical TEXT, chrom TEXT, strand TEXT, blocks TEXT,
            PRIMARY KEY(canonical,chrom,strand,blocks)) WITHOUT ROWID;
          CREATE TABLE associations (model TEXT, peptide TEXT, canonical TEXT, chrom TEXT, strand TEXT,
            blocks TEXT, aa_start INTEGER, aa_end INTEGER,
            PRIMARY KEY(model,peptide,blocks,aa_start)) WITHOUT ROWID;
          CREATE TABLE tryptic_associations (model TEXT, peptide TEXT, blocks TEXT, aa_start INTEGER,
            PRIMARY KEY(model,peptide,blocks,aa_start)) WITHOUT ROWID;
          CREATE INDEX hits_canonical ON hits(canonical);
          CREATE INDEX associations_canonical ON associations(canonical);
        ''')

        condition = 'WHERE EXISTS (SELECT 1 FROM tryptic_hits t WHERE t.canonical=h.canonical AND t.chrom=h.chrom AND t.strand=h.strand AND t.blocks=h.blocks)' if require_tryptic else ''
        self.db.execute('CREATE VIEW active_hits AS SELECT h.* FROM hits h '+condition)
        condition = 'WHERE EXISTS (SELECT 1 FROM tryptic_associations t WHERE t.model=a.model AND t.peptide=a.peptide AND t.blocks=a.blocks AND t.aa_start=a.aa_start)' if require_tryptic else ''
        self.db.execute('CREATE VIEW active_associations AS SELECT a.* FROM associations a '+condition)

    def add(self, peptide, chrom, strand, blocks, matched, evidence, tryptic=False):
        chrom, strand, blocks = locus_key(chrom,strand,blocks)
        self.db.execute('''INSERT INTO hits VALUES (?,?,?,?,?,?,?) ON CONFLICT(peptide,chrom,strand,blocks)
          DO UPDATE SET evidence=CASE WHEN instr(','||hits.evidence||',',','||excluded.evidence||',')>0
          THEN hits.evidence ELSE hits.evidence||','||excluded.evidence END''',
          (peptide,canonical(peptide),chrom,strand,blocks,matched,evidence))
        if tryptic:self.db.execute('INSERT OR IGNORE INTO tryptic_hits VALUES (?,?,?,?)',(canonical(peptide),chrom,strand,blocks))

    def import_genomic(self, path):
        for row in read_table(path):
            self.add(row['peptide'],row['chrom'],row['strand'],json.loads(row['blocks']),row['matched_sequence'],row['evidence'],tryptic=row.get('query_tryptic')=='True')
        self.db.commit()

    def associate(self, models, peptides, directory, reference, genetic_code, matcher=None):
        """Search the union of all proteins once; validate every projected hit."""
        proteins = directory/'all_model_proteins.faa'
        write_proteins(models, proteins)
        lookup = {m.id:m for m in models}
        if not any(m.protein for m in models):return
        (matcher or ProteinMatcher(peptides)).map(proteins, directory/'protein_mapping.tsv')
        for row in read_table(directory/'protein_mapping.tsv'):
            if row['protein'] == 'UNMAPPED':continue
            model = lookup[row['protein']]; peptide=row['peptide']
            start = int(row['location']) - 1; end = start + len(peptide)
            blocks = model.project(start*3,end*3)
            matched = translate(join_blocks(reference,model.chrom,model.strand,blocks),genetic_code)
            if canonical(matched) != canonical(peptide):
                raise RuntimeError(f'CDS projection verification failed: {model.id} {peptide}')
            context=protein_context(model.protein+('*' if model.complete_stop else ''),start,end,peptide)
            self.add(peptide,model.chrom,model.strand,blocks,matched,model.source,tryptic=context['query_tryptic'])
            chrom,strand,encoded=locus_key(model.chrom,model.strand,blocks)
            self.db.execute('INSERT OR IGNORE INTO associations VALUES (?,?,?,?,?,?,?,?)',
                            (model.id,peptide,canonical(peptide),chrom,strand,encoded,start,end))
            if context['query_tryptic']:
                self.db.execute('INSERT OR IGNORE INTO tryptic_associations VALUES (?,?,?,?)', (model.id,peptide,encoded,start))
        self.db.commit()

    def counts(self, all_mappings=False):
        # I/L-equivalent observed sequences constitute one evidence sequence.
        table='hits' if all_mappings else 'active_hits'
        return dict(self.db.execute(f'SELECT canonical,count(*) FROM (SELECT DISTINCT canonical,chrom,strand,blocks FROM {table}) GROUP BY canonical'))

    def model_hits(self, all_mappings=False, model_ids=None):
        result=defaultdict(list)
        table='associations' if all_mappings else 'active_associations'
        for model,peptide,sequence,chrom,strand,encoded,start,end in self.db.execute(f'SELECT * FROM {table} ORDER BY model,canonical,aa_start'):
            if model_ids is not None and model not in model_ids:continue
            result[model].append(dict(peptide=peptide,canonical=sequence,chrom=chrom,strand=strand,blocks=json.loads(encoded),aa_start=start,aa_end=end,mapping_id=stable_id('L',sequence,chrom,strand,encoded)))
        return result

    def export(self, directory, peptides, metadata):
        counts=self.counts(all_mappings=True)
        analysis_counts=self.counts()
        tryptic_counts=dict(self.db.execute('SELECT canonical,count(*) FROM tryptic_hits GROUP BY canonical'))
        fields=['mapping_id','peptide_id','peptide','chrom','start','end','strand','blocks','matched_sequence','match_type','mapping_type','evidence','loci_for_peptide','identification_quality','tryptic_compatible','tryptic_loci_for_peptide']
        mapped=set();stats=Counter()
        with (directory/'mapped_peptides.tsv').open('w') as tsv, (directory/'mapped_peptides.bed').open('w') as bed, (directory/'mapped_peptides.igv.bed').open('w') as igv, (directory/'splice_peptides.bed').open('w') as splice, (directory/'tryptic_peptides.bed').open('w') as tryptic_bed:
            writer=csv.DictWriter(tsv,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
            igv.write('track name="GAPMS_v1" itemRgb="On" description="Display only; reverse-strand labels reversed"\n')
            query="SELECT h.*,t.canonical IS NOT NULL FROM hits h LEFT JOIN tryptic_hits t ON (h.canonical=t.canonical AND h.chrom=t.chrom AND h.strand=t.strand AND h.blocks=t.blocks) ORDER BY h.chrom,CAST(json_extract(h.blocks,'$[0][0]') AS INTEGER),h.strand,h.peptide,h.blocks"
            for peptide,sequence,chrom,strand,encoded,matched,evidence,compatible in self.db.execute(query):
                blocks=json.loads(encoded);left=blocks[0][0];right=blocks[-1][1]
                unique=counts[sequence]==1;spliced=len(blocks)>1
                row=dict(mapping_id=stable_id('L',sequence,chrom,strand,encoded),peptide_id=metadata[peptide]['peptide_id'],peptide=peptide,
                         chrom=chrom,start=left,end=right,strand=strand,blocks=encoded,matched_sequence=matched,
                         match_type='exact' if peptide==matched else 'il_equivalent',mapping_type='spliced' if spliced else 'contiguous',
                         evidence=evidence,loci_for_peptide=counts[sequence],identification_quality=metadata[peptide]['quality'],tryptic_compatible=bool(compatible),tryptic_loci_for_peptide=tryptic_counts.get(sequence,0))
                writer.writerow(row)
                fields_bed=[chrom,left,right,peptide,1000 if unique else 0,strand,left,right,'210,100,35' if spliced else '45,110,180',len(blocks),','.join(str(b-a) for a,b in blocks),','.join(str(a-left) for a,b in blocks)]
                line='\t'.join(map(str,fields_bed))+'\n';bed.write(line)
                if spliced:splice.write(line)
                if compatible:
                    tryptic_fields=list(fields_bed);tryptic_fields[4]=1000 if tryptic_counts[sequence]==1 else 0
                    tryptic_bed.write('\t'.join(map(str,tryptic_fields))+'\n');stats['tryptic_mapping_rows']+=1
                fields_bed[3]=peptide[::-1] if strand=='-' else peptide
                igv.write('\t'.join(map(str,fields_bed))+'\n')
                mapped.add(peptide);stats['mapping_rows']+=1
                if spliced:stats['spliced_mapping_rows']+=1
        write_table(directory/'unmapped_peptides.tsv',['peptide_id','peptide','length','status'],
                    (dict(peptide_id=metadata[p]['peptide_id'],peptide=p,length=len(p),status='unmapped_in_searched_space') for p in peptides if p not in mapped))
        stats.update(dict(peptides=len(peptides),mapped_peptides=len(mapped),unmapped_peptides=len(peptides)-len(mapped),
                          unique_peptides=sum(counts.get(canonical(p))==1 for p in peptides),multiple_peptides=sum(counts.get(canonical(p),0)>1 for p in peptides),
                          canonical_evidence_sequences=len({canonical(p) for p in peptides}),
                          tryptic_mapped_peptides=sum(tryptic_counts.get(canonical(p),0)>0 for p in peptides),
                          tryptic_unique_peptides=sum(tryptic_counts.get(canonical(p),0)==1 for p in peptides),
                          tryptic_multiple_peptides=sum(tryptic_counts.get(canonical(p),0)>1 for p in peptides)))
        write_table(directory/'peptide_locus_summary.tsv',['peptide','length','loci','status'],
                    (dict(peptide=p,length=len(p),loci=counts.get(canonical(p),0),status='unique' if counts.get(canonical(p))==1 else 'multiple' if counts.get(canonical(p),0)>1 else 'unmapped') for p in peptides))
        stats=dict(stats)
        stats.update(analysis_evidence='tryptic' if self.require_tryptic else 'all_exact',
                     analysis_mapped_peptides=sum(analysis_counts.get(canonical(p),0)>0 for p in peptides),
                     analysis_unique_peptides=sum(analysis_counts.get(canonical(p),0)==1 for p in peptides),
                     analysis_multiple_peptides=sum(analysis_counts.get(canonical(p),0)>1 for p in peptides))
        return dict(stats)

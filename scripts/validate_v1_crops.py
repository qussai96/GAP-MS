"""Independent checks for completed GAPMS v1 crop runs; BioPython translation."""
from collections import Counter
import csv
import json
from pathlib import Path
import random
import re
from Bio import SeqIO
import sqlite3
import sys
import traceback

from Bio.Seq import Seq
import pysam


def rows(path):
 with path.open() as stream:yield from csv.DictReader(stream,delimiter='\t')


def validate(record):
 out=Path(record['output']);manifest=json.loads((out/'run_manifest.json').read_text())
 assert manifest['status']=='completed',manifest.get('error',manifest['status'])
 summary=json.loads((out/'reports/summary.json').read_text())
 statistics=summary['mapping'];assert statistics['mapped_peptides']+statistics['unmapped_peptides']==statistics['peptides']
 identification_rows=sum(1 for _ in rows(out/'mapping/peptide_identifications.tsv'))
 assert identification_rows==summary['input']['input_records'], 'original identification records were dropped'
 assert statistics['unique_peptides']+statistics['multiple_peptides']==statistics['mapped_peptides']
 rng=random.Random(42);sample=[];count=0;bed_count=0;spliced=0;lengths=Counter()
 for count,row in enumerate(rows(out/'mapping/mapped_peptides.tsv'),1):
  blocks=json.loads(row['blocks']);assert blocks and all(a<b for a,b in blocks)
  assert int(row['start'])==blocks[0][0] and int(row['end'])==blocks[-1][1]
  assert all(a[1]<=b[0] for a,b in zip(blocks,blocks[1:]))
  assert row['strand'] in '+-' and row['identification_quality']=='unknown'
  if len(blocks)>1:spliced+=1
  if len(sample)<2000:sample.append(row)
  else:
   chosen=rng.randrange(count)
   if chosen<len(sample):sample[chosen]=row
 assert count==statistics['mapping_rows']
 with (out/'mapping/mapped_peptides.bed').open() as stream:
  for line in stream:
   fields=line.rstrip('\n').split('\t');assert len(fields)==12
   assert len(fields[10].split(','))==int(fields[9])==len(fields[11].split(','))
   bed_count+=1
 assert bed_count==count
 assert sum(1 for _ in rows(out/'mapping/unmapped_peptides.tsv'))==statistics['unmapped_peptides']
 index=next((out/'work/reference').glob('*.fai'))
 with pysam.FastaFile(record['assembly'],filepath_index=str(index)) as reference:
  for row in sample:
   nt=''.join(reference.fetch(row['chrom'],a,b) for a,b in json.loads(row['blocks']))
   if row['strand']=='-':nt=str(Seq(nt).reverse_complement())
   protein=str(Seq(nt).translate(table=manifest['parameters']['genetic_code']))
   assert protein==row['matched_sequence'],row
   assert protein.replace('I','L')==row['peptide'].replace('I','L'),row
 # Independently translate whole candidate ORFs, not just peptide substrings.
 model_identity={r['model_id']:r['original_transcript_id'] for r in rows(out/'orfs/orf_peptide_evidence.tsv')}
 proteins={r.id:str(r.seq) for r in SeqIO.parse(out/'orfs/candidate_orfs.faa','fasta')}
 sampled_models=rng.sample(sorted(proteins),min(2000,len(proteins)))
 selected_original={model_identity[mid]:mid for mid in sampled_models}
 cds={tid:[] for tid in selected_original}
 with (out/'orfs/candidate_orfs.gtf').open() as stream:
  for line in stream:
   if line.startswith('#') or not line.strip():continue
   fields=line.rstrip('\n').split('\t')
   if fields[2]!='CDS':continue
   match=re.search(r'transcript_id "([^"]+)"',fields[8]);assert match
   tid=match.group(1)
   if tid in cds:cds[tid].append((fields[0],int(fields[3])-1,int(fields[4]),fields[6],fields[7]))
 with pysam.FastaFile(record['assembly'],filepath_index=str(index)) as reference:
  for tid,segments in cds.items():
   assert segments
   segments.sort(key=lambda item:item[1]);strand=segments[0][3]
   nt=''.join(reference.fetch(chrom,a,b) for chrom,a,b,_,phase in segments)
   if strand=='-':nt=str(Seq(nt).reverse_complement())
   phase=segments[0 if strand=='+' else -1][4];nt=nt[int(phase) if phase!='.' else 0:]
   protein=str(Seq(nt[:len(nt)//3*3]).translate(table=manifest['parameters']['genetic_code']))
   if protein.endswith('*'):protein=protein[:-1]
   assert protein==proteins[selected_original[tid]],tid
 for row in rows(out/'mapping/peptide_locus_summary.tsv'):
  lengths[(len(row['peptide']),row['status'])]+=1
 events=Counter();confidence=Counter();diagnostic_count=0
 connection=sqlite3.connect(f'file:{out}/work/evidence.sqlite?mode=ro',uri=True)
 for row in rows(out/'refinement/annotation_change_proposals.tsv'):
  assert row['confidence']!='high_confidence_candidate','Unknown-quality Casanovo sequences must not get high confidence'
  events[row['event_type']]+=1;confidence[row['confidence']]+=1
  sequences=json.loads(row['peptides']);assert sequences
  model=row['candidate_model_id']
  if set(proteins[model])-set('ACDEFGHIKLMNPQRSTVWY*'):
   assert 'unresolved_amino_acids_in_candidate_ORF' in json.loads(row['limitations'])
  for sequence in sequences:
   found=connection.execute('SELECT blocks,aa_start FROM associations WHERE model=? AND canonical=?',(model,sequence)).fetchall();assert found
   if row['event_type']=='alternative_start':assert any(start in [0,1] for encoded,start in found)
   if row['event_type']=='alternative_splice':
    ref_junctions=set()
    for blocks in json.loads(row['reference_cds_blocks']).values():
     ref_junctions.update((a[1],b[0]) for a,b in zip(blocks,blocks[1:]))
    assert any(any((a[1],b[0]) not in ref_junctions for a,b in zip(json.loads(encoded),json.loads(encoded)[1:])) for encoded,start in found)
   diagnostic_count+=1
 connection.close()
 assert dict(events)==summary['assessment']['annotation_events']
 from reference_comparison import compare
 reference=compare(record)
 return dict(reference_comparison=reference,crop=record['crop'],job_id=record['job_id'],status='validated',peptides=statistics['peptides'],mapped_peptides=statistics['mapped_peptides'],
             mapped_percent=round(100*statistics['mapped_peptides']/statistics['peptides'],3),unique_peptides=statistics['unique_peptides'],multiple_peptides=statistics['multiple_peptides'],
             mapping_rows=count,spliced_mapping_rows=spliced,independently_translated_sample=len(sample),independently_translated_orfs=len(sampled_models),diagnostic_event_peptide_checks=diagnostic_count,
             seconds=summary['seconds'],assessment=summary['assessment'],ambiguity_by_length=[dict(length=n,status=s,peptides=c) for (n,s),c in sorted(lengths.items())])


if __name__=='__main__':
 manifest=Path(sys.argv[1]);records=json.loads(manifest.read_text());results=[]
 for record in records:
  try:result=validate(record)
  except Exception as error:
   result=dict(crop=record['crop'],job_id=record['job_id'],status='failed_validation',error=str(error),traceback=traceback.format_exc())
  results.append(result);print(result['crop'],result['status'],flush=True)
  (manifest.parent/'validation.json').write_text(json.dumps(results,indent=2)+'\n')
 assert len(results)==9 and all(r['status']=='validated' for r in results)

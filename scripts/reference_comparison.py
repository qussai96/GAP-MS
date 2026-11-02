"""Concordance with historical protein-inferred Casanovo BED, not truth precision."""
import csv,json
from pathlib import Path


def blocks(fields):
 start=int(fields[1]);sizes=[int(n) for n in fields[10].rstrip(',').split(',')];offsets=[int(n) for n in fields[11].rstrip(',').split(',')]
 result=[]
 for size,offset in zip(sizes,offsets):
  left,right=start+offset,start+offset+size
  if result and left==result[-1][1]:result[-1]=(result[-1][0],right)
  else:result.append((left,right))
 return tuple(result)


def compare(record):
 out=Path(record['output']);crop=out.parents[2];reference=crop/'outputs/gapms/casanovo/helixer/prediction_search/peptides_mapped.bed'
 if not reference.exists():return {'available':False}
 params=json.loads((out/'run_manifest.json').read_text())['parameters']
 queries={line.strip().replace('I','L') for line in (out/'mapping/queries.txt').read_text().splitlines()}
 all_reference=set();selected=set();excluded=0;not_queried=0
 with reference.open() as stream:
  for line in stream:
   if not line.strip() or line.startswith(('track','browser','#')):continue
   fields=line.rstrip('\n').split('\t');assert len(fields)>=12
   peptide=fields[3][::-1] if fields[5]=='-' else fields[3]
   peptide=peptide.replace('I','L')
   key=(peptide,fields[0],fields[5],blocks(fields));all_reference.add(key)
   if len(peptide)<params['min_length']:excluded+=1;continue
   if peptide not in queries:not_queried+=1;continue
   selected.add(key)
 found=set();new_count=0
 with (out/'mapping/mapped_peptides.tsv').open() as stream:
  for row in csv.DictReader(stream,delimiter='\t'):
   new_count+=1
   key=(row['peptide'].replace('I','L'),row['chrom'],row['strand'],tuple(tuple(b) for b in json.loads(row['blocks'])))
   if key in selected:found.add(key)
 missing=selected-found
 result=dict(available=True,reference=str(reference),reference_distinct_loci=len(all_reference),eligible_reference_loci=len(selected),
             shared_loci=len(found),missed_reference_loci=len(missing),eligible_reference_coverage_percent=round(100*len(found)/len(selected),4) if selected else None,
             excluded_below_length_rows=excluded,reference_rows_not_queried=not_queried,missing_examples=sorted(missing)[:10],
             interpretation='Coverage/concordance with historical protein-inferred coordinates, not a truth-set precision estimate. Negative-strand display labels restored; I/L and adjacent blocks normalized.')
 return result

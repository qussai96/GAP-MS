"""Default workflow mapper builds one automaton across both branches."""
import json,sys
from contextlib import ExitStack
from pathlib import Path
from unittest.mock import patch

from gapms import main as workflow
from gapms.protein_mapper import ProteinMatcher


def test_one_automaton_shared_by_both_branches(tmp_path):
 genome=tmp_path/'genome.fa';genome.write_text('>chr1\nATG\n')
 peptides=tmp_path/'peptides.txt';peptides.write_text('PEPTIDE\n')
 prediction=tmp_path/'predictions.gtf';prediction.write_text('chr1\tx\tCDS\t1\t3\t.\t+\t0\tParent=p\n')
 protein=tmp_path/'prediction.fa';protein.write_text('>p\nMPEPTIDE\n')
 bam=tmp_path/'reads.bam';bam.write_text('placeholder')
 output=tmp_path/'results';bam_dir=output/'bam_search';bam_dir.mkdir(parents=True)
 bam_gtf=bam_dir/'models.gtf';bam_gtf.write_text(prediction.read_text())
 bam_proteins=bam_dir/'proteins.fa';bam_proteins.write_text('>b\nPEPTIDEPEPTIDE\n')
 argv=['gapms','-a',str(genome),'-p',str(peptides),'-g',str(prediction),'-f',str(protein),'-b',str(bam),'-o',str(output)]
 with ExitStack() as stack:
  stack.enter_context(patch.object(sys,'argv',argv))
  stack.enter_context(patch.object(workflow,'prepare_bam_search_inputs',return_value={'gtf':bam_gtf,'protein_fasta':bam_proteins}))
  for name in ['extract_features','filter_predictions','map_peptides_to_genome','compare_bam_support_to_input_gtf','write_merged_supported_models','plot_parent_run_summary']:
   stack.enter_context(patch.object(workflow,name))
  factory=stack.enter_context(patch.object(ProteinMatcher,'from_file',wraps=ProteinMatcher.from_file))
  workflow.main()
  assert factory.call_count==1
 first=json.loads((output/'prediction_search/protein_mapping.summary.json').read_text())
 second=json.loads((output/'bam_search/protein_mapping.summary.json').read_text())
 assert first['automaton_searches']==1 and second['automaton_searches']==2
 assert second['mapping_rows']==2

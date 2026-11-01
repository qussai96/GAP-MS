"""Reusable orchestration for post-exact one-edit rescue and benchmarks."""
from collections import Counter, deque
from concurrent.futures import ProcessPoolExecutor
import logging
import multiprocessing
from pathlib import Path
import time

import pysam

from .junctions import splice_paths
from .rescue import init_rescue_worker, make_rescue_automaton, rescue_chunk, rescue_path
from .rescue_output import RescueStore
from .search import chunk_tasks

LOG=logging.getLogger("peplocus")


def _bounded(executor,tasks,workers):
    tasks=iter(tasks);pending=deque()
    for _ in range(workers*2):
        task=next(tasks,None)
        if task is not None:pending.append(executor.submit(rescue_chunk,task))
    while pending:
        yield pending.popleft().result()
        task=next(tasks,None)
        if task is not None:pending.append(executor.submit(rescue_chunk,task))


def run_rescue(genome,index,lengths,peptides,unmapped_ids,directory,junctions=(),min_length=12,
               threads=1,chunk_size=2000000,table=1,il_equivalent=True,max_junctions=3,
               max_paths=1000,verify=True,igv_display=False):
    if min_length < 2:
        raise ValueError("one-edit rescue minimum length must be at least 2")
    started=time.monotonic();directory=Path(directory)
    unmapped_ids=tuple(unmapped_ids)
    eligible=tuple(pid for pid in unmapped_ids if len(peptides[pid])>=min_length)
    stats=Counter();timings={}
    store=RescueStore(directory/"rescue.sqlite",peptides,unmapped_ids,eligible,il_equivalent)
    try:
        if eligible:
            LOG.info("One-edit rescue: %s eligible exact-unmapped peptides (minimum length %s)",len(eligible),min_length)
            automaton=make_rescue_automaton(peptides,eligible,il_equivalent)
            max_length=max(len(peptides[pid]) for pid in eligible)
            start=time.monotonic();report=start;done=0
            tasks=chunk_tasks(lengths,chunk_size,max_length+1)
            args=(genome,index,automaton,peptides,table,il_equivalent)
            if threads==1:
                init_rescue_worker(*args)
                try:
                    for hits,counters in map(rescue_chunk,tasks):
                        store.add(hits);stats.update(counters);done+=1
                        if time.monotonic()-report>20:
                            LOG.info("One-edit rescue: %s chunks, %s candidate alignments",done,stats["candidate_alignments"]);report=time.monotonic()
                finally:
                    from . import rescue
                    rescue._WORKER[0].close();rescue._WORKER=None
            else:
                context=multiprocessing.get_context("fork" if "fork" in multiprocessing.get_all_start_methods() else "spawn")
                with ProcessPoolExecutor(threads,initializer=init_rescue_worker,initargs=args,mp_context=context) as executor:
                    for hits,counters in _bounded(executor,tasks,threads):
                        store.add(hits);stats.update(counters);done+=1
                        if time.monotonic()-report>20:
                            LOG.info("One-edit rescue: %s chunks, %s candidate alignments",done,stats["candidate_alignments"]);report=time.monotonic()
            timings["contiguous_seconds"]=round(time.monotonic()-start,3)
            with pysam.FastaFile(genome,filepath_index=index) as fasta:
                if junctions:
                    start=time.monotonic();path_stats=Counter();buffer=[];report=start
                    for path in splice_paths(junctions,lengths,3*(max_length+1),max_junctions,max_paths,path_stats):
                        buffer.extend(rescue_path(fasta,path,automaton,peptides,table,il_equivalent,stats))
                        if len(buffer)>=5000:store.add(buffer);buffer.clear()
                        if time.monotonic()-report>20:
                            LOG.info("One-edit splice rescue: %s paths",path_stats["splice_paths_searched"]);report=time.monotonic()
                    store.add(buffer);stats.update(path_stats)
                    timings["supported_splice_seconds"]=round(time.monotonic()-start,3)
                start=time.monotonic()
                results=store.finalize(directory,fasta,table,verify,igv_display)
                timings["output_verification_seconds"]=round(time.monotonic()-start,3)
        else:
            with pysam.FastaFile(genome,filepath_index=index) as fasta:
                results=store.finalize(directory,fasta,table,verify,igv_display)
    finally:
        store.close()
    result={"results":results,"counters":dict(stats),"timings":timings,
            "elapsed_seconds":round(time.monotonic()-started,3),"minimum_length":min_length,
            "scope":"exactly one amino-acid edit; exact-unmapped queries only; all six frames plus bounded supplied splice paths; no novel approximate junction discovery",
            "edit_convention":"query relative to assembly translation: I consumes query only; D consumes translation only; X substitutes; I/L equivalence costs zero when enabled",
            "limits":"no seed/candidate/locus caps; supplied splice graph path limits still apply",
            "genetic_code":table,"il_equivalent":il_equivalent,"coordinate_verification":verify}
    LOG.info("One-edit rescue finished: %s/%s eligible peptides, %s loci, %.1f s",results["rescued_peptides"],len(eligible),results["mapping_loci"],result["elapsed_seconds"])
    return result

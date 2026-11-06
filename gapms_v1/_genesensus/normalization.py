"""CDS frame validation against supplied predictor proteins, never reference genes."""

from dataclasses import replace
from Bio import SeqIO
from Bio.Seq import Seq


def normalize_phases(source_models, specs, genome):
    report = {}
    for source, models in source_models.items():
        spec = specs[source]
        proteins = {}
        if spec.get("protein_fasta"):
            for record in SeqIO.parse(spec["protein_fasta"], "fasta"):
                seq = str(record.seq).upper().rstrip("*")
                if record.id in proteins and proteins[record.id] != seq:
                    raise ValueError(
                        f"Ambiguous predictor protein ID: {source}:{record.id}"
                    )
                proteins[record.id] = seq
        elif spec.get("residue_profiles"):
            try:
                import h5py
            except ImportError as exc:
                raise RuntimeError(
                    "Residue profiles require genesensus[molecular]"
                ) from exc
            with h5py.File(spec["residue_profiles"], "r") as h:
                idkey = next((k for k in ["protein_id", "protein_ids"] if k in h), None)
                if idkey and "sequence" in h:
                    proteins = {
                        i: s.upper().rstrip("*")
                        for i, s in zip(h[idkey].asstr()[:], h["sequence"].asstr()[:])
                    }
        corrected = []
        missing = 0
        for i, m in enumerate(models):
            expected = proteins.get(m.identifier)
            if expected is None:
                missing += 1
                continue
            dna = genome.coding_sequence(m)

            def translate(phase):
                s = dna[phase:]
                return str(Seq(s[: len(s) // 3 * 3]).translate()).rstrip("*")

            if translate(m.phase) == expected:
                continue
            alternatives = [p for p in range(3) if translate(p) == expected]
            if len(alternatives) != 1:
                raise ValueError(
                    f"Predictor protein/genome sequence mismatch or ambiguous frame: {source}:{m.identifier}"
                )
            new = replace(m, phase=alternatives[0])
            models[i] = new
            corrected.append(
                dict(
                    identifier=m.identifier,
                    contig=m.contig,
                    old_phase=m.phase,
                    new_phase=new.phase,
                    reason="unique_translation_match_to_supplied_predictor_protein",
                )
            )
        report[source] = dict(
            frame_repairs=corrected,
            proteins_missing_for_frame_validation=missing,
            supplied_protein_sequences=len(proteins),
        )
    return report

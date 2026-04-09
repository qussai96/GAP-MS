from pathlib import Path
import textwrap

import pandas as pd


def write_text(path: Path, content: str) -> Path:
    path = Path(path)
    path.write_text(textwrap.dedent(content).strip() + "\n")
    return path


def create_mini_gapms_dataset(base_dir: Path):
    """Create a tiny synthetic GAP-MS dataset for unit/integration tests."""
    base_dir = Path(base_dir)
    base_dir.mkdir(parents=True, exist_ok=True)

    output_dir = base_dir / "output"
    for subdir in [
        output_dir,
        output_dir / "Figures",
        output_dir / "Txt",
        output_dir / "Compare_to_Reference",
        output_dir / "Compare_to_Reference" / "Novel",
    ]:
        subdir.mkdir(parents=True, exist_ok=True)

    prediction_gtf = base_dir / "predictions.gtf"
    protein_fasta = base_dir / "predictions.faa"
    peptides_txt = base_dir / "peptides.txt"
    mapping_file = base_dir / "mapping.tsv"
    scores_csv = base_dir / "external_scores.csv"
    reference_gtf = base_dir / "reference.gtf"
    tmap_file = base_dir / "report.tmap"
    peptides_bed = base_dir / "peptides_mapped.bed"
    all_scores_tsv = base_dir / "all_proteins_scores.tsv"

    write_text(
        prediction_gtf,
        """
        chr1\tsrc\tmRNA\t1\t33\t0.95\t+\t.\tID=protA;Parent=geneA;gene=geneA
        chr1\tsrc\tCDS\t1\t18\t0.95\t+\t0\tID=cdsA1;Parent=protA;protein_id=protA;gene=geneA
        chr1\tsrc\tCDS\t19\t33\t0.95\t+\t0\tID=cdsA2;Parent=protA;protein_id=protA;gene=geneA
        chr1\tsrc\tmRNA\t101\t130\t0.40\t+\t.\tID=protB;Parent=geneB;gene=geneB
        chr1\tsrc\tCDS\t101\t130\t0.40\t+\t0\tID=cdsB1;Parent=protB;protein_id=protB;gene=geneB
        """,
    )

    write_text(
        reference_gtf,
        """
        chr1\tsrc\tmRNA\t4\t33\t0.99\t+\t.\tID=refA;Parent=geneRefA;gene=geneRefA
        chr1\tsrc\tCDS\t4\t18\t0.99\t+\t0\tID=refA1;Parent=refA;protein_id=refA;gene=geneRefA
        chr1\tsrc\tCDS\t19\t33\t0.99\t+\t0\tID=refA2;Parent=refA;protein_id=refA;gene=geneRefA
        chr1\tsrc\tmRNA\t101\t130\t0.80\t+\t.\tID=refB;Parent=geneRefB;gene=geneRefB
        chr1\tsrc\tCDS\t101\t130\t0.80\t+\t0\tID=refB1;Parent=refB;protein_id=refB;gene=geneRefB
        """,
    )

    write_text(
        protein_fasta,
        """
        >protA
        MACDEFGHAAA
        >protB
        MNPQRSTVWYA
        """,
    )

    write_text(
        peptides_txt,
        """
        MACD
        FGHA
        RSTV
        """,
    )

    write_text(
        mapping_file,
        """
        Peptide\tProtein\tlocation
        MACD\tprotA\t1
        FGHA\tprotA\t2
        RSTV\tprotB\t3
        """,
    )

    write_text(
        scores_csv,
        """
        description,external_score
        protA,0.95
        protB,0.15
        """,
    )

    write_text(
        tmap_file,
        """
        qry_id\tref_id\tclass_code
        protA\trefA\tj
        protB\trefB\t=
        """,
    )

    write_text(
        peptides_bed,
        """
        chr1\t0\t12\tMACD\t0\t+
        """,
    )

    pd.DataFrame(
        [
            {
                "Protein": "protA",
                "mapped_peptides": 2,
                "sequence_coverage": 0.7273,
                "N_terminal_peptides": 1,
                "C_terminal_peptides": 0,
                "splice_peptides": 0,
            },
            {
                "Protein": "protB",
                "mapped_peptides": 1,
                "sequence_coverage": 0.3636,
                "N_terminal_peptides": 0,
                "C_terminal_peptides": 0,
                "splice_peptides": 0,
            },
        ]
    ).to_csv(all_scores_tsv, sep="\t", index=False)

    return {
        "output_dir": output_dir,
        "prediction_gtf": prediction_gtf,
        "protein_fasta": protein_fasta,
        "peptides_txt": peptides_txt,
        "mapping_file": mapping_file,
        "scores_csv": scores_csv,
        "reference_gtf": reference_gtf,
        "tmap_file": tmap_file,
        "peptides_bed": peptides_bed,
        "all_scores_tsv": all_scores_tsv,
    }

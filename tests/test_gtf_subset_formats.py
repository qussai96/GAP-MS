"""
Manual test: exercise save_gtf_subset_all_features against three real GTF formats.
Run with: python tests/test_gtf_subset_formats.py
"""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from gapms.gtf_utils import save_gtf_subset_all_features

BASE = Path("/home/students/q.abbas/GAPMS/0-raw_files/GCF_001625215.2_Daucus_carota/inputs")

GTF_FILES = {
    "braker":    BASE / "braker/GCF_001625215.2_Daucus_carota.braker.gtf",
    "annevo":    BASE / "annevo/GCF_001625215.2_Daucus_carota.annevo.gtf",
    "reference": BASE / "reference/GCF_001625215.2_Daucus_carota.ref.gtf",
}

# One test case per format: (feature_type_to_scan, attr_key_regex)
# We collect 5 real IDs from each file by scanning the relevant mRNA/transcript rows.
import re

def collect_ids(gtf_path, n=5):
    """Return up to n real transcript/mRNA protein identifiers from a GTF file."""
    ids = []
    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            ftype = parts[2].lower()
            attrs = parts[8]
            if ftype not in ('transcript', 'mrna'):
                continue
            # Try GFF3 ID=
            m = re.search(r'ID=([^;,\s]+)', attrs)
            if m:
                ids.append(m.group(1))
            else:
                # Try GTF transcript_id "value"
                m = re.search(r'transcript_id\s*"([^"]+)"', attrs)
                if m:
                    ids.append(m.group(1))
                elif re.match(r'^[\w.\-]+$', attrs.strip()):
                    # Bare ID (e.g. braker transcript rows)
                    ids.append(attrs.strip())
            if len(ids) >= n:
                break
    return set(ids)


def run_test(label, gtf_path, supported):
    print(f"\n{'='*60}")
    print(f"Format : {label}")
    print(f"GTF    : {gtf_path.name}")
    print(f"Subset : {sorted(supported)}")
    with tempfile.TemporaryDirectory() as tmpdir:
        out_dir = Path(tmpdir)
        save_gtf_subset_all_features(gtf_path, supported, out_dir, "out.gtf")
        out = out_dir / "out.gtf"
        lines = out.read_text().splitlines()
        print(f"Output lines : {len(lines)}")
        if lines:
            # Show feature-type breakdown
            from collections import Counter
            types = Counter(l.split('\t')[2] for l in lines if '\t' in l)
            print(f"Feature types: {dict(types)}")
            print(f"Sample line  : {lines[0][:120]}")
        else:
            print("OUTPUT IS EMPTY — no lines matched the subset!")


if __name__ == "__main__":
    ok = True
    for label, gtf_path in GTF_FILES.items():
        supported = collect_ids(gtf_path, n=5)
        run_test(label, gtf_path, supported)
        with tempfile.TemporaryDirectory() as tmpdir:
            out_dir = Path(tmpdir)
            save_gtf_subset_all_features(gtf_path, supported, out_dir, "out.gtf")
            n_lines = len((out_dir / "out.gtf").read_text().splitlines())
        if n_lines == 0:
            print(f"  FAIL: {label} produced 0 output lines")
            ok = False
        else:
            print(f"  PASS: {label} ({n_lines} lines)")
    print()
    sys.exit(0 if ok else 1)

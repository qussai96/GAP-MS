"""Build the pinned GAPMS 0.3.0 dependency without modifying either namespace."""
import argparse
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=Path('dist'))
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='gapms-legacy-build-') as directory:
        stage=Path(directory)
        shutil.copytree(root/'gapms',stage/'gapms',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        for name in ['setup.py','pyproject.toml','MANIFEST.in']:
            shutil.copy2(root/'packaging/legacy'/name,stage/name)
        for name in ['README.md','LICENSE','THIRD_PARTY.md','xgb_protein_classifier.pkl']:
            shutil.copy2(root/name,stage/name)
        subprocess.run([sys.executable,'-m','pip','wheel','--no-deps','--no-build-isolation','--wheel-dir',str(out),str(stage)],check=True)


if __name__=='__main__':main()

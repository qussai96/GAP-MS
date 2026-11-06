"""Lossless optional streaming compression for large audit tables."""
from contextlib import contextmanager
from pathlib import Path
import gzip,io

@contextmanager
def diagnostic_writer(path,compressed=False):
    path=Path(path)
    if not compressed:
        with path.open('w',encoding='utf-8',newline='') as handle:yield handle
    else:
        with path.open('wb') as raw:
            with gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0,compresslevel=1) as zipped:
                with io.TextIOWrapper(zipped,encoding='utf-8',newline='') as handle:yield handle

def diagnostic_reader(folder,name):
    path=Path(folder)/name
    if path.exists():return path.open('r',encoding='utf-8',newline='')
    return gzip.open(str(path)+'.gz','rt',encoding='utf-8',newline='')

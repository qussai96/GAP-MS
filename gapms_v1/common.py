"""Deterministic identifiers, auditable tables and safe subprocess records."""
import csv
import hashlib
import json
import logging
from pathlib import Path
import subprocess
import time

LOG = logging.getLogger('gapms_v1')


def stable_id(prefix, *values):
    encoded = json.dumps(values, separators=(',', ':'), sort_keys=True).encode()
    return prefix + hashlib.sha256(encoded).hexdigest()[:20]


def write_table(path, fields, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, separators=(',', ':')) if isinstance(v, (list, dict, tuple)) else v for k, v in row.items() if k in fields})


def read_table(path):
    with Path(path).open() as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def write_json(path, value):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + '.partial')
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    tmp.replace(path)


def run_command(command, cwd, log, commands):
    command = [str(v) for v in command]
    entry = dict(command=command, cwd=str(cwd), started=time.time())
    commands.append(entry)
    LOG.info('Running %s', ' '.join(command))
    with Path(log).open('w') as stream:
        result = subprocess.run(command, cwd=cwd, stdout=stream, stderr=subprocess.STDOUT)
    entry.update(returncode=result.returncode, seconds=time.time() - entry['started'])
    if result.returncode:
        raise RuntimeError(f'{command[0]} failed ({result.returncode}); see {log}')


def input_identity(path, checksum_limit=32 * 1024**2):
    path = Path(path).resolve()
    st = path.stat()
    identity = dict(path=str(path), bytes=st.st_size, mtime_ns=st.st_mtime_ns)
    if st.st_size <= checksum_limit:
        identity['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    else:
        with path.open('rb') as stream:
            first = stream.read(1024**2)
            stream.seek(max(0, st.st_size - 1024**2))
            last = stream.read(1024**2)
        identity['sampled_sha256'] = hashlib.sha256(first + last).hexdigest()
        identity['checksum_scope'] = 'first_and_last_1MiB_plus_size_mtime; not full checksum'
    return identity

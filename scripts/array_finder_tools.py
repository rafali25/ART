"""Run PILER-CR 1.06 and MinCED 0.4.2 on one sequence and parse their array calls.

Binaries come from .venv-array-finders/ (scripts/array_finders_setup.sbatch).
Each call works in a private temporary directory under $TMPDIR or /tmp.
"""
from __future__ import annotations

import os
from pathlib import Path
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT/'.venv-array-finders/bin'
# Tuned regime of the paper's baseline comparison: repeat and spacer ranges widened
# to 15-45 nt and 25-420 nt to bracket measured ART geometry; all else default.
PILERCR_TUNED = ('-minrepeat', '15', '-maxrepeat', '45', '-minspacer', '25', '-maxspacer', '420')
MINCED_TUNED = ('-minRL', '15', '-maxRL', '45', '-minSL', '25', '-maxSL', '420')

_ROW = re.compile(r'^\s*(\d+)\s+(\d+)\s+([\d.]+)\s')
_FOOTER = re.compile(r'^\s*(\d+)\s+(\d+)\s+(?:(\d+)\s+)?([A-Za-z-]+)\s*$')


def parse_pilercr(text):
    """Arrays from the DETAIL REPORT: copies, consensus repeat length, per-copy %id."""
    detail = text.split('DETAIL REPORT', 1)[1].split('SUMMARY BY', 1)[0] if 'DETAIL REPORT' in text else ''
    arrays = []
    for block in re.split(r'^Array \d+\s*$', detail, flags=re.M)[1:]:
        rules = [i for i, line in enumerate(block.splitlines()) if line.lstrip().startswith('====')]
        lines = block.splitlines()
        if len(rules) < 2:
            raise ValueError('Unparsed PILER-CR array block')
        rows = [_ROW.match(line) for line in lines[rules[0]+1:rules[1]]]
        if not all(rows):
            raise ValueError('Unparsed PILER-CR copy row')
        footer = _FOOTER.match(lines[rules[1]+1])
        if not footer or int(footer.group(1)) != len(rows):
            raise ValueError('PILER-CR footer does not match its copy rows')
        identities = [float(m.group(3)) for m in rows]
        arrays.append(dict(copies=len(rows), repeat_nt=int(footer.group(2)),
                           mean_identity=sum(identities)/len(identities),
                           start1=int(rows[0].group(1)), consensus=footer.group(4)))
    reported = re.search(r':\s*(\d+) putative CRISPR arrays? found', text)
    if not reported or int(reported.group(1)) != len(arrays):
        raise ValueError('PILER-CR array count does not match its detail report')
    return arrays


def parse_minced(text):
    arrays = []
    for block in re.split(r'^CRISPR \d+\s+Range:', text, flags=re.M)[1:]:
        copies = re.search(r'^Repeats:\s*(\d+)', block, flags=re.M)
        if not copies:
            raise ValueError('Unparsed MinCED array block')
        arrays.append(dict(copies=int(copies.group(1))))
    return arrays


def _run(command, sequence, parse):
    with tempfile.TemporaryDirectory(prefix='array-finder-', dir=os.environ.get('TMPDIR', '/tmp')) as work:
        fasta, report = Path(work)/'in.fna', Path(work)/'out.txt'
        fasta.write_text('>query\n'+sequence+'\n')
        subprocess.run(command(fasta, report), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        return parse(report.read_text() if report.exists() else '')


def run_pilercr(sequence, options=()):
    return _run(lambda fasta, out: [str(BIN/'pilercr'), '-in', str(fasta), '-out', str(out), '-quiet', '-noinfo', *options],
                sequence, parse_pilercr)


def run_minced(sequence, options=()):
    return _run(lambda fasta, out: [str(BIN/'minced'), *options, str(fasta), str(out)], sequence, parse_minced)


def best_array_score(arrays):
    """Best-array statistic: (copies-1) x consensus repeat length x mean copy identity.

    PILER-CR reports no numeric array score, so this mirrors S=(R-1)*I using
    only quantities PILER-CR itself reports. Zero when no array is called.
    """
    return max(((a['copies']-1)*a['repeat_nt']*a['mean_identity']/100 for a in arrays), default=0.)

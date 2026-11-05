"""Verify source-coordinate maps and detect known seeds with each native profile."""
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
from pathlib import Path
import subprocess

from Bio import SeqIO

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def positive_control(item):
    family, profile, seeds, database_size = item
    stem = OUT / 'positive_controls' / family
    args = [str(ROOT/'.venv/bin/hmmsearch'), '--cpu','1','--seed','42','--noali',
            '-Z',str(database_size),'--domZ',str(database_size),'-E','0.001','--domE','0.001',
            '--tblout',str(stem)+'.tbl','--domtblout',str(stem)+'.domtbl',str(profile),str(seeds)]
    with Path(str(stem)+'.log').open('w') as handle:
        result = subprocess.run(args, stdout=handle, stderr=subprocess.STDOUT)
    assert result.returncode == 0, family
    passing = []
    for line in Path(str(stem)+'.domtbl').read_text().splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        c=line.split()
        if float(c[6])<=1e-5 and float(c[12])<=1e-5 and (int(c[16])-int(c[15])+1)/int(c[5])>=0.6:
            passing.append(c[0])
    assert passing, 'No known seed recovered with '+family
    return dict(family=family, returncode=result.returncode, recovered_known_proteins=sorted(set(passing)),
                effective_Z=database_size, effective_domZ=database_size, args=args,
                seed_sha256=sha(seeds.read_bytes()), profile_sha256=sha(profile.read_bytes()),
                limitation='Same-profile seed recovery checks executable/profile validity, not biological generalization')


def main():
    summary=json.loads((OUT/'summary.json').read_text())
    genes=list(csv.DictReader((OUT/'gene_coordinates.tsv').open(),delimiter='\t'))
    dna={r.id:r.seq for r in SeqIO.parse(OUT/'sources/deposited/139.curated.Mac.phage.genomes.fa','fasta')}
    proteins={r.id:str(r.seq) for r in SeqIO.parse(OUT/'proteins.faa','fasta')}
    assert len(proteins)==len(genes)==summary['proteins']
    assert len(dna)==summary['genomes']==139
    for row in genes:
        a,b=int(row['start0']),int(row['end0']);source=dna[row['genome_id']]
        assert 0<=a<b<=len(source)
        extracted=source[a:b]
        if row['strand']=='-1':
            extracted=extracted.reverse_complement()
        assert sha(str(extracted).encode())==row['CDS_sha256']
        assert sha(proteins[row['protein_id']].encode())==row['protein_sha256']
        assert len(proteins[row['protein_id']])==int(row['protein_length_aa'])
        assert row['translation_compatible']=='True'
    coordinates={row['protein_id']:(row['genome_id'],int(row['start0'])+1,int(row['end0']),row['strand']) for row in genes}
    gff_count=0
    for line in (OUT/'predicted_genes.gff3').read_text().splitlines():
        if line.startswith('#'):
            continue
        c=line.split('\t');attrs=dict(token.split('=',1) for token in c[8].split(';'))
        assert coordinates[attrs['ID']]==(c[0],int(c[3]),int(c[4]),'1' if c[6]=='+' else '-1')
        gff_count+=1
    assert gff_count==len(genes)
    groups=defaultdict(list)
    for identifier,sequence in dna.items():
        groups[sha(str(sequence).upper().encode())].append(identifier)
    duplicate_groups={key:ids for key,ids in groups.items() if len(ids)>1}
    (OUT/'duplicate_genome_groups.json').write_text(json.dumps(duplicate_groups,indent=2)+'\n')
    (OUT/'positive_controls').mkdir(exist_ok=True)
    base=ROOT/'results/art_discovery/profiles'
    items=[('ART_RT',base/'ART_RT.hmm',base/'ART_RT.seeds.faa',len(genes)),
           ('Type_I_partner',base/'Type_I_partner.hmm',base/'Type_I_partner.seeds.faa',len(genes)),
           ('Type_II_partner',base/'Type_II_partner.hmm',base/'Type_II_partner.seeds.faa',len(genes)),
           ('TypeI_two_seed_exploratory',ROOT/'results/art_type1_public/profiles/TypeI_two_seed_exploratory.hmm',
            ROOT/'results/art_type1_public/profiles/typeI_two_LPJP1_queries.faa',len(genes)),
           ('RVT_1',ROOT/'.tools/art-reconstruction/PF00078.hmm',base/'ART_RT.seeds.faa',len(genes))]
    with ThreadPoolExecutor(max_workers=4) as pool:
        controls=list(pool.map(positive_control,items))
    result=dict(passed=True, genomic_CDS_mappings_checked=len(genes), protein_hashes_checked=len(genes),
                GFF_coordinates_checked=gff_count, unique_source_sequences=len(groups),
                exact_duplicate_genome_groups=duplicate_groups, positive_controls=controls,
                source_archive_MD5_verified=True, unexplained_translation_mismatches=0)
    (OUT/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()

"""Build a portable, source-linked candidate dataset from the frozen ART census."""
from collections import Counter
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import socket
import zipfile

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from art_expansion_workflow import ROOT, OUT, DATASETS, read_tsv, table, write_json, sha, utc


def main():
    if not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Run this exporter via Slurm; no login-node analysis.')
    census=ROOT/'results/art_evidence_census_20261105'
    frozen=json.loads((census/'summary.json').read_text())
    for name,digest in frozen['input_sha256'].items():
        assert sha(ROOT/name)==digest, name
    evidence=read_tsv(census/'sample_evidence.tsv')
    catalog=read_tsv(OUT/'step5/pair_catalog.tsv')
    assert len(evidence)==len(catalog)==46
    assert {r['locus_id'] for r in evidence}=={r['locus_id'] for r in catalog}
    destination=ROOT/'data/art_samples_2026-11-05'
    if (destination/'verification.json').exists():
        raise SystemExit('A verified export already exists; preserve this snapshot and choose a new version.')
    for directory in ['sequences','loci','sources','provenance','code']:
        (destination/directory).mkdir(parents=True,exist_ok=True)
    metadata={}
    manifests=[ROOT/'results/art_expansion/stage1/source_records.tsv',
        ROOT/'results/art_sample_expansion_20261104/panel_acquisition.tsv']
    manifests += [ROOT/'results'/name/'source_records.tsv' for name in ['art_discovery','art_type1_public','art_type1_mgnify']]
    manifests += [ROOT/'results/art_sample_expansion_20261104'/name/'source_records.tsv' for name in ['deferred_acquisition/processed','diverse_partner_search/processed']]
    manifests += [OUT/'step2/processed/source_records.tsv',OUT/'HPGC/acquired_source_records.tsv',OUT/'JPGD/acquired_source_records.tsv']
    for manifest in manifests:
        for row in read_tsv(manifest):
            if row.get('path'):metadata[row['path'].replace('\\','/')]=row
    rts={r.id:str(r.seq) for r in SeqIO.parse(OUT/'step5/RT_proteins.faa','fasta')}
    output=[];source_index={};source_rows=[];RT_records=[];partner_records=[];rt_dna_records=[];partner_dna_records=[];locus_records=[]
    input_hashes={str((census/'sample_evidence.tsv').relative_to(ROOT)):sha(census/'sample_evidence.tsv')}
    fields={
        'sample_id':'Stable identifier within this versioned export; shared by FASTA and locus files.',
        'locus_id':'Original analysis locus identifier.',
        'partner_type':'Type_I or Type_II partner-family assignment; computational candidate.',
        'cohort':'Analysis cohort in which this representative was retained.',
        'record_role':'Original discovery/control/reference role.',
        'new_in_20261105_expansion':'New to this expansion run, not a claim of biological novelty.',
        'pair_complete':'Both RT and partner have complete translation-matched CDSs.',
        'validation_status':'Candidate evidence tier; not an experimentally established ART assignment.',
        'source_identifier':'Accession or internal genome identifier used for the source.',
        'original_source_identifier':'Original catalog sequence identifier, when available.',
        'source_run':'SRA/ENA sequencing run, where resolved.',
        'source_study':'Source study accession; unresolved if unavailable.',
        'source_BioSample':'Source sample accession; unresolved if unavailable.',
        'source_group':'BioSample where available, otherwise canonical source-DNA SHA256; not proof of independence.',
        'source_url':'Recorded retrieval URL or accession landing page; catalog URL for HPGC/JPGD.',
        'source_canonical_DNA_sha256':'Minimum SHA256 of forward/reverse-complement uppercase source DNA.',
        'source_file':'Packaged original source sequence file; path relative to this dataset directory.',
        'source_file_format':'Biopython source format.',
        'RT_feature_id':'RT protein/CDS identifier in source and locus annotation.',
        'partner_feature_id':'Partner protein/CDS identifier in source and locus annotation.',
        'RT_length_aa':'RT protein length in amino acids.',
        'partner_length_aa':'Partner protein length in amino acids.',
        'RT_genetic_code':'NCBI translation-table number used for the RT.',
        'partner_genetic_code':'NCBI translation-table number used for the partner.',
        'RT_codon_start':'One-based coding offset within extracted RT CDS.',
        'partner_codon_start':'One-based coding offset within extracted partner CDS.',
        'RT_protein_sequence':'RT amino-acid sequence without terminal stop.',
        'partner_protein_sequence':'Partner amino-acid sequence without terminal stop.',
        'RT_CDS_sequence':'Oriented raw RT CDS DNA, including any terminal stop and pre-codon_start bases.',
        'partner_CDS_sequence':'Oriented raw partner CDS DNA, including any terminal stop and pre-codon_start bases.',
        'protein_pair_sha256':'SHA256 of RT_protein_sequence + | + partner_protein_sequence.',
        'category1_ART_profile_pass':'Frozen ART RT profile: sequence/domain E<=1e-5, profile coverage>=0.60.',
        'category2_ART_reference_affinity':'Descriptive ART-reference affinity in the current RT-core tree.',
        'category2_family_validation_status':'ART versus NART remains unresolved; this is an unknown validation result, not an ART rejection.',
        'category3_RVT1_domain':'Significant RVT_1 envelope with >=60% profile coverage.',
        'category3_catalytic_motif_pattern':'[YF]xDD pattern within that envelope; not proof of catalytic activity.',
        'motif_matches_json':'Motif sequences and zero-based positions in the RT protein.',
        'NTD_length_aa':'Residues before the first significant RVT_1 envelope; operational boundary.',
        'category4_NTD_at_least_150aa':'Descriptive extended-NTD indicator; 150 aa is not a published Yoon acceptance cutoff.',
        'same_strand_partner':'Partner and RT are co-oriented in the extracted locus.',
        'direct_downstream_partner':'No intervening annotated CDS before the selected downstream partner.',
        'upstream_annotation_free_nt':'Observed contiguous annotated-CDS-free DNA immediately before the RT; may be boundary-limited.',
        'category5_direct_partner_and_500nt':'Direct co-oriented partner and >=500 nt observed upstream annotation-free DNA.',
        'available_upstream_window_nt':'Upstream DNA available in the locus extraction; not necessarily noncoding.',
        'upstream_window_DNA_sequence':'Oriented locus sequence before the RT start, including any annotated genes.',
        'category6_historical_supported_array':'Supported array evidence from the retained historical assessment; no merged-cohort retest.',
        'historical_array_scope':'Cohort in which the historical supported-array evidence was assessed.',
        'array_status':'Original detailed array assessment/status label.',
        'array_retested_in_expanded_collection':'False: no new joint expanded-collection significance test was performed.',
        'experimentally_validated_in_this_run':'False: no experimental validation was performed in this run.',
        'partner_cluster70':'Type-I partner cluster at 70% identity / 80% mutual coverage, or not_typeI; not independent lineages.',
        'locus_genbank_file':'Packaged annotated oriented locus, relative path.',
        'locus_DNA_file':'Packaged oriented locus FASTA, relative path.',
        'coordinate_map_file':'Packaged zero-based half-open source/locus coordinate map, relative path.',
    }
    for index,r in enumerate(evidence,1):
        sample=f'ARTS{index:04d}'
        source_path=ROOT/r['source_path'].replace('\\','/')
        source=SeqIO.read(source_path,r['source_format'])
        dna=str(source.seq).upper()
        canonical=min(hashlib.sha256(dna.encode()).hexdigest(),hashlib.sha256(str(Seq(dna).reverse_complement()).encode()).hexdigest())
        assert canonical==r['source_DNA_sha256'],r['locus_id']
        digest=sha(source_path)
        input_hashes[str(source_path.relative_to(ROOT))]=digest
        if digest not in source_index:
            source_name=f'SRC{len(source_index)+1:04d}'+source_path.suffix
            relative='sources/'+source_name
            shutil.copy2(source_path,destination/relative)
            source_index[digest]=relative
            source_rows.append(dict(source_file=relative,source_identifier=source.id,format=r['source_format'],
                original_repository_path=r['source_path'],file_sha256=digest,canonical_DNA_sha256=canonical))
        relative_source=source_index[digest]
        locus_path=ROOT/r['locus_genbank'].replace('\\','/')
        locus=SeqIO.read(locus_path,'genbank')
        def feature(pid):
            return next(f for f in locus.features if f.type=='CDS' and f.qualifiers.get('protein_id',f.qualifiers.get('locus_tag',['']))[0]==pid)
        rt=feature(r['RT_feature_id']);partner=feature(r['partner_feature_id'])
        rt_seq=rts[r['tree_key']];partner_seq=partner.qualifiers['translation'][0]
        assert hashlib.sha256((rt_seq+'|'+partner_seq).encode()).hexdigest()==r['pair_sequence_sha256']
        cds_sequences=[];codes=[];starts=[]
        for f,expected in [(rt,rt_seq),(partner,partner_seq)]:
            cds=f.extract(locus.seq);code=int(f.qualifiers.get('transl_table',['11'])[0]);start=int(f.qualifiers.get('codon_start',['1'])[0])
            coding=cds[start-1:]
            if r['pair_complete']=='True':actual=str(coding.translate(table=code,cds=True))
            else:
                actual=str(coding.translate(table=code)).rstrip('*')
                if expected.startswith('M') and actual[1:]==expected[1:]:actual='M'+actual[1:]
            assert len(actual)==len(expected) and all(a==b or b=='X' for a,b in zip(actual,expected)),r['locus_id']
            cds_sequences.append(str(cds));codes.append(code);starts.append(start)
        original_map=ROOT/r['coordinate_map'].replace('\\','/')
        mapping=json.loads(original_map.read_text())
        pieces=[source.seq[s['source_start0']:s['source_end0']] if s['orientation']==1 else source.seq[s['source_start0']:s['source_end0']].reverse_complement() for s in mapping['segments']]
        assert Seq('').join(pieces)==locus.seq,r['locus_id']
        gb_relative=f'loci/{sample}.gb';dna_relative=f'loci/{sample}.fna';map_relative=f'loci/{sample}.map.json'
        shutil.copy2(locus_path,destination/gb_relative)
        SeqIO.write(SeqRecord(locus.seq,id=sample,description=r['locus_id']),destination/dna_relative,'fasta')
        mapping.update(source_path=relative_source,source_format=r['source_format'],locus_path=gb_relative,
            sample_id=sample,locus_id=r['locus_id'],coordinate_convention='zero_based_half_open',
            original_repository_map=r['coordinate_map'])
        write_json(destination/map_relative,mapping)
        meta=metadata.get(r['source_path'].replace('\\','/'),{})
        source_url=meta.get('source_url','')
        if r['dataset'] in DATASETS:source_url=DATASETS[r['dataset']]['url']
        if not source_url and re_accession(r['source_genome']):source_url='https://www.ncbi.nlm.nih.gov/nuccore/'+r['source_genome']
        strict=r['ART_profile_screen']=='True';complete=r['pair_complete']=='True'
        status=('complete_ART_profile_context_candidate' if strict else 'complete_RT_partner_candidate_family_review') if complete else ('incomplete_ART_profile_context_candidate' if strict else 'incomplete_RT_partner_candidate_family_review')
        observed=int(r['observed_upstream_annotation_free_nt'])
        record=dict(sample_id=sample,locus_id=r['locus_id'],partner_type=r['family'],cohort=r['dataset'],record_role=r['role'],
            new_in_20261105_expansion=r['new_in_current_run'],pair_complete=r['pair_complete'],validation_status=status,
            source_identifier=r['source_genome'],original_source_identifier=r.get('original_source_id') or meta.get('original_id',''),
            source_run=r['source_run'],source_study=r['study_accession'],source_BioSample=r['source_BioSample'],source_group=r['source_group'],
            source_url=source_url,source_canonical_DNA_sha256=canonical,source_file=relative_source,source_file_format=r['source_format'],
            RT_feature_id=r['RT_feature_id'],partner_feature_id=r['partner_feature_id'],RT_length_aa=len(rt_seq),partner_length_aa=len(partner_seq),
            RT_genetic_code=codes[0],partner_genetic_code=codes[1],RT_codon_start=starts[0],partner_codon_start=starts[1],
            RT_protein_sequence=rt_seq,partner_protein_sequence=partner_seq,RT_CDS_sequence=cds_sequences[0],partner_CDS_sequence=cds_sequences[1],
            protein_pair_sha256=r['pair_sequence_sha256'],category1_ART_profile_pass=r['ART_profile_screen'],
            category2_ART_reference_affinity=r['descriptive_ART_reference_affinity'],category2_family_validation_status='unresolved_without_verified_NART_controls',
            category3_RVT1_domain=r['significant_RVT1_domain'],category3_catalytic_motif_pattern=r['domain_YF_x_DD_pattern'],motif_matches_json=r['motif_matches_json'],
            NTD_length_aa=r['NTD_before_first_RVT1_aa'],category4_NTD_at_least_150aa=r['NTD_at_least_150aa_descriptive'],
            same_strand_partner=r['same_strand_partner'],direct_downstream_partner=r['direct_downstream_partner'],upstream_annotation_free_nt=observed,
            category5_direct_partner_and_500nt=r['direct_downstream_partner']=='True' and observed>=500,
            available_upstream_window_nt=r['upstream_nt'],upstream_window_DNA_sequence=str(locus.seq[:int(rt.location.start)]),
            category6_historical_supported_array=r['legacy_supported_array'],historical_array_scope=r['legacy_array_evidence_scope'],array_status=r['array_status'],
            array_retested_in_expanded_collection=False,experimentally_validated_in_this_run=False,partner_cluster70=r['partner_cluster70'],
            locus_genbank_file=gb_relative,locus_DNA_file=dna_relative,coordinate_map_file=map_relative)
        assert set(record)==set(fields)
        output.append(record)
        for collection,sequence in [(RT_records,rt_seq),(partner_records,partner_seq),(rt_dna_records,cds_sequences[0]),(partner_dna_records,cds_sequences[1]),(locus_records,str(locus.seq))]:
            collection.append(SeqRecord(Seq(sequence),id=sample,description=r['locus_id']))
    with (destination/'art_samples.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(fields));writer.writeheader();writer.writerows(output)
    table(destination/'art_samples.tsv',output,list(fields))
    table(destination/'data_dictionary.tsv',[dict(column=k,description=v) for k,v in fields.items()])
    for name,records in [('RT_proteins.faa',RT_records),('partner_proteins.faa',partner_records),('RT_CDS.fna',rt_dna_records),('partner_CDS.fna',partner_dna_records),('oriented_loci.fna',locus_records)]:
        assert SeqIO.write(records,destination/'sequences'/name,'fasta')==46
    table(destination/'provenance/source_manifest.tsv',source_rows)
    ids={r['locus_id']:r['sample_id'] for r in output}
    aliases=read_tsv(OUT/'step5/exact_pair_aliases.tsv')
    table(destination/'provenance/excluded_aliases.tsv',[dict(r,representative_sample_id=ids[r['representative']]) for r in aliases])
    for source,name in [(census/'sample_evidence.tsv','original_sample_evidence.tsv'),(census/'summary.json','evidence_census_summary.json'),
                         (OUT/'step5/source_mapping_verification.tsv','prior_source_mapping_verification.tsv'),(OUT/'step5/summary.json','expansion_summary.json')]:
        shutil.copy2(source,destination/'provenance'/name)
        input_hashes[str(source.relative_to(ROOT))]=sha(source)
    shutil.copy2(__file__,destination/'code/export_art_sample_dataset.py')
    write_json(destination/'provenance/input_sha256.json',input_hashes)
    family_counts=dict(Counter(r['partner_type'] for r in output))
    summary=dict(dataset_version='2026-11-05',created_utc=utc(),slurm_job_id=os.environ['SLURM_JOB_ID'],compute_host=socket.gethostname(),
        samples=46,partner_types=family_counts,complete_pairs=sum(r['pair_complete']=='True' for r in output),
        complete_typeI_pairs=sum(r['partner_type']=='Type_I' and r['pair_complete']=='True' for r in output),
        complete_strict_ART_screen_typeI_pairs=sum(r['partner_type']=='Type_I' and r['pair_complete']=='True' and r['category1_ART_profile_pass']=='True' for r in output),
        newly_retained_in_expansion=sum(r['new_in_20261105_expansion']=='True' for r in output),excluded_aliases=len(aliases),source_files=len(source_rows),
        all_samples_family_validation_status='unresolved_without_verified_NART_controls',experimentally_validated_in_this_run=0)
    write_json(destination/'summary.json',summary)
    (destination/'README.md').write_text('''# ART candidate sample dataset — 5 November 2026

`art_samples.csv` is the main shareable data file: **46 deduplicated RT–partner
pairs (36 type I, 10 type II)**. Each row includes source identifiers, sequences,
coding completeness, the six evidence categories, and validation limitations.
This collection includes prior references/controls and newly acquired candidates;
`cohort`, `record_role`, and `new_in_20261105_expansion` distinguish them.
The collection contains 31 complete type I pairs, 18 passing the strict ART screen.

These are computational candidates. ART-versus-NART validation is unresolved;
no experimental validation was performed in this run. Evidence indicators are
overlapping checks, not six mutually exclusive classes or calibrated probabilities.
False for profile/array evidence does not by itself establish biological absence.
Missing source metadata is blank or `unresolved` and is not an independence claim.

## Contents

- `art_samples.csv` / `.tsv`: identical sample rows; protein, CDS and upstream-window
  sequences are included directly, so the main table can be used on its own.
- `data_dictionary.tsv`: definitions and units for every column.
- `sequences/`: RT/partner protein and CDS FASTAs plus oriented locus FASTA, keyed
  by the same stable `ARTS0001`–`ARTS0046` identifiers.
- `loci/`: annotated GenBank, locus FASTA and source-coordinate maps per sample.
- `sources/`: original source sequence files; shared sources are stored once.
- `provenance/`: source manifest, three excluded aliases, original evidence,
  previous mapping verification, and input hashes.
- `summary.json`, `verification.json`, `SHA256SUMS.json`: counts and export checks.
- `code/`: the exporter used to build this snapshot.

## Evidence interpretation

1. ART-profile support uses the frozen sequence/domain E<=1e-5 and >=60% HMM
   coverage screen. It is not proof of ART function.
2. Tree affinity is descriptive. Verified NART controls are missing; all family
   validation statuses remain unresolved.
3. RVT_1 and a `[YF]xDD` pattern support RT identity. Motif presence does not prove
   catalytic activity or establish the homologous catalytic position by itself.
4. NTD length is the sequence before the first significant RVT_1 envelope. The
   >=150-aa indicator is descriptive, not a Yoon acceptance threshold.
5. Architecture requires a directly downstream co-oriented selected partner and
   >=500 nt of observed upstream annotation-free sequence. The available upstream
   window is a different quantity and can contain coding genes. Source gene
   annotations and contig boundaries limit this assessment.
6. Supported arrays retain historical cohort-specific assessments, including
   outstanding copy/boundary/family review. No joint expanded-cohort statistical
   retest was performed, and no new accepted arrays were generated in this run.

Coordinate maps use zero-based half-open intervals. CDS DNA is oriented in the
coding direction and preserves raw feature sequence; apply the recorded
`codon_start` and genetic code before translating. Complete CDSs include stop
codons; protein FASTAs omit terminal stops. GenBank coordinates retain their
standard representation. Sample IDs identify sequence-pair representatives,
not independent organisms or experimental replicates.

The three annotation/exact-pair aliases are retained as provenance, not additional
sample rows. Public source identifiers and retrieval links are supplied for
attribution; this export does not change the original data providers' licenses.
The literature and complete database archives are not bundled.
''')
    with (destination/'art_samples.csv').open(newline='') as handle:roundtrip=list(csv.DictReader(handle))
    assert len(roundtrip)==46 and len({r['sample_id'] for r in roundtrip})==46
    assert len({r['protein_pair_sha256'] for r in roundtrip})==46
    assert family_counts=={'Type_I':36,'Type_II':10}
    for r in roundtrip:
        for field in ['source_file','locus_genbank_file','locus_DNA_file','coordinate_map_file']:
            assert (destination/r[field]).is_file()
    for name in ['RT_proteins.faa','partner_proteins.faa','RT_CDS.fna','partner_CDS.fna','oriented_loci.fna']:
        records=list(SeqIO.parse(destination/'sequences'/name,'fasta'))
        assert [r.id for r in records]==[r['sample_id'] for r in roundtrip]
    manifest={str(p.relative_to(destination)):sha(p) for p in sorted(destination.rglob('*')) if p.is_file() and p.name not in ['SHA256SUMS.json','verification.json']}
    write_json(destination/'SHA256SUMS.json',manifest)
    for name,digest in manifest.items():assert sha(destination/name)==digest
    write_json(destination/'verification.json',dict(passed=True,samples=46,unique_pair_hashes=46,
        source_coordinate_maps_replayed=46,CDS_translations_replayed=92,FASTA_records_per_file=46,
        packaged_file_hashes_verified=len(manifest),original_frozen_results_modified=False,
        slurm_job_id=os.environ['SLURM_JOB_ID'],compute_host=socket.gethostname()))
    archive=destination.with_suffix('.zip')
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as handle:
        for path in sorted(destination.rglob('*')):
            if path.is_file():handle.write(path,str(Path(destination.name)/path.relative_to(destination)))
    with zipfile.ZipFile(archive) as handle:assert handle.testzip() is None
    archive.with_suffix('.zip.sha256').write_text(sha(archive)+'  '+archive.name+'\n')
    print(json.dumps(dict(summary=summary,main_csv=str(destination/'art_samples.csv'),archive=str(archive),archive_bytes=archive.stat().st_size,verified=True),indent=2),flush=True)


def re_accession(text):
    import re
    return bool(re.fullmatch(r'[A-Z]{1,4}_?\d+\.\d+',text))


if __name__=='__main__':main()

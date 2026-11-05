"""Package mapped-locus exports, independent analyses, costs and explicit limitations."""
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import shutil

from art_reconstruction_common import ROOT,sha
from extract_art_loci import read_tsv,write_json


def artifact_hashes(folder):
    """Hash deliverables, excluding MMseqs scratch databases and Cygwin links."""
    hashes={};excluded=[]
    for parent,dirs,files in os.walk(folder,followlinks=False):
        for name in list(dirs):
            if name.endswith('_tmp') or name=='latest':
                excluded.append((Path(parent)/name).relative_to(folder).as_posix());dirs.remove(name)
        for name in files:
            if name in ('output_sha256.json','latest'):continue
            path=Path(parent)/name
            hashes[path.relative_to(folder).as_posix()]=sha(path)
    return hashes,excluded


def run():
    base=ROOT/'results/art_expansion';stage=base/'stage1';family=ROOT/'results/art_family_reconstruction'
    metrics=json.loads((stage/'cost_and_yield.json').read_text());recon=json.loads((family/'summary.json').read_text())
    domain=json.loads((family/'domain_phylogeny/summary.json').read_text());verified=json.loads((base/'verification.json').read_text())
    for name in ('candidate_loci.tsv','sequence_family_placement.tsv','observed_locus_architecture.tsv'):
        shutil.copyfile(stage/name,base/name)
    rows=read_tsv(base/'candidate_loci.tsv');new=[r for r in rows if r['record_role']=='public_stage1_pilot']
    panel=read_tsv(base/'stage2/prepared_assembly_sources.tsv')
    summary=dict(stage1_public_records_searched=metrics['records_searched'],new_mapped_RT_candidates=metrics['new_RT_candidates'],
        mapped_RT_loci_total=len(rows),stage2_records_prepared=len(panel),stage2_profile_search_executed=False,stage3_search_executed=False,
        RT_domain_tree_tips=domain['tips'],sequence_component_proteins=recon['replication']['component_proteins'],
        longest_contig_cluster_representatives=recon['replication']['cluster_representatives'],
        verified_NART_reference_count=domain['NART_reference_count'],full_paper_dataset_reproduced=False,
        complete_ART_locus_claims=0,validation_gate=metrics['validation_gate'],verification_passed=verified['passed'])
    write_json(base/'summary.json',summary)
    lines=['# Staged ART search and family reconstruction','',
        '**Completed a bounded annotated-phage pilot and two separate sequence analyses.** The full validation gate remains incomplete, so all newly exported family/locus categories are provisional. Stage 2 inputs are prepared; stage 2/3 profile searches were not executed.','',
        '## Stage 1: additional public annotated phages','',
        '| Genome | Accession/version | Result |','| --- | --- | --- |',
        '| Staphylococcus phage Machias | MW349128.1 | One mapped ART-like RT candidate, QQO92593.1 |',
        '| Bacillus phage AR9 | KU878088.1 | No RT candidate at the available ART/partner thresholds; not a verified NART control |',
        '| Yersinia phage phiR1-37 | NC_016163.1 | No candidate at the available thresholds |','',
        f"The three new genomes contain {metrics['nucleotide_nt']:,} nt and {metrics['new_proteins']:,} annotated translations. With eight baseline genomes, the search covers {metrics['proteins_searched']:,} proteins. Five ambiguous baseline proteins were retained and flagged, rather than discarded. ART RT, type I, type II and generic Pfam RVT_1 profiles were searched over the complete mapped protein set; the two type III families remain unavailable.", '',
        'Both RT-first and partner-first searches recover the Machias candidate. Partner-first records all first-three-upstream genes, including opposing orientations; a generic GNAT label never substitutes for partner-family homology. Partner windows, orphan reporting and approximately 6-kb flanks are supported independently of a strong ART-profile hit.', '']
    lines += ['Partner evidence remains limited: the type I profile has one seed and no positive holdout; type II has only one conservative lineage group. These profiles are retrieval aids rather than well-calibrated universal family classifiers. Type III first and second genes remain separate unavailable families.', '']
    for r in new:
        c=json.loads(r['rt_coordinates_and_strand']);p=json.loads(r['art_profile_score_and_coverage']);ntd=json.loads(r['ntd_boundary'])
        lines += [f"Machias RT: **{r['RT_protein_id']}**, source interval **[{c['start0']}, {c['end0']})**, strand {c['strand']:+d}; 1-based inclusive **{c['start0']+1}–{c['end0']}**. ART profile: {p['bits']} domain bits, {p['coverage']:.3f} coverage. The 517-aa source CDS translates correctly. Its N-terminal interval is [0, {ntd['end0']}) aa, ending before the first significant {ntd['profile_accession']} domain envelope. Downstream QQO92594.1 matches the type II partner profile.", '',
        'The exploratory upstream array scan yields **no supported array call** for Machias. This is not a biological absence claim: the scanner still has validation failures and this new locus has no calibrated cohort array p-value. The RT-domain tree shows ART-reference affinity. That sequence result is retained independently of its unresolved array architecture. This is an additional candidate in our panel, not a claim of a newly discovered enzyme family.', '']
    lines += ['## Two analyses, separate files','',
        '- `sequence_family_placement.tsv`: protein/profile evidence, generic RT-domain/NTD boundaries and sequence-tree affinity.',
        '- `observed_locus_architecture.tsv`: partners, non-coding distance, repeats/spacers and array status.',
        '- `candidate_loci.tsv`: a joined review/export table containing every requested field, with stable source and coordinate-map links.', '',
        'No array or partner feature supplies the sequence-placement label. Protein-only myRT/curated references can enter sequence trees but cannot become genomic locus rows. The nine baseline ART reference labels are inherited from the curated publication-based panel; they are not independent validation of the historical family definition. Baseline queries have identical reference counterparts in the tree, so their nearest-reference agreement is not a generalization estimate. Machias is also a related Staphylococcus phage, not a new unrelated lineage.', '',
        '## Separate published-procedure reconstruction','',
        f"The executable replication branch uses DIAMOND all-against-all (`--very-sensitive`, E≤1e-5, no top-hit cap), undirected links at **DIAMOND score ≥140 bits**, the MarsHill-connected component, optional mapped database-cluster expansion at ≥400 aa, MMseqs2 90% identity/80% mutual coverage, longest-contig representatives, FAMSA and FastTree LG. The component has **{recon['replication']['component_proteins']} RT proteins**, no myRT reference protein, and **{recon['replication']['cluster_representatives']} representatives**. This 140-bit graph threshold is not an HMM profile cutoff.", '',
        'The original paper’s 12-phage/2,019-reference membership, representative database and 706 extra cluster members were not supplied. This run uses nine baseline RTs, one additional mapped candidate and 2,328 canonical myRT references from the available 2,339-entry collection. It does **not** recreate the reported 823 proteins, 230 clusters or 95 ART representatives. The 350–900-aa range is restricted to the replication pool; out-of-range candidates remain in a sensitivity FASTA and the independent domain analysis. Additional database-cluster relations may be supplied through `--cluster-members`; every expanded member must have a real nucleotide mapping.', '',
        'The architecture-defined “smallest clade containing type I + ≥500-nt non-coding loci” is recorded **only** under `replication/architecture_defined_clade.json`. There is only one eligible representative anchor in this panel, so a meaningful multi-lineage architecture-defined clade cannot be reconstructed. Its output is not used as the sequence classifier.', '',
        'On Windows, the first MMseqs run exceeded its 1-GB prefilter allowance. The bounded run uses `-k 5 --split-memory-limit 1G`; all requested identity/coverage settings remain unchanged. This computational setting and the missing source dataset prevent a claim of exact numerical paper replication. Commands, failures, completed logs and checksums are retained. The upstream wrapper’s UAC branch is not used.', '',
        '## Independent RT-domain phylogeny','',
        f"The domain tree has **{domain['tips']} tips** and **{domain['columns']} alignment columns**, including ART references, retrons, DGRs, CRISPR-associated RTs, group II introns, group-II-like RTs, UG and Abi proteins. Reference selection first clusters at 50% identity/80% mutual coverage, then samples by labelled class and lineage; it does not use array/partner annotations.", '',
        'Domains are the first significant Pfam RVT_1 envelope ±40 aa, aligned with MAFFT L-INS-i, trimmed at 50% column occupancy and inferred with FastTree LG+Gamma. This is an approximate evaluation tree, not the paper’s trimAl/IQ-TREE 774-tip analysis. FastTree local support is not the paper’s UFBoot/SH-aLRT statistic. Unplaced reference proteins are reported rather than called non-RT.', '',
        '**No verified NART reference is available.** ART-like queries therefore receive `ART_reference_affinity_NART_unresolved`; no definitive ART/NART boundary or assignment is claimed. Additional class-labelled FASTA references can be added with `--extra-reference-fasta` and an evidence/source manifest. A long N-terminus, RNA-polymerase neighbor or absence of arrays is not used to invent a NART label.', '',
        '## Stage 2 panel prepared, not searched','',
        '| Source | Fixed accession | Preparation |','| --- | --- | --- |']
    for r in panel:
        method='Retain deposited CDS annotations' if r['annotation_method']=='source_original' else '191 Pyrodigal-gv candidate CDSs; original unannotated record retained'
        lines.append(f"| {r['source_record_status']} | {r['accession']} | {method} |")
    lines += ['',
        'PX848769.1 is a wastewater-metagenome-derived Klebsiella phage assembly with 279 deposited CDSs. PV020493.1 is an UNVERIFIED partial Staphylococcus viral contig with no original CDS calls. Its preparation uses Pyrodigal-gv 0.3.2/Pyrodigal 3.7.1, meta mode; the selected model is recorded. Source nucleotide identity, predicted coordinates, strand, partial flags and genetic code are preserved. No gene call is treated as experimental verification, and the unverified/partial source status remains visible.', '',
        'The exact accessions, sequences, original annotations, selection query and file hashes are frozen. The gate prevents stage 2 profile/array expansion until the unresolved copy-accuracy and reference-validation requirements are addressed. Stage 3 additionally requires measured stage 2 cost/yield, an explicit collection manifest and a record budget. There is no automatic whole-database download.', '',
        '## Measured cost and yield','',
        f"Latest cached-input stage 1 run: **{metrics['search_seconds']:.3f} s** for the four HMM searches over {metrics['proteins_searched']:,} proteins; **{metrics['elapsed_seconds']:.3f} s** for source ingestion, searches, locus extraction and exploratory architecture export. New source files occupy {metrics['input_file_bytes']:,} bytes. The first run including public retrieval took 16.719 s (HMM search 0.703 s). These are wall times on this workstation, not cluster-wide throughput estimates.", '',
        f"The latest complete reconstruction took **{recon['elapsed_seconds']:.3f} s**, including the full DIAMOND reference-pool comparison and domain analysis. Its all-against-all TSV is {(family/'replication/all_against_all.tsv').stat().st_size:,} bytes. Reference-pool all-against-all work is a separate cost from searching phage proteomes. Larger-run estimates must account for pool size, candidate count, mapping retrieval and array-null tests.", '',
        'Yield is one mapped RT candidate from three purposively chosen genomes. This panel was chosen for a near-relative challenge and jumbo-phage comparisons; it does not estimate ART prevalence, discovery sensitivity or a population hit rate. Only two partner-family profiles are available, and annotation/gene-calling errors can still hide RTs.', '',
        '## Verification and rerun','',
        f"Independent verification checks {verified['mapped_candidate_loci']} nucleotide-linked candidate loci, oriented RT CDSs, {verified['repeat_source_mappings_checked']} original-genome repeat mappings, DIAMOND connectivity, longest-contig representative choices, domain extraction and alignment-column maps, and both stage gates. The repeat detector’s earlier biological benchmark remains **27/30**; this work does not relabel its failures as successful validation.", '',
        '```powershell', '.venv/Scripts/python.exe scripts/setup_art_reconstruction_tools.py', '.venv/Scripts/python.exe scripts/art_staged_search.py --stage 1',
        '.venv/Scripts/python.exe scripts/reconstruct_art_family.py', '.venv/Scripts/python.exe scripts/art_staged_search.py --stage 2 --prepare-only',
        '.venv/Scripts/python.exe scripts/art_staged_search.py --stage 2', '# The preceding search reports the unmet gate and does not run.',
        '.venv/Scripts/python.exe scripts/verify_art_expansion.py', '.venv/Scripts/python.exe scripts/report_art_expansion.py',
        ".venv/Scripts/python.exe -m unittest discover -s scripts -p 'test_*art*.py'", '```', '',
        'Canonical outputs: [candidate loci](candidate_loci.tsv), [sequence placement](sequence_family_placement.tsv), [locus architecture](observed_locus_architecture.tsv), [full-length replication tree](../art_family_reconstruction/replication/full_length.FastTree.nwk), [RT-domain tree](../art_family_reconstruction/domain_phylogeny/RT_domain.FastTree.nwk), [stage 2 manifest](../../config/art_stage2_assemblies.tsv), [verification](verification.json).', '',
        'Method sources: supplied Yoon PDF, family/phylogeny methods (PDF pp.31–33); [ART preprint](https://doi.org/10.64898/2026.10.05.755809); [DIAMOND](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options); [MMseqs2](https://github.com/soedinglab/MMseqs2/wiki); [FAMSA](https://github.com/refresh-bio/FAMSA); [FastTree](https://github.com/morgannprice/fasttree); [Pyrodigal-gv](https://github.com/althonos/pyrodigal-gv); [Machias source](https://www.ebi.ac.uk/ena/browser/view/MW349128.1).', '']
    (base/'REPORT.md').write_text('\n'.join(lines),encoding='utf8')
    for folder in (base,family):
        code=folder/'code';code.mkdir(exist_ok=True)
        for name in ('art_reconstruction_common.py','art_staged_search.py','reconstruct_art_family.py','verify_art_expansion.py','report_art_expansion.py','setup_art_reconstruction_tools.py','test_art_reconstruction.py'):
            shutil.copyfile(ROOT/'scripts'/name,code/name)
        write_json(folder/'frozen_inputs.json',{str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'config/art_expansion.json',ROOT/'config/art_stage2_assemblies.tsv',ROOT/'requirements-art-reconstruction.txt',ROOT/'.tools/art-reconstruction/provenance.json',ROOT/'results/art_validation/summary.json']})
        hashes,excluded=artifact_hashes(folder)
        write_json(folder/'checksum_scope.json',dict(excluded_MMseqs_scratch_directories=excluded,scientific_exports_and_inputs_hashed=True))
        hashes['checksum_scope.json']=sha(folder/'checksum_scope.json')
        write_json(folder/'output_sha256.json',hashes)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':run()

"""Re-run the bounded archived-lead pilot; remote queues are managed separately."""
import json
import subprocess
import sys
from art_reconstruction_common import ROOT, Runner
from art_type1_pipeline import profile_setup, process, null_assessment, diversity
from art_type1_archive_leads import main as leads, download_current_prefixes, extend_two_prefixes, reconcile_current
from art_type1_tree import main as tree
from finalize_art_type1_public import main as finalize
from verify_art_type1_public import main as verify


def main():
    cfg=json.loads((ROOT/"config/art_type1_public.json").read_text())
    out=ROOT/cfg["output"];out.mkdir(parents=True,exist_ok=True)
    print("Build exploratory retrieval profile",flush=True)
    profile_setup(out,cfg,Runner(out))
    print("Retrieve fixed four-target panel and current assemblies",flush=True)
    leads();download_current_prefixes();extend_two_prefixes();leads();reconcile_current()
    print("Process mapped loci and test array backgrounds",flush=True)
    process(cfg,out);null_assessment(cfg,out)
    print("Sequence-only tree and primary-partner diversity",flush=True)
    tree();diversity(cfg,out);finalize()
    with (out/"software_tests.log").open("w",encoding="utf8") as stream:
        subprocess.run([sys.executable,"-m","unittest","discover","-s","scripts","-p","test_art*.py"],
            cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT,check=True)
    verify()


if __name__=="__main__":main()

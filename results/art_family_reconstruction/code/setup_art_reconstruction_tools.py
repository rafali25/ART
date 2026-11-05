"""Install pinned official reconstruction binaries locally; record URLs and hashes."""
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tarfile
import urllib.request
import zipfile

ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/'.tools/art-reconstruction'
SOURCES={
 'diamond':('https://github.com/bbuchfink/diamond/releases/download/v2.1.15/diamond-windows.zip','diamond.zip'),
 'mmseqs':('https://github.com/soedinglab/MMseqs2/releases/download/17-b804f/mmseqs-win64.zip','mmseqs.zip'),
 'famsa':('https://github.com/refresh-bio/FAMSA/releases/download/v2.5.2/famsa-v2.5.2-x64_windows.tar.gz','famsa.tar.gz'),
 'fasttree':('https://raw.githubusercontent.com/morgannprice/fasttree/a780aa4f6a18791ea76c04afaf37aafe43e51d59/FastTree.c','FastTree.c'),
 'RVT_1':('https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00078/?annotation=hmm','PF00078.hmm.gz'),
}

def download(item):
 name,(url,file)=item; path=DEST/file
 if not path.exists():
  with urllib.request.urlopen(url,timeout=60) as response:path.write_bytes(response.read())
 return name,dict(url=url,file=str(path.relative_to(ROOT)),sha256=hashlib.sha256(path.read_bytes()).hexdigest())

def run():
 DEST.mkdir(parents=True,exist_ok=True)
 with ThreadPoolExecutor(max_workers=4) as ex:provenance=dict(ex.map(download,SOURCES.items()))
 for name in ('diamond','mmseqs','famsa'):
  folder=DEST/name;folder.mkdir(exist_ok=True);path=ROOT/provenance[name]['file']
  if path.suffix=='.zip':
   with zipfile.ZipFile(path) as archive:
    for member in archive.infolist():
     if not (folder/member.filename).resolve().is_relative_to(folder.resolve()):raise ValueError('Unsafe archive member')
    archive.extractall(folder)
  else:
   with tarfile.open(path) as archive:archive.extractall(folder,filter='data')
 text=gzip.decompress((DEST/'PF00078.hmm.gz').read_bytes())
 if b'ACC   PF00078.' not in text or b'NAME  RVT_1' not in text:raise ValueError('Not a Pfam RVT_1 HMM')
 (DEST/'PF00078.hmm').write_bytes(text)
 env=os.environ.copy();env['PATH']=str(ROOT/'.tools/cygwin/bin')+os.pathsep+env['PATH']
 command=[str(ROOT/'.tools/cygwin/bin/gcc.exe'),'-O3','-o',str(DEST/'FastTree.exe'),str(DEST/'FastTree.c'),'-lm']
 with (DEST/'FastTree.build.log').open('w') as handle:subprocess.run(command,env=env,stdout=handle,stderr=subprocess.STDOUT,check=True)
 binaries={name:str(next((DEST/name).rglob(pattern)).relative_to(ROOT)) for name,pattern in [('diamond','diamond.exe'),('mmseqs','mmseqs.bat'),('famsa','famsa.exe')]}
 binaries['fasttree']=str((DEST/'FastTree.exe').relative_to(ROOT));binaries['RVT_1']=str((DEST/'PF00078.hmm').relative_to(ROOT))
 # Use local hardlinks for BusyBox helpers; never invoke the wrapper's optional
 # UAC/symlink installation branch. Reconstruction calls mmseqs.exe directly.
 helper=Path(binaries['mmseqs']).parent/'bin/busybox.exe'
 subprocess.run([str(ROOT/helper),'--install',str((ROOT/helper).parent)],check=True)
 provenance['binaries']=binaries
 provenance['binary_sha256']={name:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for name,path in binaries.items()}
 (DEST/'provenance.json').write_text(json.dumps(provenance,indent=2))
 print(json.dumps(binaries,indent=2))

if __name__=='__main__':run()

"""Install project-local MAFFT and compile official HMMER using a local Cygwin tree.

No WSL repair, registry edits, PATH changes, or system installation is performed.
Cygwin package hashes are checked against its downloaded package index.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib
import io
import json
import lzma
import re
import shutil
import subprocess
import tarfile
import urllib.request
import zipfile
import zstandard

ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / '.tools'
CACHE = TOOLS / 'downloads'
CYG = TOOLS / 'cygwin'
MIRROR = 'https://mirrors.kernel.org/sourceware/cygwin/'


def download(url, path):
    if not path.exists():
        urllib.request.urlretrieve(url, path)
    return path


def extract_tar(path, dest, cygwin=False):
    raw = path.read_bytes()
    if path.suffix == '.zst':
        raw = zstandard.ZstdDecompressor().stream_reader(io.BytesIO(raw)).read()
    links = []
    with tarfile.open(fileobj=io.BytesIO(raw), mode='r:*') as archive:
        for item in archive:
            name = item.name.removeprefix('./')
            # POSIX manual-page names may contain characters forbidden by NTFS.
            # They are not executable build dependencies.
            if any(re.search(r'[<>:"|?*]', part) for part in name.split('/')):
                continue
            if cygwin:
                for prefix in ('usr/bin/', 'usr/lib/'):
                    if name.startswith(prefix):
                        name = name[4:]
                        break
            target = (dest / name).resolve()
            if not target.is_relative_to(dest.resolve()):
                raise ValueError(f'Unsafe archive path: {name}')
            if item.isdir():
                target.mkdir(parents=True, exist_ok=True)
            elif item.isfile():
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(archive.extractfile(item).read())
            elif item.issym() or item.islnk():
                links.append((target, item.linkname, item.islnk()))
    # Required shell alias is supplied explicitly below. Other documentation
    # symlinks are not necessary for compiling HMMER.
    return links


def materialize_links(links):
    # NTFS-safe copies for executable aliases and compiler-library hardlinks.
    for _ in range(5):
        for target, link, hardlink in links:
            if not target.is_relative_to(CYG) or target.exists():
                continue
            if not target.is_relative_to(CYG/'bin') and not target.is_relative_to(CYG/'lib'):
                continue
            if link.startswith('/') or hardlink:
                name=link.lstrip('/').removeprefix('./')
                if name.startswith(('usr/bin/','usr/lib/')):name=name[4:]
                src=(CYG/name).resolve()
            else:
                src=(target.parent/link).resolve()
            if not src.exists() and src.suffix!='.exe':
                src=Path(str(src)+'.exe')
            if src.is_file() and src.is_relative_to(CYG):
                if src.suffix=='.exe' and target.suffix!='.exe':target=Path(str(target)+'.exe')
                target.parent.mkdir(parents=True,exist_ok=True)
                shutil.copy2(src,target)


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    if not (CACHE / 'cygwin_setup.ini').exists():
        compressed = download(MIRROR + 'x86_64/setup.xz', CACHE / 'cygwin_setup.xz')
        (CACHE / 'cygwin_setup.ini').write_bytes(lzma.decompress(compressed.read_bytes()))
    index = (CACHE / 'cygwin_setup.ini').read_text(encoding='utf-8')
    packages = {}
    for block in index.split('\n@ ')[1:]:
        current = re.split(r'\n\[(?:prev|test)\]', block)[0]
        name = current.splitlines()[0]
        match = re.search(r'^install: (\S+) (\d+) (\w+)$', current, re.M)
        deps = re.search(r'^depends2: (.*)$', current, re.M)
        if match:
            packages[name] = dict(name=name, path=match[1], size=int(match[2]),
                                  sha512=match[3], dependencies=re.findall(
                                      r'(?:^|,)\s*([\w+.-]+)', deps[1]) if deps else [])
    selected = {}

    def include(name):
        if name.startswith('_') or name in selected:
            return
        p = packages[name]
        selected[name] = p
        for dep in p['dependencies']:
            include(dep)

    for name in ('gcc-core', 'make', 'bash', 'coreutils', 'sed', 'grep', 'gawk', 'diffutils'):
        include(name)
    print(f'Local Cygwin build dependencies: {len(selected)} packages, '
          f'{sum(p["size"] for p in selected.values()) / 1e6:.1f} MB', flush=True)

    def fetch(p):
        path = download(MIRROR + p['path'], CACHE / Path(p['path']).name)
        if hashlib.sha512(path.read_bytes()).hexdigest() != p['sha512']:
            raise ValueError(f'Package checksum mismatch: {p["name"]}')
        return p, path

    CYG.mkdir(parents=True, exist_ok=True)
    links=[]
    with ThreadPoolExecutor(max_workers=6) as pool:
        for p, path in pool.map(fetch, selected.values()):
            links.extend(extract_tar(path, CYG, cygwin=True))
            print(f'Extracted {p["name"]}', flush=True)
    materialize_links(links)
    shutil.copy2(CYG / 'bin/bash.exe', CYG / 'bin/sh.exe')
    shutil.copy2(CYG / 'bin/ld.bfd.exe', CYG / 'bin/ld.exe')
    (CYG / 'tmp').mkdir(exist_ok=True)
    (CYG / 'etc').mkdir(exist_ok=True)
    (CYG / 'etc/fstab').write_text('', encoding='ascii')
    (TOOLS / 'cygwin_packages.json').write_text(
        json.dumps(list(selected.values()), indent=2), encoding='utf-8')

    url = 'https://mafft.ddbj.nig.ac.jp/alignment/software/mafft-7.526-win64-signed.zip'
    mafft = download(url, CACHE / 'mafft-7.526-win64-signed.zip')
    if not (TOOLS / 'mafft-win').exists():
        with zipfile.ZipFile(mafft) as z:
            z.extractall(TOOLS)
    source_url = 'https://eddylab.org/software/hmmer/hmmer-3.4.tar.gz'
    source = download(source_url, CACHE / 'hmmer-3.4.tar.gz')
    if not (TOOLS / 'hmmer-3.4/configure').exists():
        extract_tar(source, TOOLS)
    manifest = dict(mafft_url=url, mafft_sha256=hashlib.sha256(mafft.read_bytes()).hexdigest(),
                    hmmer_url=source_url, hmmer_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
    (TOOLS / 'source_downloads.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    posix = '/cygdrive/' + str(ROOT.drive)[0].lower() + ROOT.as_posix()[2:]
    command = (f'export PATH=/usr/bin:/bin; cd {posix}/.tools/hmmer-3.4 && '
               f'./configure --prefix={posix}/.tools/hmmer-install && make -j4 && make install')
    with (TOOLS / 'hmmer_build.log').open('w', encoding='utf-8') as log:
        subprocess.run([str(CYG / 'bin/bash.exe'), '-c', command], check=True,
                       stdout=log, stderr=subprocess.STDOUT)
    print('MAFFT downloaded; HMMER compiled and installed locally.', flush=True)


if __name__ == '__main__':
    main()

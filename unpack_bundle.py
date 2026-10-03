from pathlib import Path
import hashlib,tarfile,shutil
root=Path(__file__).resolve().parent
archive=root/'verified_bundle_v1.tar.xz'
expected=(root/'BUNDLE_SHA256.txt').read_text().split()[0]
assert hashlib.sha256(archive.read_bytes()).hexdigest()==expected,'Archive SHA-256 mismatch'
dest=root/'bundle'
if dest.exists(): raise RuntimeError('bundle already exists; extract in a fresh directory')
with tarfile.open(archive,'r:xz') as stream:
    for m in stream.getmembers():
        target=(root/m.name).resolve()
        target.relative_to(dest.resolve())
        if m.isdir(): target.mkdir(parents=True,exist_ok=True)
        elif m.isfile():
            target.parent.mkdir(parents=True,exist_ok=True)
            with stream.extractfile(m) as src,target.open('wb') as out: shutil.copyfileobj(src,out)
        else: raise RuntimeError('Unexpected link or special archive member: '+m.name)
print('Verified bundle extracted. Run: python -B bundle/verify_release.py')

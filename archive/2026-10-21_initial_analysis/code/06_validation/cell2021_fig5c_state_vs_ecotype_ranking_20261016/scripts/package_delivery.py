"""Copy verified new artifacts to the two requested delivery locations; no source mutation."""
from pathlib import Path
import hashlib,json,shutil,zipfile
import pandas as pd

OUT=Path(__file__).resolve().parents[1]
CHAT=Path('C:/Users/马泽凯/Documents/Codex/2026-10-16/referenced-chatgpt-conversation-this-is-an/outputs')
TARGETS=[Path('D:/IBD_ecotype/初步结果汇报/相关结果/20261015_BayesPrism_CORE12_53states_7ecotypes/14_Cell2021_Fig5C_rankings_20261016'),CHAT/'IBD_Fig5C_20261016']
summary=json.loads((OUT/'audit/delivery_QA_summary.json').read_text())
assert summary['Failed']==0 and summary['Figures']==11 and summary['InputsUnchanged']
for target in TARGETS:
    assert not target.exists(), 'Refuse to overwrite existing delivery folder: '+str(target)

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def files():return sorted(p for p in OUT.rglob('*') if p.is_file() and '__pycache__' not in p.parts)
def manifest():
    rows=[dict(RelativePath=str(p.relative_to(OUT)).replace('\\','/'),Bytes=p.stat().st_size,SHA256=sha(p)) for p in files() if p.name!='artifact_SHA256.tsv']
    pd.DataFrame(rows).to_csv(OUT/'audit/artifact_SHA256.tsv',sep='\t',index=False)
manifest()
for target in TARGETS:
    target.parent.mkdir(exist_ok=True,parents=True)
    shutil.copytree(OUT,target,ignore=shutil.ignore_patterns('__pycache__'))
verifications=[]
for target in TARGETS:
    for p in files():
        other=target/p.relative_to(OUT)
        equal=other.exists() and sha(p)==sha(other)
        verifications.append(dict(Destination=str(target),RelativePath=str(p.relative_to(OUT)),ExactCopy=equal))
        assert equal,str(other)
pd.DataFrame(verifications).to_csv(OUT/'audit/delivery_copy_verification.tsv',sep='\t',index=False)
manifest()
for target in TARGETS:
    for rel in ['audit/delivery_copy_verification.tsv','audit/artifact_SHA256.tsv']:
        shutil.copy2(OUT/rel,target/rel)
archive=CHAT/'IBD_Fig5C_20261016_完整交付.zip'
with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in sorted(TARGETS[1].rglob('*')):
        if p.is_file():z.write(p,p.relative_to(TARGETS[1]))
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    archived_files=len(z.namelist())
    assert archived_files==len([p for p in TARGETS[1].rglob('*') if p.is_file()])
checksum=dict(Archive=str(archive),SHA256=sha(archive),Bytes=archive.stat().st_size,ArchivedFiles=archived_files,CopyVerificationRows=len(verifications),AllCopiesExact=True,Directories=list(map(str,TARGETS)))
(CHAT/'IBD_Fig5C_20261016_交付校验.json').write_text(json.dumps(checksum,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(checksum,ensure_ascii=False,indent=2))

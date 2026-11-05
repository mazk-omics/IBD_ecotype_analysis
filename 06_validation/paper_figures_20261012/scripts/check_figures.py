# 检查最终图件的字体、碰撞、对齐及冻结输入SHA；不进行新的生物学计算。
# 服务器ibd_refaudit运行；输出逐图QA和下载清单，允许同一SSH隧道逐文件回传。
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import pandas as pd
out=Path('/home/mazekai/IBD_EcoTyper/06_validation/paper_figures_20261012')
qa=out/'scripts/figure_qa'
rows=[]
for pdf in sorted((out/'figures').glob('Fig*.pdf')):
    if '.collision-audit' in pdf.name: continue
    for script,extra in [('audit_pdf_text.py',['--min-pt','5','--json']),('audit_figure_collisions.py',['--json-out',str(pdf.with_suffix('.collision-audit.json')),'--overlay-pdf',str(pdf.with_suffix('.collision-audit.pdf'))])]:
        result=subprocess.run([sys.executable,str(qa/script),str(pdf)]+extra,capture_output=True,text=True)
        report=pdf.with_suffix('.text-audit.json' if script=='audit_pdf_text.py' else '.collision-audit.log')
        report.write_text(result.stdout+'\n'+result.stderr)
        rows.append({'Figure':pdf.name,'Check':script,'ExitCode':result.returncode})
        print(pdf.name,script,result.returncode,flush=True)
pd.DataFrame(rows).to_csv(out/'reports/figure_QA_status.tsv',sep='\t',index=False)
manifest=pd.read_csv(out/'source_data/frozen_inputs_SHA256.tsv',sep='\t')
checks=[]
for p in manifest.itertuples():
    actual=hashlib.sha256(Path(p.Path).read_bytes()).hexdigest()
    checks.append({'Path':p.Path,'Unchanged':actual==p.SHA256})
assert all(p['Unchanged'] for p in checks)
pd.DataFrame(checks).to_csv(out/'reports/frozen_inputs_after_check.tsv',sep='\t',index=False)
print('FROZEN_INPUTS_UNCHANGED',len(checks),flush=True)

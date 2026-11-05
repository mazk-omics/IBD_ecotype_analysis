# 盘点并直接复用官方已生成的Fig2A对应分类型热图及Fig4A生态型热图，不重画这些数据。
# 服务器ibd_refaudit运行；原图逐字节复制并核对SHA，12类热图仅拼为多页PDF便于阅读。
from pathlib import Path
import shutil
import hashlib
import pandas as pd
import pymupdf
root=Path('/home/mazekai/IBD_EcoTyper')
out=root/'06_validation/paper_figures_20261012'
original=root/'05_EcoTyper/output/MSCCR_CORE12/ctspec_region_ibd'
reused=out/'reused_originals'
reused.mkdir(exist_ok=True)
states=pd.read_csv(out/'source_data/frozen_states.tsv',sep='\t')
rows=[]
atlas=pymupdf.open()
for ct in states.CellType.unique():
    for suffix in ['pdf','png']:
        p=original/ct/('state_assignment_heatmap.'+suffix)
        target=reused/('Fig2A_'+ct+'_official.'+suffix)
        shutil.copyfile(p,target)
        sha=hashlib.sha256(p.read_bytes()).hexdigest()
        assert hashlib.sha256(target.read_bytes()).hexdigest()==sha
        rows.append({'PaperPanel':'Fig2A','CellType':ct,'Source':str(p),'Deliverable':str(target.relative_to(out)),'SHA256':sha,'Action':'REUSED_BYTE_IDENTICAL'})
    doc=pymupdf.open(original/ct/'state_assignment_heatmap.pdf')
    atlas.insert_pdf(doc)
    doc.close()
atlas.save(reused/'Fig2A_official_profiles_12pages.pdf')
atlas.close()
for suffix in ['pdf','png']:
    p=original/'Ecotypes'/('heatmap_assigned_samples_viridis.'+suffix)
    target=reused/('Fig4A_official_ecotype_abundance.'+suffix)
    shutil.copyfile(p,target)
    sha=hashlib.sha256(p.read_bytes()).hexdigest()
    assert hashlib.sha256(target.read_bytes()).hexdigest()==sha
    rows.append({'PaperPanel':'Fig4A','CellType':'All','Source':str(p),'Deliverable':str(target.relative_to(out)),'SHA256':sha,'Action':'REUSED_BYTE_IDENTICAL'})
pd.DataFrame(rows).to_csv(out/'reports/reused_figure_inventory.tsv',sep='\t',index=False)
print('REUSED',len(rows),'files; assembled 12-page atlas without rerendering original plots')

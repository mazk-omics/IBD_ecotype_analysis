# 核对 SeqOut 新命中的肠道单细胞研究，下载小型元数据而不下载表达文件。
# 在本机 SeqOut Python 3.13 环境运行；输出在 recovery_candidates_20261007/scRNA_extension。
# 样本记录、供者及实验通道可能不同，下面输出用于人工判断而不是自动给队列排名。
from pathlib import Path
import json
import os
import subprocess
import pandas as pd

# Step1. 从实际搜索命中中选取黏膜全谱系、上皮与 CD 间质候选进行深入核对。
project = Path(__file__).resolve().parents[3]
outdir = project / 'IBD_ecotype_analysis/06_validation/recovery_candidates_20261007/scRNA_extension'
seqout = project / '.tools/seqout/environment/Scripts/seqout.exe'
environment = os.environ.copy()
environment['PYTHONUTF8'] = '1'
accessions = ['GSE243625', 'GSE250487', 'GSE214695', 'GSE190564', 'GSE150115', 'GSE134881', 'GSE185224', 'GSE164985', 'GSE116222', 'GSE290695', 'GSE266546', 'GSE275145', 'GSE255722']
record = []
study_rows = []
sample_rows = []
file_rows = []

# Step2. 保存原始元数据、链接及逐样本属性；重复出现的文件链接不意味着需要重复下载。
for accession in accessions:
    path = outdir / (accession + '_metadata.json')
    command = [str(seqout), 'download', accession, '-o', str(path)]
    result = subprocess.run(command, capture_output=True, text=True, encoding='utf-8', env=environment, timeout=180)
    (outdir / (accession + '_console.txt')).write_text(result.stdout + result.stderr, encoding='utf-8')
    record.append({'accession': accession, 'command': command, 'returncode': result.returncode, 'output_exists': path.exists()})
    (outdir / 'metadata_execution_record.json').write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding='utf-8')
    if result.returncode != 0 or not path.exists():
        print(accession, '读取失败，见 console', flush=True)
        continue
    data = json.loads(path.read_text(encoding='utf-8'))
    study = data['project']
    study_rows.append({'accession': accession, 'title': study.get('title'), 'summary': study.get('summary'), 'design': study.get('overall_design'), 'geo_sample_records': len(data.get('samples', [])), 'pubmed': json.dumps(study.get('pubmed_ids', [])), 'study_supplementary': json.dumps(study.get('supplementary_data', []), ensure_ascii=False)})
    for item in study.get('supplementary_data', []):
        url = item[0] if isinstance(item, list) else item  # SeqOut 的系列链接带有 [URL, 文件类型] 两列。
        if isinstance(url, str) and url not in ['', 'NONE']:
            file_rows.append({'dataset': accession, 'sample': '', 'url': url, 'downloaded': False})
    for sample in data.get('samples', []):
        channels = sample.get('channels') or []
        attributes = channels[0].get('characteristics', {}) if channels else {}
        sample_rows.append({'dataset': accession, 'sample': sample['accession'], 'title': sample.get('title'), 'source': channels[0].get('source') if channels else None, 'attributes': json.dumps(attributes, ensure_ascii=False), 'description': sample.get('description')})
        for url in sample.get('supplementary_data', []):
            if isinstance(url, str) and url not in ['', 'NONE']:
                file_rows.append({'dataset': accession, 'sample': sample['accession'], 'url': url, 'downloaded': False})
    print(accession, '样本记录:', len(data.get('samples', [])), '标题:', study.get('title'), flush=True)

# Step3. 保存便于审阅的中间表，供者数和细胞注释是否可用由原始资料确认。
pd.DataFrame(study_rows).to_csv(outdir / '重点研究设计.tsv', sep='\t', index=False)
pd.DataFrame(sample_rows).to_csv(outdir / '重点研究样本属性.tsv', sep='\t', index=False)
pd.DataFrame(file_rows).to_csv(outdir / '重点研究文件链接_尚未下载.tsv', sep='\t', index=False)
print('研究:', len(study_rows), '样本记录:', len(sample_rows), '文件链接:', len(file_rows), flush=True)

# recovery完成并检查通过后生成回传清单；原始计数和置换缓存保留在服务器。
# 在服务器ibd_refaudit环境运行；本次目录内的小表、官方输出、日志、执行脚本通过既有SSH隧道回传。
from pathlib import Path
import hashlib
import json

# Step1. 限定本轮输出和已完成状态，不发送原始大矩阵或其他项目文件
root = Path('/home/mazekai/IBD_EcoTyper/06_validation/GSE214695_scRNA_recovery_20261009')
assert (root / 'audit/final_data_check.json').exists()
assert (root / 'formal_exit_code.txt').read_text().strip() == '0'
paths = [p for p in (root / 'audit').rglob('*') if p.is_file()]
paths += [p for p in (root / 'formal').rglob('*') if p.is_file()]
paths += [p for p in root.glob('*') if p.is_file() and p.suffix in ['.txt', '.tsv', '.json', '.log', '.sha256', '.py', '.sh']]
paths += [p for p in (root / 'input').glob('*') if p.is_file() and p.name not in ['GSE214695_Pilot18_20261009_CPM.txt', 'GSE214695_Pilot18_20261009_Smoke_CPM.txt']]
paths = sorted(set(paths) - {root / 'download_manifest.json'})
assert all(p.resolve().is_relative_to(root.resolve()) for p in paths)
assert max(p.stat().st_size for p in paths) < 100_000_000, '异常大文件应先核对，不能把随机缓存当审阅附件'

# Step2. 文件逐个记录字节数与哈希，供本地接收逐项检查
manifest = []
for path in paths:
    with path.open('rb') as handle:
        digest = hashlib.file_digest(handle, 'sha256').hexdigest()
    manifest.append({'path': path.relative_to(root).as_posix(), 'bytes': path.stat().st_size, 'sha256': digest})
(root / 'download_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
print('回传文件：', len(manifest), '；总MB：', round(sum(x['bytes'] for x in manifest) / 1e6, 2), flush=True)

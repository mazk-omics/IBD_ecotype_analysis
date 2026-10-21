# 解包通过同一SSH会话取回的小型审计文件，不传输原始矩阵或推断大对象。
# 本机运行；只允许写入本次独立目录，防止压缩包路径越界。
from pathlib import Path
import sys
import tarfile

out = Path(__file__).resolve().parents[1]
archive = Path(sys.argv[1])
with tarfile.open(archive, 'r:gz') as handle:
    members = handle.getmembers()
    for member in members:
        target = (out / member.name).resolve()
        assert target.is_relative_to(out.resolve()) and not member.issym() and not member.islnk()
    handle.extractall(out, members=members, filter='data')
print('AUDIT_DELIVERY_EXTRACTED', len(members), 'files to', out)

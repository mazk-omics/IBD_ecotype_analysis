# 核查公开仓库的固定提交与文件清单，区分教程、EcoTyper接口与CIBERSORTx数值内核。
# 本机Python运行，只读GitHub公开API，不克隆或执行第三方程序；保存来源与检索范围。
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError
import json

# Step1. 只核查本轮网页已识别的三个仓库；不把有限检索说成整个GitHub都没有源码
project = Path(__file__).resolve().parents[2]
outdir = project / 'CIBERSORTx_server_plan_20261011/public_repository_audit'
outdir.mkdir(exist_ok=True)
repositories = ['digitalcytometry/ecotyper', 'ysuzukilab/Cibersortx', 'cmasotti-Lab/CIBERSORTx_Guide']
records = []
for repository in repositories:
    record = {'repository': repository, 'url': f'https://github.com/{repository}', 'scope': 'Public metadata, default-branch commit and recursive file tree; no third-party code executed'}
    try:
        with urlopen(Request(f'https://api.github.com/repos/{repository}', headers={'User-Agent': 'IBD-scientific-method-audit'}), timeout=20) as response:
            metadata = json.load(response)
        branch = metadata['default_branch']
        with urlopen(Request(f'https://api.github.com/repos/{repository}/commits/{branch}', headers={'User-Agent': 'IBD-scientific-method-audit'}), timeout=20) as response:
            commit = json.load(response)
        sha = commit['sha']
        with urlopen(Request(f'https://api.github.com/repos/{repository}/git/trees/{sha}?recursive=1', headers={'User-Agent': 'IBD-scientific-method-audit'}), timeout=20) as response:
            tree = json.load(response)
        files = [item['path'] for item in tree['tree'] if item['type'] == 'blob']
        record.update({'status': 'read_success', 'default_branch': branch, 'commit': sha, 'tree_truncated': tree.get('truncated', False), 'license_metadata': metadata.get('license'), 'file_count': len(files), 'relevant_paths': [name for name in files if any(term in name.lower() for term in ['cibersort', 'hires', 'signature', 'fraction', 'readme', 'license'])]})
        (outdir / f"{repository.replace('/', '_')}_file_tree.json").write_text(json.dumps({'commit': sha, 'files': files}, ensure_ascii=False, indent=2), encoding='utf-8')
    except (HTTPError, URLError, TimeoutError) as error:
        record.update({'status': 'read_failed_no_retry', 'error': str(error)})
    records.append(record)
    print(repository, record['status'], record.get('file_count', 'unavailable'))
(outdir / 'repository_inventory.json').write_text(json.dumps(records, ensure_ascii=False, indent=2), encoding='utf-8')

# 为同一SSH隧道提供临时文件传输，仅绑定服务器127.0.0.1，文件限制在独立传输目录。
# Python运行：python 本脚本 ROOT PORT；任务结束显式停止服务，不新增SSH/SCP连接。
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import unquote, urlparse
import sys

root = Path(sys.argv[1]).resolve()
root.mkdir(parents=True, exist_ok=True)
port = int(sys.argv[2])

# HTTP库要求处理器类；它只做传输，不参与研究计算或改变原始文件。
class TransferHandler(SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(root), **kwargs)

    def do_PUT(self):
        relative = unquote(urlparse(self.path).path).lstrip('/')
        target = (root / relative).resolve()
        if root not in target.parents or target.exists():
            self.send_error(409, 'Outside root or file already exists')
            return
        target.parent.mkdir(parents=True, exist_ok=True)
        remaining = int(self.headers['Content-Length'])
        with target.open('xb') as stream:
            while remaining:
                block = self.rfile.read(min(1024 * 1024, remaining))
                if not block:
                    raise IOError('Incomplete transfer')
                stream.write(block)
                remaining -= len(block)
        self.send_response(201)
        self.end_headers()

ThreadingHTTPServer(('127.0.0.1', port), TransferHandler).serve_forever()

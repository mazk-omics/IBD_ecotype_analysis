# Xiyou 单连接运行规则

按用户要求，采用一个持续 SSH 会话；后续命令与文件传输复用这一会话。
禁止每条命令重新 SSH、额外 SCP/SFTP、失败自动重连。断开后停止并报告。

本机配置位于 C:/Users/马泽凯/.ssh/config；备份 config.before_single_session_20261006_172037。
仅 Xiyou 块添加指定密钥、publickey、BatchMode、单次连接尝试、严格主机校验及保活。
这些参数只限制单个 SSH 进程；保活本身不能避免多个客户端创建新连接。
启动器 Open-XiyouPersistentSession.ps1 用互斥锁避免它自身被同时启动两次；不控制 VS Code 或其他客户端。

Windows 原生 OpenSSH 不支持 Unix ControlMaster 连接共享，因此使用同一 exec session ID 持续写入。
小脚本经该会话传输；较大结果可通过同一 SSH 的本地转发下载，仍仅一条 SSH 连接。
不修改服务器安全策略，不保证连接频率不会受到服务商其他规则限制。

2026-10-06：用户先撤回自主登录授权，随后明确授权继续抽样并回传文件；本次按单连接执行。

参考：
- https://github.com/PowerShell/Win32-OpenSSH/wiki/Project-Scope
- https://man.openbsd.org/ssh_config

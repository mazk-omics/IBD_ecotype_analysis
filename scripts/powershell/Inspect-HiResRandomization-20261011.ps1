# 本机静态核对三个科学数值函数；复用既有工具和只读项目，不连接服务器。
# 导出列随机化、取列及 Spearman 返回类型的解释，不修改二进制或读取认证代码。
$taskProject = Split-Path (Split-Path $PSScriptRoot -Parent) -Parent
$taskTools = Join-Path (Split-Path $taskProject -Parent) '.tools\cibersortx_static'
$taskAudit = Join-Path $taskProject 'CIBERSORTx_server_plan_20261011'

# Step1. 固定已有 Ghidra/JDK，限制本机内存与线程
$env:JAVA_HOME = Join-Path $taskTools 'jdk-25.0.4.1+1'
$env:GHIDRA_HEADLESS_MAXMEM = '2G'
$env:GHIDRA_HEADLESS_JAVA_OPTIONS = "-Dapplication.settingsdir=$taskTools\settings -Dapplication.cachedir=$taskTools\cache -Dapplication.tempdir=$taskTools\temp"
$taskHeadless = Join-Path $taskTools 'ghidra_12.1.4_PUBLIC\support\analyzeHeadless.bat'

# Step2. 禁止保存类型修改，科学伪代码另存到本机排除目录
& $taskHeadless (Join-Path $taskAudit 'decompiler_workspace') CibersortxNumericAudit -process CIBERSORTxHiRes -noanalysis -readOnly -scriptPath (Join-Path $taskProject 'scripts\java') -postScript InspectHiResRandomization.java $taskAudit -max-cpu 2 *> (Join-Path $taskAudit 'ghidra_hires_randomization.log')
if ($LASTEXITCODE -ne 0) { throw 'Numerical inspection failed; inspect local log.' }
Write-Output 'Read-only numerical randomization audit completed.'

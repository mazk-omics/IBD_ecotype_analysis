# 本机给 HiRes 分析项目临时补充参数类型，检查 double 与整数寄存器分配；不运行供应商 ELF。
# 工具沿用 .tools/cibersortx_static；另存类型证据，-readOnly 不保存 Ghidra 项目改动。
$taskProject = Split-Path (Split-Path $PSScriptRoot -Parent) -Parent
$taskRoot = Split-Path $taskProject -Parent
$taskTools = Join-Path $taskRoot '.tools\cibersortx_static'
$taskAudit = Join-Path $taskProject 'CIBERSORTx_server_plan_20261011'

# Step1. Java 与缓存沿用已审计的本机工具版本
$env:JAVA_HOME = Join-Path $taskTools 'jdk-25.0.4.1+1'
$env:GHIDRA_HEADLESS_MAXMEM = '2G'
$env:GHIDRA_HEADLESS_JAVA_OPTIONS = "-Dapplication.settingsdir=$taskTools\settings -Dapplication.cachedir=$taskTools\cache -Dapplication.tempdir=$taskTools\temp"
$taskHeadless = Join-Path $taskTools 'ghidra_12.1.4_PUBLIC\support\analyzeHeadless.bat'
$taskWorkspace = Join-Path $taskAudit 'decompiler_workspace'
$taskScripts = Join-Path $taskProject 'scripts\java'

# Step2. 复用单个项目；禁用分析与保存，避免覆盖最初的自动反编译证据
& $taskHeadless $taskWorkspace CibersortxNumericAudit -process CIBERSORTxHiRes -noanalysis -readOnly -scriptPath $taskScripts -postScript InspectCibersortxHiResParameters.java $taskAudit -max-cpu 2 *> (Join-Path $taskAudit 'ghidra_hires_typed_parameters.log')
if ($LASTEXITCODE -ne 0) { throw 'Typed parameter audit failed; inspect local log.' }
Write-Output 'HiRes read-only typed parameter inspection completed.'

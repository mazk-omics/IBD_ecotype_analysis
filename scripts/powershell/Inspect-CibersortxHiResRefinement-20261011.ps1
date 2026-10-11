# 在本机复用固定 HiRes 静态分析项目，补齐 vector 返回值后核对论文的收缩公式。
# 使用既有 Ghidra/JDK；新输出另存，-readOnly 丢弃临时类型修改，不运行供应商程序。
$taskProject = Split-Path (Split-Path $PSScriptRoot -Parent) -Parent
$taskRoot = Split-Path $taskProject -Parent
$taskTools = Join-Path $taskRoot '.tools\cibersortx_static'
$taskAudit = Join-Path $taskProject 'CIBERSORTx_server_plan_20261011'

# Step1. 沿用已审计的工具版本与缓存
$env:JAVA_HOME = Join-Path $taskTools 'jdk-25.0.4.1+1'
$env:GHIDRA_HEADLESS_MAXMEM = '2G'
$env:GHIDRA_HEADLESS_JAVA_OPTIONS = "-Dapplication.settingsdir=$taskTools\settings -Dapplication.cachedir=$taskTools\cache -Dapplication.tempdir=$taskTools\temp"
$taskHeadless = Join-Path $taskTools 'ghidra_12.1.4_PUBLIC\support\analyzeHeadless.bat'
$taskWorkspace = Join-Path $taskAudit 'decompiler_workspace'
$taskScripts = Join-Path $taskProject 'scripts\java'

# Step2. 禁用重新分析与保存，保留早期自动导出证据
& $taskHeadless $taskWorkspace CibersortxNumericAudit -process CIBERSORTxHiRes -noanalysis -readOnly -scriptPath $taskScripts -postScript InspectCibersortxHiResRefinement.java $taskAudit -max-cpu 2 *> (Join-Path $taskAudit 'ghidra_hires_refinement_parameters.log')
if ($LASTEXITCODE -ne 0) { throw 'Refinement audit failed; inspect local log.' }
Write-Output 'HiRes read-only refinement inspection completed.'

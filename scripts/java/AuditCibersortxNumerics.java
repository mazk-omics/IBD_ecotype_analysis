// 本机 Ghidra headless 脚本：导出已命名数值函数的静态关系和反编译伪代码。
// 只分析 Linux ELF 文件，不运行程序；认证相关函数或含认证标记的结果不导出。
// @category IBD.SourceAudit
import ghidra.app.script.GhidraScript;
import ghidra.app.decompiler.DecompInterface;
import ghidra.app.decompiler.DecompileResults;
import ghidra.program.model.listing.Function;
import ghidra.program.model.symbol.Symbol;
import java.nio.file.Path;
import java.nio.file.Files;
import java.nio.charset.StandardCharsets;
import java.io.BufferedWriter;
import java.util.Set;
import java.util.regex.Pattern;

public class AuditCibersortxNumerics extends GhidraScript {
    public void run() throws Exception {
        // Step1. 输出目录由命令行明确传入；类与回调是 Ghidra 接口的必要结构
        String label = getScriptArgs().length > 1 ? getScriptArgs()[1] : currentProgram.getName();
        Path out = Path.of(getScriptArgs()[0], label);
        Files.createDirectories(out);
        Pattern numerical = Pattern.compile("CIBERSORTxHiRes::|CIBERSORTxFractions::|GETable::|LinearRegression::|CellStatTable::|SVRWrapper::|Table::|^svm_|^Solver::|^Solver_NU::|^Kernel::|^Cache::|^SVC_Q::|^ONE_CLASS_Q::|^SVR_Q::|^solve_|^sigmoid_train$|^sigmoid_predict$|^multiclass_probability$|^runHiRes$|^CIBERSORTxHiRes$|^CIBERSORTxFractions$|^CoreAlg$|^doPerm$|^buildSignatureMatrix$|^makeRefandClassFiles$|^FDR$|^adjustMixtures$|^adjustSigMatrix$|^runrmbatchSmode$|^runrmbatchBmode$|^doBootstrap$|^nnls$|^nnlsFTWrapper$|^subsetNNLS$|^lh_|^transfEigenMatrix$|^transfSTDMatrix$|^res$|^resid$|^subsetResid$|^cholInvWT$|^findBestRSSThres$|^colSums$|^rowSums$|^pmax$|^pmin$|^log10$|^stats::|^runR::doBootstrap$|^runR::qtmpAddRnorm$|^runR::randTable$|^runR::CIBERSORTxGEP$|^runR::makeHiResHeatmap$|^runR::coreSVR$|^runR::combat$|^runR::normalizeQuantiles$|^runR::adjustPseudo$");
        Pattern sensitive = Pattern.compile("(?i)auth|token|password|credential|username|login|license|https?://|http::|[a-z0-9._%+-]+@[a-z0-9.-]+\\.[a-z]{2,}");
        // 抽样、排序、矩阵算术与 NNLS 的低层依赖需同时保留，才可追踪算法；不选择系统命令封装。
        Pattern helpers = Pattern.compile("^utility::(sequence|sample<|order<|countEqThresh<|countNA<|isUnique<|subsetVector<|split<|uniqueUnion<|which_min<|which<|whichAboveThresh<|whichBelowThresh<|copySamplingState$|resetSamplingState$|mCol$|writeVec<|writeMatrix<)|^operator[-+*/]$|^(abs|diag|pow|sqrt|transpose|vecSum)<|^(diff_|g1_|h12_)$");
        numerical = Pattern.compile("(?:" + numerical.pattern() + ")|(?:" + helpers.pattern() + ")");
        numerical = Pattern.compile(numerical.pattern() + "|^(utility::|GETable::)?rmTokensAndIG$");
        int inspected = 0;
        int exported = 0;
        int withheld = 0;
        DecompInterface decompiler = new DecompInterface();
        decompiler.openProgram(currentProgram);

        // Step2. 先保存算法函数目录，再逐个恢复伪代码；保留结果状态而不掩盖失败
        try (BufferedWriter table = Files.newBufferedWriter(out.resolve("numerical_functions.tsv"), StandardCharsets.UTF_8);
             BufferedWriter graph = Files.newBufferedWriter(out.resolve("numerical_call_edges.tsv"), StandardCharsets.UTF_8);
             BufferedWriter dependencies = Files.newBufferedWriter(out.resolve("numerical_dependency_edges.tsv"), StandardCharsets.UTF_8);
             BufferedWriter screening = Files.newBufferedWriter(out.resolve("withheld_marker_categories.tsv"), StandardCharsets.UTF_8);
             BufferedWriter candidates = Files.newBufferedWriter(out.resolve("application_function_candidates.tsv"), StandardCharsets.UTF_8)) {
            table.write("Address\tName\tBodyBytes\tStatus\tPseudocodeFile\n");
            graph.write("Caller\tCallee\n");
            dependencies.write("Caller\tCallee\tAddress\tSelected\n");
            candidates.write("Address\tName\tBodyBytes\n");
            screening.write("Address\tName\tMarker\tCount\tBoostTokenizerPresent\n");
            var functions = currentProgram.getFunctionManager().getFunctions(true);
            while (functions.hasNext() && !monitor.isCancelled()) {
                Function function = functions.next();
                String name = function.getSymbol().getName(true);
                // 经独立分类确认的唯一标识符误判；仅筛选输入替换此名称，导出代码保持原文。
                String screenName = name.replace("rmTokensAndIG", "rmGeneLabelsAndIG");
                boolean applicationCandidate = !name.contains("std::") && !name.contains("boost::") && !name.contains("Eigen::") && !name.contains("__gnu") && !name.startsWith("FUN_") && !name.startsWith("_") && !function.isExternal() && !function.isThunk() && !sensitive.matcher(screenName).find();
                if (applicationCandidate) candidates.write(function.getEntryPoint().toString() + "\t" + name.replace('\t', ' ') + "\t" + function.getBody().getNumAddresses() + "\n");
                if (!numerical.matcher(name).find() || sensitive.matcher(screenName).find() || function.isExternal() || function.isThunk()) continue;
                inspected++;
                String address = function.getEntryPoint().toString();
                Set<Function> called = function.getCalledFunctions(monitor);
                for (Function target : called) {
                    String targetName = target.getSymbol().getName(true);
                    String screenTargetName = targetName.replace("rmTokensAndIG", "rmGeneLabelsAndIG");
                    if (!sensitive.matcher(screenTargetName).find()) dependencies.write(name + "\t" + targetName.replace('\t', ' ') + "\t" + target.getEntryPoint().toString() + "\t" + numerical.matcher(targetName).find() + "\n");
                    if (numerical.matcher(targetName).find() && !sensitive.matcher(screenTargetName).find()) graph.write(name + "\t" + targetName + "\n");
                }
                DecompileResults result = decompiler.decompileFunction(function, 60, monitor);
                String status = "decompilation_failed";
                String filename = "";
                if (result.decompileCompleted() && result.getDecompiledFunction() != null) {
                    String code = result.getDecompiledFunction().getC();
                    if (sensitive.matcher(code.replace("rmTokensAndIG", "rmGeneLabelsAndIG")).find()) {
                        status = "withheld_sensitive_markers";
                        withheld++;
                        // 只记录命中类别与数量，避免输出字符串上下文或潜在凭据。
                        for (String marker : new String[]{"auth", "token", "password", "credential", "username", "login", "license", "https?://", "http::", "[a-z0-9._%+-]+@[a-z0-9.-]+\\.[a-z]{2,}"}) {
                            var matches = Pattern.compile(marker, Pattern.CASE_INSENSITIVE).matcher(code);
                            int count = 0;
                            while (matches.find()) count++;
                            if (count > 0) screening.write(address + "\t" + name + "\t" + marker + "\t" + count + "\t" + code.contains("boost::tokenizer") + "\n");
                        }
                    } else {
                        filename = "numeric_" + address + ".c";
                        Files.writeString(out.resolve(filename), "/* Ghidra pseudocode, not original source; not compiled or validated. */\n" + code, StandardCharsets.UTF_8);
                        status = "pseudocode_exported";
                        exported++;
                    }
                }
                table.write(address + "\t" + name.replace('\t', ' ') + "\t" + function.getBody().getNumAddresses() + "\t" + status + "\t" + filename + "\n");
                table.flush();
                graph.flush();
            }
        } finally {
            decompiler.dispose();
        }
        // Step3. 仅打印统计，不将原始字符串或整个反编译结果送入日志
        println("Numerical audit: " + currentProgram.getName() + " inspected=" + inspected + " exported=" + exported + " withheld=" + withheld);
    }
}

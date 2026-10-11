// 仅核对 HiRes 列随机化和 Spearman 数值返回类型，不执行供应商程序或读取认证逻辑。
// 在既有 Ghidra -readOnly 项目运行，临时 ABI 类型解释另存；不是原始源码。
// @category IBD.SourceAudit
import ghidra.app.script.GhidraScript;
import ghidra.app.decompiler.DecompInterface;
import ghidra.program.model.listing.Function;
import ghidra.program.model.listing.ParameterImpl;
import ghidra.program.model.data.*;
import ghidra.program.model.symbol.SourceType;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.charset.StandardCharsets;
import java.util.regex.Pattern;

public class InspectHiResRandomization extends GhidraScript {
    public void run() throws Exception {
        // Step1. 白名单只含三个科学数值函数，输出不得覆盖早期导出
        Path out = Path.of(getScriptArgs()[0], "decompiled_numeric", "CIBERSORTxHiRes_randomization_typed");
        Files.createDirectories(out);
        DataType pointer = new PointerDataType(currentProgram.getDataTypeManager());
        DataType integer = IntegerDataType.dataType;
        Function randomizer = getFunctionAt(toAddr("001fc958"));
        Function getCol = getFunctionAt(toAddr("001fc71a"));
        Function tail = getFunctionAt(toAddr("00129f15"));
        if (randomizer == null || getCol == null || tail == null || !randomizer.getName().equals("randomizeInCols") || !getCol.getName().equals("getCol") || !tail.getName().equals("pspearman")) throw new IllegalStateException("Numerical whitelist mismatch.");

        // Step2. 成员 this 与返回 vector 的隐藏指针分开；RNG 位于第五个整数寄存器
        randomizer.replaceParameters(Function.FunctionUpdateType.DYNAMIC_STORAGE_ALL_PARAMS, true, SourceType.USER_DEFINED, new ParameterImpl("table", pointer, currentProgram), new ParameterImpl("option_a", pointer, currentProgram), new ParameterImpl("option_b", pointer, currentProgram), new ParameterImpl("option_c", integer, currentProgram), new ParameterImpl("rng", pointer, currentProgram));
        randomizer.setReturnType(VoidDataType.dataType, SourceType.USER_DEFINED);
        getCol.replaceParameters(Function.FunctionUpdateType.DYNAMIC_STORAGE_ALL_PARAMS, true, SourceType.USER_DEFINED, new ParameterImpl("output_vector", pointer, currentProgram), new ParameterImpl("table", pointer, currentProgram), new ParameterImpl("column_index", integer, currentProgram));
        tail.setReturnType(DoubleDataType.dataType, SourceType.USER_DEFINED);

        // Step3. 类型只在分析器内修改；输出前继续沿用凭据筛选
        Pattern sensitive = Pattern.compile("(?i)auth|token|password|credential|username|login|license|https?://|http::|[a-z0-9._%+-]+@[a-z0-9.-]+\\.[a-z]{2,}");
        DecompInterface decompiler = new DecompInterface();
        decompiler.openProgram(currentProgram);
        try {
            for (Function function : new Function[]{randomizer, getCol, tail}) {
                var result = decompiler.decompileFunction(function, 60, monitor);
                if (!result.decompileCompleted() || result.getDecompiledFunction() == null) throw new IllegalStateException("Numerical decompilation failed.");
                String code = result.getDecompiledFunction().getC();
                if (sensitive.matcher(code).find()) throw new IllegalStateException("Output withheld by screening.");
                Files.writeString(out.resolve("typed_" + function.getEntryPoint() + ".c"), "/* Analyst-supplied numerical ABI interpretation; not original source or validated executable. */\n" + code, StandardCharsets.UTF_8);
            }
        } finally {
            decompiler.dispose();
        }
        println("Numerical randomization inspection completed; discard project changes with -readOnly.");
    }
}

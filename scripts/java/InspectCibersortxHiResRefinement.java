// 本机补齐 runHiRes 的参数及 vector 返回值，核对原论文收缩公式的算术次序，不执行供应商程序。
// 在 Ghidra -readOnly 项目中运行；另存带类型假设的伪代码与参数存储表，不覆盖早期导出或二进制。
// @category IBD.SourceAudit
import ghidra.app.script.GhidraScript;
import ghidra.app.decompiler.DecompInterface;
import ghidra.program.model.listing.Function;
import ghidra.program.model.listing.Parameter;
import ghidra.program.model.listing.ParameterImpl;
import ghidra.program.model.data.DataType;
import ghidra.program.model.data.IntegerDataType;
import ghidra.program.model.data.DoubleDataType;
import ghidra.program.model.data.PointerDataType;
import ghidra.program.model.symbol.SourceType;
import java.nio.file.Path;
import java.nio.file.Files;
import java.nio.charset.StandardCharsets;
import java.util.regex.Pattern;

public class InspectCibersortxHiResRefinement extends GhidraScript {
    public void run() throws Exception {
        // Step1. 仅选择已导出的数值函数；C++ 非平凡对象按 invisible reference 传入
        Function function = getFunctionAt(toAddr("0010cc9b"));
        if (function == null || !function.getName().equals("runHiRes")) throw new IllegalStateException("Pinned numerical function not found.");
        Path out = Path.of(getScriptArgs()[0], "decompiled_numeric", "CIBERSORTxHiRes_refinement_typed");
        Files.createDirectories(out);
        DataType integer = IntegerDataType.dataType;
        DataType floating = DoubleDataType.dataType;
        DataType pointer = new PointerDataType(currentProgram.getDataTypeManager());
        String[] names = {"worker_id", "fractions", "mixtures", "group_adjusted_pvalues", "cell_type_count", "sample_count", "window_size", "partition_size_hint", "bootstrap_primary", "bootstrap_secondary", "sample_metric", "gene_partitions", "output_prefix", "result_cube"};
        DataType[] types = {integer, pointer, pointer, pointer, integer, integer, floating, floating, integer, integer, pointer, pointer, pointer, pointer};
        Parameter[] parameters = new Parameter[names.length];
        for (int i = 0; i < names.length; i++) parameters[i] = new ParameterImpl(names[i], types[i], currentProgram);

        // Step2. 根据 Linux x86-64 的整数/浮点独立寄存器分配，保留返回类型未定这一限制
        function.replaceParameters(Function.FunctionUpdateType.DYNAMIC_STORAGE_ALL_PARAMS, true, SourceType.USER_DEFINED, parameters);
        StringBuilder table = new StringBuilder("Ordinal\tName\tType\tStorage\n");
        for (Parameter parameter : function.getParameters()) table.append(parameter.getOrdinal() + "\t" + parameter.getName() + "\t" + parameter.getDataType().getName() + "\t" + parameter.getVariableStorage() + "\n");
        Files.writeString(out.resolve("runHiRes_parameter_storage.tsv"), table.toString(), StandardCharsets.UTF_8);

        // 返回非平凡 vector 时 RDI 是隐式结果指针；补齐两个矩阵算术函数丢失的 RDX 操数
        Function matrixAddition = getFunctionAt(toAddr("00138cf8"));
        Function matrixDivision = getFunctionAt(toAddr("0013903a"));
        if (matrixAddition == null || matrixDivision == null || !matrixAddition.getName().equals("operator+") || !matrixDivision.getName().equals("operator/")) throw new IllegalStateException("Pinned matrix operations not found.");
        matrixAddition.replaceParameters(Function.FunctionUpdateType.DYNAMIC_STORAGE_ALL_PARAMS, true, SourceType.USER_DEFINED, new ParameterImpl("output_matrix", pointer, currentProgram), new ParameterImpl("left_matrix", pointer, currentProgram), new ParameterImpl("right_matrix", pointer, currentProgram));
        matrixDivision.replaceParameters(Function.FunctionUpdateType.DYNAMIC_STORAGE_ALL_PARAMS, true, SourceType.USER_DEFINED, new ParameterImpl("output_matrix", pointer, currentProgram), new ParameterImpl("input_matrix", pointer, currentProgram), new ParameterImpl("divisor_reference", pointer, currentProgram));

        // Step3. 补齐 vector 的隐藏结果指针，核对论文收缩公式的实际运算次序
        String[] binaryHelpers = {"00125228", "001253c3", "00134906", "00134a0f", "001372e1", "00137409", "00138036", "001381b2", "001382bd", "001392b0", "00139d0e", "00132ebb"};
        for (String address : binaryHelpers) {
            Function helper = getFunctionAt(toAddr(address));
            if (helper == null) throw new IllegalStateException("Numerical helper missing.");
            helper.replaceParameters(Function.FunctionUpdateType.DYNAMIC_STORAGE_ALL_PARAMS, true, SourceType.USER_DEFINED, new ParameterImpl("output_vector", pointer, currentProgram), new ParameterImpl("left_reference", pointer, currentProgram), new ParameterImpl("right_reference", pointer, currentProgram));
        }
        String[] unaryHelpers = {"00125128", "001371e1", "001fa870", "001fab5c", "001fc33c"};
        for (String address : unaryHelpers) {
            Function helper = getFunctionAt(toAddr(address));
            if (helper == null) throw new IllegalStateException("Numerical reduction missing.");
            helper.replaceParameters(Function.FunctionUpdateType.DYNAMIC_STORAGE_ALL_PARAMS, true, SourceType.USER_DEFINED, new ParameterImpl("output_vector", pointer, currentProgram), new ParameterImpl("input_reference", pointer, currentProgram));
        }

        // Step4. 新输出明确标为分析器类型假设；凭据筛选与原有数值导出相同
        DecompInterface decompiler = new DecompInterface();
        decompiler.openProgram(currentProgram);
        try {
            var result = decompiler.decompileFunction(function, 90, monitor);
            if (!result.decompileCompleted() || result.getDecompiledFunction() == null) throw new IllegalStateException("Typed decompilation failed.");
            String code = result.getDecompiledFunction().getC();
            Pattern sensitive = Pattern.compile("(?i)auth|token|password|credential|username|login|license|https?://|http::|[a-z0-9._%+-]+@[a-z0-9.-]+\\.[a-z]{2,}");
            if (sensitive.matcher(code).find()) throw new IllegalStateException("Output withheld by credential screening.");
            Files.writeString(out.resolve("runHiRes_typed.c"), "/* Static pseudocode with analyst-supplied parameter types; not original source or validated executable. Return and object layouts remain incomplete. */\n" + code, StandardCharsets.UTF_8);
            // 两个小函数另存，供核对隐藏返回指针与实际逐元素运算；不增加无关函数导出
            for (Function helper : new Function[]{matrixAddition, matrixDivision}) {
                var helperResult = decompiler.decompileFunction(helper, 60, monitor);
                if (!helperResult.decompileCompleted() || helperResult.getDecompiledFunction() == null) throw new IllegalStateException("Typed matrix helper failed.");
                String helperCode = helperResult.getDecompiledFunction().getC();
                if (sensitive.matcher(helperCode).find()) throw new IllegalStateException("Matrix output withheld by screening.");
                Files.writeString(out.resolve("typed_" + helper.getEntryPoint() + ".c"), "/* Analyst ABI interpretation of a numerical matrix helper; not original or validated executable source. */\n" + helperCode, StandardCharsets.UTF_8);
            }
        } finally {
            decompiler.dispose();
        }
        println("Typed numerical refinement audit completed; project changes must be discarded with -readOnly.");
    }
}

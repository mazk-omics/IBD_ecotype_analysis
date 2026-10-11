// 只复核已定位 signature 构建函数的筛选命中类型，不输出源码、字符串或凭据上下文。
// 本机 Ghidra headless 运行；读固定 Fractions 静态项目，写一个筛选诊断 TSV。
// @category IBD.SourceAudit
import ghidra.app.script.GhidraScript;
import ghidra.app.decompiler.DecompInterface;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.regex.Pattern;

public class InspectCibersortxScreening extends GhidraScript {
    public void run() throws Exception {
        // Step1. 仅定位一个既有函数；类别判断不生成可运行程序
        var function = currentProgram.getFunctionManager().getFunctionAt(toAddr(0x10dd29));
        if (function == null || !function.getName().equals("buildSignatureMatrix")) throw new IllegalStateException("Fixed binary/function mismatch");
        var decompiler = new DecompInterface();
        decompiler.openProgram(currentProgram);
        try {
            var result = decompiler.decompileFunction(function, 60, monitor);
            if (!result.decompileCompleted()) throw new IllegalStateException("Decompilation incomplete");
            String code = result.getDecompiledFunction().getC();
            // Step2. 先标记注释和字面字符串范围，再判断命中是否属于代码标识符；不保存其内容
            var literals = Pattern.compile("(?s)/\\*.*?\\*/|//[^\\n]*|\"(?:\\\\.|[^\"\\\\])*\"|'(?:\\\\.|[^'\\\\])*'").matcher(code);
            boolean[] literalPositions = new boolean[code.length()];
            while (literals.find()) for (int i = literals.start(); i < literals.end(); i++) literalPositions[i] = true;
            var markers = Pattern.compile("(?i)auth|token|password|credential|username|login|license|https?://|http::").matcher(code);
            String table = "Function\tMarker\tKind\tRecognizedTextSplittingIdentifier\tCodeIdentifier\n";
            while (markers.find()) {
                int left = markers.start();
                int right = markers.end();
                while (left > 0 && Character.isJavaIdentifierPart(code.charAt(left - 1))) left--;
                while (right < code.length() && Character.isJavaIdentifierPart(code.charAt(right))) right++;
                String identifier = code.substring(left, right);
                boolean known = !literalPositions[markers.start()] && java.util.Set.of("tokenize", "tokenizer", "token_compress_on", "token_compress_off", "token_type", "token_iterator").contains(identifier);
                String safeName = !literalPositions[markers.start()] && identifier.matches("[A-Za-z_][A-Za-z_0-9]{0,79}") ? identifier : "not_exported";
                table += "buildSignatureMatrix\t" + markers.group().toLowerCase() + "\t" + (literalPositions[markers.start()] ? "literal_or_comment" : "code_identifier") + "\t" + known + "\t" + safeName + "\n";
            }
            // Step3. 输出类别、布尔值和代码标识符，字符串内容不导出；不能据此自动放宽筛选
            Files.writeString(Path.of(getScriptArgs()[0]), table);
            println("Signature construction screening categories saved; source content not exported.");
        } finally {
            decompiler.dispose();
        }
    }
}

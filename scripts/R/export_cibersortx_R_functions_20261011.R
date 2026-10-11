# 按 R 语法树读取官方应用源码，只导出不涉及认证的函数片段；不执行这些函数。
# 本机 R 4.4.1 运行，输入 vendor_reference/parse_inputs.json，输出函数目录与只读片段。
library(jsonlite)

# Step1. 读取源码文本的临时审计输入，不 source() 供应商程序
args = commandArgs(trailingOnly=TRUE)
inputs = fromJSON(args[1], simplifyVector=FALSE)
outdir = args[2]
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
records = list()
record_index = 0L
parse_status = list()

for (input in inputs) {  # 多镜像、多个源码版本必须分别记录
    parse_text = gsub("\r\n", "\n", input$text, fixed=TRUE)  # parse(text=) 先统一 CRLF；不修改原始文件或算法，行号不变
    parse_error = ""
    parsed = tryCatch(parse(text=parse_text, keep.source=TRUE), error=function(e) { parse_error <<- strsplit(conditionMessage(e), "\n", fixed=TRUE)[[1]][1]; NULL })  # 只保留首行位置，不输出原始源码上下文
    safe_error = gsub("([\"']).*?\\1", "[literal withheld]", parse_error)
    parse_status[[length(parse_status) + 1L]] = list(image=input$image, path=input$path, source_sha256=input$sha256, status=if (is.null(parsed)) "parse_failed" else "parsed", parse_input_line_endings="CRLF_to_LF", error=safe_error)
    if (is.null(parsed)) {
        cat("Parse failed:", input$image, input$path, safe_error, "\n")
        next
    }
    ranges = attr(parsed, "srcref")
    source_lines = strsplit(parse_text, "\n", fixed=TRUE)[[1]]

    # Step2. 用语法树区分函数定义与顶层驱动，保留原始行号方便复核
    for (j in seq_along(parsed)) {
        expression = parsed[[j]]
        if (!is.call(expression) || length(expression) < 3L || !(as.character(expression[[1]]) %in% c("<-", "="))) next
        if (!is.symbol(expression[[2]]) || !is.call(expression[[3]]) || as.character(expression[[3]][[1]]) != "function") next
        name = as.character(expression[[2]])
        identifiers = all.names(expression[[3]], functions=TRUE, unique=TRUE)
        start = as.integer(ranges[[j]][1])
        end = as.integer(ranges[[j]][3])
        fragment = paste(source_lines[start:end], collapse="\n")
        sensitive = grepl("token|username|auth|password|credential|login|license|https?://|[a-z0-9._%+-]+@[a-z0-9.-]+\\.[a-z]{2,}", fragment, ignore.case=TRUE)
        status = if (sensitive) "withheld_sensitive_markers" else "original_function_fragment_exported"
        filename = ""
        if (!sensitive) {
            filename = paste0(input$image, "_", substr(input$sha256, 1, 12), "_", name, ".R")
            writeLines(c("# 原始函数片段，保留原始算法；仅用于静态审阅，未执行或验证，不是可独立运行的软件。", fragment), file.path(outdir, filename), useBytes=TRUE)
        }
        safe_identifiers = identifiers[!grepl("token|username|auth|password|credential|login|license", identifiers, ignore.case=TRUE)]
        record_index = record_index + 1L
        records[[record_index]] = list(image=input$image, path=input$path, source_sha256=input$sha256, function_name=if (grepl("token|username|auth|password|credential|login|license", name, ignore.case=TRUE)) "[withheld]" else name, start_line=start, end_line=end, status=status, file=filename, identifiers=safe_identifiers)
    }
}

# Step3. 保存状态和函数位置；完整原始文本不写入结构报告
write_json(records, file.path(dirname(dirname(outdir)), "R_function_expression_audit.json"), pretty=TRUE, auto_unbox=TRUE)
write_json(parse_status, file.path(dirname(dirname(outdir)), "R_parse_status.json"), pretty=TRUE, auto_unbox=TRUE)
cat("Parsed source versions:", length(inputs), "function definitions:", length(records), "exported:", sum(vapply(records, function(x) x$status == "original_function_fragment_exported", logical(1))), "\n")

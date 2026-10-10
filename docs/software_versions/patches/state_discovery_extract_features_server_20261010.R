suppressPackageStartupMessages({
library(data.table)
source("lib/misc.R")
})

args = c("Codefacs", "All_genes", "T.cell.Group2") 
args = commandArgs(T)
dataset = args[1]
fractions = args[2]
cell_type = args[3]
scaling_column = args[4] 
filter_genes = args[5]

if(is.na(filter_genes))
{
	filter_genes = F
}

specific_genes_dir = file.path("../EcoTyper", dataset, fractions, "Analysis", "Cell_type_specific_genes")
input_dir = file.path("../CIBERSORTx/hires", dataset, fractions) 

raw_input = read.delim(file.path(input_dir, paste0(cell_type, ".txt")), row.names = 1)
if(nrow(raw_input) < 50)
{
	cat(paste0("Warning: Only ", nrow(raw_input), " genes are available for '", cell_type, "'. At least 50 genes are required for cell state discovery. Skipping cell state discovery for this cell type!\n"))
	quit()
}

output_dir = file.path("../EcoTyper", dataset, fractions, "Cell_States", "discovery", cell_type) 
dir.create(output_dir, recursive = T, showWarning = F)

if(min(raw_input, na.rm = T) >= 1)
{
	log_data = log2(raw_input)	
}else{
	if(min(raw_input, na.rm = T) >= 0)
	{
		log_data = log2(raw_input + 1)
	}else{
		warning(paste0("Negative values found in the expression matrix of '", cell_type, "'. Not performing log-transformation...\n"))
		log_data = raw_input
	}
}

if(filter_genes)
{
	cat(paste0("Filtering '", cell_type, "' profiles for cell type specific genes...\n")) 
	if(!file.exists(file.path(specific_genes_dir, paste0(cell_type, "_cell_type_specific_genes.txt"))))
	{
		stop(paste0("List of cell type specific genes not found for cell type '", cell_type, "'. Please make sure that step 1 finished sucessfully!\n"))
	}
	gene_info = read.delim(file.path(specific_genes_dir, paste0(cell_type, "_cell_type_specific_genes.txt")))
	gene_info = gene_info[gene_info$CellType == cell_type,]
	log_data = log_data[rownames(log_data) %in% gene_info$Gene,]
	if(nrow(log_data) < 50)
	{
		cat(paste0("Warning: Only ", nrow(log_data), " genes are available for '", cell_type, "'. At least 50 genes are required for cell state discovery. Skipping cell state discovery for this cell type!\n"))
		quit()
	}
}else{
	cat(paste0("Not filtering '", cell_type, "' profiles for cell type specific genes...\n")) 
} 

write.table(log_data, file.path(output_dir, "expression_full_matrix_log2.txt"), sep = "\t", row.names = T)

data = get_variable_genes(log_data)
write.table(data, file.path(output_dir, "expression_top_genes_log2.txt"), sep = "\t", row.names = T)

clinical = read_clinical(
    colnames(data),
    dataset = dataset,
    dataset_type = "discovery"
)

# Verify that annotation has been reordered exactly to match
# the expression matrix sample order.
if (!identical(
    as.character(clinical$ID),
    as.character(colnames(data))
)) {
    stop(
        "Annotation sample order does not exactly match expression matrix."
    )
}

if (
    is.null(scaling_column) ||
    is.na(scaling_column) ||
    scaling_column == "NULL"
) {

    by = NULL

    cat(
        paste0(
            "Performing univariance normalization of CIBERSORTx ",
            "high resolution genes imputed for '",
            cell_type,
            "', across all samples...\n"
        )
    )

} else {

    if (!(scaling_column %in% colnames(clinical))) {

        stop(
            paste0(
                "Scaling column '",
                scaling_column,
                "' not found in annotation. Available columns: ",
                paste(
                    colnames(clinical),
                    collapse = ", "
                )
            )
        )
    }

    # IMPORTANT:
    # use the actual sample-level annotation vector,
    # not the literal column name string.
    by = as.character(
        clinical[[scaling_column]]
    )

    if (length(by) != ncol(data)) {

        stop(
            paste0(
                "Scaling vector length mismatch: ",
                length(by),
                " annotation values for ",
                ncol(data),
                " expression samples."
            )
        )
    }

    if (anyNA(by) || any(by == "")) {

        stop(
            paste0(
                "Missing values detected in scaling column '",
                scaling_column,
                "'."
            )
        )
    }

    group_sizes = table(by)

    if (any(group_sizes < 2)) {

        stop(
            paste0(
                "At least one scaling group contains fewer than 2 samples: ",
                paste(
                    names(group_sizes),
                    group_sizes,
                    sep = "=",
                    collapse = ", "
                )
            )
        )
    }

    cat(
        paste0(
            "Performing univariance normalization of CIBERSORTx ",
            "high resolution genes imputed for '",
            cell_type,
            "', within groups defined by column '",
            scaling_column,
            "': ",
            paste(
                names(group_sizes),
                group_sizes,
                sep = "=",
                collapse = ", "
            ),
            "\n"
        )
    )
}

scaled_data = scale_data(
    data,
    by = by
)

write.table(
    scaled_data,
    file.path(
        output_dir,
        "expression_top_genes_scaled.txt"
    ),
    sep = "\t",
    row.names = T
)

scaled_data = scale_data(
    log_data,
    by = by
)

write.table(
    scaled_data,
    file.path(
        output_dir,
        "expression_full_matrix_scaled.txt"
    ),
    sep = "\t",
    row.names = T
)

#' Predict Molecular Subtypes
#'
#' This function predicts molecular subtypes based on the input expression matrix.
#' It supports both RNA-seq and PCR data types.
#'
#' @param expression_matrix A matrix where rows are genes and columns are samples.
#' @param type A string specifying the type of input data: "RNA-seq" (default) or "PCR".
#' @return A dataframe with two columns: `SampleName` and `PredictedSubtype`.
#' @export
predict_molecular_subtype <- function(expression_matrix, type = "RNA-seq") {
  # Load goal_gene_pairs from package data
  data("goal_gene_pairs", package = "pandora")

  if (type == "RNA-seq") {
    data("unique_genes_df", package = "pandora")

    # Remove ; only pass expression_matrix and goal_gene_pairs to pair_genes
    paired_genesT <- pair_genes(expression_matrix)

    predicted_subtypes <- predict_subtype(paired_genesT)
    expression_matrix2<-expression_matrix[,-1]

    result <- data.frame(
      SampleName = colnames(expression_matrix2),
      PredictedSubtype = predicted_subtypes
    )

  } else if (type == "PCR") {
    expression_matrix2 <- expression_matrix[!duplicated(expression_matrix[, 1]), ]

    # Assign the first column as row names (gene names)
    rownames(expression_matrix2) <- expression_matrix2[, 1]

    # Remove the first column (now it's used as row names)
    expression_matrix3 <- expression_matrix2[, -1]
    expression_matrix_t <- t(expression_matrix3)

    # Load the PCR-specific model
    data("rf_model_0907", package = "pandora")

    gene_pairs <- c("AP1S3.CD79B", "MCM2.PARD6B", "AP1S3.PI16",
                    "CD93.LGALS2", "ADAMTS2.LGALS2", "AP1S3.JSRP1",
                    "C5AR1.LGALS2", "ADAMTS2.AGR3", "ADAM9.GATM", "BCAS1.FPR1")

    result_df <- data.frame(matrix(NA, nrow = length(gene_pairs), ncol = nrow(expression_matrix_t)))
    rownames(result_df) <- gene_pairs
    colnames(result_df) <- rownames(expression_matrix_t)

    for (pair in gene_pairs) {
      genes <- strsplit(pair, "\\.")[[1]]
      geneA <- genes[1]
      geneB <- genes[2]

      if (geneA %in% colnames(expression_matrix_t) & geneB %in% colnames(expression_matrix_t)) {
        result_df[pair, ] <- ifelse(expression_matrix_t[, geneA] > expression_matrix_t[, geneB], 0, 1)
      } else {
        warning(paste("One or both genes", geneA, "and", geneB, "are not present in the expression matrix"))
      }
    }

    cz_pair_qpcr <- t(result_df)
    predicted_subtypes <- predict(rf_model_0907, newdata = cz_pair_qpcr)

    result <- data.frame(
      SampleName = rownames(expression_matrix_t),
      PredictedSubtype = predicted_subtypes
    )

  } else {
    stop("Invalid input type. Please choose either 'RNA-seq' or 'PCR'.")
  }

  return(result)
}

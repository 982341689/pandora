#' Pair Genes for Subtype Prediction
#'
#' This function pairs genes from the expression matrix based on the goal gene pairs.
#'
#' @param expression_matrix A matrix where rows are genes and columns are samples.
#' @param goal_gene_pairs A vector of gene pairs in the format "geneA.geneB".
#' @return A dataframe of paired genes with values 0 or 1, where 1 means geneA > geneB.
#' @export
pair_genes <- function(expression_matrix) {


  # Remove rows with duplicate values in the first column (which would be used as row names)
  expression_matrix <- expression_matrix[!duplicated(expression_matrix[, 1]), ]

  # Assign the first column as row names
  row.names(expression_matrix) <- expression_matrix[, 1]

  # Remove the first column (since it is now the row names)
  expression_matrix <- expression_matrix[, -1]
  unique_genes<-unique_genes_df$genes



  # Step 2: Create a mapping for specific gene names
  gene_mapping <- c("AC027601.1" = "TMEM105", "AC022498.2" = "FLJ42393")

  # Step 3: Update the rownames in expression_matrix based on the mapping
  row.names(expression_matrix) <- sapply(row.names(expression_matrix), function(x) ifelse(x %in% names(gene_mapping), gene_mapping[x], x))

  # Step 1: Find genes that are in unique_genes_df but not in expression_matrix rownames
  missing_genes <- setdiff(unique_genes_df$genes, row.names(expression_matrix))

  # Step 2: Check if there are missing genes and print a warning or success message
  if (length(missing_genes) > 0) {
    warning(paste("The following genes are not present in expression_matrix:", paste(missing_genes, collapse = ", ")))
  } else {
    message("All genes in unique_genes_df are present in expression_matrix.")
  }



  matching_genes <- intersect(row.names(expression_matrix), unique_genes_df$genes)
  expression_matrix_filtered <- expression_matrix[rownames(expression_matrix) %in% matching_genes, ]



  #pairing_240907



  # Step 2: Split each gene pair into geneA and geneB
  split_pairs <- strsplit(goal_gene_pairs, "\\|")

  # Step 3: Create an empty list to store results (to avoid row name issues)
  cz_pair_results <- list()

  # Get the sample names (column names) from expression_matrix_filtered
  sample_names <- colnames(expression_matrix_filtered)

  # Step 4: Loop through the goal gene pairs and calculate relative expression
  for (pair in split_pairs) {
    geneA <- pair[[1]]  # First gene in the pair
    geneB <- pair[[2]]  # Second gene in the pair

    # Ensure both genes exist in expression_matrix_filtered
    if (geneA %in% rownames(expression_matrix_filtered) & geneB %in% rownames(expression_matrix_filtered)) {
      # Calculate relative expression: 1 if geneA > geneB, otherwise 0
      comparison <- ifelse(expression_matrix_filtered[geneA, ] > expression_matrix_filtered[geneB, ], 1, 0)

      # Extract the vector from the comparison (since it's a single row data frame)
      comparison_vector <- as.vector(comparison)

      # Debugging: Print the structure of comparison result
      print(paste("Comparison result for", geneA, "and", geneB, ":"))
      print(comparison_vector)
      print(paste("Length of comparison:", length(comparison_vector)))
      print(paste("Number of samples (sample_names):", length(sample_names)))

      # Ensure the length of the comparison result matches the number of samples
      if (length(comparison_vector) == length(sample_names)) {
        # Create a dataframe for this gene pair
        pair_df <- data.frame(t(comparison_vector))

        # Assign column names to match the sample names
        colnames(pair_df) <- sample_names

        # Assign the row name in the format "geneA|geneB"
        rowname <- paste0(geneA, "|", geneB)
        cz_pair_results[[rowname]] <- pair_df  # Add result to list with appropriate name
      } else {
        warning(paste("Mismatch in the number of samples for gene pair:", geneA, "|", geneB))
      }
    } else {
      warning(paste("One or both genes", geneA, "and", geneB, "are not present in expression_matrix_filtered"))
    }
  }

  # Step 5: Combine the results into a single dataframe
  cz_pair_final <- do.call(rbind, cz_pair_results)

  paired_genesT<-t(cz_pair_final)
  colnames(paired_genesT) <- gsub("[-|]", ".", colnames(paired_genesT))

  # Return the transposed and cleaned dataframe
  return(paired_genesT)
}


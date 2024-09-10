#' Predict Molecular Subtypes Using a Pre-trained Model
#'
#' This function uses a pre-trained model to predict molecular subtypes based on the paired gene data.
#'
#' @param paired_genesT A transposed dataframe of paired genes with values 0 or 1.
#' @param model A pre-trained machine learning model. If NULL, the model will be loaded from the package.
#' @return A vector of predicted subtypes.
#' @export
predict_subtype <- function(paired_genesT, model = NULL) {
  # Load the model from the package if not provided
  if (is.null(model)) {
    data("large_model", package = "pandora")
  }

  # Predict subtypes
  predicted_subtypes <- predict(large_model, newdata = paired_genesT)

  return(predicted_subtypes)
}

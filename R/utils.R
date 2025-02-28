#' Check if an object is a phyloseq object
#'
#' @param x Object to check
#'
#' @return Boolean indicating whether x is a phyloseq object
#' @keywords internal
is_phyloseq <- function(x) {
  # Check if phyloseq package is available
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package 'phyloseq' needed for this function. Please install it.")
  }
  
  # Use the class to determine if it's a phyloseq object
  # (is.phyloseq isn't exported from the phyloseq namespace)
  return(inherits(x, "phyloseq"))
}

#' Check if a data frame can be used as metadata
#'
#' Checks if the data frame has row names that can be used as sample identifiers
#'
#' @param x Data frame to check
#'
#' @return Boolean indicating whether x can be used as metadata
#' @keywords internal
is_valid_metadata <- function(x) {
  # Must be a data frame
  if (!is.data.frame(x)) {
    return(FALSE)
  }
  
  # Must have row names
  if (is.null(rownames(x)) || length(rownames(x)) == 0) {
    return(FALSE)
  }
  
  # Row names should be unique identifiers
  if (length(unique(rownames(x))) != length(rownames(x))) {
    return(FALSE)
  }
  
  # Row names shouldn't be just the default 1, 2, 3, etc.
  if (identical(rownames(x), as.character(1:nrow(x)))) {
    return(FALSE)
  }
  
  return(TRUE)
}

#' Standardize sample names using common transformations
#'
#' @param names Character vector of sample names
#' @param transformations List of transformations to apply (default set includes common patterns)
#'
#' @return Character vector of standardized sample names
#' @keywords internal
standardize_names <- function(names, transformations = NULL) {
  if (is.null(transformations)) {
    transformations <- list(
      # Replace underscores with dots
      function(x) gsub("_", ".", x),
      # Replace hyphens with dots
      function(x) gsub("-", ".", x),
      # Remove leading/trailing whitespace
      function(x) trimws(x),
      # Convert to lowercase
      function(x) tolower(x),
      # Remove special characters
      function(x) gsub("[^[:alnum:]\\._-]", "", x)
    )
  }
  
  # Apply each transformation
  result <- names
  for (transform in transformations) {
    result <- transform(result)
  }
  
  return(result)
}
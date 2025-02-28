#' Render a phyloseq metadata check report
#'
#' @param physeq A phyloseq object
#' @param metadata A data frame containing metadata with samples as row names
#' @param output_file Output file name (default: "phylometacheck_report.html")
#' @param output_dir Directory to save the report
#' @param ... Additional parameters passed to rmarkdown::render
#'
#' @return Path to the generated HTML report
#' @export
#'
#' @examples
#' \dontrun{
#' data(GlobalPatterns, package = "phyloseq")
#' meta <- data.frame(sample_data(GlobalPatterns))
#' # Modify some names to create mismatches
#' rownames(meta)[1:3] <- paste0(rownames(meta)[1:3], "_modified")
#' render_report(GlobalPatterns, meta)
#' }
render_report <- function(physeq, metadata, 
                          output_file = "phylometacheck_report.html",
                          output_dir = getwd(), ...) {
  
  # Make sure required packages are installed
  required_pkgs <- c("rmarkdown", "phyloseq", "dplyr", "ggplot2", 
                    "DT", "stringdist", "shiny", "htmltools")
  
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  
  if (length(missing_pkgs) > 0) {
    stop("The following required packages are missing: ", 
         paste(missing_pkgs, collapse = ", "), 
         ". Please install them using install.packages().")
  }
  
  # Get the path to the template
  template_path <- get_template_path()
  if (!file.exists(template_path)) {
    stop("Could not find the phyloseq-meta-check template.")
  }
  
  # Create a temporary file for the report
  temp_file <- tempfile(fileext = ".Rmd")
  on.exit(unlink(temp_file), add = TRUE)
  
  # Copy the template to the temporary file
  file.copy(template_path, temp_file)
  
  # Set up the output file path
  output_path <- file.path(output_dir, output_file)
  
  # Render the report
  rmarkdown::render(
    input = temp_file,
    output_file = output_path,
    params = list(
      physeq_obj = physeq,
      metadata_obj = metadata
    ),
    ...
  )
  
  # Return the path to the generated report
  return(output_path)
}
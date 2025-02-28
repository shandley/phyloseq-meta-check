#' Get the path to the phyloseq-meta-check template
#'
#' This function returns the path to the R Markdown template for phyloseq-meta-check
#'
#' @return Path to the template
#' @keywords internal
get_template_path <- function() {
  system.file("rmarkdown/templates/phyloseq-meta-check/skeleton/skeleton.Rmd",
              package = "phylometacheck")
}

#' Use phyloseq-meta-check template to create a new R Markdown document
#'
#' @param file Path where the new file will be created
#'
#' @return Path to the created file (invisibly)
#' @export
#'
#' @examples
#' \dontrun{
#' use_template("my_report.Rmd")
#' }
use_template <- function(file) {
  template_path <- get_template_path()
  if (!file.exists(template_path)) {
    stop("Could not find the phyloseq-meta-check template.")
  }
  file.copy(template_path, file)
  invisible(file)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_name_functions.R
\name{suggest_name_fixes}
\alias{suggest_name_fixes}
\title{Suggest fixes for mismatched sample names based on fuzzy matching}
\usage{
suggest_name_fixes(physeq, metadata, max_distance = 0.2, method = "jw")
}
\arguments{
\item{physeq}{A phyloseq object}

\item{metadata}{A data frame containing metadata with samples as row names}

\item{max_distance}{Maximum string distance to consider for matches}

\item{method}{String distance method to use (default: "jw" for Jaro-Winkler)}
}
\value{
A data frame with suggested matches and confidence scores
}
\description{
Suggest fixes for mismatched sample names based on fuzzy matching
}
\examples{
\dontrun{
data(GlobalPatterns, package = "phyloseq")
meta <- data.frame(sample_data(GlobalPatterns))
# Modify some names to create mismatches
rownames(meta)[1:3] <- paste0(rownames(meta)[1:3], "_modified")
suggestions <- suggest_name_fixes(GlobalPatterns, meta)
}
}

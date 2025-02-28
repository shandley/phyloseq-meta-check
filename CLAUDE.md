# Commands & Guidelines for phyloseq-meta-check

## R Commands
- **Run Rmarkdown**: `Rscript -e "rmarkdown::render('your_file.Rmd', output_format='html_document')"`
- **Lint R code**: `Rscript -e "lintr::lint('your_file.R')"`
- **Test package**: `Rscript -e "devtools::test()"`
- **Run single test**: `Rscript -e "devtools::test_file('tests/testthat/test_file.R')"`

## Code Style Guidelines
- **Formatting**: Follow tidyverse style guide (2-space indentation, 80-char line length)
- **Naming**: snake_case for variables/functions, TitleCase for classes
- **R packages**: Load with library(), prefer tidyverse ecosystem (dplyr, ggplot2, etc.)
- **Use phyloseq syntax**: Consistent with [phyloseq package](https://joey711.github.io/phyloseq/)
- **Documentation**: Roxygen2 for functions (params, return values, examples)
- **Error handling**: Use tryCatch() for robust error handling in R
- **Vectorize operations**: Avoid for-loops when possible in R

## Project Overview
This tool helps check and fix sample names to ensure they are compatible for merging into a phyloseq object.
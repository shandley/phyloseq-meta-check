# phylometacheck

A R package for checking and fixing sample names to make phyloseq objects and metadata mergeable.

## Overview

`phylometacheck` helps microbiome researchers identify and resolve sample name mismatches between phyloseq objects and metadata files. Using fuzzy matching algorithms, the package can detect similar names and suggest corrections, making it easier to merge data from different sources.

## Features

* Identify exact matches and mismatches between sample names
* Apply fuzzy matching to find similar names with configurable matching criteria
* Generate interactive HTML reports showing matching status and suggestions
* Visualize sample name relationships
* Apply corrections to either phyloseq objects or metadata
* Export correction code for reproducibility

## Installation

```r
# Install from GitHub
devtools::install_github("ScottHandley/phylometacheck")
```

## Usage

### Basic sample name checks

```r
library(phyloseq)
library(phylometacheck)
library(dplyr)

# Load your data
data(GlobalPatterns, package = "phyloseq")
meta <- data.frame(sample_data(GlobalPatterns))
# Create some mismatches for demonstration
rownames(meta)[1:3] <- paste0(rownames(meta)[1:3], "_modified")

# Find mismatches
mismatches <- find_sample_mismatches(GlobalPatterns, meta)
print(mismatches)

# Get suggestions for fixes
suggestions <- suggest_name_fixes(GlobalPatterns, meta)
print(suggestions)

# Apply fixes to metadata
fixed_metadata <- apply_name_fixes(
  suggestions[1:2,],  # Select which fixes to apply
  metadata = meta,
  target = "metadata"
)
```

### Generate an interactive report

```r
# Create a full interactive HTML report
report_path <- render_report(
  physeq = GlobalPatterns,
  metadata = meta,
  output_file = "sample_name_check.html"
)

# Open the report
browseURL(report_path)
```

### Use as an R Markdown template

```r
# Create a new report from template
rmarkdown::draft("my_report.Rmd", 
                 template = "phyloseq-meta-check", 
                 package = "phylometacheck")
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

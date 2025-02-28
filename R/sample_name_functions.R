#' Find mismatches between phyloseq sample names and metadata
#'
#' @param physeq A phyloseq object
#' @param metadata A data frame containing metadata with samples as row names
#' @param exact_match If TRUE, only find exact matches/mismatches; if FALSE (default),
#'   include fuzzy matching results
#'
#' @return A list containing exact matches, samples only in phyloseq,
#'   samples only in metadata, and if exact_match=FALSE, potential fuzzy matches
#' @export
#'
#' @examples
#' \dontrun{
#' data(GlobalPatterns, package = "phyloseq")
#' meta <- data.frame(sample_data(GlobalPatterns))
#' # Modify some names to create mismatches
#' rownames(meta)[1:3] <- paste0(rownames(meta)[1:3], "_modified")
#' mismatches <- find_sample_mismatches(GlobalPatterns, meta)
#' }
find_sample_mismatches <- function(physeq, metadata, exact_match = FALSE) {
  # Check inputs
  if (!is_phyloseq(physeq)) {
    stop("physeq must be a phyloseq object")
  }
  
  if (!is_valid_metadata(metadata)) {
    stop("metadata must be a data frame with row names")
  }
  
  # Get sample names
  physeq_samples <- phyloseq::sample_names(physeq)
  meta_samples <- rownames(metadata)
  
  # Find exact matches and mismatches
  exact_matches <- intersect(physeq_samples, meta_samples)
  physeq_only <- setdiff(physeq_samples, meta_samples)
  meta_only <- setdiff(meta_samples, physeq_samples)
  
  # Return results
  result <- list(
    exact_matches = exact_matches,
    physeq_only = physeq_only,
    meta_only = meta_only
  )
  
  # If not in exact_match mode and there are mismatches, add fuzzy matches
  if (!exact_match && (length(physeq_only) > 0 || length(meta_only) > 0)) {
    result$fuzzy_matches <- suggest_name_fixes(physeq, metadata)
  }
  
  return(result)
}

#' Suggest fixes for mismatched sample names based on fuzzy matching
#'
#' @param physeq A phyloseq object
#' @param metadata A data frame containing metadata with samples as row names
#' @param max_distance Maximum string distance to consider for matches
#' @param method String distance method to use (default: "jw" for Jaro-Winkler)
#'
#' @return A data frame with suggested matches and confidence scores
#' @export
#'
#' @examples
#' \dontrun{
#' data(GlobalPatterns, package = "phyloseq")
#' meta <- data.frame(sample_data(GlobalPatterns))
#' # Modify some names to create mismatches
#' rownames(meta)[1:3] <- paste0(rownames(meta)[1:3], "_modified")
#' suggestions <- suggest_name_fixes(GlobalPatterns, meta)
#' }
suggest_name_fixes <- function(physeq, metadata, max_distance = 0.2, 
                               method = "jw") {
  
  # Check inputs
  if (!is_phyloseq(physeq)) {
    stop("physeq must be a phyloseq object")
  }
  
  if (!is_valid_metadata(metadata)) {
    stop("metadata must be a data frame with row names")
  }
  
  # Check if stringdist package is available
  if (!requireNamespace("stringdist", quietly = TRUE)) {
    stop("Package 'stringdist' needed for fuzzy matching. Please install it.")
  }
  
  # Get sample names
  physeq_samples <- phyloseq::sample_names(physeq)
  meta_samples <- rownames(metadata)
  
  # Find mismatches
  physeq_only <- setdiff(physeq_samples, meta_samples)
  meta_only <- setdiff(meta_samples, physeq_samples)
  
  # If no mismatches, return empty data frame
  if (length(physeq_only) == 0 || length(meta_only) == 0) {
    return(data.frame(
      physeq_name = character(0),
      meta_name = character(0),
      match_quality = numeric(0),
      patterns = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # Create standardized versions of sample names for better matching
  physeq_std <- standardize_names(physeq_only)
  meta_std <- standardize_names(meta_only)
  
  # Compute distance matrix between all physeq and metadata samples
  dist_matrix <- stringdist::stringdistmatrix(physeq_std, meta_std, method = method)
  rownames(dist_matrix) <- physeq_only
  colnames(dist_matrix) <- meta_only
  
  # Convert distance to similarity (1 = identical, 0 = completely different)
  max_dist <- max(dist_matrix)
  if (max_dist > 0) {
    sim_matrix <- 1 - (dist_matrix / max_dist)
  } else {
    sim_matrix <- matrix(1, nrow = nrow(dist_matrix), ncol = ncol(dist_matrix))
  }
  
  # Find best matches for each physeq sample
  suggestions <- data.frame(
    physeq_name = character(0),
    meta_name = character(0),
    match_quality = numeric(0),
    patterns = character(0),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:length(physeq_only)) {
    # Find best match
    best_idx <- which.max(sim_matrix[i, ])
    best_score <- sim_matrix[i, best_idx]
    
    # Only include if it meets the minimum quality
    if (best_score >= (1 - max_distance)) {
      # Determine the pattern of difference
      physeq_name <- physeq_only[i]
      meta_name <- meta_only[best_idx]
      
      # Detect common patterns
      patterns <- c()
      
      # Case difference
      if (tolower(physeq_name) == tolower(meta_name)) {
        patterns <- c(patterns, "case")
      }
      
      # Special character replacements
      if (gsub("_", ".", physeq_name) == meta_name || 
          physeq_name == gsub("_", ".", meta_name)) {
        patterns <- c(patterns, "underscore/dot")
      }
      
      # Prefix/suffix
      if (startsWith(meta_name, physeq_name) || startsWith(physeq_name, meta_name)) {
        patterns <- c(patterns, "prefix")
      }
      
      if (endsWith(meta_name, physeq_name) || endsWith(physeq_name, meta_name)) {
        patterns <- c(patterns, "suffix")
      }
      
      # If no specific patterns detected, it's a more general string similarity
      if (length(patterns) == 0) {
        patterns <- c("similar")
      }
      
      # Add to suggestions
      suggestions <- rbind(suggestions, data.frame(
        physeq_name = physeq_name,
        meta_name = meta_name,
        match_quality = best_score,
        patterns = paste(patterns, collapse = ", "),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Sort by match quality
  suggestions <- suggestions[order(-suggestions$match_quality), ]
  
  return(suggestions)
}

#' Apply suggested name fixes to either phyloseq object or metadata
#'
#' @param fixes A data frame of fixes, typically from suggest_name_fixes()
#' @param physeq A phyloseq object to rename samples in (if target="physeq")
#' @param metadata A data frame to rename samples in (if target="metadata")
#' @param target Whether to rename samples in "physeq" or "metadata"
#'
#' @return Updated phyloseq object or metadata data frame with corrected names
#' @export
#'
#' @examples
#' \dontrun{
#' data(GlobalPatterns, package = "phyloseq")
#' meta <- data.frame(sample_data(GlobalPatterns))
#' # Modify some names to create mismatches
#' rownames(meta)[1:3] <- paste0(rownames(meta)[1:3], "_modified")
#' suggestions <- suggest_name_fixes(GlobalPatterns, meta)
#' # Select which suggestions to apply
#' fixes_to_apply <- suggestions[1:2, ]
#' # Apply to metadata
#' fixed_meta <- apply_name_fixes(fixes_to_apply, metadata = meta, target = "metadata")
#' }
apply_name_fixes <- function(fixes, physeq = NULL, metadata = NULL, 
                             target = c("physeq", "metadata")) {
  
  # Check inputs
  if (!inherits(fixes, "data.frame")) {
    stop("fixes must be a data frame")
  }
  
  if (!all(c("physeq_name", "meta_name") %in% colnames(fixes))) {
    stop("fixes must have columns 'physeq_name' and 'meta_name'")
  }
  
  target <- match.arg(target)
  
  if (target == "physeq") {
    if (is.null(physeq)) {
      stop("physeq must be provided when target='physeq'")
    }
    
    if (!is_phyloseq(physeq)) {
      stop("physeq must be a phyloseq object")
    }
    
    # Create a mapping of old to new names
    old_names <- phyloseq::sample_names(physeq)
    new_names <- old_names
    
    # Update names based on fixes
    for (i in 1:nrow(fixes)) {
      idx <- which(old_names == fixes$physeq_name[i])
      if (length(idx) > 0) {
        new_names[idx] <- fixes$meta_name[i]
      }
    }
    
    # Update the sample names in the phyloseq object
    phyloseq::sample_names(physeq) <- new_names
    
    return(physeq)
    
  } else if (target == "metadata") {
    if (is.null(metadata)) {
      stop("metadata must be provided when target='metadata'")
    }
    
    if (!is_valid_metadata(metadata)) {
      stop("metadata must be a data frame with row names")
    }
    
    # Create a new data frame with updated row names
    fixed_metadata <- metadata
    
    # Create a mapping from old to new names
    name_map <- setNames(fixes$physeq_name, fixes$meta_name)
    
    # Update row names
    meta_names <- rownames(metadata)
    for (old_name in names(name_map)) {
      idx <- which(meta_names == old_name)
      if (length(idx) > 0) {
        meta_names[idx] <- name_map[old_name]
      }
    }
    
    # Apply the new row names
    rownames(fixed_metadata) <- meta_names
    
    return(fixed_metadata)
  }
}

#' Plot sample name matches and mismatches
#'
#' @param matches Results from find_sample_mismatches()
#' @param suggestions Optional results from suggest_name_fixes()
#'
#' @return A ggplot2 visualization
#' @export
#'
#' @examples
#' \dontrun{
#' data(GlobalPatterns, package = "phyloseq")
#' meta <- data.frame(sample_data(GlobalPatterns))
#' # Modify some names to create mismatches
#' rownames(meta)[1:3] <- paste0(rownames(meta)[1:3], "_modified")
#' mismatches <- find_sample_mismatches(GlobalPatterns, meta)
#' suggestions <- suggest_name_fixes(GlobalPatterns, meta)
#' plot_sample_matches(mismatches, suggestions)
#' }
plot_sample_matches <- function(matches, suggestions = NULL) {
  # Check if ggplot2 is available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for plotting. Please install it.")
  }
  
  # Create data frames for visualization
  exact_count <- length(matches$exact_matches)
  physeq_only_count <- length(matches$physeq_only)
  meta_only_count <- length(matches$meta_only)
  
  # Create data frame for matches/mismatches counts
  df <- data.frame(
    category = c("Exact Matches", "Phyloseq Only", "Metadata Only"),
    count = c(exact_count, physeq_only_count, meta_only_count)
  )
  
  # Basic bar plot for counts
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = category, y = count, fill = category)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(title = "Sample Name Matches and Mismatches",
         x = "Category", y = "Count") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
  
  # If suggestions are provided, add more detailed visualization
  if (!is.null(suggestions) && nrow(suggestions) > 0) {
    # Add match quality visualization
    suggestions$physeq_name <- factor(suggestions$physeq_name, 
                                   levels = suggestions$physeq_name)
    
    p2 <- ggplot2::ggplot(suggestions, 
                 ggplot2::aes(x = physeq_name, y = meta_name, fill = match_quality)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient(low = "#f8d7da", high = "#d4edda") +
      ggplot2::labs(title = "Suggested Name Matches",
           x = "Phyloseq Sample", y = "Metadata Sample",
           fill = "Match Quality") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    
    # Combine plots with patchwork if available
    if (requireNamespace("patchwork", quietly = TRUE)) {
      p <- p1 / p2
      return(p)
    } else {
      # Just return the match quality plot if patchwork is not available
      return(p2)
    }
  } else {
    # Return the basic counts plot if no suggestions
    return(p1)
  }
}

#' Create interactive report for phyloseq and metadata sample name checking
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
#' create_phylometacheck_report(GlobalPatterns, meta)
#' }
create_phylometacheck_report <- function(physeq, metadata, 
                                         output_file = "phylometacheck_report.html",
                                         output_dir = getwd(), ...) {
  # Call the render_report function from report_functions.R
  return(render_report(physeq, metadata, output_file, output_dir, ...))
}
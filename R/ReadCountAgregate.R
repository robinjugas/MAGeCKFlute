#' Aggregate read counts and collapse into gene clusters
#'
#' Go over sgRNA counts and filter & mark the unsatisfying ones
#'
#' @docType methods
#' @name ReadCountAgregate
#' @rdname ReadCountAgregate
#' @aliases ReadCountAgregate
#'
#' @param normalised_count table of normalized counts
#' @param Squareview_table Squareview table
#'
#' @return DF
#'
#' @author Matej Jasik
#' @import openxlsx
#' @export

ReadCountAgregate <- function(normalised_count, Squareview_table){
  
  # requireNamespace("openxlsx", quietly = TRUE) || stop("need openxlsx package")
  message(Sys.time(), "  Collapsing sgRNA counts into gene cluster")
  
  ##############################################################################
  # Read data from the file
  df <- read.table(normalised_count, header = TRUE, sep = "\t")
  # Filter out rows containing the string "Non-Target-Control"
  df <- df[!grepl("Non-Targeting-Control", df$Gene),]
  
  ##############################################################################
  # Function to round numeric values to zero decimal places
  round_values <- function(df) {
    # Identify numeric columns excluding "Gene"
    numeric_cols <- sapply(df, is.numeric) & names(df) != "Gene"
    
    # Round numeric columns to zero decimal places
    df[numeric_cols] <- lapply(df[numeric_cols], function(x) round(x, 0))
    
    return(df)
  }
  
  # Function to aggregate values for each gene across columns
  aggregate_values <- function(df) {
    # Initialize an empty data frame to store aggregated results
    result <- data.frame(Gene = character(), stringsAsFactors = FALSE)
    
    # Get unique gene names
    genes <- unique(df$Gene)
    
    # Loop through each gene
    for (gene in genes) {
      # Subset the data for the current gene
      gene_data <- df[df$Gene == gene, ]
      
      # Create a row for the current gene
      gene_row <- data.frame(Gene = gene, stringsAsFactors = FALSE)
      
      # Aggregate values for each subsequent column
      for (i in 3:ncol(df)) { # Start from the second column assuming Gene is the first column
        # Concatenate values for the current column
        condition_values <- paste(gene_data[[i]], collapse = "|")
        # Assign to the gene_row using original column name
        gene_row[[names(df)[i]]] <- condition_values
      }
      
      # Append gene_row to result data frame
      result <- rbind(result, gene_row)
    }
    
    return(result)
  }
  
  ##############################################################################
  # Round numeric values in df
  df <- round_values(df)
  
  # Get aggregated result
  aggregated_df <- aggregate_values(df)
  
  # merge with Squareview_table
  Squareview_data <- merge(Squareview_table, aggregated_df, by = "Gene")

  ##############################################################################
  # Write aggregated data to Excel file
  # write.xlsx(aggregated_df, file = "aggregated_data_round.xlsx")
  # # Read the original aggregated data from Excel
  # original_data <- read.xlsx("aggregated_data_round.xlsx")
  # # Read the data from the Read_file
  # new_data <- read.table(normalised_count, header = TRUE, sep = "\t")
  # # Filter out rows containing "Non-Targeting-Control"
  # filtered_data <- new_data[grepl("Non-Targeting-Control", new_data$Gene), ]
  # # Round numeric values in filtered_data
  # filtered_data <- round_values(filtered_data)
  # # Remove the second column (assuming it's named "Condition")
  # filtered_data <- filtered_data[, -2]
  # # Ensure that filtered_data has the same columns as original_data
  # if (ncol(filtered_data) != ncol(original_data)) {
  #   stop("Number of columns in filtered_data does not match original_data.")
  # }
  # # Match column names of filtered_data with original_data
  # colnames(filtered_data) <- colnames(original_data)
  # # Append the filtered rows to the original data
  # updated_data <- rbind(original_data, filtered_data)
  # # Write the updated data back to Excel
  # write.xlsx(updated_data, file = "aggregated_data_round.xlsx")
  message(Sys.time(), "  Done")
  
  return(Squareview_data)
}
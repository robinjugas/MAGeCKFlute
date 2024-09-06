#' Filter bad sgRNA guides
#'
#' Go over sgRNA counts and filter & mark the unsatisfying ones
#'
#' @docType methods
#' @name filterBADsgRNA
#' @rdname filterBADsgRNA
#' @aliases filterBADsgRNA
#'
#' @param fileCountNormalized A order ranked numeric vector with geneid as names
#' @param Squareview_table "Entrez", "Ensembl", or "Symbol"
#'
#'
#' @return DF
#'
#' @author Robin Jugas, Matej Jasik
#'
#' @export


filterBADsgRNA <- function(fileCountNormalized, Squareview_table)
{
  
  ################################################################################
  # CALCULATE GUIDE NUMBERS
  
  CountNormalized <- fread(fileCountNormalized, header = TRUE)
  T0col <- grep("_T0", names(CountNormalized)) #T0 sloupec
  otherCol <- setdiff((3:ncol(CountNormalized)),T0col) #dalsi COUNT sloupce
  GenesUnique <- unique(CountNormalized$Gene) #nazvy genu unikatne
  T0threshold <- 0 #prah pro T0 filtr
  otherColthreshold <- 0 #prah pro dalsi COUNT sloupce
  # Define your thresholds
  T0threshold_sum <- 500
  otherColthreshold_sum <- 500
  
  TO_issue_genes <- c() #nazvy genu s problemem v T0
  low_count_issue_genes <- c() #nazvy genu s problemem v countech obecne
  
  
  finalDF <- data.frame(Gene = GenesUnique, T0_warning = character(length(GenesUnique)), warning = character(length(GenesUnique)))
  
  test <- function(x){
    c(mean = mean(x),sd = sd(x))
  }
  
  stats <- CountNormalized[ ,unlist(lapply(.SD, test)), .SDcols = 3:ncol(CountNormalized)]
  
  # 
  # which(finalDF$gene=="TAF2")
  # tempgene <- 16830
  
  for (tempgene in 1:nrow(finalDF)) {
    temp <- CountNormalized[Gene == finalDF[tempgene, "Gene"]]
    
    # Instance summary check
    if (sum(temp[, ..T0col] <= T0threshold, na.rm = TRUE) > 0) {
      finalDF[tempgene, "T0_warning_instances"] <- paste0(sum(temp[, ..T0col] <= T0threshold, na.rm = TRUE), "_T0_guides_ZERO")
    } else {
      finalDF[tempgene, "T0_warning_instances"] <- ""
    }
    
    if (sum(temp[, ..otherCol] <= otherColthreshold, na.rm = TRUE) > 0) {
      finalDF[tempgene, "warning_instances"] <- paste0(sum(temp[, ..otherCol] <= otherColthreshold, na.rm = TRUE), "_guides_ZERO")
    } else {
      finalDF[tempgene, "warning_instances"] <- ""
    }
    
    # Total sum warning check
    sum_T0 <- sum(temp[, ..T0col], na.rm = TRUE)
    if (sum_T0 < T0threshold_sum) {
      finalDF[tempgene, "T0_warning_sum"] <- paste0("Sum of T0 values (", sum_T0, ") is below ", T0threshold_sum)
    } else {
      finalDF[tempgene, "T0_warning_sum"] <- ""
    }
    
    sum_other <- sum(temp[, ..otherCol], na.rm = TRUE)
    if (sum_other < otherColthreshold_sum) {
      finalDF[tempgene, "warning_sum"] <- paste0("Sum of otherCol values (", sum_other, ") is below ", otherColthreshold_sum)
    } else {
      finalDF[tempgene, "warning_sum"] <- ""
    }
  }
  
  setDT(finalDF)
  ################################################################################
  ## merge with Squareview_table
  
  # Squareview_data <- fread(fileSquareview_data, header = TRUE)
  
  Squareview_data <- merge(Squareview_table, finalDF, by = "Gene")
  return(Squareview_data)
  
  
  # fwrite(Squareview_data, file = gsub(".txt","_guide_warnings.txt",fileSquareview_data), sep = "\t", col.names = TRUE, quote = FALSE, append = FALSE)
  # EXCEL
  # wb <- openxlsx::createWorkbook()
  # openxlsx::addWorksheet(wb, sheetName = "List1")
  # openxlsx::writeData(wb, sheet = "List1", Squareview_data)
  # openxlsx::saveWorkbook(wb, gsub(".txt","_guide_warnings.xlsx",fileSquareview_data), overwrite = TRUE, returnValue = FALSE)
  
  
}

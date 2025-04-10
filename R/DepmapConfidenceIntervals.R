#' Load processed Depmap data
#'
#' @docType methods
#' @name DepmapConfidenceIntervals
#' @rdname DepmapConfidenceIntervals
#'
#' @author Robin Jugas
#'
#' @return data.frame with confidence intervals added
#' @import depmap
#' @export

#test
dd <- gdata_cc
##############################################################################

DepmapConfidenceIntervals <- function(dd){

  # load depmap
  crispr <- depmap::depmap_crispr()
  crispr <- matrix(crispr$dependency, nrow = length(unique(crispr$depmap_id)),
                   ncol = length(unique(crispr$gene_name)),
                   dimnames = list(rownames = unique(crispr$depmap_id),
                                   colnames = unique(crispr$gene_name)))
  crispr <- as.data.frame(t(crispr))
  crisprX <- data.frame(gene=rownames(crispr), minConf=rep(0,nrow(crispr)), maxConf=rep(0,nrow(crispr)),
                        estimate=rep(0,nrow(crispr)), p.value=rep(0,nrow(crispr))
  )
  
  
  # do t.test
  i=1
  for(i in 1:nrow(crispr)){
    x <- crispr[i,]
    res <- t.test(x)
    
    crisprX[i,"minConf"] <- res$conf.int[1]
    crisprX[i,"maxConf"] <- res$conf.int[2]
    crisprX[i,"estimate"] <- res$estimate
    crisprX[i,"p.value"] <- res$p.value
    
  }
  
  dd <- merge(dd, crisprX , by="gene")
  
  
  
  #
  summary(crispr[which(rownames(crispr)=="AAMP"),])
  
  
  
  return(dd)
}

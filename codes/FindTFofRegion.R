###
#   File name : FindTFofRegion.R
#   Author    : Hyunjin Kim
#   Date      : May 1, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Find TF candidates that may be bound in regions of interest
#
#   Instruction
#               1. Source("FindTFofRegion.R")
#               2. Run the function "findTF" - specify the input file (diffbind result) and output path
#               3. TF candidates and their info will be saved as tab-separated file in the output path
#
#   Example
#               > source("The_directory_of_FindTFofRegion.R/FindTFofRegion.R")
#               > findTF(regionPath="./data/ATAC_Seq/Mouse_AA_cd4_cd8/01.AACD4vsCtrlCD4_atfc1.5p0.05_1180_peaks_annotation.csv",
#                        outputPath="./results/CD4_TFs.txt")
###

findTF <- function(regionPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Christiano/Eddy/2018/ATAC_Seq/Mouse_AA_cd4_cd8/data/01.AACD4vsCtrlCD4_atfc1.5p0.05_1180_peaks_annotation.csv",
                   outputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Christiano/Eddy/2018/ATAC_Seq/Mouse_AA_cd4_cd8/TF_footprint/CD4_TFs.txt") {
  
  ### load library
  if(!require(biomaRt)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("biomaRt")
    library(biomaRt)
  }
  
  ### set biomart parameters
  ensembl <- useMart("ENSEMBL_MART_FUNCGEN", dataset = "mmusculus_motif_feature")
  
  ### get TF info
  tfInfo <- getBM(attributes = listAttributes(ensembl)$name, mart = ensembl)
  
  ### load DB regions
  dbRegions <- read.csv(file = regionPath, stringsAsFactors = FALSE, check.names = FALSE)
  
  ### clean & order tfInfo
  tfInfo <- tfInfo[,-c(5,6,7,10)]
  tfInfo <- tfInfo[order(tfInfo$chromosome_name, tfInfo$chromosome_start, tfInfo$chromosome_end),]
  
  ### find TFs of the regions
  cnt <- 1
  result <- data.frame(matrix(NA, 1, (ncol(dbRegions)+ncol(tfInfo))))
  colnames(result) <- c(colnames(dbRegions), colnames(tfInfo))
  for(i in 1:nrow(dbRegions)) {
    chr <- substr(dbRegions$seqnames[i], 4, nchar(dbRegions$seqnames[i]))
    start <- dbRegions$start[i]
    end <- dbRegions$end[i]
    
    tfInfoTemp <- tfInfo[which(tfInfo$chromosome_name == chr),]
    
    idx <- intersect(which(tfInfoTemp$chromosome_start >= start), which(tfInfoTemp$chromosome_end <= end))
    
    if(length(idx) > 0) {
      result <- rbind(result, cbind(dbRegions[i,], tfInfoTemp[idx,], row.names=NULL))
      cnt <- cnt+1
    }
    
  }
  result <- result[-1,]
  
  ### rename some of the colnames of the result
  colnames(result)[23:28] <- c("TF_binding_matrix_id", "TFBS_chr", "TFBS_start", "TFBS_end", "TF_name", "TFBS_so_accession")
  
  ### print some info
  print(paste("Total number of TFs in the result:", nrow(result)))
  print(paste("Total number of the unique DB regions in the result:",
              length(which(duplicated(result[,c(1,2,3)]) == FALSE)), "/", nrow(dbRegions)))
  
  ### save the result
  write.table(result, file = outputPath, sep = "\t", row.names = FALSE)
  
}

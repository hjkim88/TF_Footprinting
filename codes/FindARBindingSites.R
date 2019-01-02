###
#   File name : FindARBindingSites.R
#   Author    : Hyunjin Kim
#   Date      : May 31, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Find binding sites of AR from region of interest
#
#   Instruction
#               1. Source("FindARBindingSites.R")
#               2. Run the function "findARBS" - specify the input file (fasta), MA_ID, and the output file path
#               3. AR binding site candidates and their info will be saved as tab-separated file in the output path
#
#   Example
#               > source("The_directory_of_FindARBindingSites.R/FindARBindingSites.R")
#               > findARBS(fastaPath="./data/Another_Project/CXCL8.fasta",
#                          maID="MA0007.2",
#                          outputPath="./results/AR_MA0007_2.txt")
###

findARBS <- function(fastaPath="./data/Another_Project/CXCL8.fasta",
                     maID="MA0007.2",
                     outputPath="./results/AR_MA0007_2.txt") {
  
  ### load library
  if(!require(TFBSTools)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("TFBSTools")
    library(TFBSTools)
  }
  if(!require(JASPAR2018)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("JASPAR2018")
    library(JASPAR2018)
  }
  if(!require(Biostrings)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("Biostrings")
    library(Biostrings)
  }
  
  
  ### load fasta
  target_sequence <- readDNAStringSet(filepath = fastaPath, format = "fasta")[[1]]
  
  
  ### load PWM
  pfm <- getMatrixByID(JASPAR2018, ID=maID)
  pwm <- toPWM(pfm, type = "prob")
  
  
  ### search TFBS
  set.seed(1234)
  siteset <- searchSeq(pwm, target_sequence, seqname="target", strand="*", min.score="80%")
  
  
  ### make an empty result object
  result <- data.frame(matrix(NA, 1, 9))
  colnames(result) <- c("TF_Name", "TFBS_Start", "TFBS_End",
                        "Motif_ID", "Sequence", "Score", "Sequence_Strand",
                        "PValue_Exact", "PValue_Sampling")
  
  
  ### iteratively add candidates to the result object
  if(length(siteset) > 0) {
    for(k in 1:length(siteset)) {
      result <- rbind(result,
                      cbind(TF_Name=siteset@pattern@name,
                            TFBS_Start=start(siteset)[k], TFBS_End=end(siteset)[k],
                            Motif_ID=siteset@pattern@ID, Sequence=as.character(views(siteset))[k],
                            Score=score(siteset)[k], Sequence_Strand=strand(siteset)[k],
                            PValue_Exact=pvalues(siteset, type="TFMPvalue")[k],
                            PValue_Sampling=pvalues(siteset, type="sampling")[k]))
    }
  } else {
    writeLines("No matched sequences")
  }
  
  
  ### remove the first NA row of the result object
  result <- result[-1,]
  
  
  ### save the result
  write.table(result, file = outputPath, sep = "\t", row.names = FALSE)
  
}



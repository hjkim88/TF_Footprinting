###
#   File name : FindTFofRegion2.R
#   Author    : Hyunjin Kim
#   Date      : May 10, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Find TF candidates that may be bound in regions of interest
#
#   Instruction
#               1. Source("FindTFofRegion2.R")
#               2. Run the function "findTF2" - specify the input files (diffbind result & PWM) and the output file path
#               3. TF candidates and their info will be saved as tab-separated file in the output path
#
#   Example
#               > source("The_directory_of_FindTFofRegion2.R/FindTFofRegion2.R")
#               > findTF2(regionPath="./data/ATAC_Seq/Mouse_AA_cd4_cd8/01.AACD4vsCtrlCD4_atfc1.5p0.05_1180_peaks_annotation.csv",
#                         pwmRDAPath="./data/CISBP_Mus_musculus/mus_musculus_pwms.rda",
#                         outputPath="./results/CD4_TFs2.txt")
###

###
#   Details
#   
#   We have differentially bound regions from two ATAC-Seq (treatment vs control).
#   Now we would like to know which TFs have TFBSs in those differentially bound regions.
#   It is possible to have sequences of those regions on the mm10 reference genome.
#   And then, All the TFs' valdiated PWMs are matched to the sequences.
#   PWMs can be downloaded from CIS-BP database.
#   With the matched PWMs, we could know whith TFs would bound to the differentially bound regions.
###

findTF2 <- function(regionPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Christiano/Eddy/2018/ATAC_Seq/Mouse_AA_cd4_cd8/data/01.AACD4vsCtrlCD4_atfc1.5p0.05_1180_peaks_annotation.csv",
                    pwmRDAPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Christiano/Eddy/2018/ATAC_Seq/Mouse_AA_cd4_cd8/data/mus_musculus_pwms.rda",
                    outputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Christiano/Eddy/2018/ATAC_Seq/Mouse_AA_cd4_cd8/TF_footprint/CD4_TFs2.txt") {
  
  ### load library
  if(!require(Biostrings)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("Biostrings")
    library(Biostrings)
  }
  if(!require(BSgenome.Mmusculus.UCSC.mm10)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("BSgenome.Mmusculus.UCSC.mm10")
    library(BSgenome.Mmusculus.UCSC.mm10)
  }
  
  
  ### load datasets
  load(pwmRDAPath)
  genome_mm10 <- BSgenome.Mmusculus.UCSC.mm10
  dbRegions <- read.csv(file = regionPath, stringsAsFactors = FALSE, check.names = FALSE)
  
  
  ### get indicies of non-direct motifs 
  indirect_idx <- which(is.na(PWM_TFs_Direct))
  
  
  ### using direct interactions only
  PWMs <- PWMs[-indirect_idx]
  PWM_TFs_Direct <- PWM_TFs_Direct[-indirect_idx]
  
  
  ### get chromosome names
  chr_names <- names(genome_mm10)
  
  
  ### make an empty result object
  result <- data.frame(matrix(NA, 1, ncol(dbRegions)+6))
  colnames(result) <- c(colnames(dbRegions), "TF_Name", "TF_Start", "TF_End", "Motif_ID", "Score", "Sequence")
  
  
  ### iteratively get TF candidates for each differentially bound region
  for(i in 1:nrow(dbRegions)) {
    ### get the sequence of a differentially bound region
    sequence <- subseq(genome_mm10[[dbRegions$seqnames[i]]],
                       start = dbRegions$start[i], end = dbRegions$end[i])
    
    ### iteratively add TF info whenever a given TF's motif matched to the sequence 
    for(j in 1:length(PWMs)) {
      ### get hit regions
      hits <- matchPWM(PWMs[[j]], sequence, min.score="90%", with.score=TRUE)
      
      ## add info
      if(length(hits) > 0) {
        ### motif can be matched to the sequence more than once
        for(k in 1:length(hits)) {
          ### motif can have many TFs
          for(l in 1:length(PWM_TFs_Direct[[j]])) {
            result <- rbind(result, cbind(dbRegions[i,], TF_Name=PWM_TFs_Direct[[j]][l],
                                          TF_Start=start(hits)[k], TF_End=end(hits)[k], Motif_ID=names(PWMs)[j],
                                          Score=mcols(hits)$score[k], Sequence=as.character(hits)[k]))
          }
        }
      }
    }
    
    ### garbage collection
    gc()
    
    ### progress bar
    writeLines(paste(i, nrow(dbRegions), sep = " / "))
  }
  
  
  ### remove the first NA row of the result object
  result <- result[-1,]
  
  
  ### print some info
  print(paste("Total number of TFs in the result:", nrow(result)))
  print(paste("Total number of the unique DB regions in the result:",
              length(which(duplicated(result[,c(1,2,3)]) == FALSE)), "/", nrow(dbRegions)))
  
  
  ### save the result
  write.table(result, file = outputPath, sep = "\t", row.names = FALSE)
  
}


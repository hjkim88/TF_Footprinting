###
#   File name : FindTFofRegion3.R
#   Author    : Hyunjin Kim
#   Date      : May 11, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Find TF candidates that may be bound in regions of interest
#
#   Instruction
#               1. Source("FindTFofRegion3.R")
#               2. Run the function "findTF3" - specify the input files (diffbind result & PWM) and the output file path
#               3. TF candidates and their info will be saved as tab-separated file in the output path
#
#   Example
#               > source("The_directory_of_FindTFofRegion3.R/FindTFofRegion3.R")
#               > findTF3(regionPath="./data/ATAC_Seq/Mouse_AA_cd4_cd8/01.AACD4vsCtrlCD4_atfc1.5p0.05_1180_peaks_annotation.csv",
#                         pwmRDAPath="./data/CISBP_Mus_musculus/mus_musculus_pwms.rda",
#                         outputPath="./results/CD4_TFs3.txt")
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
#
#   This code is different from FindTFofRegion2.R since this uses TFBSTools
#   and the FindTFofRegion2.R uses Biostrings.
#   This code additionally generates p-values of the PWM matchings.
###

findTF3 <- function(regionPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Christiano/Eddy/2018/ATAC_Seq/Mouse_AA_cd4_cd8/data/01.AACD4vsCtrlCD4_atfc1.5p0.05_1180_peaks_annotation.csv",
                    pwmRDAPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Christiano/Eddy/2018/ATAC_Seq/Mouse_AA_cd4_cd8/data/mus_musculus_pwms.rda",
                    outputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Christiano/Eddy/2018/ATAC_Seq/Mouse_AA_cd4_cd8/TF_footprint/CD4_TFs3.txt") {
  
  ### load library
  if(!require(TFBSTools)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("TFBSTools")
    library(TFBSTools)
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
  
  
  ### matrix -> PWMatrix
  PWMs2 <- list()
  for(i in 1:length(PWMs)) {
    PWMs2[[i]] <- PWMatrix(ID = names(PWMs)[i], name = names(PWMs)[i], profileMatrix = PWMs[[i]])
  }
  names(PWMs2) <- names(PWMs)
  
  
  ### make an empty result object
  result <- data.frame(matrix(NA, 1, ncol(dbRegions)+9))
  colnames(result) <- c(colnames(dbRegions), "TF_Name", "TFBS_Start", "TFBS_End",
                        "Motif_ID", "Sequence", "Score", "Sequence_Strand",
                        "PValue_Exact", "PValue_Sampling")
  
  
  ### iteratively get TF candidates for each differentially bound region
  set.seed(1234)
  for(i in 1:nrow(dbRegions)) {
    ### get the sequence of a differentially bound region
    sequence <- subseq(genome_mm10[[dbRegions$seqnames[i]]],
                       start = dbRegions$start[i], end = dbRegions$end[i])
    
    ### iteratively add TF info whenever a given TF's motif matched to the sequence 
    for(j in 1:length(PWMs2)) {
      ### get hit regions
      siteset <- searchSeq(PWMs2[[j]], sequence, seqname=paste0("seq", i), strand="*", min.score="90%")
      
      ## add info
      if(length(siteset) > 0) {
        ### motif can be matched to the sequence more than once
        for(k in 1:length(siteset)) {
          ### motif can have many TFs
          for(l in 1:length(PWM_TFs_Direct[[j]])) {
            result <- rbind(result, cbind(dbRegions[i,], TF_Name=PWM_TFs_Direct[[j]][l],
                                          TFBS_Start=start(siteset)[k], TFBS_End=end(siteset)[k],
                                          Motif_ID=names(PWMs2)[j], Sequence=as.character(views(siteset))[k],
                                          Score=score(siteset)[k], Sequence_Strand=strand(siteset)[k],
                                          PValue_Exact=pvalues(siteset, type="TFMPvalue")[k],
                                          PValue_Sampling=pvalues(siteset, type="sampling")[k]))
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


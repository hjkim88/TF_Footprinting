###
#   File name : GetTFBSRegions.R
#   Author    : Hyunjin Kim
#   Date      : May 8, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Find TFBS regions of all known TFs with their motif info on mouce mm10 genome
#
#   Instruction
#               1. Source("GetTFBSRegions.R")
#               2. Run the function "getRegions" - specify the input file (PWM RDA) and the output file path
#               3. TFBS regions on mm10 will be generated in the output path
#
#   Example
#               > source("The_directory_of_GetTFBSRegions.R/GetTFBSRegions.R")
#               > getRegions(pwmRDAPath="./data/CISBP_Mus_musculus/mus_musculus_pwms.rda",
#                            outputPath="./results/TF_TFBS_mm10.txt")
###

getRegions <- function(pwmRDAPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Christiano/Eddy/2018/ATAC_Seq/Mouse_AA_cd4_cd8/data/mus_musculus_pwms.rda",
                       outputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Christiano/Eddy/2018/ATAC_Seq/Mouse_AA_cd4_cd8/TF_footprint/TF_TFBS_mm10.txt") {
  
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
  
  
  ### get chromosome names
  chr_names <- names(genome_mm10)
  
  
  ### get indicies of non-direct motifs 
  indirect_idx <- which(is.na(PWM_TFs_Direct))
  
  
  ### using direct interactions only
  PWMs <- PWMs[-indirect_idx]
  PWM_TFs_Direct <- PWM_TFs_Direct[-indirect_idx]
  
  
  # ### get start, end, and score from DNA String object
  # getInfoFromDNAStr <- function(tf_name, pwm_id, hits) {
  #   
  #   info <- cbind(rep(tf_name, length(hits[[1]])), rep(pwm_id, length(hits[[1]])),
  #                 rep(names(hits)[1], length(hits[[1]])),
  #                 start(hits[[1]]), end(hits[[1]]), mcols(hits[[1]])$score)
  #   
  #   for(i in 2:length(hits)) {
  #     info <- rbind(info, cbind(rep(tf_name, length(hits[[i]])), rep(pwm_id, length(hits[[i]])),
  #                               rep(names(hits)[i], length(hits[[i]])),
  #                               start(hits[[i]]), end(hits[[i]]), mcols(hits[[i]])$score))
  #   }
  #   
  #   colnames(info) <- c("TF_Name", "PWM_ID", "Chr", "Start", "End", "Score")
  #   
  #   return(info)
  #   
  # }
  
  
  ### create temp directory
  tempDirName <- paste0("./", as.numeric(Sys.time()))
  dir.create(tempDirName)
  
  
  ### a function to get PWM-Genome hits
  getHits <- function(iteration=1) {
    for(i in iteration:length(PWMs)) {
      tryCatch(
        {
          result <- matrix(NA, 1, 6)
          colnames(result) <- paste0(rep("V", ncol(result)), seq(ncol(result)))
          
          # j in 1:22 because 1 - chr1 to 22 - chrM; only they are meaningful chromosomes
          for(j in 1:22) {
             hits <- suppressWarnings(matchPWM(PWMs[[i]], genome_mm10[[j]], min.score="90%", with.score=TRUE))
             result <- rbind(result, cbind(rep(PWM_TFs_Direct[[i]], length(hits)),
                                           rep(names(PWMs)[i], length(hits)),
                                           rep(chr_names[j], length(hits)),
                                           start(hits), end(hits), mcols(hits)$score))
          }
          
          ### set column names
          colnames(result) <- c("TF_Name", "PWM_ID", "Chr", "Start", "End", "Score")
          
          ### remove the first NAs
          result <- result[-1,]
          
          ### save the partial result
          write.table(result, file = paste0(tempDirName, "/temp_", i, ".txt"),
                      sep = "\t", row.names = FALSE)
          
          ### garbage collection
          gc()
          
          ### progress bar
          writeLines(paste(i, length(PWMs), sep = " / "))
        },
        error = function(e) {
          writeLines(paste("ERROR:", e))
          gc()
          getHits(iteration = i)
        }
      )
    }
  }
  
  
  ### PWM-Genome hits
  getHits()
  
  
  ### load the hits and combine them
  # collect result files from the temp dir
  f <- list.files(paste0(tempDirName, "/"))
  # read and combine them
  result <- read.table(file = paste0(tempDirName, "/", f[1]), header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, check.names = FALSE)
  for(i in 2:length(f)) {
    temp <- read.table(file = paste0(tempDirName, "/", f[i]), header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, check.names = FALSE)
    result <- rbind(result, temp)
  }
  
  
  ### save the result
  save(list = c("result"), file = outputPath)
  
  
  ### remove the temp directory
  unlink(tempDirName, recursive = TRUE)
  
}




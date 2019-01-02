###
#   File name : LoadnSavePWMs.R
#   Author    : Hyunjin Kim
#   Date      : May 8, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : load PWM files and save them in R Biostrings object
#
#   Instruction
#               1. Source("LoadnSavePWMs.R")
#               2. Run the function "loadnSavePWMs" - specify the input directory (PWMs & TF info) and output file path
#               3. The result R object will be saved in the output path
#
#   Example
#               > source("The_directory_of_LoadnSavePWMs.R/LoadnSavePWMs.R")
#               > loadnSavePWMs(pwmPath="./data/CISBP_Mus_musculus/pwms_all_motifs/",
#                               tfInfoPath="./data/CISBP_Mus_musculus/TF_Information_all_motifs_plus.txt",
#                               outputPath="./data/CISBP_Mus_musculus/mus_musculus_pwms.rda")
###

loadnSavePWMs <- function(pwmPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Christiano/Eddy/2018/ATAC_Seq/Mouse_AA_cd4_cd8/data/pwms_all_motifs/",
                          tfInfoPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Christiano/Eddy/2018/ATAC_Seq/Mouse_AA_cd4_cd8/data/TF_Information_all_motifs_plus.txt",
                          outputPath="//isilon.c2b2.columbia.edu/ifs/archive/shares/af_lab/Christiano/Eddy/2018/ATAC_Seq/Mouse_AA_cd4_cd8/data/mus_musculus_pwms.rda") {
  
  ### collect PWMs from pwmPath
  f <- list.files(pwmPath)
  f <- f[which(endsWith(f, ".txt"))]
  
  
  ### iteratively load PWMs and make a list
  PWMs <- list()
  for(i in 1:length(f)) {
    PWMs[[i]] <- as.matrix(t(read.table(file = paste0(pwmPath, f[i]), header = TRUE, sep = "\t",
                              row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)))
  }
  names(PWMs) <- sapply(f, function(x) substr(x, 1, nchar(x)-4))
  
  
  ### load TF info
  tfInfo <- read.table(file = tfInfoPath, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, check.names = FALSE, comment.char = "")
  tfInfo_D <- tfInfo[which(tfInfo$TF_Status == "D"),]
  
  
  ### get TF names for every PWM
  # Direct
  PWM_TFs_Direct <- list()
  PWM_TFs_Direct <- lapply(names(PWMs),
                      function(x) {
                        y <- tfInfo_D[which(tfInfo_D$Motif_ID == x), "TF_Name"]
                        if(length(y) > 0) {
                          return(y)
                        } else {
                          return(NA)
                        }
                      }
                    )
  # All
  PWM_TFs_All <- list()
  PWM_TFs_All <- lapply(names(PWMs),
                         function(x) {
                           y <- tfInfo[which(tfInfo$Motif_ID == x), "TF_Name"]
                           if(length(y) > 0) {
                             return(y)
                           } else {
                             return(NA)
                           }
                         }
                  )
  
  
  ### remove the empty PWMs
  # get indicies
  idx <- NULL
  for(i in 1:length(PWMs)) {
    if(length(PWMs[[i]]) == 0) {
      idx <- c(idx, i)
    }
  }
  # remove those
  PWMs <- PWMs[-idx]
  PWM_TFs_Direct <- PWM_TFs_Direct[-idx]
  PWM_TFs_All <- PWM_TFs_All[-idx]
  
  
  ### save the results
  save(list = c("PWMs", "PWM_TFs_Direct", "PWM_TFs_All"), file = outputPath)
  
}


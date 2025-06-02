# Generate and write LCAM per sample to file

pb_sample <- read.csv("src_output/pb_sample.csv",r=1,h=1,stringsAsFactors = F, check.names = FALSE)

source("scripts/get_LCAM_scores.R")
scores <- get_LCAM_scores(as.matrix(pb_sample))

# Write to file
write.csv(
  scores,
  file = "output/1. data preprocessing/lcam/lcam_sample.csv"
)

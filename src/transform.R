# This script converts the leader_et_al dataset from scRNAseq to bulk
# Additionally it maps and filters on metabolic genes using KEGG

if (!exists("lung_ldm")) {
  message("Loading Mt. Sinai scRNAseq data into R")
  load("data/lung_ldm.rd")
}



rm("lung_ldm")

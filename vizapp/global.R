library(shiny)
library(tidyverse)

setwd("../")

# Load Data
annot_vcf <- read_table(Sys.glob("workflow/output/*/*.annot.vcf.gz"))
annot_vcf <- annot_vcf[-grep("^#", annot_vcf$`##fileformat=VCFv4.2`),] %>%
  separate(`##fileformat=VCFv4.2`, into = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "GM12878"), sep = "\t")

#Initialize Table
Final <- data.frame(matrix(ncol = 12, nrow = nrow(annot_vcf)))
colnames(Final) <- c("Chromosome", "Position", "Gene_Variant", "Amino_Acid_Change", "Exon_Number", "Effects", "Reference_Allele", "Alternate_Allele", "Total_Depth",
                     "Alt_Percent", "Allele_Frequency(PopMax)", "ID")

#Add Table Data
Final$Chromosome <- annot_vcf$CHROM
Final$Position <- annot_vcf$POS
Final$ID <- annot_vcf$ID
Final$Reference_Allele <- annot_vcf$REF
Final$Alternate_Allele <- annot_vcf$ALT
for (i in 1:length(annot_vcf$INFO)) {
  Final$Gene_Variant[i] <- ((gsub(";.*", "", substring(annot_vcf$INFO[i], regexpr("ANN", annot_vcf$INFO[i]) + 4)) %>% strsplit(",") %>% unlist())[1] %>%
                             strsplit("\\|") %>% unlist())[4]
  Final$Amino_Acid_Change[i] <- ((gsub(";.*", "", substring(annot_vcf$INFO[i], regexpr("ANN", annot_vcf$INFO[i]) + 4)) %>% strsplit(",") %>% unlist())[1] %>%
                                  strsplit("\\|") %>% unlist())[11]
  Final$Exon_Number[i] <- (((gsub(";.*", "", substring(annot_vcf$INFO[i], regexpr("ANN", annot_vcf$INFO[i]) + 4)) %>% strsplit(",") %>% unlist())[1] %>%
                            strsplit("\\|") %>% unlist())[9] %>% strsplit("/") %>% unlist())[1]
  Final$Effects[i] <- ((gsub(";.*", "", substring(annot_vcf$INFO[i], regexpr("ANN", annot_vcf$INFO[i]) + 4)) %>% strsplit(",") %>% unlist())[1] %>%
                        strsplit("\\|") %>% unlist())[2]
  Final$Alt_Percent[i] <- ((strsplit(annot_vcf$GM12878[i], ":") %>% unlist())[4] %>% strtoi())/((strsplit(annot_vcf$GM12878[i], ":") %>% unlist())[2] %>% strtoi()) *100
  if (grepl("DP", annot_vcf$INFO[i])) {
    Final$Total_Depth[i] <- gsub(";.*", "", substring(annot_vcf$INFO[i], regexpr("DP", annot_vcf$INFO[i]) + 3))
  }
  if (grepl("AF_POPMAX", annot_vcf$INFO[i])) {
    Final$`Allele_Frequency(PopMax)`[i] <- gsub(";.*", "", substring(annot_vcf$INFO[i], regexpr("AF_POPMAX", annot_vcf$INFO[i]) + 10))
  }
}
Final$Alt_Percent <- round(Final$Alt_Percent, digits = 2)
library(readxl)
library(maftools)
library(readr)
library("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)

Patho_MAF_File <- read_delim("H:/My Drive/Pathogenic_Landscape/data/absolute/clinical_research_filtered_combined/Absolute_Patho_MAF.maf", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
Indig_Patho_MAF <- read_delim("H:/My Drive/Pathogenic_Landscape/data/indigene/Indig_Patho_MAF.maf", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
combined_patho <- read_delim("H:/My Drive/Pathogenic_Landscape/assets/absolute_indie_patho/Combined_Patho_MAF.maf", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

patho = read.maf(maf = Patho_MAF_File)
patho = read.maf(maf = Indig_Patho_MAF)
patho = read.maf(maf = combined_patho)

plotmafSummary(maf = patho, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE, showBarcodes = TRUE, textSize = 0.5)
# too many sample, variants per sample crowded with sample names, so disabling barcodes
plotmafSummary(maf = patho, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE, showBarcodes = FALSE, textSize = 0.5)

# way too many samples, plot broke
oncoplot(maf = patho, top = 10, draw_titv = TRUE)

patho.titv = titv(maf = patho, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = patho.titv)

# Mutex/CoOcc
somaticInteractions(maf = patho, top = 30, pvalue = c(0.05, 0.1), fontSize = 0.8, showCounts = FALSE, leftMar = 6, topMar = 6, countStats = "all")
comut <- somaticInteractions(maf = patho, top = 30, pvalue = c(0.05, 0.1), fontSize = 0.8, showCounts = FALSE, leftMar = 7, topMar = 7, countStats = "all")
write.csv(comut, "indig_comut.csv")
write.csv(comut, "absolute_comut.csv")
write.csv(comut, "combined_comut.csv")


# Example: Extract mutation data for a specific gene pair (e.g., TP53 and KRAS)
genes_of_interest <- c("APC", "KRAS")

# Subset the MAF to only include these two genes
subset_maf <- subsetMaf(maf = patho, genes = genes_of_interest)

# Extract unique samples with mutations in the selected genes
samples_with_mutations <- unique(subset_maf@data$Tumor_Sample_Barcode)

# View the involved samples
print(samples_with_mutations)






#-------------------------------------------------------------------------------
# didnt work, messed up aachange col i think, and oncodrive issue
patho.sig = oncodrive(maf = patho, AACol = 'aaChnage', minMut = 5, pvalMethod = 'zscore')

# Drug-Gene Interactions
dgi = drugInteractions(maf = patho, fontSize = 0.7)

#
pws = pathways(maf = patho, plotType = 'treemap', panelWidths = c(100,200,200))
plotPathways(maf = patho, pathlist = pws)

#
patho.tnm = trinucleotideMatrix(maf = patho, add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
plotApobecDiff(tnm = patho.tnm, maf = patho, pVal = 0.2)
library('NMF')
patho.sign = estimateSignatures(mat = patho.tnm, nTry = 6) # error

#
patho.mutsig.corrected = prepareMutSig(maf = patho)
write.table(patho.mutsig.corrected, "patho_mutsig.maf", sep = "\t", quote = FALSE, row.names = FALSE)



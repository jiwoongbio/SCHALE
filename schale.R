args <- commandArgs(TRUE)
functionFile <- args[1]
outputPrefix <- args[2]
variantAlleleDepthFile <- args[3]
chromosomeFile <- args[4]
centromereFile <- args[5]
normal.sample <- args[6]
tumor.sample <- args[7]

source(functionFile)
variant.table <- read.table(variantAlleleDepthFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

chromosome.table <- read.table(chromosomeFile, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
chromosome.lengths <- chromosome.table[, 2]
names(chromosome.lengths) <- chromosome.table[, 1]

centromere.table <- read.table(centromereFile, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(centromere.table) <- c("chromosome", "start", "end")

variant.table <- schale.variant.table(variant.table, normal.sample = normal.sample)
segment.table <- schale.segment.table(variant.table, centromere.table = centromere.table)

CNA.offset <- schale.calculate.CNA.offset(segment.table)
variant.table <- schale.normalize.CNA.offset(variant.table, CNA.offset)
segment.table <- schale.normalize.CNA.offset(segment.table, CNA.offset)

pdf(file = paste(outputPrefix, "schale.genome.pdf", sep = "."), width = 24, height = 6)
schale.genome.plot(variant.table, segment.table, tumor.sample, chromosome.lengths = chromosome.lengths)
dev.off()

pdf(file = paste(outputPrefix, "schale.segment.pdf", sep = "."), width = 8, height = 6)
schale.segment.plot(segment.table, tumor.sample, col = "gray")
dev.off()

getwd()
setwd('C:\\Users\\chaud\\Pictures\\edgene\\try-1\\raw_seq')
library(dada2)
list.files()
path = 'C:\\Users\\chaud\\Pictures\\edgene\\try-1\\raw_seq'
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
# Extracting sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
fnFs

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, # Remove any reads with ambiguous bases
                     compress=TRUE, 
                     multithread=FALSE)

head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, minOverlap = 12)
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

getwd()
#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/chaud/Pictures/edgene/try-1/raw_seq/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "~/tax/silva_species_assignment_v132.fa.gz")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

write.csv(track,"track.csv")
saveRDS(seqtab.nochim,"seqtab.nochim.rds")
saveRDS(taxa, "taxa.rds")

library(phyloseq)
otu.tab <- otu_table(object = seqtab.nochim, taxa_are_rows = FALSE)
tax.tab <- tax_table(object = taxa)

ps <- phyloseq(otu.tab, tax.tab) #one phyloseq object
ps

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

otu.tab <- as.data.frame(otu_table(ps))
otu.tab <- t(otu.tab)
tax.tab <- as.data.frame(tax_table(ps))

write.table(x = otu.tab, file = './otu_tab.txt', quote = F, sep = '\t', col.names = NA)
write.table(x = tax.tab, file = './tax_tab.txt', quote = F, sep = '\t', col.names = NA)

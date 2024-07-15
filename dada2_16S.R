require(dada2)
require(Biostrings)
require(DECIPHER)
require(phyloseq)
require(seqinr)
require(data.table)
require(metagMisc)
require(tibble)
library(phyloseq)
library(tidyverse)
library(vegan)

mult = 10
mlt = 10
path_trein_set = "~/storage/somebases/silva_nr99_v138.1_train_set.fa.gz"
path_trein_set_species = "~/storage/somebases/silva_species_assignment_v138.1.fa.gz"
name_Brief = "nameBrief.txt"
path <- "/home/gladkov/storage/cellulolit/paper/16/raw"
active_dir <- "/home/gladkov/storage/cellulolit/paper/16"
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
on.exit()
filtFs <- file.path(active_dir, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(active_dir, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,180), trimLeft=c(19,20), maxN=0, truncQ=2, maxEE=c(2,5), rm.phix=TRUE, compress=TRUE, multithread=mult)
errF <- learnErrors(filtFs, multithread=mult)
errR <- learnErrors(filtRs, multithread=mult)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err=errF, multithread=mult)
dadaRs <- dada(derepRs, err=errR, multithread=mult)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
getN <- function(x) sum(getUniques(x))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=mult, verbose=TRUE)

briefToSeq <- colnames(seqtab.nochim)
names(briefToSeq) <- paste0("Seq", seq(ncol(seqtab.nochim))) 
st.brief <- seqtab.nochim
colnames(st.brief) <- names(briefToSeq)

write.table(briefToSeq, file = name_Brief, sep = "\t")

taxa.dada2 <- assignTaxonomy(briefToSeq, path_trein_set , multithread=mult)
taxa.dada2.species <- assignSpecies(briefToSeq, path_trein_set_species) # maybe use not seqs but brief 
rownames(taxa.dada2.species) <- rownames(briefToSeq)
briefToSeq.df <- data.frame(briefToSeq)
rownames(taxa.dada2.species) <- rownames(briefToSeq.df)
rownames(taxa.dada2) <- rownames(taxa.dada2.species)
taxa <- cbind2(taxa.dada2, taxa.dada2.species[,2])
colnames(taxa)[7] <- "Species"

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.table(track, file = "track.tsv", sep= "\t", col.names = NA, quote=FALSE)

st.brief.t.df <- data.frame(t(st.brief))
write.table(st.brief.t.df, file = "otu_table.txt", sep= "\t", col.names = NA, quote=FALSE)

briefToSeq.ls <- as.list(briefToSeq.df[,c("briefToSeq")])
briefToSeq.names <- as.list(rownames(briefToSeq.df))
seqinr::write.fasta(briefToSeq.ls, briefToSeq.names , "rep_seq.fasta", as.string = FALSE)

write.table(taxa, file = "taxa.txt", sep= "\t", col.names = NA, quote=FALSE)

otu_table <- read_tsv("16/otu_table.txt") %>% 
  column_to_rownames("...1")

hue <- otu_table %>% 
  colnames()

colnames(otu_table) <- unlist(purrr::map(stringr::str_split(hue, "\\.", 3), function(x) paste0(x[[2]], '_' ,x[[3]]))) 
hue_mod <- colnames(otu_table)

map <- data.frame(
  ID = hue_mod,
  Substrate = unlist(purrr::map(stringr::str_split(hue_mod, "\\_", 3), function(x) x[[2]] %>% substring(1, 1) %>% as.factor())),
  Repeat = unlist(purrr::map(stringr::str_split(hue_mod, "\\_", 3), function(x) x[[2]] %>% substring(1, 3) %>% as.factor())),
  Time = unlist(purrr::map(stringr::str_split(hue_mod, "\\_", 3), function(x) x[[2]] %>% substring(2, 3) %>% as.double()))
) %>% column_to_rownames("ID")

otus <- otu_table %>% 
  as.matrix() %>% 
  t()

taxa <- read_tsv("16/taxa.txt") %>% 
  column_to_rownames("...1") %>% 
  as.matrix()


taxa[taxa == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "Allorhizobium"

ps <- phyloseq(otu_table(otus, taxa_are_rows=FALSE), 
               sample_data(map),
               tax_table(taxa))

ps <- merge_phyloseq(ps, refseq(Biostrings::readDNAStringSet("16/rep_seq.fasta")))

ps.f <- filter_chrmitnas(ps)

Biostrings::writeXStringSet(ps.f@refseq, file="16/ref_filt.fasta")

tree_new <- read_tree(treefile="16/tree.nwk")

ps.f@phy_tree <- tree_new

ps.f  <- prune_taxa(taxa_sums(ps.f) > 0, ps.f)
ps.f@sam_data <- read_tsv('map_resp.tsv') %>% 
  column_to_rownames('ID') %>% 
  sam_data()

ps.f@tax_table <- ps.f@tax_table@.Data %>% 
  data.frame() %>% 
  mutate(Phylum = str_replace_all(Phylum, c("Firmicutes" = "Bacillota",
                                            "Chloroflexi" = "Chloroflexota",
                                            "Proteobacteria" = "Pseudomonadota",
                                            "Desulfobacterota" = "Thermodesulfobacteriota",
                                            "Crenarchaeota" = "Thermoproteota",
                                            "Cyanobacteria" = "Cyanobacteria"))) %>% 
  as.matrix() %>% 
  tax_table()


write_rds(ps.f, "ps.f")


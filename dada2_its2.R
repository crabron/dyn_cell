library(dada2)
library(ShortRead)
library(Biostrings)
library(seqinr)
library(tidyverse)
library(fuzzyjoin)
require(dada2)
require(Biostrings)
require(DECIPHER)
require(seqinr)
require(data.table)
require(metagMisc)
require(tibble)
library(tidyverse)


mult = 10
mlt = 10
path_trein_set = "~/storage/somebases/unite_short_2023_ref.fa"
path_trein_set_2 = "~/storage/somebases/unite_2023_ref.fa"
name_Brief = "nameBrief.txt"
path <- "/home/gladkov/storage/cellulolit/paper/its/raw"
active_dir <- "/home/gladkov/storage/cellulolit/paper/its"
FWD <- "GCATCGATGAAGAACGCAGC"  ## CHANGE ME to your forward primer sequence
REV <- "TCCTCCGCTTATTGATATGC"  ## CHANGE ME...

fnFs <- sort(list.files(path, pattern = "R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2_001.fastq.gz", full.names = TRUE))


allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = reverse(dna), 
        RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

cutadapt <- "/home/gladkov/.conda/envs/pandas/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))


FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
R3.flags <- paste("--minimum-lengt 10")
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags,R3.flags, "-n", 2,  # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files             
}

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 3), 
    truncQ = 3, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
errF <- learnErrors(filtFs, multithread=mult)
errR <- learnErrors(filtRs, multithread=mult)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err=errF, multithread=mult, pool=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=mult, pool=TRUE)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)

getN <- function(x) sum(getUniques(x))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=mult, verbose=TRUE)

briefToSeq <- colnames(seqtab.nochim)
names(briefToSeq) <- paste0("Seq", seq(ncol(seqtab.nochim))) 
st.brief <- seqtab.nochim
colnames(st.brief) <- names(briefToSeq)

write.table(briefToSeq, file = name_Brief, sep = "\t")


taxa.dada2 <- assignTaxonomy(briefToSeq, path_trein_set , multithread=mult, tryRC=TRUE)
rownames(taxa.dada2) <- rownames(briefToSeq)
briefToSeq.df <- data.frame(briefToSeq)
rownames(taxa.dada2) <- rownames(briefToSeq.df)


write.table(taxa.dada2, file = "taxa_1.txt", sep= "\t", col.names = NA, quote=FALSE)

tx.fail <- taxa.dada2 %>% 
  as.data.frame() %>% 
  filter(is.na(Phylum)) %>% 
  rownames()
  
tx.nice <- taxa.dada2 %>% 
  as.data.frame() %>% 
  filter(!is.na(Phylum))

brSeq.fail <- briefToSeq[tx.fail] 
taxa.dada3 <- assignTaxonomy(brSeq.fail, path_trein_set_2 , multithread=mult, tryRC=TRUE)


rownames(taxa.dada3) <- tx.fail
write.table(taxa.dada3, file = "taxa_2.txt", sep= "\t", col.names = NA, quote=FALSE)

taxa <- rbind(tx.nice, taxa.dada3)
taxa <- taxa %>% 
  arrange(naturalsort::naturalsort(row.names(.), decreasing = TRUE))
write.table(taxa, file = "taxa.txt", sep= "\t", col.names = NA, quote=FALSE)

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


saveRDS(errF, 'errF')
saveRDS(errR, 'errR')

otu_table <- read_tsv("its/otu_table.txt") %>% 
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

taxa <- read_tsv("its/taxa.txt") %>% 
  column_to_rownames("...1") %>% 
  as.matrix()


taxa[taxa == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "Allorhizobium"


psits <- phyloseq(otu_table(otus, taxa_are_rows=FALSE), 
                  sample_data(map),
                  tax_table(taxa))

psits <- merge_phyloseq(psits, refseq(Biostrings::readDNAStringSet("its/rep_seq.fasta")))

psits.f <- filter_chrmitnas(psits)

Biostrings::writeXStringSet(psits.f@refseq, file="its/ref_filt_its.fasta")

plot_rich_reads_samlenames_lm(psits.f, 'Repeat')
psits.f <- prune_samples(sampleNames(psits.f) != 'its2_C14.4', psits.f)
plot_rich_reads_samlenames_lm(psits.f, 'Repeat')

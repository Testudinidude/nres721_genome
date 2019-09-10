# navigate to directory
setwd("/Volumes/G-DRIVE/Dropbox/Teaching/2019_NRES_721/nres721_genome/mapped")

# read in BAM files
Brid_R1R2_BAM <- readLines("Ubrucei.txt")
Cala_R1R2_BAM <- readLines("../mapped/Cala_R1R2_mapped.txt")
Cala_male_R1R2_BAM <- readLines("../mapped/Cala_male_R1R2_mapped.txt")

# mine useful information
Brid_R1R2_seqs <- matrix(nrow=length(Brid_R1R2_BAM),ncol=4)
for(i in 1:length(Brid_R1R2_BAM)){
	X <- strsplit(Brid_R1R2_BAM[i],"\t")[[1]]
	Xn <- nchar(X[10])
	Brid_R1R2_seqs[i,] <- c(X[c(1,3,10)],Xn)
}

Cala_R1R2_seqs <- matrix(nrow=length(Cala_R1R2_BAM),ncol=4)
for(i in 1:length(Cala_R1R2_BAM)){
  X <- strsplit(Cala_R1R2_BAM[i],"\t")[[1]]
  Xn <- nchar(X[10])
  Cala_R1R2_seqs[i,] <- c(X[c(1,3,10)],Xn)
}

Cala_male_R1R2_seqs <- matrix(nrow=length(Cala_male_R1R2_BAM),ncol=4)
for(i in 1:length(Cala_male_R1R2_BAM)){
  X <- strsplit(Cala_male_R1R2_BAM[i],"\t")[[1]]
  Xn <- nchar(X[10])
  Cala_male_R1R2_seqs[i,] <- c(X[c(1,3,10)],Xn)
}

# keep only loci that mapped
Brid_R1R2_seqMatching <- Brid_R1R2_seqs[Brid_R1R2_seqs[,2]!="*",]
Xn <- as.numeric(Brid_R1R2_seqMatching[,4])
Brid_R1R2_seqMatching <- as.data.frame(Brid_R1R2_seqMatching[,-4])
Brid_R1R2_seqMatching <- data.frame(Brid_R1R2_seqMatching,Xn)
names(Brid_R1R2_seqMatching) <- c("locus","genome","sequence","length")

Cala_R1R2_seqMatching <- Cala_R1R2_seqs[Cala_R1R2_seqs[,2]!="*",]
Xn <- as.numeric(Cala_R1R2_seqMatching[,4])
Cala_R1R2_seqMatching <- as.data.frame(Cala_R1R2_seqMatching[,-4])
Cala_R1R2_seqMatching <- data.frame(Cala_R1R2_seqMatching,Xn)
names(Cala_R1R2_seqMatching) <- c("locus","genome","sequence","length")

Cala_male_R1R2_seqMatching <- Cala_male_R1R2_seqs[Cala_male_R1R2_seqs[,2]!="*",]
Xn <- as.numeric(Cala_male_R1R2_seqMatching[,4])
Cala_male_R1R2_seqMatching <- as.data.frame(Cala_male_R1R2_seqMatching[,-4])
Cala_male_R1R2_seqMatching <- data.frame(Cala_male_R1R2_seqMatching,Xn)
names(Cala_male_R1R2_seqMatching) <- c("locus","genome","sequence","length")

dim(Brid_R1R2_seqMatching)[1]/dim(Brid_R1R2_seqs)[1] # >98% (of matches kept; not of loci)
dim(Cala_R1R2_seqMatching)[1]/dim(Cala_R1R2_seqs)[1] # >99% (of matches kept; not of loci)
dim(Cala_male_R1R2_seqMatching)[1]/dim(Cala_male_R1R2_seqs)[1] # >99% (of matches kept; not of loci)

# remove loci that map to more than one contig
# first, remove duplicate rows where a locus maps to the same contig twice
Brid_R1R2_seqMatching_nodups <- Brid_R1R2_seqMatching[!duplicated(Brid_R1R2_seqMatching[,1:2]),]
dim(Brid_R1R2_seqMatching_nodups)[1]/dim(Brid_R1R2_seqMatching)[1] # got rid of 62% of rows

# then, remove rows with duplicate locus names
# (because now, duplicate locus names occur only when a locus mapped to more than one contig)
Brid_R1R2_seqMatching_nodupsnodups <- Brid_R1R2_seqMatching_nodups[!(duplicated(Brid_R1R2_seqMatching_nodups[,1]) | duplicated(Brid_R1R2_seqMatching_nodups[,1], fromLast = TRUE)), ]
dim(Brid_R1R2_seqMatching_nodupsnodups)[1]/dim(Brid_R1R2_seqMatching_nodups)[1] # got rid of ~15% of loci

# remove loci that map to more than one contig
# first, remove duplicate rows where a locus maps to the same contig twice
Cala_R1R2_seqMatching_nodups <- Cala_R1R2_seqMatching[!duplicated(Cala_R1R2_seqMatching[,1:2]),]
dim(Cala_R1R2_seqMatching_nodups)[1]/dim(Cala_R1R2_seqMatching)[1] # got rid of 58% of rows
# then, remove rows with duplicate locus names
# (because now, duplicate locus names occur only when a locus mapped to more than one contig)
Cala_R1R2_seqMatching_nodupsnodups <- Cala_R1R2_seqMatching_nodups[!(duplicated(Cala_R1R2_seqMatching_nodups[,1]) | duplicated(Cala_R1R2_seqMatching_nodups[,1], fromLast = TRUE)), ]
dim(Cala_R1R2_seqMatching_nodupsnodups)[1]/dim(Cala_R1R2_seqMatching_nodups)[1] # got rid of ~1% of loci

# remove loci that map to more than one contig
# first, remove duplicate rows where a locus maps to the same contig twice
Cala_male_R1R2_seqMatching_nodups <- Cala_male_R1R2_seqMatching[!duplicated(Cala_male_R1R2_seqMatching[,1:2]),]
dim(Cala_male_R1R2_seqMatching_nodups)[1]/dim(Cala_male_R1R2_seqMatching)[1] # got rid of 58% of rows
# then, remove rows with duplicate locus names
# (because now, duplicate locus names occur only when a locus mapped to more than one contig)
Cala_male_R1R2_seqMatching_nodupsnodups <- Cala_male_R1R2_seqMatching_nodups[!(duplicated(Cala_male_R1R2_seqMatching_nodups[,1]) | duplicated(Cala_male_R1R2_seqMatching_nodups[,1], fromLast = TRUE)), ]
dim(Cala_male_R1R2_seqMatching_nodupsnodups)[1]/dim(Cala_male_R1R2_seqMatching_nodups)[1] # got rid of ~1% of loci

# get number of unique contigs matched
numUniqContigs_Brid <- length(unique(Brid_R1R2_seqMatching_nodupsnodups$genome)) #766
UniqContigs_Brid <- unique(Brid_R1R2_seqMatching_nodupsnodups$genome)
# plot all
barplot(table(Brid_R1R2_seqMatching_nodupsnodups$genome))
# plot only chromosomes
barplot(table(as.character(Brid_R1R2_seqMatching_nodupsnodups$genome[grep("NC_", Brid_R1R2_seqMatching_nodupsnodups$genome)])))

numUniqContigs_Cala <- length(unique(Cala_R1R2_seqMatching_nodupsnodups$genome)) #182
UniqContigs_Cala <- unique(Cala_R1R2_seqMatching_nodupsnodups$genome)
# plot all
barplot(table(Cala_R1R2_seqMatching_nodupsnodups$genome))
# plot only chromosomes
barplot(table(as.character(Cala_R1R2_seqMatching_nodupsnodups$genome[grep("NC_", Cala_R1R2_seqMatching_nodupsnodups$genome)])))
# how many map to x chromosome?
length(Cala_R1R2_seqMatching_nodupsnodups$genome[Cala_R1R2_seqMatching_nodupsnodups$genome=="NC_006621.3"])
# how many map to a chromosome?
length(Cala_male_R1R2_seqMatching_nodupsnodups$locus[grep("CM", Cala_male_R1R2_seqMatching_nodupsnodups$genome)])

numUniqContigs_Cala_male <- length(unique(Cala_male_R1R2_seqMatching_nodupsnodups$genome)) #182
UniqContigs_Cala_male <- unique(Cala_male_R1R2_seqMatching_nodupsnodups$genome)
# plot all
barplot(table(Cala_male_R1R2_seqMatching_nodupsnodups$genome))
# plot only chromosomes
barplot(table(as.character(Cala_male_R1R2_seqMatching_nodupsnodups$genome[grep("CM", Cala_male_R1R2_seqMatching_nodupsnodups$genome)])))
# how many map to x chromosome?
length(Cala_male_R1R2_seqMatching_nodupsnodups$genome[Cala_male_R1R2_seqMatching_nodupsnodups$genome=="CM016469.1"]) # 873
# how many map to y chromosome?
length(Cala_male_R1R2_seqMatching_nodupsnodups$genome[Cala_male_R1R2_seqMatching_nodupsnodups$genome=="CM016470.1"]) # 49

# export subset of loci for Cala
Cala_X_Chrom <- as.character(Cala_male_R1R2_seqMatching_nodupsnodups$locus[grep("CM016469.1", Cala_male_R1R2_seqMatching_nodupsnodups$genome)])
Cala_Y_Chrom <- as.character(Cala_male_R1R2_seqMatching_nodupsnodups$locus[grep("CM016470.1", Cala_male_R1R2_seqMatching_nodupsnodups$genome)])
Cala_All_Chrom <- as.character(Cala_male_R1R2_seqMatching_nodupsnodups$locus[grep("CM", Cala_male_R1R2_seqMatching_nodupsnodups$genome)])
Cala_Non_Sex_Chrom <- setdiff(Cala_All_Chrom, c(Cala_X_Chrom,Cala_Y_Chrom))

Cala_Non_Sex_Chrom_Subset <- sample(Cala_Non_Sex_Chrom, size = 20000 - length(c(Cala_X_Chrom,Cala_Y_Chrom)), replace = FALSE)
Cala_Subset <- c(Cala_X_Chrom, Cala_Y_Chrom, Cala_Non_Sex_Chrom_Subset)
length(Cala_Subset)

# export subset of loci for Brid
Brid_All_Chrom <- Brid_R1R2_seqMatching_nodupsnodups$locus[grep("NC_", Brid_R1R2_seqMatching_nodupsnodups$genome)]
# Brid_Subset <- sample(Brid_All_Chrom, size = 20000, replace = FALSE) # too few

# write loci kept
write.table(Cala_Subset, "Cala_20K_Loci.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Brid_All_Chrom, "Brid_15K_Loci.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# from: https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/
# "NC" = Complete genomic molecule, usually reference assembly
# "NW" = Contig or scaffold, primarily WGSa

pdf(file = "3RAD_loci_per_chromosome.pdf", height = 12)
par(mfrow = c(3,1), mar = c(7,4,4,2) + 0.1)
barplot(table(as.character(Cala_R1R2_seqMatching_nodupsnodups$genome[grep("NC_", Cala_R1R2_seqMatching_nodupsnodups$genome)])),
        col = "cornflowerblue", xlab = "", ylab = "",
        main = "coyotes", xaxt = 'n', space = 0, xlim = c(0,40))
axis (side = 1, labels = unique(Cala_R1R2_seqMatching_nodupsnodups$genome[grep("NC_", Cala_R1R2_seqMatching_nodupsnodups$genome)]),
      las = 2, at = 0.5+c(0:c(length(unique(Cala_R1R2_seqMatching_nodupsnodups$genome[grep("NC_", Cala_R1R2_seqMatching_nodupsnodups$genome)]))-1)),
      cex.axis = 0.5)
mtext(side = 1, text = "Chromosome", line = 5)
mtext(side = 2, text = "Number of Loci", line = 2.5)

barplot(table(as.character(Cala_male_R1R2_seqMatching_nodupsnodups$genome[grep("CM", Cala_male_R1R2_seqMatching_nodupsnodups$genome)])),
        col = "cornflowerblue", xlab = "", ylab = "",
        main = "coyotes", xaxt = 'n', space = 0)
axis (side = 1, labels = unique(Cala_male_R1R2_seqMatching_nodupsnodups$genome[grep("CM", Cala_male_R1R2_seqMatching_nodupsnodups$genome)]),
      las = 2, at = 0.5+c(0:c(length(unique(Cala_male_R1R2_seqMatching_nodupsnodups$genome[grep("CM", Cala_male_R1R2_seqMatching_nodupsnodups$genome)]))-1)),
      cex.axis = 0.5)
mtext(side = 1, text = "Chromosome", line = 5)
mtext(side = 2, text = "Number of Loci", line = 2.5)

barplot(table(as.character(Brid_R1R2_seqMatching_nodupsnodups$genome[grep("NC_", Brid_R1R2_seqMatching_nodupsnodups$genome)])),
        col = "orange", xlab = "", ylab = "",
        main = "pygmy rabbits", xaxt = 'n', space = 0)
axis (side = 1, labels = unique(Brid_R1R2_seqMatching_nodupsnodups$genome[grep("NC_", Brid_R1R2_seqMatching_nodupsnodups$genome)]),
      las = 2, at = 0.5+c(0:c(length(unique(Brid_R1R2_seqMatching_nodupsnodups$genome[grep("NC_", Brid_R1R2_seqMatching_nodupsnodups$genome)]))-1)),
      cex.axis = 0.5)
mtext(side = 1, text = "Chromosome", line = 5)
mtext(side = 2, text = "Number of Loci", line = 2.5)
dev.off()

pdf(file = "3RAD_loci_per_contig.pdf")
par(mfrow = c(2,1))
barplot(table(as.character(Cala_R1R2_seqMatching_nodupsnodups$genome)),
        col = "cornflowerblue", xlab = "", ylab = "",
        main = "coyotes", xaxt = 'n', space = 0)
axis (side = 1, labels = c("",""),
      las = 2, at = 0.5+c(0,c(length(unique(Cala_R1R2_seqMatching_nodupsnodups$genome[grep("NC_", Cala_R1R2_seqMatching_nodupsnodups$genome)]))-1)),
      cex.axis = 0.5, tick = TRUE)
axis (side = 1, labels = c("chrom."),
      las = 1, at = mean(c(0,c(length(unique(Cala_R1R2_seqMatching_nodupsnodups$genome[grep("NC_", Cala_R1R2_seqMatching_nodupsnodups$genome)]))-1))),
      cex.axis = 1, tick = FALSE)

mtext(side = 1, text = "Contig/Scaffold", line = 2)
mtext(side = 2, text = "Number of Loci", line = 2.5)

barplot(table(as.character(Brid_R1R2_seqMatching_nodupsnodups$genome)),
        col = "orange", xlab = "", ylab = "",
        main = "pygmy rabbits", xaxt = 'n', space = 0)
axis (side = 1, labels = c("",""),
      las = 2, at = 0.5+c(0,c(length(unique(Brid_R1R2_seqMatching_nodupsnodups$genome[grep("NC_", Brid_R1R2_seqMatching_nodupsnodups$genome)]))-1)),
      cex.axis = 0.5, tick = TRUE)
axis (side = 1, labels = c("chrom."),
      las = 1, at = mean(c(0,c(length(unique(Brid_R1R2_seqMatching_nodupsnodups$genome[grep("NC_", Brid_R1R2_seqMatching_nodupsnodups$genome)]))-1))),
      cex.axis = 1, tick = FALSE)
mtext(side = 1, text = "Contig/Scaffold", line = 2)
mtext(side = 2, text = "Number of Loci", line = 2.5)
dev.off()
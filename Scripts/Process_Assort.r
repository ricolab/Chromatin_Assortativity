# Script to calculate abundances of ChIP-seq features in PHiC fragments 
# Divide network into promoter - promoter and promoter - other end networks
# calculate assortativity of features in networks

library("igraph")
library("GenomicRanges")
library("tidyverse")
#file.edit("/data/Projects/kat/Projects/Assortativity/ChAs_collated.R")

# Load mESC PHiC interaction map
PCHiC_map <- read.table(file = "Data/PCHiC_interaction_map.txt", header = TRUE, sep = "\t")

# Concatenate "chr" onto chromosome numbers only if not present already
PCHiC_map$baitChr <- paste("chr", PCHiC_map$baitChr, sep = "")
PCHiC_map$oeChr <- paste("chr", PCHiC_map$oeChr, sep = "")

# Filter interactions with > 5 in wildtype cells 
PCHiC_wt <- dplyr::filter(PCHiC_map, mESC_wt > 5)
PCHiC_wt$baits <- paste(PCHiC_wt$baitChr, PCHiC_wt$baitStart, sep = "_")
PCHiC_wt$OEs <- paste(PCHiC_wt$oeChr, PCHiC_wt$oeStart, sep = "_")

# Get unique bait regions and promoter OE regions and non-promoter OE regions
PCHiC_baits <- unique(PCHiC_wt$baits)
PCHiC_OE_all <- unique(PCHiC_wt$OEs)

PCHiC_PromOE <- PCHiC_OE_all[which(PCHiC_OE_all %in% PCHiC_baits)]
PCHiC_nonPromOE <- PCHiC_OE_all[-which(PCHiC_OE_all %in% PCHiC_baits)]

# Create igraph networks
PCHiC_wholenet <- graph_from_data_frame(PCHiC_wt[,c(15:16)], directed = FALSE)
PCHiC_PP <- delete.vertices(PCHiC_wholenet, V(PCHiC_wholenet)[which(V(PCHiC_wholenet)$name %in% PCHiC_nonPromOE)])
PCHiC_POE <- delete.vertices(PCHiC_wholenet, V(PCHiC_wholenet)[which(V(PCHiC_wholenet)$name %in% PCHiC_PromOE)])

# Put all bait and all OE regions into a BED file format and then GRanges object

colnames(PCHiC_wt)[c(1:3, 6:8)] <- rep(c("chr", "start", "end"), 2)
PCHiC_bed <- unique(rbind(PCHiC_wt[,c(1:3)], PCHiC_wt[,c(6:8)]))
PCHiC_GRange <- with(PCHiC_bed, GRanges(chr, IRanges(start, end)))
PCHiC_GRange$ID <- paste(PCHiC_bed$chr, PCHiC_bed$start, sep = "_")

# Load in peak/binarised matrix of ChIP-seq features
binarised_all <- read.table("Data/original_chromFeatures.txt", header = TRUE, sep = "\t")

# Put binarised peaks into GRanges object 
bedepi <- with(binarised_all, GRanges(chr, ranges=IRanges(start, end),strand=Rle(strand(rep("*", nrow(binarised_all))))))
mcols(bedepi) <- binarised_all[,-c(1,2,3)]
names(bedepi) <- paste(binarised_all[,1], binarised_all[,2], sep="_")

# Overlap peaks with PHiC fragments
overlaps <- findOverlaps(PCHiC_GRange, bedepi)

## matching overlapping bedchicmore IDs (start) with ChIP-seq features
match_hit <- data.frame(PCHiC_GRange$ID[queryHits(overlaps)],as.data.frame(mcols(bedepi)[subjectHits(overlaps),]),stringsAsFactors=T)
colnames(match_hit)[1] <- "fragment"

# aggregate windows in fragments, collapse fragments
data.dt <- data.table(match_hit)
setkey(data.dt, fragment) #sorts ascending by fragment
agchic <- data.frame(data.dt[, lapply(.SD, mean), by = fragment]) #mean of ChIP-seq features by fragment
agchic <- data.frame(agchic[,-1], row.names = agchic[,1])


### Function to calculate assortativity of chromatin features in a network

## G = igraph network
## data = feature abundance matrix

## original function

calc_assort <- function(G, data){
  names <- colnames(data)
  ass <- list()
  G_epi <- list()
  for (i in c(1:ncol(data))){
    G <- set.vertex.attribute(G, names[i],value = data[V(G)$name, names[i]])
    attsel <- which(names(vertex.attributes(G)) == names[i])
    G_epi[[names[i]]] <- delete.vertices(G, V(G)[is.na(vertex.attributes(G)[[attsel]])])
    ass[[names[i]]] <- assortativity(G_epi[[names[[i]]]], types1 = vertex.attributes(G_epi[[names[i]]])[[attsel]], directed = F)
  }
  
  return(ass)
}


# Calculate average abundance of ChIP-seq feature in fragment
# Calculate assortativity

ab_PCHiC_all <- colMeans(agchic[which(rownames(agchic) %in% V(PCHiC_wholenet)$name),])
ab_PCHiC_PP <- colMeans(agchic[which(rownames(agchic) %in% V(PCHiC_PP)$name),])
ab_PCHiC_POE <- colMeans(agchic[which(rownames(agchic) %in% V(PCHiC_POE)$name),])

ass_PCHiC_all <- calc_assort(PCHiC_wholenet, agchic)
ass_PCHiC_PP <- calc_assort(PCHiC_PP, agchic)
ass_PCHiC_POE <- calc_assort(PCHiC_POE, agchic)


# Colours for plots
nodecats <- read.table("/data/Projects/kat/Projects/Assortativity/epinet_nodescats.txt", sep="\t")

colnames(agchic) <- gsub("X5", "5", colnames(agchic))
colnames(agchic) <- gsub("Ezh2", "EZH2", colnames(agchic))
namesvec <- colnames(agchic)
namesvec <- gsub('\\.', '_', namesvec)

rownames(nodecats) <- nodecats[,1]
cats <- as.vector(unique(nodecats[,2]))
cols10 <- rainbow(11)
names(cols10) <- cats

cols10['Other'] <- 'grey'
cols4plot <- rep('grey', 78)
names(cols4plot) <- namesvec
for (c in cats){
  rel <- which(nodecats[,2]==c)
  cols4plot[rownames(nodecats)[rel]] <- cols10[c]
}


# Plots

# PCHiC general network ChAs vs abundance

labs <- c(1:ncol(agchic))
labs2 <- rep("", ncol(agchic))
sel <- which(ass_PCHiC_all > 0.05 | ass_PCHiC_all < 0 | ab_PCHiC_all > 0.03)
labs2[sel] <- namesvec[sel]
plot(ab_PCHiC_all, unlist(ass_PCHiC_all), pch = 20, xlim = c(0,0.25), col = cols4plot, main = "PCHi-C network", xlab = "Abundance", ylab = "ChAs")
text((ab_PCHiC_all),jitter(unlist(ass_PCHiC_all), 2), pos = 2, offset = 0.2, labels = labs, cex = 0.7)
text((ab_PCHiC_all),jitter(unlist(ass_PCHiC_all), 2), pos = 4, offset = 0.2, labels = labs2, cex = 0.7)
abline(h = 0)
legend("bottomright", legend = as.vector(unique(nodecats[,2])), bg = "white", pch = 20, col = cols10[as.vector(unique(nodecats[,2]))], cex = 0.7)


# PCHiC Promoter-Promoter network ChAs vs abundance

labs <- c(1:ncol(agchic))
labs2 <- rep("", ncol(agchic))
sel <- which(ass_PCHiC_PP > 0.05 | ass_PCHiC_PP < 0 | ab_PCHiC_PP > 0.03)
labs2[sel] <- namesvec[sel]
plot(ab_PCHiC_PP, unlist(ass_PCHiC_PP), pch = 20, xlim = c(0,0.3), col = cols4plot, main = "PCHi-C P-P network", xlab = "Abundance", ylab = "ChAs")
text((ab_PCHiC_PP),jitter(unlist(ass_PCHiC_PP), 2), pos = 2, offset = 0.2, labels = labs, cex = 0.7)
text((ab_PCHiC_PP),jitter(unlist(ass_PCHiC_PP), 2), pos = 4, offset = 0.2, labels = labs2, cex = 0.7)
abline(h = 0)
legend("bottomright", legend = as.vector(unique(nodecats[,2])), bg = "white", pch = 20, col = cols10[as.vector(unique(nodecats[,2]))], cex = 0.7)


# PCHiC Promoter-nonPromoter network ChAs vs abundance

labs <- c(1:ncol(agchic))
labs2 <- rep("", ncol(agchic))
sel <- which(ass_PCHiC_POE > 0.04 | ass_PCHiC_POE < 0 | ab_PCHiC_POE > 0.03)
labs2[sel] <- namesvec[sel]
plot(ab_PCHiC_POE, unlist(ass_PCHiC_POE), pch = 20, xlim = c(0,0.21), col = cols4plot, main = "PCHi-C P-OE network", xlab = "Abundance", ylab = "ChAs")
text((ab_PCHiC_POE),jitter(unlist(ass_PCHiC_POE), 2), pos = 2, offset = 0.2, labels = labs, cex = 0.7)
text((ab_PCHiC_POE),jitter(unlist(ass_PCHiC_POE), 2), pos = 4, offset = 0.2, labels = labs2, cex = 0.7)
abline(h = 0)
legend("bottomright", legend = as.vector(unique(nodecats[,2])), bg = "white", pch = 20, col = cols10[as.vector(unique(nodecats[,2]))], cex = 0.7)


# PCHiC PP vs POE ChAs

labs <- c(1:ncol(agchic))
labs2 <- rep("", ncol(agchic))
sel <- which(ass_PCHiC_POE > 0.04 | ass_PCHiC_POE < 0 | ass_PCHiC_PP > 0.03)
labs2[sel] <- namesvec[sel]
plot(unlist(ass_PCHiC_PP), unlist(ass_PCHiC_POE), pch = 20, xlim = c(0,0.45), col = cols4plot, main = "PCHi-C P-P vs PO ChAs", xlab = "P-P", ylab = "P-OE")
text(unlist(ass_PCHiC_PP),jitter(unlist(ass_PCHiC_POE), 2), pos = 2, offset = 0.2, labels = labs, cex = 0.7)
text(unlist(ass_PCHiC_PP),jitter(unlist(ass_PCHiC_POE), 2), pos = 4, offset = 0.2, labels = labs2, cex = 0.7)
abline(a = 0, b = 1, h = 0, v = 0)


legend("topleft", legend = as.vector(unique(nodecats[,2])), bg = "white", pch = 20, col = cols10[as.vector(unique(nodecats[,2]))], cex = 0.7)






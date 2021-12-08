rm(list = ls())

library(phytools)
library(geiger)
library(vioplot)
library(viridis)

## This function will calculate the tip ages of one or more trees and make a
## matrix with the species ages for each tree, the species name, the genus, and ## the class that they belong to
tip_ages <- function(tree, class) {

	if (class(tree) == "multiPhylo") {
		
		ages <- list()
		for (i in 1:length(tree)) {

			N <- length(tree[[i]]$tip.label)
			mat <- cbind(tree[[i]]$edge, tree[[i]]$edge.length)
			res <- mat[which(mat[, 2] < N+1), ]
			row.names(res) <- tree[[i]]$tip.label
			res <- sort(res[, 3], decreasing = T)
			ages[[i]] <- res

		}

		spp <- names(table(do.call(c, lapply(tree, function(x) x$tip.label))))
		ages_tab <- matrix(nrow = length(spp), 
		                   ncol = 4+length(tree))
		ages_tab[, 1] <- class
		ages_tab[, 2] <- NA
		ages_tab[, 4] <- spp

		for (i in 1:nrow(ages_tab)) {
			ages_tab[i, 3] <- strsplit(ages_tab[i, 4], "_")[[1]][1]

			for (j in 1:length(ages)) {
				ages_tab[i, 4+j] <- ifelse(length(ages[[j]][names(ages[[j]]) == 
														ages_tab[i, 4]]) == 0,
											NA, ages[[j]][names(ages[[j]]) == 
														ages_tab[i, 4]])
			}
		}

		return(ages_tab)

	} else {

		N <- length(tree$tip.label)
		mat <- cbind(tree$edge, tree$edge.length)
		ages <- mat[which(mat[, 2] < N+1), ]
		row.names(ages) <- tree$tip.label
		ages <- sort(ages[, 3], decreasing = T)

		ages_tab <- matrix(nrow = length(tree$tip.label), ncol = 5)
		ages_tab[, 1] <- class
		ages_tab[, 2] <- NA
		ages_tab[, 4] <- tree$tip.label

		for (i in 1:nrow(ages_tab)) {
			ages_tab[i, 3] <- strsplit(ages_tab[i, 4], "_")[[1]][1]
			ages_tab[i, 5] <- ifelse(length(ages[names(ages) == 
														ages_tab[i, 4]]) == 0,
											NA, ages[names(ages) == 
														ages_tab[i, 4]])
		}

		return(ages_tab)

	}

}

## Loading the phylogenies and taxonomies that were used to produce these 
## phylogenies
tr_am <- read.nexus("amphibia27JUL2020.nex")
for (i in 1:length(tr_am)) tr_am[[i]] <- drop.tip(tr_am[[i]], "Homo_sapiens")
tr_sq <- read.nexus("squamata27JUL20.nex")
tr_av <- read.nexus("aves_Ericson_27_JUL_20.nex")
tr_ma <- read.nexus("mammalia_node_dated_27JUL2020.nex")

dat_am <- read.csv("amph_shl_new_Classification.csv", stringsAsFactors = F)
dat_am$Scientific.Name <- gsub(" ", "_", dat_am$Scientific.Name)
dat_am <- dat_am[! grepl("Homo_sapiens", dat_am$Scientific.Name), ]

dat_sq <- read.csv("squam_shl_new_Classification.csv", 
                   stringsAsFactors = F)
dat_sq$Species <- gsub(" ", "_", dat_sq$Species)

dat_av <- read.csv("BLIOCPhyloMasterTax_birds.csv", stringsAsFactors = F)
dat_av$Scientific <- gsub(" ", "_", dat_av$Scientific)

dat_ma <- read.csv("taxonomy_mamPhy_5911species.csv", stringsAsFactors = F)


## Species ages for Amphibia

ages_amphibia <- tip_ages(tr_am, "Amphibia")
colnames(ages_amphibia) <- c("Class", "Taxon", "Genus", "Species", 
                             1:length(tr_am))

## Getting the orders of each species
for (i in 1:nrow(ages_amphibia)) {
	ages_amphibia[i, 2] <- dat_am$Taxon[dat_am$Scientific.Name == 
											ages_amphibia[i, 4]]
}

## Species ages for Squamata

Anguimorpha <- c("Anguidae", "Anniellidae", "Diploglossidae", "Helodermatidae",
                 "Lanthanotidae", "Shinisauridae", "Varanidae", "Xenosauridae")
Lacertoidea <- c("Amphisbaenidae", "Bipedidae", "Blanidae", "Cadeidae", 
                 "Gymnophthalmidae", "Lacertidae", "Rhineuridae", "Teiidae", 
                 "Trogonophiidae")
Gekkota <- c("Carphodactylidae", "Diplodactylidae", "Eublepharidae", 
             "Gekkonidae", "Phyllodactylidae", "Pygopodidae", 
             "Sphaerodactylidae")
Scincoidea <- c("Cordylidae", "Gerrhosauridae", "Scincidae", "Xantusiidae")
Iguania <- c("Agamidae", "Chamaeleonidae", "Corytophanidae", "Crotaphytidae", 
             "Dactyloidae", "Hoplocercidae", "Iguanidae", "Leiocephalidae", 
             "Leiosauridae", "Liolaemidae", "Opluridae", "Phrynosomatidae", 
             "Polychrotidae", "Tropiduridae")

ages_squamata <- tip_ages(tr_sq, "Reptilia")
colnames(ages_squamata) <- c("Class", "Taxon", "Genus", "Species", 
                             1:length(tr_sq))

## Getting the subclades of each species

for (i in 1:nrow(ages_squamata)) {
	ages_squamata[i, 2] <- ifelse(dat_sq$Family[dat_sq$Species == 
											ages_squamata[i, 4]] %in%
												Anguimorpha, 
												"Anguimorpha",
						   ifelse(dat_sq$Family[dat_sq$Species == 
											ages_squamata[i, 4]] %in%
												Lacertoidea, 
												"Lacertoidea",
						   ifelse(dat_sq$Family[dat_sq$Species == 
											ages_squamata[i, 4]] %in%
												Gekkota, 
												"Gekkota",	
						   ifelse(dat_sq$Family[dat_sq$Species == 
											ages_squamata[i, 4]] %in%
												Scincoidea, 
												"Scincoidea",
						   ifelse(dat_sq$Family[dat_sq$Species == 
											ages_squamata[i, 4]] %in%
												Iguania, 
												"Iguania",		
						   ifelse(dat_sq$Taxon[dat_sq$Species == 
											ages_squamata[i, 4]] ==
												"Serpentes", "Serpentes", NA
												))))))
}

## Species ages for Aves
ages_aves <- tip_ages(tr_av, "Aves")
colnames(ages_aves) <- c("Class", "Taxon", "Genus", "Species", 1:length(tr_av))

## Getting the orders of each species

firstup <- function(x) {
	x <- tolower(x)
	substr(x, 1, 1) <- toupper(substr(x, 1, 1))
	return(x)
}

for (i in 1:nrow(ages_aves)) {
	ages_aves[i, 2] <- firstup(dat_av$IOCOrder[dat_av$Scientific == 
															ages_aves[i, 4]])
}

## Species ages for Mammalia
ages_mammalia <- tip_ages(tr_ma, "Mammalia")
colnames(ages_mammalia) <- c("Class", "Taxon", "Genus", "Species", 
                             1:length(tr_ma))

## Getting the orders of each species

for (i in 1:nrow(ages_mammalia)) {
	ages_mammalia[i, 2] <- firstup(dat_ma$ord[dat_ma$Species_Name == 
														ages_mammalia[i, 4]])
}

ages <- rbind(ages_amphibia, ages_squamata, ages_aves, ages_mammalia)

write.csv(ages, "ages.csv", row.names = FALSE)

###############

## Repeting analyzes randomly pruning 1.25, 2.5, 5% of the tips of the
## phylogenies to consider incomplete sampling 

## Amphibia

## Getting 1.25, 2.5, 5% of the tips 
p1_am <- round(0.0125*length(tr_am[[1]]$tip))
p2_am <- round(0.025*length(tr_am[[1]]$tip))
p3_am <- round(0.05*length(tr_am[[1]]$tip))

## Dropping the tips that were not sampled
tr1_am <- tr2_am <- tr3_am <- list()
for (i in 1:length(tr_am)) {
	tr1_am[[i]] <- drop.tip(tr_am[[i]], sample(tr_am[[i]]$tip.label)[1:p1_am])
	tr2_am[[i]] <- drop.tip(tr_am[[i]], sample(tr_am[[i]]$tip.label)[1:p2_am])
	tr3_am[[i]] <- drop.tip(tr_am[[i]], sample(tr_am[[i]]$tip.label)[1:p3_am])
}
class(tr1_am) <- class(tr2_am) <- class(tr3_am) <-"multiPhylo"

## Calculating species ages for each pruned topology
ages1_amphibia <- tip_ages(tr1_am, "Amphibia")
colnames(ages1_amphibia) <- c("Class", "Taxon", "Genus", "Species", 
                              1:length(tr_am))

ages2_amphibia <- tip_ages(tr2_am, "Amphibia")
colnames(ages2_amphibia) <- c("Class", "Taxon", "Genus", "Species", 
                              1:length(tr_am))

ages3_amphibia <- tip_ages(tr3_am, "Amphibia")
colnames(ages3_amphibia) <- c("Class", "Taxon", "Genus", "Species", 
                              1:length(tr_am))


## Repeting for Squamata

p1_sq <- round(0.0125*length(tr_sq[[1]]$tip))
p2_sq <- round(0.025*length(tr_sq[[1]]$tip))
p3_sq <- round(0.05*length(tr_sq[[1]]$tip))

tr1_sq <- tr2_sq <- tr3_sq <- list()
for (i in 1:length(tr_sq)) {
	tr1_sq[[i]] <- drop.tip(tr_sq[[i]], sample(tr_sq[[i]]$tip.label)[1:p1_sq])
	tr2_sq[[i]] <- drop.tip(tr_sq[[i]], sample(tr_sq[[i]]$tip.label)[1:p2_sq])
	tr3_sq[[i]] <- drop.tip(tr_sq[[i]], sample(tr_sq[[i]]$tip.label)[1:p3_sq])
}
class(tr1_sq) <- class(tr2_sq) <- class(tr3_sq) <-"multiPhylo"

ages1_squamata <- tip_ages(tr1_sq, "Reptilia")
colnames(ages1_squamata) <- c("Class", "Taxon", "Genus", "Species", 
                              1:length(tr_sq))

ages2_squamata <- tip_ages(tr2_sq, "Reptilia")
colnames(ages2_squamata) <- c("Class", "Taxon", "Genus", "Species", 
                              1:length(tr_sq))

ages3_squamata <- tip_ages(tr3_sq, "Reptilia")
colnames(ages3_squamata) <- c("Class", "Taxon", "Genus", "Species", 
                              1:length(tr_sq))

## Aves

p1_av <- round(0.0125*length(tr_av[[1]]$tip))
p2_av <- round(0.025*length(tr_av[[1]]$tip))
p3_av <- round(0.05*length(tr_av[[1]]$tip))

tr1_av <- tr2_av <- tr3_av <- list()
for (i in 1:length(tr_av)) {
	tr1_av[[i]] <- drop.tip(tr_av[[i]], sample(tr_av[[i]]$tip.label)[1:p1_av])
	tr2_av[[i]] <- drop.tip(tr_av[[i]], sample(tr_av[[i]]$tip.label)[1:p2_av])
	tr3_av[[i]] <- drop.tip(tr_av[[i]], sample(tr_av[[i]]$tip.label)[1:p3_av])
}
class(tr1_av) <- class(tr2_av) <- class(tr3_av) <-"multiPhylo"

ages1_aves <- tip_ages(tr1_av, "Aves")
colnames(ages1_aves) <- c("Class", "Taxon", "Genus", "Species", 1:length(tr_av))

ages2_aves <- tip_ages(tr2_av, "Aves")
colnames(ages2_aves) <- c("Class", "Taxon", "Genus", "Species", 1:length(tr_av))

ages3_aves <- tip_ages(tr3_av, "Aves")
colnames(ages3_aves) <- c("Class", "Taxon", "Genus", "Species", 1:length(tr_av))

## Mammalia

p1_ma <- round(0.0125*length(tr_ma[[1]]$tip))
p2_ma <- round(0.025*length(tr_ma[[1]]$tip))
p3_ma <- round(0.05*length(tr_ma[[1]]$tip))

tr1_ma <- tr2_ma <- tr3_ma <- list()
for (i in 1:length(tr_ma)) {
	tr1_ma[[i]] <- drop.tip(tr_ma[[i]], sample(tr_ma[[i]]$tip.label)[1:p1_ma])
	tr2_ma[[i]] <- drop.tip(tr_ma[[i]], sample(tr_ma[[i]]$tip.label)[1:p2_ma])
	tr3_ma[[i]] <- drop.tip(tr_ma[[i]], sample(tr_ma[[i]]$tip.label)[1:p3_ma])
}
class(tr1_ma) <- class(tr2_ma) <- class(tr3_ma) <-"multiPhylo"

ages1_mammalia <- tip_ages(tr1_ma, "Mammalia")
colnames(ages1_mammalia) <- c("Class", "Taxon", "Genus", "Species", 
                              1:length(tr_ma))

ages2_mammalia <- tip_ages(tr2_ma, "Mammalia")
colnames(ages2_mammalia) <- c("Class", "Taxon", "Genus", "Species", 
                              1:length(tr_ma))

ages3_mammalia <- tip_ages(tr3_ma, "Mammalia")
colnames(ages3_mammalia) <- c("Class", "Taxon", "Genus", "Species", 
                              1:length(tr_ma))


ages1 <- rbind(ages1_amphibia, ages1_squamata, ages1_aves, ages1_mammalia)
write.csv(ages1, "ages_0.0125.csv", row.names = FALSE)

ages2 <- rbind(ages2_amphibia, ages2_squamata, ages2_aves, ages2_mammalia)
write.csv(ages2, "ages_0.025.csv", row.names = FALSE)

ages3 <- rbind(ages3_amphibia, ages3_squamata, ages3_aves, ages3_mammalia)
write.csv(ages3, "ages_0.05.csv", row.names = FALSE)

######################

## Figure 1 

ages <- read.csv("ages.csv")
ages1 <- read.csv("ages_0.0125.csv")
ages2 <- read.csv("ages_0.025.csv")
ages3 <- read.csv("ages_0.05.csv")

## This function will calculate the average species ages for each topology, 
## using log10 to transform the data or not. In this case, the average is taken
## for each class.
ageCalc_sim <- function(dataset, name, log = TRUE) {
	ntree <- 4+length(tr_am)
	data <- as.data.frame(subset(dataset, dataset[, 1] == name, c(1, 5:ntree)))
	data[, 2:ncol(data)] <- apply(data[, 2:ncol(data)], 2, as.numeric)

	if (log == TRUE) {
		meanAge <- numeric()
		taxa <- character()
		for (i in 1:length(tr_am)) {
			meanAge[i] <- mean(log10(t(data[, 1+i])), na.rm = T)
			taxa[i] <- name
		}
		res <- cbind(meanAge, taxa)
		res
	} else {
		meanAge <- numeric()
		taxa <- character()
		for (i in 1:length(tr_am)) {
			meanAge[i] <- mean(t(data[, 1+i]), na.rm = T)
			taxa[i] <- name
		}
		res <- cbind(meanAge, taxa)
		res
	}
}

## Average species ages for the phylogenies not pruned
res0 <- rbind(
	ageCalc_sim(ages, "Amphibia", log = FALSE),
	ageCalc_sim(ages, "Reptilia", log = FALSE),
	ageCalc_sim(ages, "Aves", log = FALSE),
	ageCalc_sim(ages, "Mammalia", log = FALSE)
)

## Average species ages for the pruned phylogenies
res1 <- rbind(
	ageCalc_sim(ages1, "Amphibia", log = FALSE),
	ageCalc_sim(ages1, "Reptilia", log = FALSE),
	ageCalc_sim(ages1, "Aves", log = FALSE),
	ageCalc_sim(ages1, "Mammalia", log = FALSE)
)

res2 <- rbind(
	ageCalc_sim(ages2, "Amphibia", log = FALSE),
	ageCalc_sim(ages2, "Reptilia", log = FALSE),
	ageCalc_sim(ages2, "Aves", log = FALSE),
	ageCalc_sim(ages2, "Mammalia", log = FALSE)
)

res3 <- rbind(
	ageCalc_sim(ages3, "Amphibia", log = FALSE),
	ageCalc_sim(ages3, "Reptilia", log = FALSE),
	ageCalc_sim(ages3, "Aves", log = FALSE),
	ageCalc_sim(ages3, "Mammalia", log = FALSE)
)

res0 <- as.data.frame(res0)
res1 <- as.data.frame(res1)
res2 <- as.data.frame(res2)
res3 <- as.data.frame(res3)

## Plotting the violin plots 

pdf("Figure1.pdf", height = 5)

cols_alpha <- c(rgb(48/255, 18/255, 59/255, 0.4),
                rgb(26/255, 228/255, 182/255, 0.4),
                rgb(250/255, 186/255, 57/255, 0.4),
                rgb(122/255, 4/255, 3/255, 0.4))

cols <- turbo(4)


par(mar = c(4, 8, 2, 3))

vioplot(as.numeric(res0[res0$taxa=="Amphibia", 1]),
	   as.numeric(res0[res0$taxa=="Reptilia", 1]),
	   as.numeric(res0[res0$taxa=="Aves", 1]),
	   as.numeric(res0[res0$taxa=="Mammalia", 1]),
	   col = cols_alpha, horizontal = T, xaxt = "n", border = cols, ylog = TRUE)

vioplot(as.numeric(res1[res1$taxa=="Amphibia", 1]),
	   as.numeric(res1[res1$taxa=="Reptilia", 1]),
	   as.numeric(res1[res1$taxa=="Aves", 1]),
	   as.numeric(res1[res1$taxa=="Mammalia", 1]), 
	   add = T, col = cols_alpha, horizontal = T, border = cols, ylog = TRUE)

vioplot(as.numeric(res2[res2$taxa=="Amphibia", 1]),
	   as.numeric(res2[res2$taxa=="Reptilia", 1]),
	   as.numeric(res2[res2$taxa=="Aves", 1]),
	   as.numeric(res2[res2$taxa=="Mammalia", 1]), 
	   add = T, col = cols_alpha, horizontal = T, border = cols, ylog = TRUE)

vioplot(as.numeric(res3[res3$taxa=="Amphibia", 1]),
	   as.numeric(res3[res3$taxa=="Reptilia", 1]),
	   as.numeric(res3[res3$taxa=="Aves", 1]),
	   as.numeric(res3[res3$taxa=="Mammalia", 1]), 
	   add = T, col = cols_alpha, horizontal = T, border = cols, ylog = TRUE)
axis(2, at = c(1:4), labels = c("Amphibia", "Squamata", "Aves",
                                "Mammalia"), las = 2)

dev.off()


## Figure 2

## This function will also calculate the average species ages, but now for 
## different orders/subclades
ageCalc <- function(dataset, name, log = TRUE) {
	data <- dataset[dataset[, 2] == name, ]
	data <- as.data.frame(data[complete.cases(data), ])
	ntree <- 4+length(tr_am)
	data[, 5:ntree] <- apply(data[, 5:ntree], 2, as.numeric)
	if (log == TRUE) {
		meanAge <- numeric()
		taxa <- character()
		for (i in 1:length(tr_am)) {
			meanAge[i] <- mean(log10(t(data[, 4+i])))
			taxa[i] <- name
		}
		res <- cbind(meanAge, taxa)
		res
	} else {
		meanAge <- numeric()
		taxa <- character()
		for (i in 1:length(tr_am)) {
			meanAge[i] <- mean(t(data[, 4+i]))
			taxa[i] <- name
		}
		res <- cbind(meanAge, taxa)
		res
	}
}

res <- rbind(
	ageCalc(ages, "Carnivora", log = FALSE),
	ageCalc(ages, "Cetartiodactyla", log = FALSE),
	ageCalc(ages, "Chiroptera", log = FALSE),
	ageCalc(ages, "Diprotodontia", log = FALSE),
	ageCalc(ages, "Primates", log = FALSE),
	ageCalc(ages, "Columbiformes", log = FALSE),
	ageCalc(ages, "Passeriformes", log = FALSE),
	ageCalc(ages, "Piciformes", log = FALSE),
	ageCalc(ages, "Psittaciformes", log = FALSE),
	ageCalc(ages, "Anguimorpha", log = FALSE),
	ageCalc(ages, "Gekkota", log = FALSE),
	ageCalc(ages, "Iguania", log = FALSE),
	ageCalc(ages, "Lacertoidea", log = FALSE),
	ageCalc(ages, "Scincoidea", log = FALSE),
	ageCalc(ages, "Serpentes", log = FALSE),
	ageCalc(ages, "Anura", log = FALSE),
	ageCalc(ages, "Caudata", log = FALSE),
	ageCalc(ages, "Gymnophiona", log = FALSE)
)

res <- as.data.frame(res)

pdf("Figure2.pdf", height = 6, width=7.5)

par(mar = c(4, 8, 2, 3))

cols_alpha <- c(rgb(48/255, 18/255, 59/255, 0.8),
                rgb(26/255, 228/255, 182/255, 0.8),
                rgb(250/255, 186/255, 57/255, 0.8),
                rgb(122/255, 4/255, 3/255, 0.8))

cols <- turbo(4)

vioplot(as.numeric(res[res$taxa == "Gymnophiona", 1]),
	   as.numeric(res[res$taxa == "Caudata", 1]),
	   as.numeric(res[res$taxa == "Anura", 1]), 
	   as.numeric(res[res$taxa == "Serpentes", 1]),
	   as.numeric(res[res$taxa == "Scincoidea", 1]),
	   as.numeric(res[res$taxa == "Lacertoidea", 1]),
	   as.numeric(res[res$taxa == "Iguania", 1]),
	   as.numeric(res[res$taxa == "Gekkota", 1]),
	   as.numeric(res[res$taxa == "Anguimorpha", 1]),
	   as.numeric(res[res$taxa == "Psittaciformes", 1]),
	   as.numeric(res[res$taxa ==	"Piciformes", 1]),
	   as.numeric(res[res$taxa == "Passeriformes", 1]),
	   as.numeric(res[res$taxa == "Columbiformes", 1]),
	   as.numeric(res[res$taxa == "Primates", 1]),
	   as.numeric(res[res$taxa == "Diprotodontia", 1]),
	   as.numeric(res[res$taxa == "Chiroptera", 1]),
	   as.numeric(res[res$taxa == "Cetartiodactyla", 1]),
	   as.numeric(res[res$taxa == "Carnivora", 1]),
	   col = c(cols_alpha[1], cols_alpha[1], cols_alpha[1],
	   		   cols_alpha[2], cols_alpha[2], cols_alpha[2], cols_alpha[2], 
	   		   cols_alpha[2], cols_alpha[2],
	           cols_alpha[3], cols_alpha[3], cols_alpha[3], cols_alpha[3],
	           cols_alpha[4], cols_alpha[4], cols_alpha[4], cols_alpha[4], 
	           cols_alpha[4]), 
	   border = c(cols[1], cols[1], cols[1],
	   		   	  cols[2], cols[2], cols[2], cols[2], cols[2], cols[2],
	           	  cols[3], cols[3], cols[3], cols[3],
	           	  cols[4], cols[4], cols[4], cols[4], cols[4]),
	   horizontal = T, xaxt = "n", ylog = TRUE)
axis(2, at = c(1:18), labels = c("Gymnophiona", "Caudata", "Anura", 
                                 "Serpentes", "Scincoidea", "Lacertoidea",
                                 "Iguania", "Gekkota", "Anguimorpha",
                                 "Psittaciformes", "Piciformes",
                                 "Passeriformes", "Columbiformes", 
                                 "Primates", "Diprotodontia", "Chiroptera",
                                 "Cetartiodactyla", "Carnivora"), las = 2)

dev.off()

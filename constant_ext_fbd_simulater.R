library(devtools)
library(TreeSim)
library(phytools)
library(phangorn)
###starting params###
#####################
morpho_rate <- 0.01
molec_rate <- 0.006
seq_length <- 2000
seq_length_molec <- 1000
fossilisation_rate <- 0.05
speciation <- 1
extinction <- 0.9

n_extant_tips <- 75
threshold <- 5

sim_fossils_on_tree_known <- function(fossilisation_rate, speciation, extinction, n_extant_tips, n_reps, morpho_rate, molec_rate, seq_length, seq_length_molec){
a=1
###RUN THE SIMULATION###
########################
number <- vector(mode="numeric", length=0)
while(length(number) < threshold){

###simulate tree###
tree <- sim.bd.taxa(n_extant_tips, 1, speciation, extinction, frac = 1, complete = TRUE, stochsampling = FALSE)[[1]] # the overall tree simulated from the origin
crown_group_tree <- extract.clade(tree, findMRCA(tree, drop.tip(tree, getExtinct(tree))$tip.label, "node")) # the tree starting from the mrca of extant species, the node we are interested in how fossils might potentially calibrate
fossil_addition_tree <- crown_group_tree # the tree on which fossil tips are sequentially added
crown_group_tree_extant <- drop.tip(crown_group_tree, getExtinct(crown_group_tree)) # get extant tips in the crown group tree
fossil_addition_tree_extant <- crown_group_tree_extant

###sim fossils on full tree###
##############################

###simulate fossil times - we do with the memoryless distribution as per the explanation set out###
table <- cbind(crown_group_tree[[1]], crown_group_tree$edge.length) # node information in crown group tree with edge lengths
fossil_vector <- vector("list", nrow(table))
fossil_names <- vector("list", nrow(table))
for (b in 1:nrow(table)){ # fossils simulated for each branch in the tree
continue = "yes" 
fossil_vector[[b]] <- append(fossil_vector[[b]], 0) # initiate the vector
while (continue == "yes"){
fossil_vector[[b]] <- append(fossil_vector[[b]], fossil_vector[[b]][[length(fossil_vector[[b]])]] + rexp(1, fossilisation_rate)) # add next fossil sampling time according to a given rate
if (fossil_vector[[b]][[length(fossil_vector[[b]])]] >= table[,3][[b]]){ # if the new fossil sampling time is greater than the length of the relevent branch you stop adding more fossils
continue = "no"
}
}
fossil_vector[[b]] <- fossil_vector[[b]][-which(fossil_vector[[b]] == 0)] # get rid of vector initiator
fossil_vector[[b]] <- fossil_vector[[b]][-which(fossil_vector[[b]] > table[,3][[b]])] # get rid of over lengthed fossils

###if there are fossils on a branch, create names for them and then add them
if (length(fossil_vector[[b]]) > 0){ 
fossil_names[[b]] <- paste(rep(paste("fossil_", b, sep=""), length(fossil_vector[[b]])), seq(1, length(fossil_vector[[b]]), 1), sep="_")
if (table[,2][[b]] > length(crown_group_tree$tip.label)){ # if its not a terminal branch you are adding the fossil tip to
for (c in 1:length(fossil_vector[[b]])){
fossil_addition_tree <- bind.tip(fossil_addition_tree, fossil_names[[b]][[c]], edge.length=0, where=findMRCA(fossil_addition_tree, extract.clade(crown_group_tree, table[,2][[b]])$tip.label, "node"), position=table[,3][[b]] - fossil_vector[[b]][[c]])
}
}
if (table[,2][[b]] <= length(crown_group_tree$tip.label)){
for (c in 1:length(fossil_vector[[b]])){
fossil_addition_tree <- bind.tip(fossil_addition_tree, fossil_names[[b]][[c]], edge.length=0, where=which(fossil_addition_tree$tip.label == crown_group_tree$tip.label[[table[,2][[b]]]]), position=table[,3][[b]] - fossil_vector[[b]][[c]])
}
}
}
}
###write the outputs###
#######################

fossil_tips <- unlist(fossil_names)
if (length(fossil_tips) > 15){
dir.create(as.character(paste("input_", a, sep="")))
###write full tree
write.tree(fossil_addition_tree, paste("input_", a, "/", "full_", a, ".tre", sep=""))
###write fbd tree
fossil_addition_tree_no_ext <- drop.tip(fossil_addition_tree, which(fossil_addition_tree$tip.label %in% c(crown_group_tree_extant$tip.label, fossil_tips) == FALSE))
write.tree(fossil_addition_tree_no_ext, paste("input_", a, "/", "fbd_", a, ".tre", sep=""))
###write sequences molecular
fossil_addition_tree_mol_rate_var <- fossil_addition_tree
for (i in 1:length(fossil_addition_tree_mol_rate_var$edge.length)){
fossil_addition_tree_mol_rate_var$edge.length[[i]] <- fossil_addition_tree_mol_rate_var$edge.length[[i]]*rlnorm(1, log(1) - ((0.2937025^2)/2), 0.2937025)
}
mol_sequences <- simSeq(fossil_addition_tree_mol_rate_var, seq_length_molec, type = "DNA", rate = molec_rate)
extant_mol_sequences <- subset(mol_sequences, drop.tip(fossil_addition_tree, getExtinct(fossil_addition_tree))$tip.label)
write.phyDat(extant_mol_sequences, paste("input_", a, "/", a, "_molec_sequences.fasta", sep=""), format="nexus")
write.tree(fossil_addition_tree_mol_rate_var, paste("input_", a, "/", a, "_rate.tre", sep=""))
###write sequences morpho
###all
sequences <- simSeq(fossil_addition_tree, seq_length, type = "USER", levels = c(0,1), rate = morpho_rate) # simulate morphological characters on the fossil addition tree
extant_and_foss_sequences <- subset(sequences, fossil_addition_tree_no_ext$tip.label)
write.phyDat(extant_and_foss_sequences, paste("input_", a, "/", a, "_morpho_sequences.fasta", sep=""), format="nexus")
###40 chars
forty_sites <- subset(extant_and_foss_sequences, select=sample(seq(1,2000,1), 40, replace=FALSE))
write.phyDat(forty_sites, paste("input_", a, "/", a, "_morpho_sequences_forty.fasta", sep=""), format="nexus")
###8 chars
eight_sites <- subset(extant_and_foss_sequences, select=sample(seq(1,2000,1), 8, replace=FALSE))
write.phyDat(eight_sites, paste("input_", a, "/", a, "_morpho_sequences_eight.fasta", sep=""), format="nexus")
###write fossil ages
tip_times <- round(max(node.depth.edgelength(fossil_addition_tree_no_ext)) - node.depth.edgelength(fossil_addition_tree_no_ext)[seq(1,length(fossil_addition_tree_no_ext$tip.label),1)],6)
taxon <- fossil_addition_tree_no_ext$tip.label
min_age <- tip_times
max_age <- tip_times
tip_ages <- data.frame(taxon, min_age, max_age)
write.table(tip_ages, paste("input_", a, "/", a, "_sampling_time.csv", sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
message(paste("success", "on rep ", a, sep=""))
###oldest fossil
oldest_fossil_age <- max(node.depth.edgelength(fossil_addition_tree_no_ext)) - min(node.depth.edgelength(fossil_addition_tree_no_ext)[seq(1, length(fossil_addition_tree_no_ext$tip.label), 1)])
write(paste("oldest_fossil <- ", oldest_fossil_age, sep=""), paste("input_", a, "/oldest_fossil_", a, ".txt", sep=""))
number <- append(number , 1)
}
a=a+1
}
}


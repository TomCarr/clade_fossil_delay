library(devtools)
library(TreeSim)
library(phytools)

###starting params###
#####################

fossilisation_rate <- 0.0001
speciation <- 1
extinction <- 0.9

n_extant_tips <- 10
threshold <- 2000

sim_fossils_on_tree_known <- function(fossilisation_rate, speciation, extinction, n_extant_tips, n_reps){

###storage_variables_for_plotting###
####################################

a=1
absolute_fossil_ages_storage <<- vector(mode="numeric", length=0)
absolute_node_ages_storage <<- vector(mode="numeric", length=0)
absolute_fossil_ages_storage_extant <<- vector(mode="numeric", length=0)
absolute_node_ages_storage_extant <<- vector(mode="numeric", length=0)

###RUN THE SIMULATION###
########################

while((length(absolute_fossil_ages_storage_extant[which(absolute_fossil_ages_storage_extant == "break")]) < threshold) || (length(absolute_fossil_ages_storage[which(absolute_fossil_ages_storage == "break")]) < threshold)){

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

fossil_addition_tree <<- fossil_addition_tree

###get fossils in each clade

fossil_tips <- unlist(fossil_names)
if (length(fossil_tips) > 0){ # if we've actually simulated fossils
for (b in 1:(length(crown_group_tree_extant$tip.label)-1)){ # go through each node
calibrated_clade_extant <- extract.clade(crown_group_tree_extant, length(crown_group_tree_extant$tip.label)+b)
calibrated_clade <- extract.clade(fossil_addition_tree, findMRCA(fossil_addition_tree, calibrated_clade_extant$tip.label, "node"))
calibrated_clade_fossil_tips <- calibrated_clade$tip.label[grep("fossil_", calibrated_clade$tip.label)]
if (length(calibrated_clade_fossil_tips) > 0){
subset_clade_fossils <- vector(mode="numeric", length=0)
if (length(calibrated_clade_extant$tip.label) > 2){
subset_clade_one <- extract.clade(crown_group_tree_extant, max(crown_group_tree_extant[[1]][,2][which(crown_group_tree_extant[[1]][,1] == length(crown_group_tree_extant$tip.label)+b)]))$tip.label
subset_clade_one_all <- extract.clade(fossil_addition_tree, findMRCA(fossil_addition_tree, subset_clade_one, "node"))$tip.label
subset_clade_fossils <- subset_clade_one_all[grep("fossil_", subset_clade_one_all)]
if (length(subset_clade_one) < (length(calibrated_clade_extant$tip.label)-1)){ 
subset_clade_two <- calibrated_clade_extant$tip.label[-which(calibrated_clade_extant$tip.label %in% subset_clade_one)]
subset_clade_two_all <- extract.clade(fossil_addition_tree, findMRCA(fossil_addition_tree, subset_clade_two, "node"))$tip.label
subset_clade_fossils <- append(subset_clade_fossils, subset_clade_two_all[grep("fossil_", subset_clade_two_all)])
}
}
if (length(subset_clade_fossils) > 0){
calibrated_clade_fossil_tips <- calibrated_clade_fossil_tips[-which(calibrated_clade_fossil_tips %in% subset_clade_fossils)]
}
###final processing of remaining fossils
if (length(calibrated_clade_fossil_tips) > 0){
all_fossil_ages <- max(node.depth.edgelength(calibrated_clade)) - node.depth.edgelength(calibrated_clade)[which(calibrated_clade$tip.label %in% calibrated_clade$tip.label[grep("fossil_", calibrated_clade$tip.label)])] ### basically refedine original calibrated clade fossil tips
main_clade_fossil_ages <- max(node.depth.edgelength(calibrated_clade)) - node.depth.edgelength(calibrated_clade)[which(calibrated_clade$tip.label %in% calibrated_clade_fossil_tips)]
oldest_fossil_all <- max(all_fossil_ages)
oldest_fossil_main <- max(main_clade_fossil_ages)
if (oldest_fossil_all == oldest_fossil_main){
absolute_fossil_ages_storage <<- append(absolute_fossil_ages_storage, oldest_fossil_main)
absolute_node_ages_storage <<- append(absolute_node_ages_storage, max(node.depth.edgelength(calibrated_clade)))
}
}
}
}
write.tree(fossil_addition_tree, paste(a, ".tre", sep=""))
absolute_fossil_ages_storage <<- append(absolute_fossil_ages_storage, "break")
absolute_node_ages_storage <<- append(absolute_node_ages_storage, "break")
message(paste("up to ",  length(absolute_fossil_ages_storage[which(absolute_fossil_ages_storage == "break")]), "on full sim", collapse=""))
}

###sim fossils on extant tree###
################################

###simulate fossil times - we do with the memoryless distribution as per the explanation set out###
table <- cbind(crown_group_tree_extant[[1]], crown_group_tree_extant$edge.length) # node information in crown group tree with edge lengths
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
if (table[,2][[b]] > length(crown_group_tree_extant$tip.label)){ # if its not a terminal branch you are adding the fossil tip to
for (c in 1:length(fossil_vector[[b]])){
fossil_addition_tree_extant <- bind.tip(fossil_addition_tree_extant, fossil_names[[b]][[c]], edge.length=0, where=findMRCA(fossil_addition_tree_extant, extract.clade(crown_group_tree_extant, table[,2][[b]])$tip.label, "node"), position=table[,3][[b]] - fossil_vector[[b]][[c]])
}
}
if (table[,2][[b]] <= length(crown_group_tree_extant$tip.label)){
for (c in 1:length(fossil_vector[[b]])){
fossil_addition_tree_extant <- bind.tip(fossil_addition_tree_extant, fossil_names[[b]][[c]], edge.length=0, where=which(fossil_addition_tree_extant$tip.label == crown_group_tree_extant$tip.label[[table[,2][[b]]]]), position=table[,3][[b]] - fossil_vector[[b]][[c]])
}
}
}
}

fossil_addition_tree_extant <<- fossil_addition_tree_extant

###get fossils in each clade

fossil_tips <- unlist(fossil_names)
if (length(fossil_tips) > 0){
for (b in 1:(length(crown_group_tree_extant$tip.label)-1)){
calibrated_clade_extant <- extract.clade(crown_group_tree_extant, length(crown_group_tree_extant$tip.label)+b)
calibrated_clade <- extract.clade(fossil_addition_tree_extant, findMRCA(fossil_addition_tree_extant, calibrated_clade_extant$tip.label, "node"))
calibrated_clade_fossil_tips <- calibrated_clade$tip.label[grep("fossil_", calibrated_clade$tip.label)]
if (length(calibrated_clade_fossil_tips) > 0){
subset_clade_fossils <- vector(mode="numeric", length=0)
if (length(calibrated_clade_extant$tip.label) > 2){
subset_clade_one <- extract.clade(crown_group_tree_extant, max(crown_group_tree_extant[[1]][,2][which(crown_group_tree_extant[[1]][,1] == length(crown_group_tree_extant$tip.label)+b)]))$tip.label
subset_clade_one_all <- extract.clade(fossil_addition_tree_extant, findMRCA(fossil_addition_tree_extant, subset_clade_one, "node"))$tip.label
subset_clade_fossils <- subset_clade_one_all[grep("fossil_", subset_clade_one_all)]
if (length(subset_clade_one) < (length(calibrated_clade_extant$tip.label)-1)){ 
subset_clade_two <- calibrated_clade_extant$tip.label[-which(calibrated_clade_extant$tip.label %in% subset_clade_one)]
subset_clade_two_all <- extract.clade(fossil_addition_tree_extant, findMRCA(fossil_addition_tree_extant, subset_clade_two, "node"))$tip.label
subset_clade_fossils <- append(subset_clade_fossils, subset_clade_two_all[grep("fossil_", subset_clade_two_all)])
}
}
if (length(subset_clade_fossils) > 0){
calibrated_clade_fossil_tips <- calibrated_clade_fossil_tips[-which(calibrated_clade_fossil_tips %in% subset_clade_fossils)]
}
###final processing of remaining fossils
if (length(calibrated_clade_fossil_tips) > 0){
all_fossil_ages <- max(node.depth.edgelength(calibrated_clade)) - node.depth.edgelength(calibrated_clade)[which(calibrated_clade$tip.label %in% calibrated_clade$tip.label[grep("fossil_", calibrated_clade$tip.label)])] ### basically refedine original calibrated clade fossil tips
main_clade_fossil_ages <- max(node.depth.edgelength(calibrated_clade)) - node.depth.edgelength(calibrated_clade)[which(calibrated_clade$tip.label %in% calibrated_clade_fossil_tips)]
oldest_fossil_all <- max(all_fossil_ages)
oldest_fossil_main <- max(main_clade_fossil_ages)
if (oldest_fossil_all == oldest_fossil_main){
absolute_fossil_ages_storage_extant <<- append(absolute_fossil_ages_storage_extant, oldest_fossil_main)
absolute_node_ages_storage_extant <<- append(absolute_node_ages_storage_extant, max(node.depth.edgelength(calibrated_clade)))
}
}
}
}
write.tree(fossil_addition_tree_extant, paste(a, "extant.tre", sep=""))
absolute_fossil_ages_storage_extant <<- append(absolute_fossil_ages_storage_extant, "break")
absolute_node_ages_storage_extant <<- append(absolute_node_ages_storage_extant, "break")
message(paste("up to ",  length(absolute_fossil_ages_storage_extant[which(absolute_fossil_ages_storage_extant == "break")]), "on extant sim", collapse=""))
}
a=a+1
}
write(absolute_fossil_ages_storage, "fossil_ages.txt", ncolumns=1)
write(absolute_node_ages_storage, "node_ages.txt", ncolumns=1)
###
write(absolute_fossil_ages_storage_extant, "fossil_ages_extant.txt", ncolumns=1)
write(absolute_node_ages_storage_extant, "node_ages_extant.txt", ncolumns=1)
write(a, "how_many_sims.txt", ncolumns=1)
###
}

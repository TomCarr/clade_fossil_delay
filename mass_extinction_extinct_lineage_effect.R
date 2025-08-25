library(devtools)
library(TreeSim)
library(phytools)

###starting params###
#####################

fossilisation_rate <- 0.0001   
speciation <- 0.1
extinction <- 0

n_extant_tips <- 75
threshold <- 2500

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

while((length(absolute_fossil_ages_storage_extant) < threshold) || (length(absolute_fossil_ages_storage) < threshold)){
#message(paste("up to ",  length(absolute_fossil_ages_storage_extant), " actual fossils on extant", collapse=""))
#message(paste("simulation", a, sep=" "))

###simulate tree###
tree <- sim.bd.taxa(n_extant_tips, 1, speciation, extinction, frac = 1, complete = TRUE, stochsampling = FALSE)[[1]] # the overall tree simulated from the origin
write.tree(tree, "doing_this.tre")

###do mass extinction###
########################

###select the branches to extinct
table <- cbind(tree[[1]], tree$edge.length)
table <- cbind(table[order(table[,2]),], c(rep(0, length(tree$tip.label)), branching.times(tree)[-1])) 
table <- cbind(table, table[,4] + table[,3])
mass_extinction_time <- 0.5*max(node.depth.edgelength(tree))
mass_extinction_probability <- 0.90
extinction_affected_branch <- vector(mode="numeric", length=0)
ext_names <- vector(mode="numeric", length=0)
for (b in 1:nrow(table)){
if ((table[,4][[b]] < mass_extinction_time) & (table[,5][[b]] > mass_extinction_time)){
if (runif(1, 0, 1) <= mass_extinction_probability){
extinction_affected_branch <- append(extinction_affected_branch, b)
ext_names <- append(ext_names, paste(b, "ext", sep=""))
}
}
}

###add on the branches that will remain and remove others

ext_tree <- tree

name_usage <- vector(mode="numeric", length=0)
for (b in extinction_affected_branch){
name_usage <- append(name_usage, 1)
if (table[,2][[b]] > length(tree$tip.label)){ # if its not a terminal branch you are adding the fossil tip to
ext_tree <- bind.tip(ext_tree, ext_names[[length(name_usage)]], edge.length=0, where=findMRCA(ext_tree, extract.clade(tree, table[,2][[b]])$tip.label, "node"), position=mass_extinction_time - max(node.depth.edgelength(extract.clade(tree, table[,2][[b]]))))
ext_tree <- drop.tip(ext_tree, extract.clade(tree, table[,2][[b]])$tip.label)
}
if (table[,2][[b]] <= length(tree$tip.label)){
ext_tree <- bind.tip(ext_tree, ext_names[[length(name_usage)]], edge.length=0, where=which(ext_tree$tip.label == tree$tip.label[[table[,2][[b]]]]), position=mass_extinction_time)
ext_tree <- drop.tip(ext_tree, which(ext_tree$tip.label == tree$tip.label[[table[,2][[b]]]]))
}
}

tree <- ext_tree
extant_tree <- drop.tip(tree, grep("ext", tree$tip.label))
if (length(tree) > 1){
if (length(grep("ext", tree$tip.label)) < length(tree$tip.label) - 2){
if (max(node.depth.edgelength(extant_tree)) > 0.5*max(node.depth.edgelength(tree))){



######################################

crown_group_tree <- extract.clade(tree, findMRCA(tree, drop.tip(tree, grep("ext", tree$tip.label))$tip.label, "node")) # the tree starting from the mrca of extant species, the node we are interested in how fossils might potentially calibrate
fossil_addition_tree <- crown_group_tree # the tree on which fossil tips are sequentially added
crown_group_tree_extant <- drop.tip(crown_group_tree, grep("ext", crown_group_tree$tip.label)) # get extant tips in the crown group tree
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
#message(paste("adding fossils called", fossil_names[[b]], " to tree", collapse="   "))
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
#message(paste("simulated a ", n_extant_tips, " tip tree with ", length(unlist(fossil_names)), " fossils", sep=""))

fossil_addition_tree <<- fossil_addition_tree

###get fossils just in overall clade

fossil_tips <- unlist(fossil_names)

subset_clade_one <- extract.clade(crown_group_tree_extant, (length(crown_group_tree_extant$tip.label)+2))$tip.label
subset_clade_one_all <- extract.clade(fossil_addition_tree, findMRCA(fossil_addition_tree, subset_clade_one, "node"))$tip.label
subset_clade_fossils <- subset_clade_one_all[grep("fossil_", subset_clade_one_all)]
if (length(subset_clade_one) < (length(crown_group_tree_extant$tip.label)-1)){ 
subset_clade_two <- crown_group_tree_extant$tip.label[-which(crown_group_tree_extant$tip.label %in% subset_clade_one)]
subset_clade_two_all <- extract.clade(fossil_addition_tree, findMRCA(fossil_addition_tree, subset_clade_two, "node"))$tip.label
subset_clade_fossils <- append(subset_clade_fossils, subset_clade_two_all[grep("fossil_", subset_clade_two_all)])
}

if (length(fossil_tips) > 0){
if (length(subset_clade_fossils) > 0){
fossil_tips <- fossil_tips[-which(fossil_tips %in% subset_clade_fossils)]
}
}

###final processing of remaining fossils
if (length(fossil_tips) > 0){
all_fossil_ages <- max(node.depth.edgelength(fossil_addition_tree)) - node.depth.edgelength(fossil_addition_tree)[seq(1, length(fossil_addition_tree$tip.label), 1)] 
main_clade_fossil_ages <- max(node.depth.edgelength(fossil_addition_tree)) - node.depth.edgelength(fossil_addition_tree)[which(fossil_addition_tree$tip.label %in% fossil_tips)] ###is the ordering correct
oldest_fossil_all <- max(all_fossil_ages)
oldest_fossil_main <- max(main_clade_fossil_ages)
if (oldest_fossil_all == oldest_fossil_main){
absolute_fossil_ages_storage <<- append(absolute_fossil_ages_storage, oldest_fossil_main)
absolute_node_ages_storage <<- append(absolute_node_ages_storage, max(node.depth.edgelength(crown_group_tree)))
write.tree(fossil_addition_tree, paste(a, ".tre", sep=""))
message(paste("success", "on rep ", a, sep=""))
} else {
#absolute_fossil_ages_storage <<- append(absolute_fossil_ages_storage, "blank")
#absolute_node_ages_storage <<- append(absolute_node_ages_storage, max(node.depth.edgelength(crown_group_tree)))
}
} else {
#absolute_fossil_ages_storage <<- append(absolute_fossil_ages_storage, "blank")
#absolute_node_ages_storage <<- append(absolute_node_ages_storage, max(node.depth.edgelength(crown_group_tree)))
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
#message(paste("adding fossils called", fossil_names[[b]], " to tree", collapse="   "))
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
#message(paste("simulated a ", n_extant_tips, "extant only tip tree with ", length(unlist(fossil_names)), " fossils", sep=""))

fossil_addition_tree_extant <<- fossil_addition_tree_extant

###get fossils just in overall clade

fossil_tips <- unlist(fossil_names)

subset_clade_one <- extract.clade(crown_group_tree_extant, (length(crown_group_tree_extant$tip.label)+2))$tip.label
subset_clade_one_all <- extract.clade(fossil_addition_tree_extant, findMRCA(fossil_addition_tree_extant, subset_clade_one, "node"))$tip.label
subset_clade_fossils <- subset_clade_one_all[grep("fossil_", subset_clade_one_all)]
if (length(subset_clade_one) < (length(crown_group_tree_extant$tip.label)-1)){ 
subset_clade_two <- crown_group_tree_extant$tip.label[-which(crown_group_tree_extant$tip.label %in% subset_clade_one)]
subset_clade_two_all <- extract.clade(fossil_addition_tree_extant, findMRCA(fossil_addition_tree_extant, subset_clade_two, "node"))$tip.label
subset_clade_fossils <- append(subset_clade_fossils, subset_clade_two_all[grep("fossil_", subset_clade_two_all)])
}

if (length(fossil_tips) > 0){
if (length(subset_clade_fossils) > 0){
fossil_tips <- fossil_tips[-which(fossil_tips %in% subset_clade_fossils)]
}
}

###final processing of remaining fossils
if (length(fossil_tips) > 0){
all_fossil_ages <- max(node.depth.edgelength(fossil_addition_tree_extant)) - node.depth.edgelength(fossil_addition_tree_extant)[seq(1, length(fossil_addition_tree_extant$tip.label), 1)] 
main_clade_fossil_ages <- max(node.depth.edgelength(fossil_addition_tree_extant)) - node.depth.edgelength(fossil_addition_tree_extant)[which(fossil_addition_tree_extant$tip.label %in% fossil_tips)] ###is the ordering correct
oldest_fossil_all <- max(all_fossil_ages)
oldest_fossil_main <- max(main_clade_fossil_ages)
if (oldest_fossil_all == oldest_fossil_main){
absolute_fossil_ages_storage_extant <<- append(absolute_fossil_ages_storage_extant, oldest_fossil_main)
absolute_node_ages_storage_extant <<- append(absolute_node_ages_storage_extant, max(node.depth.edgelength(crown_group_tree)))
write.tree(fossil_addition_tree_extant, paste(a, "extant.tre", sep=""))
} else {
#absolute_fossil_ages_storage <<- append(absolute_fossil_ages_storage, "blank")
#absolute_node_ages_storage <<- append(absolute_node_ages_storage, max(node.depth.edgelength(crown_group_tree)))
}
} else {
#absolute_fossil_ages_storage <<- append(absolute_fossil_ages_storage, "blank")
#absolute_node_ages_storage <<- append(absolute_node_ages_storage, max(node.depth.edgelength(crown_group_tree)))
}
a=a+1
}
}
}
}
write(absolute_fossil_ages_storage, "fossil_ages.txt", ncolumns=1)
write(absolute_node_ages_storage, "node_ages.txt", ncolumns=1)
###
write(absolute_fossil_ages_storage_extant, "fossil_ages_extant.txt", ncolumns=1)
write(absolute_node_ages_storage_extant, "node_ages_extant.txt", ncolumns=1)
write(a, "how_many_sims.txt", ncolumns=1)
###
}


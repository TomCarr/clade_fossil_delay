library(TreeSim)
library(phytools)
library(phangorn)

###starting params###
#####################

fossilisation_rate <- 0.01
speciation <- 0.1
extinction <- 0
morpho_rate <- 0.002

n_extant_tips <- 25
seq_length <- 40

threshold <- 500

sim_fossils_on_tree <- function(fossilisation_rate, speciation, extinction, seq_length, morpho_rate, n_extant_tips){
counter=1
a=1
###RUN THE SIMULATION###
########################

while(a < threshold){
clade_age_for_things_with_synaps_anc <- vector(mode="numeric", length=0)
clade_size_for_things_with_synaps_anc <- vector(mode="numeric", length=0)
stem_branch_storer_anc <- vector(mode="numeric", length=0)
clade_synapomorphy_storer_number_anc <- vector(mode="numeric", length=0)
number_of_clades_with_synapomorphies_storer_anc <- vector(mode="numeric", length=0)
###
clade_age_for_things_with_synaps <- vector(mode="numeric", length=0)
clade_size_for_things_with_synaps <- vector(mode="numeric", length=0)
stem_branch_storer <- vector(mode="numeric", length=0)
clade_synapomorphy_storer_number <- vector(mode="numeric", length=0)
number_of_clades_with_synapomorphies_storer <- vector(mode="numeric", length=0)
###
absolute_fossil_ages_storage <- vector(mode="numeric", length=0)
absolute_node_ages_storage <- vector(mode="numeric", length=0)
absolute_fossil_ages_storage_filtered <- vector(mode="numeric", length=0)
absolute_node_ages_storage_filtered <- vector(mode="numeric", length=0)
root_age_storage <- vector(mode="numeric", length=0)

message(paste("simulation", a, sep=" "))

###simulate tree
tree <- sim.bd.taxa(n_extant_tips, 1, speciation, extinction, frac = 1, complete = TRUE, stochsampling = FALSE)[[1]] # the overall tree simulated from the origin
crown_group_tree <- extract.clade(tree, findMRCA(tree, drop.tip(tree, getExtinct(tree))$tip.label, "node")) # the tree starting from the mrca of extant species, the node we are interested in how fossils might potentially calibrate
fossil_addition_tree <- crown_group_tree # the tree on which fossil tips are sequentially added

###simulate fossil times
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
if (length(fossil_vector[[b]]) > 0){ # if there are fossils on a branch, create names for them and then add them
fossil_names[[b]] <- paste(rep(paste("fossil_", b, sep=""), length(fossil_vector[[b]])), seq(1, length(fossil_vector[[b]]), 1), sep="_")

###add fossils to tree
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
message(paste("simulated a ", n_extant_tips, " tip tree with ", length(unlist(fossil_names)), " fossils", sep=""))
###simulate sequences###
########################
message("now simulating sequences")
sequences <- simSeq(fossil_addition_tree, seq_length, type = "USER", levels = c(0,1), rate = morpho_rate, ancestral=TRUE) # simulate morphological characters on the fossil addition tree

extant_species <- drop.tip(crown_group_tree, getExtinct(crown_group_tree))$tip.label # get extant tips in the crown group tree
extant_tree <- keep.tip(crown_group_tree, extant_species) # get the extant only crown group tree 

extant_sequences <- subset(sequences, extant_tree$tip.label) # just the sequences of extant tips in crown group tree
fossil_sequences <- subset(sequences, unlist(fossil_names)) # just the sequences of the fossils
ancestral_sequence <- subset(sequences, as.character(length(fossil_addition_tree$tip.label)+1))

###search for synapomorphies in extant_tree###
##############################################
if (length(fossil_sequences) > 0){

message("now getting clade synapomorphies in extant tree")

###initial check with ancestral sequences###
############################################
clade_synapomorphy_counter_anc <- vector(mode="numeric", length=0)
for (b in which(extant_tree[[1]][,2] > length(extant_tree$tip.label))){ # going through all nodes, excluding terminals, and the root node
clade_tips <- extract.clade(extant_tree, extant_tree[[1]][,2][[b]])$tip.label
clade_nodes <- getDescendants(fossil_addition_tree, findMRCA(fossil_addition_tree, extract.clade(extant_tree, extant_tree[[1]][,2][[b]])$tip.label, "node"), curr=TRUE)
clade_nodes <- as.character(clade_nodes[-which(clade_nodes <= length(fossil_addition_tree$tip.label))]) 
clade_tips <- c(clade_tips, clade_nodes)

non_clade_tips <- extant_tree$tip.label[-which(extant_tree$tip.label %in% clade_tips)]# get tips of everything else
non_clade_nodes <- as.character(seq(length(fossil_addition_tree$tip.label) + 1, (length(fossil_addition_tree$tip.label)*2)-1, 1))
non_clade_nodes <- non_clade_nodes[-which(non_clade_nodes %in% clade_nodes)]
non_clade_tips <- c(non_clade_tips, non_clade_nodes)

clade_sequences <- subset(sequences, clade_tips) # get sequences of the clade
non_clade_sequences <- subset(sequences, non_clade_tips) # get sequences of everything else
clade_synapomorphies <- vector("list", length=0)

for (c in 1:seq_length){ # loop progressing along sites, going to search for synapomorphies for clade at each site 
if (length(unique(as.numeric(unlist(subset(extant_sequences, select = c, site.pattern=FALSE))))) != 1){ #if its not an invariant site with respect to all extant sequences
clade_site_states <- as.numeric(unlist(subset(clade_sequences, select = c, site.pattern=FALSE))) # get states at that site within the clade of interest
most_common_at_site <- sort(table(clade_site_states),decreasing=TRUE)[[1]] # get the frequency of the most common state at that site within the clade of interest
most_common_at_site_state <- as.numeric(sort(table(clade_site_states),decreasing=TRUE)[[1]])
if (most_common_at_site > 0.9*length(clade_tips)){ #if the highest frequency state is in 90% of the tips...
if (most_common_at_site_state != as.numeric(unlist(subset(ancestral_sequence, select = c, site.pattern=FALSE)))){
pre_synapomorphy <- as.numeric(names(sort(table(clade_site_states),decreasing=TRUE)))[[1]] # assign the state to pre-synapomorphy status
non_clade_site_states <- as.numeric(unlist(subset(non_clade_sequences, select = c, site.pattern=FALSE))) # get states at the site outside the clade of interest 
if ((pre_synapomorphy %in% non_clade_site_states) == FALSE){ # see if its outside the clade too, if not we call it a synapomorphy
clade_synapomorphies[[length(clade_synapomorphies)+1]] <- c(c, pre_synapomorphy) # get a lits of all the synapomorphies for the relevent clade
}
}
}
}
}
if (length(clade_synapomorphies) > 0){
clade_age_for_things_with_synaps_anc <- append(clade_age_for_things_with_synaps_anc, max(node.depth.edgelength(extract.clade(extant_tree, extant_tree[[1]][,2][[b]]))))
clade_size_for_things_with_synaps_anc <- append(clade_size_for_things_with_synaps_anc, length(extract.clade(extant_tree, extant_tree[[1]][,2][[b]])$tip.label))
stem_branch_storer_anc <- append(stem_branch_storer_anc, max(node.depth.edgelength(extract.clade(extant_tree, extant_tree[[1]][,1][[b]]))) - max(node.depth.edgelength(extract.clade(extant_tree, extant_tree[[1]][,2][[b]]))))
clade_synapomorphy_storer_number_anc <- append(clade_synapomorphy_storer_number_anc, length(clade_synapomorphies))
clade_synapomorphy_counter_anc <- append(clade_synapomorphy_counter_anc, 1)
} else {
clade_age_for_things_with_synaps_anc <- append(clade_age_for_things_with_synaps_anc, "-")
clade_size_for_things_with_synaps_anc <- append(clade_size_for_things_with_synaps_anc, "-")
stem_branch_storer_anc <- append(stem_branch_storer_anc, "-")
clade_synapomorphy_storer_number_anc <- append(clade_synapomorphy_storer_number_anc, length(clade_synapomorphies))
}
}
number_of_clades_with_synapomorphies_storer_anc <- append(number_of_clades_with_synapomorphies_storer_anc, length(clade_synapomorphy_counter_anc))

###proper###
############
clade_synapomorphy_counter <- vector(mode="numeric", length=0)
fossil_in_clade <- vector("list", length(fossil_sequences))
for (b in which(extant_tree[[1]][,2] > length(extant_tree$tip.label))){ # going through all nodes, excluding terminals, and the root node
clade_tips <- extract.clade(extant_tree, extant_tree[[1]][,2][[b]])$tip.label # get tips of the relevent clade
non_clade_tips <- extant_tree$tip.label[-which(extant_tree$tip.label %in% clade_tips)] # get tips of everything else
clade_sequences <- subset(extant_sequences, clade_tips) # get sequences of the clade
non_clade_sequences <- subset(extant_sequences, non_clade_tips) # get sequences of everything else
clade_synapomorphies <- vector("list", length=0)

#####
for (c in 1:seq_length){ # loop progressing along sites, going to search for synapomorphies for clade at each site 
if (length(unique(as.numeric(unlist(subset(extant_sequences, select = c, site.pattern=FALSE))))) != 1){ #if its not an invariant site with respect to all extant sequences
clade_site_states <- as.numeric(unlist(subset(clade_sequences, select = c, site.pattern=FALSE))) # get states at that site within the clade of interest
most_common_at_site <- sort(table(clade_site_states),decreasing=TRUE)[[1]] # get the frequency of the most common state at that site within the clade of interest
most_common_at_site_state <- as.numeric(sort(table(clade_site_states),decreasing=TRUE)[[1]])
if (most_common_at_site > 0.9*length(clade_tips)){ #if the highest frequency state is in 90% of the tips...
if (most_common_at_site_state != as.numeric(unlist(subset(ancestral_sequence, select = c, site.pattern=FALSE)))){
pre_synapomorphy <- as.numeric(names(sort(table(clade_site_states),decreasing=TRUE)))[[1]] # assign the state to pre-synapomorphy status
non_clade_site_states <- as.numeric(unlist(subset(non_clade_sequences, select = c, site.pattern=FALSE))) # get states at the site outside the clade of interest 
if ((pre_synapomorphy %in% non_clade_site_states) == FALSE){ # see if its outside the clade too, if not we call it a synapomorphy
clade_synapomorphies[[length(clade_synapomorphies)+1]] <- c(c, pre_synapomorphy) # get a lits of all the synapomorphies for the relevent clade
}
}
}
}
}

###
if (length(clade_synapomorphies) > 0){
clade_age_for_things_with_synaps <- append(clade_age_for_things_with_synaps, max(node.depth.edgelength(extract.clade(extant_tree, extant_tree[[1]][,2][[b]]))))
clade_size_for_things_with_synaps <- append(clade_size_for_things_with_synaps, length(extract.clade(extant_tree, extant_tree[[1]][,2][[b]])$tip.label))
stem_branch_storer <- append(stem_branch_storer, max(node.depth.edgelength(extract.clade(extant_tree, extant_tree[[1]][,1][[b]]))) - max(node.depth.edgelength(extract.clade(extant_tree, extant_tree[[1]][,2][[b]]))))
clade_synapomorphy_storer_number <- append(clade_synapomorphy_storer_number, length(clade_synapomorphies))
clade_synapomorphy_counter <- append(clade_synapomorphy_counter, 1)
} else {
clade_age_for_things_with_synaps <- append(clade_age_for_things_with_synaps, "-")
clade_size_for_things_with_synaps <- append(clade_size_for_things_with_synaps, "-")
stem_branch_storer <- append(stem_branch_storer, "-")
clade_synapomorphy_storer_number <- append(clade_synapomorphy_storer_number, length(clade_synapomorphies))
}
number_of_clades_with_synapomorphies_storer <- append(number_of_clades_with_synapomorphies_storer, length(clade_synapomorphy_counter))

###find which fossils can go to the clade###
############################################

if (length(clade_synapomorphies) > 0){

message(paste("assigining fossils to clade", extant_tree[[1]][,2][[b]], sep=""))

for (c in 1:length(unlist(fossil_names))){
test_sequence <- subset(sequences, unlist(fossil_names)[[c]]) # get the sequence for the fossil we are testing
fossil_site_vector <- vector(mode="numeric", length=0)
for (d in 1:length(clade_synapomorphies)){# have to do with this loop because of site pattern problems
fossil_site_vector <- append(fossil_site_vector, as.numeric(unlist(subset(test_sequence, select = clade_synapomorphies[[d]][[1]], site.pattern=FALSE))))
}
if (sum(fossil_site_vector == sapply(clade_synapomorphies, "[[", 2)) > 0.9*length(clade_synapomorphies)){ # if the fossil has 90% of the synapomorphies
fossil_in_clade[[c]] <- append(fossil_in_clade[[c]], extant_tree[[1]][,2][[b]]) #then the fossil belongs to the clade we are investigating
}
}
}
}
if (length(unlist(fossil_in_clade)) > 0){

###now get oldest fossil in each clade###							###fossil in clade is a list thats the lengths of all fossils
#########################################

message("getting oldest fossil in each clade")

removal <- which(lengths(fossil_in_clade) == 0)
if (length(removal) > 0){
fossil_in_clade <- sapply(fossil_in_clade[-removal], max)			
final_fossil_names <- unlist(fossil_names)[-removal]
} else {															###nb confused about names - resolved - dont remove max because fossil in clade refers to all the clades that it belongs to - max would keep the most specific clade
fossil_in_clade <- sapply(fossil_in_clade, max)
final_fossil_names <- unlist(fossil_names)
}

absolute_fossil_ages <- vector(mode="numeric", length=0) 			### get the absolute age of each fossil that has been assigned to a clade in the tree
for (b in 1:length(final_fossil_names)){
absolute_fossil_ages <- append(absolute_fossil_ages, max(node.depth.edgelength(fossil_addition_tree)) - node.depth.edgelength(fossil_addition_tree)[which(fossil_addition_tree$tip.label == final_fossil_names[[b]])]) #get the absolute age of each of the fossils
}

fossil_in_clade_removal <- vector(mode="numeric", length=0) 		### this removes fossils if there is an older fossil at the same node or in descendants
for (b in 1:length(fossil_in_clade)){
if (length(fossil_in_clade) > 1){
for (c in (1:length(fossil_in_clade))[-b]){
if ((fossil_in_clade[[c]] == fossil_in_clade[[b]]) & (absolute_fossil_ages[[c]] > absolute_fossil_ages[[b]])){
fossil_in_clade_removal <- append(fossil_in_clade_removal, b)
} else if ((fossil_in_clade[[c]] %in% getDescendants(extant_tree, fossil_in_clade[[b]])) & (absolute_fossil_ages[[c]] >= absolute_fossil_ages[[b]])){
fossil_in_clade_removal <- append(fossil_in_clade_removal, b)
}
}
}
}
final_fossil_names <- final_fossil_names[-unique(fossil_in_clade_removal)]
fossil_in_clade <- fossil_in_clade[-unique(fossil_in_clade_removal)]
absolute_fossil_ages <- absolute_fossil_ages[-unique(fossil_in_clade_removal)]

node_ages <- vector(mode="numeric", length=0)
node_numbers <- vector(mode="numeric", length=0)
if (length(absolute_fossil_ages) > 0){
for (b in 1:length(absolute_fossil_ages)){
node_ages <- append(node_ages, max(node.depth.edgelength(extract.clade(extant_tree, Ancestors(extant_tree, fossil_in_clade[[b]])[[1]])))) 
node_numbers <- append(node_numbers, Ancestors(extant_tree, fossil_in_clade[[b]])[[1]])
}
#################################################
#################################################
absolute_fossil_ages_storage <- append(absolute_fossil_ages_storage, c(absolute_fossil_ages, "break"))
absolute_node_ages_storage <- append(absolute_node_ages_storage, c(node_ages, "break"))
root_age_storage <- append(root_age_storage, max(node.depth.edgelength(crown_group_tree))) 
write.tree(fossil_addition_tree, paste(counter, ".tre", sep=""))
write(absolute_fossil_ages_storage, "fossil_ages.txt", ncolumns=1,append=TRUE)
write(absolute_node_ages_storage, "node_ages.txt", ncolumns=1,append=TRUE)
write(root_age_storage, "root_ages.txt", ncolumns=1,append=TRUE)
###filter because things are put on stem node###
################################################
if (length(unique(node_numbers)) != length(node_numbers)){
for (b in 1:length(unique(node_numbers))){
if (length(which(node_numbers == unique(node_numbers)[[b]])) > 1){
fossil_ages_for_node <- absolute_fossil_ages[which(node_numbers == unique(node_numbers[[b]]))]
oldest_fossil_for_node <- max(fossil_ages_for_node)
absolute_fossil_ages <- absolute_fossil_ages[-which(node_numbers == unique(node_numbers)[[b]])[-which(fossil_ages_for_node == oldest_fossil_for_node)]]
node_ages <- node_ages[-which(node_numbers == unique(node_numbers)[[b]])[-which(fossil_ages_for_node == oldest_fossil_for_node)]]
}
}
}
################################################
################################################
absolute_fossil_ages_storage_filtered <- append(absolute_fossil_ages_storage_filtered, c(absolute_fossil_ages, "break"))
absolute_node_ages_storage_filtered <- append(absolute_node_ages_storage_filtered, c(node_ages, "break"))
write(absolute_fossil_ages_storage_filtered, "fossil_ages_filtered.txt", ncolumns=1,append=TRUE)
write(absolute_node_ages_storage_filtered, "node_ages_filtered.txt", ncolumns=1,append=TRUE)
a=a+1
}
}
###write###
###########
number_of_clades_with_synapomorphies_storer_anc <- as.data.frame(as.list(number_of_clades_with_synapomorphies_storer_anc))
write.table(number_of_clades_with_synapomorphies_storer_anc, "number_of_clades_with_synapomorphies_anc.txt", row.names=FALSE, col.names=FALSE, append=TRUE)
###
clade_age_for_things_with_synaps_anc <- as.data.frame(as.list(clade_age_for_things_with_synaps_anc))
write.table(clade_age_for_things_with_synaps_anc, "clade_age_for_things_with_synaps_anc.txt", row.names=FALSE, col.names=FALSE, append=TRUE)
###
clade_size_for_things_with_synaps_anc <- as.data.frame(as.list(clade_size_for_things_with_synaps_anc))
write.table(clade_size_for_things_with_synaps_anc, "clade_size_for_things_with_synaps_anc.txt", row.names=FALSE, col.names=FALSE, append=TRUE)
###
stem_branch_storer_anc <- as.data.frame(as.list(stem_branch_storer_anc))
write.table(stem_branch_storer_anc, "stem_branch_storer_anc.txt", row.names=FALSE, col.names=FALSE, append=TRUE)
###
clade_synapomorphy_storer_number_anc <- as.data.frame(as.list(clade_synapomorphy_storer_number_anc))
write.table(clade_synapomorphy_storer_number_anc, "clade_synapomorphy_storer_number_anc.txt", row.names=FALSE, col.names=FALSE, append=TRUE)
#########
#########
number_of_clades_with_synapomorphies_storer <- as.data.frame(as.list(number_of_clades_with_synapomorphies_storer))
write.table(number_of_clades_with_synapomorphies_storer, "number_of_clades_with_synapomorphies.txt", row.names=FALSE, col.names=FALSE, append=TRUE)
###
clade_age_for_things_with_synaps <- as.data.frame(as.list(clade_age_for_things_with_synaps))
write.table(clade_age_for_things_with_synaps, "clade_age_for_things_with_synaps.txt", row.names=FALSE, col.names=FALSE, append=TRUE)
###
clade_size_for_things_with_synaps <- as.data.frame(as.list(clade_size_for_things_with_synaps))
write.table(clade_size_for_things_with_synaps, "clade_size_for_things_with_synaps.txt", row.names=FALSE, col.names=FALSE, append=TRUE)
###
stem_branch_storer <- as.data.frame(as.list(stem_branch_storer))
write.table(stem_branch_storer, "stem_branch_storer.txt", row.names=FALSE, col.names=FALSE, append=TRUE)
###
clade_synapomorphy_storer_number <- as.data.frame(as.list(clade_synapomorphy_storer_number))
write.table(clade_synapomorphy_storer_number, "clade_synapomorphy_storer_number.txt", row.names=FALSE, col.names=FALSE, append=TRUE)
}
counter <- counter+1
}
}



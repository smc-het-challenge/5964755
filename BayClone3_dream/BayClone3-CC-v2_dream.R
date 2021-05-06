rm(list=ls(all=TRUE))
library(combinat)
library(DPpackage)
set.seed(79861)

args <- commandArgs(TRUE); 

#### reading from the input file
CHR_NAME_INDX = 1
CHR_POS_INDX = 2
TOTAL_RD_INDX = 3
VARIANT_RD_INDX = 4
CN_INDX = 5


input.file <- args[[1]]

out1A_file <- "subchallenge1A.txt"
out1B_file <- "subchallenge1B.txt"
out1C_file <- "subchallenge1C.txt"
#out3A_file <- "subchallenge3A.txt"

#out3_file <- args[[5]]

#curr_dir = cat(getwd())

#### for BayClone3
R_fn_BayClone3 = sprintf("/opt/BayClone3/fn_BayClone3-CC-v2_dream.R")
source(R_fn_BayClone3)
so_file_BayClone3 = sprintf("/opt/BayClone3/BayClone3-CC-v2_dream.so")
dyn.load(so_file_BayClone3)
#### for trees
tree_generate_file = sprintf("/opt/BayClone3/RobustFunctionsV5_ForICGC_v4_dream.R")
source(tree_generate_file)

msg = sprintf("R code started with %s !!\n",input.file);
print(paste0(msg))

dat <- as.matrix(read.table(input.file, header=FALSE))
###################################
chrName <- as.matrix(dat[,CHR_NAME_INDX])  
chrPos <- as.matrix(dat[,CHR_POS_INDX])
N <- as.matrix(as.numeric(dat[,TOTAL_RD_INDX]))  
n <- as.matrix(as.numeric(dat[,VARIANT_RD_INDX]))
M_B <- as.matrix(as.numeric(dat[,CN_INDX]))

S <- nrow(N) # THE NUMBER OF SNP (I.E. NUMBER OF ROWS OF Z)
T <- 1 #THE NUMBER OF TISSUES

total_SNV <- S
##  read in data...  Here!!!!!
### we need N, n, M_B, S, T=1

#######################################

####################################
###  pre-processing
## remove extreme values
tmp_dat <- fn.preprocessing(N, M_B, n, S)
N <- as.matrix(tmp_dat$N)
n <- tmp_dat$n
M_B <- tmp_dat$M_B
S <- tmp_dat$S
removed.ind <- tmp_dat$removed_ind ### binary indicator of length original S, if removed, then 1. not, 0

####################################
###  Create indicator for M_B
## if (ind_B[i_s] == 1) {  // if BB output is not integer
Ind_B <- rep(1, S) #array(1, dim=c(S, T))
ave_M_B <- apply(M_B, 1, mean)  ## ave(M_B) over t is integer
Ind_B[ave_M_B==ceiling(M_B)] <- 0

###################################
#HYPER-PARAMETER
ind.LM <- 0 ###  if we want to use LM to estimate w_star
ind.DP <- 1 ###  if we want to use DP to estimate w_star
###  if both are 1, then use the average of the two estimates.
hyper <- fn.elicit.hyper.B(ceiling(max(M_B)), N, n, S, T, M_B, ind.LM, ind.DP)



################################################################################
#FIT THE MODEL
set.seed(1415)
#### MCMC parameters
#n.sam <- 4000; burn.in <- 3000

n.sam <- 3000; burn.in <- 2000

Min_C <- 2
Max_C <- 11

## the last element is the indicator;
ind.fixed <- 0  ## random w0 and p0 ==> ind.fixed == 0 vs fixed w0 and p0 ==> ind.fixed == 1
MCMC.sam <- BayClone2_B(Ind_B, M_B, n, N, S, T, hyper, 0.9, Min_C, Max_C, burn.in, n.sam, ind.fixed)


### remove subclones with negligible weights
#fn.BayClone2.post.process <- function(MCMC_sam, SS, TT, n_sam, ind_w0_fixed, w0_fixed, p0_fixed)
MCMC.sam.post <- fn.BayClone2.post.process.B(MCMC.sam, S, T, n.sam, 0.03, ind.fixed, hyper$w0, hyper$p0_z)

set.seed(79861)
####SUMMARY FOR FIGURES######################
###  finding a point estimate of L, Z, w
point.est <- fn.BayClone2.point.B(MCMC.sam.post, S, T, Min_C, Max_C, n.sam, ind.fixed)
if(point.est$C_star > 1)
{
    #### subchallenge 1A for dream-smc-het
    cellularity = 1.0-mean(point.est$w_star_dist)
    fid1A <- file(out1A_file, "w")
    cat(cellularity,file=fid1A)
    close(fid1A)
    
    #### For this we need to run tree algorithm
    treeResult = ConstructTree(SM = point.est$Z, CM = point.est$L, w = point.est$w, prunTh = 0.1)
    treeResult$tree = calibrateWandIndex(treeResult$tree)
    
    save(treeResult,file="./treeResult.rda")
    
    #### subchallenge 1B for dream-smc-het
    nNodes = length(treeResult$tree$Nodes)-1
    fid1B <- file(out1B_file, "w")
    cat(nNodes,file=fid1B)
    close(fid1B)
    
    #### subchallenge 1C for dream-smc-het
        
    child_cluster_str_set = sapply(strsplit(names(treeResult$tree$Edges),split = "->"),function(x)x[2])
    child_cluster_id_set = as.integer(sapply(strsplit(child_cluster_str_set,split = "_"),function(x)x[2]))
         
    total_SNV_from_tree = 0
    #     
    #     for (cluster_indx in 1:length(child_cluster_id_set))
    #     {
    #         cluster_tree_id = which(child_cluster_id_set == cluster_indx)
    #         n_ssms = length(treeResult$tree$Edges[[cluster_tree_id]][grep("SNV", treeResult$tree$Edges[[cluster_tree_id]])])
    #         total_SNV_from_tree  = total_SNV_from_tree + n_ssms
    #         w = treeResult$tree$w[cluster_tree_id+1] #### +1 because number of nodes = number of edges+1
    #         cat(cluster_indx,"\t",n_ssms,"\t",w,"\n", file=fid1C)
    #     }    
    
    
    nNodes = length(treeResult$tree$Nodes)-1
    nEdges = length(treeResult$tree$Edges)
    #### we need nEdges number of lists
    SNV_list_in_edge = vector(mode = "list", length = nEdges)
    
    #### take care of all Z entries 0 part
    tree_input_SNV_set = 1:length(which(removed.ind==0))
    tree_output_SNV_set = NULL
    for(cluster_indx in 1:nEdges)
    {
        cluster_tree_id = which(child_cluster_id_set == cluster_indx)
        SNV_CNV_str = sapply(strsplit(treeResult$tree$Edges[[cluster_tree_id]], split = ":"), function(x)x[1])
        
        SNV_str = SNV_CNV_str[grep("SNV",SNV_CNV_str)]
        SNV_id = sapply(strsplit(SNV_str, split = "_"), function(x)x[2])
        SNV_list_in_edge[cluster_indx] = list(as.integer(SNV_id))
        tree_output_SNV_set = c(tree_output_SNV_set,as.integer(SNV_id))
    }
    tree_output_SNV_set = unique(tree_output_SNV_set)
    input_output_set_diff = setdiff(tree_input_SNV_set,tree_output_SNV_set)
    SS_or <- length(removed.ind)  ## SS before removing
    #cat("SS_or = ",SS_or,"\n") #### debug
    #aux_indx = vector(mode="numeric",length=SS_or)
    aux_indx = rep(0,length=SS_or)
    aux_indx[which(removed.ind == 0)] = 1:length(which(removed.ind == 0))
    
    mut_assign <- rep(NA, SS_or)
    
    total_SNV_from_tree = 0
    
    #### cluster_nssms computes the cluster size which sums to total SNVs in VCF
    #### initialize to 0 
    cluster_nssms = vector(mode='integer',length=nEdges)
    mult_pr = c(rep(1,times=nEdges))/nEdges
    
    for(j in 1:nNodes)
    {    
        total_SNV_from_tree = total_SNV_from_tree+length(SNV_list_in_edge[[j]])
        #cat("cluster ",j," size = ",length(SNV_list_in_edge[[j]]),"\n")
    }
    removed_SNV_earlier = 0
    SNV_with_all_Z_zero = 0

    normal_clone_nssms = 0

    total_SNV_in_tree = 0
    for (i in 1:SS_or)
    {
        ### either it was removed earlier 
        if( aux_indx[i] == 0) #### write NA
        {
            #cat(chrName[i],"\t",chrPos[i],"\t","NA\n")
            removed_SNV_earlier = removed_SNV_earlier+1
            k  = which(rmultinom(1,1,mult_pr) == 1)
            cluster_nssms[k] = cluster_nssms[k] + 1
        } else {
                if (aux_indx[i] %in% input_output_set_diff)
                    ## or it was not included in tree output
                {
                    SNV_with_all_Z_zero = SNV_with_all_Z_zero + 1
                    #k  = which(rmultinom(1,1,mult_pr) == 1)
                    #cluster_nssms[k] = cluster_nssms[k] + 1
                    normal_clone_nssms = normal_clone_nssms + 1
                } else { 
                        #### get the id that tree algorithm is using
                        id = aux_indx[i]
                        nFound = 0 ### number of times it is found
                        cluster_id_list= vector(mode='integer')
                        total_SNV_in_tree = total_SNV_in_tree + 1
                        for(j in 1:nNodes)
                        {
                            if(id %in% SNV_list_in_edge[[j]])  #### this SNV is in that cluster
                            {
                                nFound = nFound+1
                                cluster_id_list[nFound] = j
                            }    
                        }
                        if(nFound == 1)
                        {
                            rand_cluster_id = cluster_id_list[nFound]
                            cluster_nssms[rand_cluster_id] = cluster_nssms[rand_cluster_id] + 1  
                        
                        }
                        else ## it is > 1 and that means repetition
                        {
                            pr_found = c(rep(1,times=nFound))/nFound
                            #cat("nFound = ",nFound,"\n")
                            cat(pr_found,"\n")
                            k = which(rmultinom(1,1,pr_found) == 1)
                            rand_cluster_id = cluster_id_list[k]
                            cluster_nssms[rand_cluster_id] = cluster_nssms[rand_cluster_id]  + 1
                        }
                }
        }    
        
    }
    #### write the result to subchallenge 1C
    fid1C <- file(out1C_file, "w")
    
    for (cluster_indx in 1:length(child_cluster_id_set))
    {
        cluster_tree_id = which(child_cluster_id_set == cluster_indx)
        n_ssms = cluster_nssms[cluster_indx]    
        w = treeResult$tree$w[cluster_tree_id+1]
        cat(cluster_indx,"\t",n_ssms,"\t",w,"\n", file=fid1C)
    }
    if(normal_clone_nssms > 0)
    {
         cat(length(child_cluster_id_set)+1,"\t",normal_clone_nssms,"\t",0.000000,"\n", file=fid1C)
    }
    close(fid1C)
    
    ##### subchallenge 3A
    #fid3 <- file(out3A_file, "w")
    
    #parent_list = sapply(strsplit(names(treeResult$tree$Edges),split = "->"),function(x)x[1])
    #child_list  = sapply(strsplit(names(treeResult$tree$Edges),split = "->"),function(x)x[2])
    #parent_id = sapply(strsplit(parent_list, split="_"),function(x)x[2])
    #child_id  = sapply(strsplit(child_list, split="_"),function(x)x[2])
    #for (i in 1:length(parent_id))
    #{
    #    cat(child_id[i],"\t",parent_id[i],"\n",file=fid3)
    #}
    #close(fid3)

} else { 
    print(paste0("Number of tumor-subclones is 0 with no result file\n"))
}

#print(paste0("total_SNV = ",total_SNV,"\n"))
#print(paste0("total_SNV_from_tree = ", total_SNV_from_tree,"\t","removed_SNV_earlier = ", removed_SNV_earlier,"\t","SNV_with_all_Z_zero = ",SNV_with_all_Z_zero,"\n"))
#print(paste0("whole total = ",total_SNV_from_tree+removed_SNV_earlier+SNV_with_all_Z_zero,"\n"))
### done

#print(paste0("total_SNV_in_tree = ",total_SNV_in_tree,"\n"))

L = length(cluster_nssms)
total_size = 0
for(j in 1:L)
{    
    cat("final cluster ",j," size = ",cluster_nssms[j],"\n")
    total_size = total_size + cluster_nssms[j]
}
print(paste0("original n_snv = ",total_SNV,"\t","final n_snv = ",total_size,"\n"))
             
print(paste0("BayClone3 on given sample completed\n"))


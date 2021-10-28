## Consensus Phylogenetic tree 
## 
## Camila Pereira Perico (camilapperico) - 2021-sep-17
## UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/
   
## import libraries
    library(ape)
    library(tictoc)

## input files
    sweepfilename = 'SARSCoV2_sweep_matrix.rds'   # vectorized sequences

## define principal parameters
    nBP = 100        # number of replica for bootstrap
    ncores = 2      # number of CPU cores available

    # choose a sample to root the phylogenetic tree (generally the Wuhan Reference Sequence)
    # just for exemplifying the use of this script, the sample below was choosen:
    outgroupname = "EPI_ISL_2274979" 

## provide the parameters above, than copy and paste the commands below .....


# ----------------------------------------------------------------------------- #
## import the sweep matrix
    sweep_matrix = readRDS(sweepfilename)
    dim(sweep_matrix) # the output must be  Nsequences  600


## create distance matrix and the neighbor-joining tree
    estimate_tr <- function(m) root(nj(dist(m, method="euclidean")), outgroup = outgroupname, resolve.root = TRUE)
    mdist <- estimate_tr(sweep_matrix)
    saveRDS(mdist,'DistanceMatrix.rds')

    print('Distance Matrix saved')

    print('starting bootstrap')

    tic()
    bs <- boot.phylo(mdist, sweep_matrix, estimate_tr, trees=TRUE,mc.cores=ncores ,B = nBP)  # B = number of replica
    toc()

    consensustree <- consensus(bs$trees, p=0.5)
    consensustree2 <- ladderize(consensustree) #  reorganizes the internal structure of the tree

    saveRDS(bs,'NJ_bootphylo.rds')
    saveRDS(consensustree2,'NJ_consensus_tree.rds')
    write.tree(consensustree2, file = "NJ_consensus_tree.tree")
    saveRDS(bs$BP/nBP,'NJ_consensus_BPvalues.rds')

## simplified visualization of the tree with corresponding bootstraps
# plot(consensustree2)
# nodelabels(text = bs$BP/nBP)

    print('results saved')


    print('output file names:')
    print(' distance matrix        : DistanceMatrix.rds')
    print(' all generated trees    : NJ_bootphylo.rds')
    print(' consensus tree (rds)   : NJ_consensus_tree.rds')
    print(' consensus tree (newick): NJ_consensus_tree.tree')
    print(' bootstrap values (%)   : NJ_consensus_BPvalues.rds')
    print('process finished')

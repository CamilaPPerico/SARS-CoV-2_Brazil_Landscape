## Consensus Clustering step
## 
## Camila Pereira Perico (camilapperico) - 2021-sep-17
## UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/
   
## import libraries
    library(ConsensusClusterPlus)  # a bioconductor package
    library(tictoc)

## inputfile - vectorized sequences
    sweepfilename = 'SARSCoV2_sweep_matrix.rds'

## define principal parameters
## (other parameters can be set in the ConsensusClusterPlus function)
    maxK = 15                           # maximum number of clusters
    reps = 1000                         # resampling number
    outputfoldername = 'consensusCLU'   # create a folder for the output files

## provide the parameters above, than copy and paste the commands below .....

# ----------------------------------------------------------------------------- #

## import the sweep matrix
    sweep_matrix = readRDS(sweepfilename)
    dim(sweep_matrix) # the output must be  Nsequences  600
    print('sweep_matrix loaded')

      data = t(sweep_matrix)

      print('Starting clustering')

      tic()
      results = ConsensusClusterPlus( d=data, 
                                      maxK = maxK,              # maximum number of clusters
                                      reps=reps,                # resampling number
                                      distance="euclidean",     # distance metric 
                                      clusterAlg="pam",         # algorithm
                                      title=outputfoldername,   # create a folder
                                      plot="pdf"                # export a pdf report (recommended)
                                    ) 
      toc()

      # save the results in previously created folder
    saveRDS(results, file = paste(outputfoldername,"/ConsensusResult.rds",sep=''))

## Extract the consensus class of each sample. 
 # Each column correspond to the number of clusters of consensus.
    ccp_class = matrix(0,dim(sweep_matrix)[1],maxK)
    for (k in 2:maxK){
            ccp_class[,k] = results[[k]][["consensusClass"]]
    }
    # save in previously created folder
    saveRDS(ccp_class, file= paste(outputfoldername,"/ConsensusClass.rds",sep=''))

    print('results and class saved')

    print('outputfoldername:')
    print(paste('    ',outputfoldername,'/',sep=''))
    print('process finished')

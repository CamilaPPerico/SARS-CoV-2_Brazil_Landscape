## Visualization of clusters
## 
## Camila Pereira Perico (camilapperico) - 2021-sep-17
## UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/
   
## import libraries
    library(Rtsne)

## input files 
    outputfoldername = 'consensusCLU'             # folder of consensus clustering
    sweepfilename = 'SARSCoV2_sweep_matrix.rds'   # vectorized sequences

## define principal parameters
# After the visual analysis of the CDF curve, suppose that the optimal number of clusters is 4
    nclu = 4
    
## provide the parameters above, than copy and paste the commands below .....


# ----------------------------------------------------------------------------- #
## import the sweep matrix
    sweep_matrix = readRDS(sweepfilename)
    dim(sweep_matrix) # the output must be  Nsequences  600
    
## import clusters 
    ccp_class = readRDS(paste(outputfoldername,"/ConsensusClass.rds",sep=''))
    clusters = ccp_class[,nclu]

## PCA - Principal Component Analysis
    # do the PCA - Principal Component Analysis
    pca_output <- prcomp(sweep_matrix, scale = FALSE)
    
    saveRDS(pca_output,'PCA_output.rds')

    
## Visualization of PCA componentes 1 and 2 colored by cluster and save
    par(oma=c(0,0,2,0))
      plot( pca_output$x[,1],pca_output$x[,2],
            xlab='PC-1',ylab='PC-2',col=clusters,pch=20)
      legend( min(pca_output$x[,1]),max(pca_output$x[,2]),
              c(1:nclu),col=as.character(c(1:nclu)),pch=20)
    dev.print(pdf, 'PCA.pdf',width=10, height=6)
    dev.off()


## t-SNE t-Distributed Stochastic Neighbor Embedding
    
    # the t-SNE require removing duplicates, 
    # then change the row names to indices
    therownames = rownames(sweep_matrix) # save if necessary
    rownames(sweep_matrix) <- 1:dim(sweep_matrix)[1]
    sweep_matrix_uniques = unique(sweep_matrix)
   
    # extract the index of unique sequences (vectors)
    idxunique = as.numeric(rownames(sweep_matrix_uniques))
    
    # do the t-SNE and save
    # the Barnes-Hut algorithm is the default 
    tsne_output <- Rtsne(sweep_matrix_uniques, 
                          dims=2, 
                          pca=FALSE, 
                          perplexity=20, 
                          dist=euclidean) 

    saveRDS(tsne_output,'tSNE_output.rds')

## Visualization of t-SNE colored by cluster and save
    cl_unique = clusters[idxunique] # use only the information of the unique vectors

    par(oma=c(0,0,2,0))
      plot(tsne_output$Y,col=cl_unique, asp=1,xlab='t-SNE 1',ylab='t-SNE 2',pch=20)
      legend( min(tsne_output$Y[,1]),max(tsne_output$Y[,2]),
              legend=as.character(c(1:nclu)),pch=20,col=as.character(c(1:nclu)))
    dev.print(pdf, 't-SNE.pdf',width=10, height=6)
    dev.off()
    
    print('resultssaved')


    print('output file names:')
    print('data matrices:')
    print('  - PCA_output.rds')
    print('  - tSNE_output.rds')
    print('pdf graphs:')
    print('  - PCA.pdf')
    print('  - t-SNE.pdf')
    print('process finished')

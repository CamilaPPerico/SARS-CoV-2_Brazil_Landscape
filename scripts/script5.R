## SWeeP projection of sequences in vectors
## 
## Camila Pereira Perico (camilapperico) - 2021-sep-17
## UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/

## import libraries
  library(rSWeeP)  # a bioconductor package
  library(pracma)
  library(stats)
  library(Biostrings)  # a bioconductor package
  library(tictoc)

## filenames
  fastafile = 'MultifastaSARSCoV2_Filtered.fasta'     # input filename - fasta sequences filtered
  epicodefile = 'MultifastaSARSCoV2_epi_Filtered.txt' # input filename - file with epi code filtered
  sweepfilename = 'SARSCoV2_sweep_matrix.rds'         # output filename - vectorized sequences
  orthbasefilename = "orthBase.rds"                   # output filename - orthonormal matrix (if not preexisting)

# If desirable, to perform multiple projections with different data, 
# it is necessary to always use the same orthbase to ensure comparability between projections.  
# Different orthonormal matrices generate different results that are not comparable with each other. 

## provide the parameters above, than copy and paste the commands below .....

# ----------------------------------------------------------------------------- #
    

## generate and save the orthonormal matrix for first usage
 baseMatrix <- orthBase(160000,600) 
 saveRDS(baseMatrix, file = orthbasefilename) 
 print('orthbase created')

## uncomment bellow for reuse the same baseMatrix for more data
   # baseMatrix = readRDS(orthbasefilename)
   # print('orthbase loaded')

## do projection using rSWeeP
  tic()
  X <- sWeeP(fastafile,baseMatrix)
  toc()

## rename the row of data with EPI code
  epi = unlist(read.table(epicodefile))
  rownames(X) <- epi

## save the results
  # save the results in RDS format
  saveRDS(X, file=sweepfilename)

  # uncomment bellow for save in plain text, if you desire to use in other platforms
  # write.table(X, file='SARSCoV2_sweep_matrix.txt')

  print('sweep saved')
  
  print('outputfilenames:')
  print(paste('vectorized seqs: ',sweepfilename,sep=""))
  print(paste('orthbase file  : ',orthbasefilename,sep=""))


  print('process finished')

## Selection of quality sequences
## 
## Camila Pereira Perico (camilapperico) - 2021-sep-17
## UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/

 
## import libraries
  library(tidyr)
  library(stringr)

  inputfilepath = 'MultifastaSARSCoV2.fasta'  # provide the multifasta filename 
  listqualityname = 'indexoffiltered.txt'     # index list of filtered sequences
## provide the parameters above, than copy and paste the commands below .....

# ----------------------------------------------------------------------------- #

## reading the multifasta
  file_multi = file(inputfilepath, "r")
  k=1
  
  good_list = list()

  while ( TRUE ) {
    line1 = readLines(file_multi, n = 1) # header
    if ( length(line1) == 0 ) { # end of file
      break
    }

    line2 = readLines(file_multi, n = 1) # protein sequence

    epiM = data.frame(h=sub(">", "", line1)) # remove '>' symbol
    # count of X
    ctX = str_count(line2, "X") # requires stringr
    # count of X proteins - how many '*****'
    ctpr = str_count(line2, "\\*\\*\\*\\*\\*") # requires stringr

    # A good string is one that has no X and and has 26 or more proteins
    if (ctpr>=25 & ctX==0){
      good_list[k] = epiM
      print(k)
      k=k+1     
    }

  }

  close(file_multi) # close the opened multifasta

  good_list = unlist(good_list)
  write.table(good_list,listqualityname,row.names=FALSE,col.names=FALSE)
  


  print('outputfilename:')
  print('index list of filtered sequences:')
  print('    indexoffiltered.txt')

  print('process finished')

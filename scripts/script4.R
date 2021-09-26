## Filtering quality sequences
## 
## Camila Pereira Perico (camilapperico) - 2021-sep-17
## UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/

## import libraries
    library(stringi)

## filenames
    listqualityname        = 'indexoffiltered.txt'                  # index list of filtered sequences
    input_csvdataname      = 'Brazildataframe.csv'                  # csv with basics information
    output_csvdataname     = 'Brazildataframe_Filtered.csv'         # csv with basics information (filtered)
    input_multifastaname   = 'MultifastaSARSCoV2.fasta'             # multifasta
    output_multifastaname  = 'MultifastaSARSCoV2_Filtered.fasta'    # multifasta (filtered)
    input_headerfilename   = 'MultifastaSARSCoV2_h.txt'             # headers file
    output_headerfilename  = 'MultifastaSARSCoV2_h_Filtered.txt'    # headers file (file)
    input_epicodefilename  = 'MultifastaSARSCoV2_epi.txt'           # epi code file 
    output_epicodefilename = 'MultifastaSARSCoV2_epi_Filtered.txt'  # epi code file (filtered)

## provide the parameters above, than copy and paste the commands below .....

# ----------------------------------------------------------------------------- #


## import the list of filtered
  epiok = unlist(read.table(listqualityname))
  
## import the csv data and filter
  dt = read.csv(csvdataname)
  write.table('epi,date,code',file=output_csvdataname,
                quote = FALSE, row.names=FALSE,col.names=FALSE)
  for (k in 1:length(epiok)){
    write.table(dt[which(dt$epi==epiok[k]),],file=output_csvdataname,
                  sep=',',append=TRUE , quote = FALSE, 
                  row.names=FALSE, col.names=FALSE)
  }
  

## import the headers and filter
  h = unlist(read.csv(input_headerfilename))
  for (k in 1:length(epiok)){
    idx = which(stri_detect_fixed(h,epiok[k]))
    write.table(h[idx],file=output_headerfilename,
                    sep=',',append=TRUE , quote = FALSE, 
                    row.names=FALSE, col.names=FALSE)
  }

## import the EPI list and filter (maintaining original order)
  h = unlist(read.table(input_epicodefilename))
  for (k in 1:length(h)){
    if (sum(stri_detect_fixed(h[k],epiok))){  # is there a valid epi in the header?
      write.table(h[k],file=output_epicodefilename,
                    sep=',',append=TRUE , quote = FALSE, 
                    row.names=FALSE, col.names=FALSE)  
    }
    
  }


## Define the print function
  write2file <- function(line1, line2,outputname){
    write.table(line1,file=outputname, 
                sep="", append=TRUE , quote = FALSE, 
                row.names=FALSE, col.names=FALSE)
    write.table(line2,file=outputname, 
                sep="", append=TRUE , quote = FALSE, 
                row.names=FALSE, col.names=FALSE)
  }


## import the Multifasta and filter
  file_multi = file(input_multifastaname, "r")
  k=1  
  while ( TRUE ) {
    line1 = readLines(file_multi, n = 1) # header
    if ( length(line1) == 0 ) { # end of file
      break
    }
    line2 = readLines(file_multi, n = 1) # sequence

    epiM = sub(">", "", line1) # remove '>' symbol

    if ( sum(epiM == epiok)!=0 ){
      # create another file with filtered sequences
      write2file(line1, line2,output_multifastaname)
    }
    print(k)
    k=k+1
  }
  close(file_multi) # close the opened data


    print('outputfilenames:')
    print(paste('basics info: ',output_csvdataname,sep=""))
    print(paste('multifasta : ',output_multifastaname,sep=""))
    print(paste('headers    : ',output_headerfilename,sep=""))
    print(paste('EPI code   : ',output_epicodefilename,sep=""))

    print('process finished')

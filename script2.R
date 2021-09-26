## Concatenation of sequences in multifasta file
## listing of headers and EPI codes in own files 
## 
## Camila Pereira Perico (camilapperico) - 2021-sep-17
## UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/


## import libraries
    library(Biostrings) # a bioconductor package
    library(tidyr)
    library(stringi)

## provide the path to the folder containing the FASTA files. 
#  Each FASTA corresponds to an organism
    path = 'seqs/'
    outputfilename = 'MultifastaSARSCoV2' # provide the prefix filename of outputs

## provide the extension of the files in path, 
#  or leave it blank between the quotes ("") if there is no pattern
    extension = '.faa' 

## provide the parameters above, than copy and paste the commands below .....

# ----------------------------------------------------------------------------- #


## The function - copy and paste the function below  before use
    concatenafasta <- function(path,outputfilename,extension) {
        
    ## Create the lists of files and names
        fastalist = list.files(path,pattern=extension,
                               full.names=TRUE,recursive=TRUE)
        namesfasta = list.files(path,pattern=extension,
                                full.names=FALSE,recursive=TRUE)
        namesfasta = sub(".faa", "", namesfasta)
        namesorg = list()
        EPIs = list()
        allfastas = list()
        toprint=list()

        for (k in 1:length(fastalist)){
            fastaFile <- readBStringSet(fastalist[k]) # require Biostrings
            seq_name = names(fastaFile)
            sequence = paste(fastaFile)
            df <- data.frame(seq_name, sequence)

            allfastas[k] = stri_join_list(list(sequence), 
                                          sep = "*****", 
                                          collapse = NULL
                                         ) # require stringi
                        aux = data.frame(h=unlist(seq_name[1]))
            aux2 = separate( data=aux,
                      col=h,
                      into = c("pr","code","date","EPI","other"),
                      sep='\\|'
                    )

            ## The header of each concatenated sequence in multifasta is 
            #  the corresponding filename
            toprint[1] = paste('>',namesfasta[k],sep="")
            toprint[2] = allfastas[k]

            ## Save the final multifasta and the first header of each FASTA 
            #  file (in the text file)

            write.table(unlist(toprint),
                    paste(outputfilename,'.fasta',sep=""), append=TRUE , 
                    col.names=FALSE,row.names=FALSE,quote = FALSE)
            write.table(unlist(seq_name[1]),
                        paste(outputfilename,'_h.txt',sep=""), append=TRUE ,
                        col.names=FALSE,row.names=FALSE,quote = FALSE)
            write.table(unlist(aux2$EPI),
                        paste(outputfilename,'_epi.txt',sep=""), append=TRUE ,
                        col.names=FALSE,row.names=FALSE,quote = FALSE)

	    print(paste(as.character(k),'of', as.character(length(fastalist))))
        }

    print('outputfilenames:')
    print(paste('multifasta: ',outputfilename,'.fasta',sep=""))
    print(paste('headers   : ',outputfilename,'_h.txt',sep=""))
    print(paste('EPI code  : ',outputfilename,'_epi.txt',sep=""))

    print('process finished')
    
}   
## ==============================

# execute
concatenafasta(path,outputfilename,extension)

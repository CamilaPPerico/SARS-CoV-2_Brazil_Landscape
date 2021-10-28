## Extraction of proteomes of interest from GISAID's allprotRELEASE.fasta file
## 
## Camila Pereira Perico (camilapperico) - 2021-sep-17
## UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/


## import libraries
  library(tidyr)
  library(stringi)

## define the parameters here
  inputfilepath = 'Sample_GISAID_Release0609.fasta' # the GISAID file 
  outputfoldername = 'seqs'                         # folder name for outputs fasta
  onlythiscountry = TRUE                            # if you want all the data, change to FALSE
  country = 'Brazil'                                # specify the country that you desire to study
  dataframename = 'Brazildataframe'                 # name of the csv with basics information


## provide the parameters above, than copy and paste the commands below .....

# ----------------------------------------------------------------------------- #

## writing function
  write2file <- function(line1, line2,h, outputfoldername, oldEPI){
    write.table(line1,file=paste(outputfoldername,'/',h$EPI,'.faa',sep=""), 
                sep="", append=TRUE , quote = FALSE, 
                row.names=FALSE, col.names=FALSE)
    write.table(line2,file=paste(outputfoldername,'/',h$EPI,'.faa',sep=""), 
                sep="", append=TRUE , quote = FALSE, 
                row.names=FALSE, col.names=FALSE)
  
    if (oldEPI != h$EPI){
      dt = data.frame(epi = h$EPI, date=h$date, code=h$code)
      write.table(dt,file=paste(dataframename,'.csv',sep=""),
                  sep=',',append=TRUE , quote = FALSE, 
                  row.names=FALSE, col.names=FALSE)
          }
  }
  
  print('starting reading...')

  #create the folder of the outputs
  dir.create(outputfoldername)

  con = file(inputfilepath, "r")
  k=1
  oldEPI="" # start empty variable
  write.table('epi,date,code',file=paste(dataframename,'.csv',sep=""),
                quote = FALSE, row.names=FALSE,col.names=FALSE)

  while ( TRUE ) {
    line1 = readLines(con, n = 1)     # header
    if ( length(line1) == 0 ) {       # end of file
      break
    }

    line2 = readLines(con, n = 1) # protein sequence

    line1_ = data.frame(h=sub(">", "", line1)) # remove '>' symbol
    h = separate( data=line1_,
                  col=h,
                  into = c("pr","code","date","EPI","other"),
                  sep='\\|'
                )

    if (onlythiscountry){
           if( stri_detect_fixed(h$code,country)){ 
           write2file(line1, line2,h, outputfoldername, oldEPI)
           print(k)
           k=k+1
      	}
    }
    else {
        write2file(line1, line2,h, outputfoldername, oldEPI)
        print(k)
        k=k+1
    }
      
      oldEPI = h$EPI

  }

  close(con) # close the opened GISAID file

  print('outputfilenames:')
  print(paste('basics information: ',dataframename,'.csv',sep=""))
  print(paste('output folder     : ',outputfoldername,'/',sep=""))

  print('process finished')

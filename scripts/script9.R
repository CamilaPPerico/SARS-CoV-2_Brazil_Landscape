## Consensus Phylogenetic tree visualization
## 
## Camila Pereira Perico (camilapperico) - 2021-sep-17
## UFPR Bioinformatics team - http://www.bioinfo.ufpr.br/
   
## import libraries
  library(wrapr)    
  library(phytools)
  library(ggtree)
  library(ggplot2)
  library(ggnewscale)

## input files
    consensusTreefilename = 'NJ_consensus_tree.tree'     # consensus tree (newick format)
    BPfilename            = 'NJ_consensus_BPvalues.rds'  # bootstrap values (%)
    clustersfilename      = 'clusters.rds'               # consensus clustering of each sample
    
## provide the parameters above, than copy and paste the commands below .....


# ----------------------------------------------------------------------------- #



## import neighbor-joining tree
  mytree   = read.tree(consensusTreefilename)
  BP       = readRDS(BPfilename)
  clusters = readRDS('clusters.rds')

  # sort the tree data the order of the branches in the 
  # tree using the EPI as a reference
  p <- match_order(clusters$sample, mytree$tip.label) 
  clusters = clusters[p,]

  groupInfo <- split(mytree$tip.label, clusters$cl)    
  mytree <- groupOTU(mytree, groupInfo)

  my_colors <- c( "#ff0000ff", "#0000ffff", "#28ff00ff", "#8a2be2ff", "#0fc0fcff", "#fdee00ff", 
                  "#eb5100ff", "#008000ff", "#ff007fff", "#999999ff", "#fc89acff", "#704214ff", 
                  "#caab00ff", "#750025ff", "#66fd79ff", "#ff6e4aff", "#000000ff")

  mytree$node.label<-round(rev(BP),digits=2)

  ## to get a tree with the branches length in the image, 
  #  remove the  option: branch.length="none"
  #  here the tree core is created and colored by cluster
  plotoftree<-ggtree(mytree, aes(color=group), layout='circular',branch.length="none") +
          scale_colour_manual(values=c("black", my_colors)) + 
          geom_tiplab(size=1, aes(angle=angle))+ 
          guides(shape = guide_legend(override.aes = list(size = 10)))+ 
          guides(color = guide_legend(override.aes = list(size = 10)))+
          theme(legend.title = element_text(size = 10), 
                  legend.text = element_text(size = 10))+
          geom_nodelab(aes(label=label),hjust=0.2, size=1)
  
  # plot and save the tree in pdf format
  par(oma=c(0,0,2,0))
  plotoftree
  dev.print(pdf,'NJrooted_tree.pdf',width=10, height=10)
  dev.off()

  print('results saved')


  print('output file names:')
  print(' plot of Phylogenetic tree : NJrooted_tree.pdf')
  print('process finished')

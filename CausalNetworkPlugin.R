#rm(list = ls())

library(bnlearn)
#setwd("D:\\FIU\\Causal Network Inference in Cancer") # This is input folder, input and output file will be in this directory

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
    pfix = prefix()
  if (length(pfix) != 0) {
     pfix <- paste(pfix, "/", sep="")
  }

#bootfile <- "boot_strength.csv"
#data_file_name <- "asv.transpose.norm.csv" #input file
#graph_id <-  1001
#filename_dir <- "directed.csv"
#filename_undir <- "undirected.csv"

bootfile <<- paste(pfix, parameters["bootfile", 2], sep="/")
data_file_name <<- paste(pfix, parameters["datafile", 2], sep="/")
graph_id <<- as.integer(parameters["graphid", 2])
filename_dir <<- paste(pfix, parameters["directed", 2], sep="/")
filename_undir <<- paste(pfix, parameters["undirected", 2], sep="/")
filename_group1 <<- paste(pfix, parameters["group1", 2], sep="/")
filename_group2 <<- paste(pfix, parameters["group2", 2], sep="/")


#group2 <- c("Selen1","Selen2","Porph1","Porph2","Trepo2")
#group1 <- c("Veill1","Actin2","Strep1","Strep2","Actin1","Strep3","Actin3")
group1 <<- c()
cont1 <- readLines(filename_group1)
group1 <<- c(group1, cont1)
group2 <<- c()
cont2 <- readLines(filename_group2)
group2 <<- c(group2, cont2)

}



run <- function() {
##########################################################################################
#### Drop out of redundant edges
print("REMOVING REDUNDANT UNDIR")
data <- read.csv(filename_undir,header = TRUE, colClasses=c("from"="character","to"="character"))
bn_df <- data.frame(data)


redundant_undi_edge <- c()

for (i in 1:((nrow(bn_df))))
{
  for (j in 1:nrow(bn_df))
  {
    if(i!=j & !(i %in% redundant_undi_edge))
    {
      if(bn_df[i,1]==bn_df[j,2] & bn_df[i,2]==bn_df[j,1])
      {
        redundant_undi_edge <- append(redundant_undi_edge,j)
      }
      
    }
  }
  
}

#### Boot strength

print("ASSEMBLING UNDIR WEIGHTS")
boot_data <- read.csv(bootfile,header = TRUE, colClasses=c("from"="character","to"="character","strength"="double","direction"="double"))


#### Undirected Edge
undirected_edge <<- bn_df [c(redundant_undi_edge),]
undirected_edge$directed <<- rep(FALSE,nrow(undirected_edge))

undirected_weights <- c()

print(nrow(undirected_edge))
print(nrow(boot_data))
for (i in 1:nrow(undirected_edge))
{
	print(i)
  for(j in 1:nrow(boot_data))
  {
    if (undirected_edge[i,1]==boot_data[j,1] & undirected_edge [i,2] == boot_data [j,2])
    {
      undirected_weights <- append(undirected_weights,boot_data[j,3])
    }
    
  }
  
}

undirected_edge$weight <<- undirected_weights

#### Directed Edge

print("ADDING DIRECTED")
directed_edge <<- read.csv(filename_dir,header = TRUE, colClasses=c("from"="character","to"="character"))
directed_edge <<- data.frame(directed_edge)

directed_edge_weight <- c()

print(nrow(directed_edge))

for (i in 1:nrow(directed_edge))
  
{
	print(i)
  for( j in 1:nrow(boot_data))
  {
    if(directed_edge[i,1]== boot_data[j,1] & directed_edge[i,2]== boot_data[j,2])
    {
      directed_edge_weight <- append(directed_edge_weight,boot_data[j,4])
    }
  }
  
}

directed_edge$directed <<- rep(TRUE,nrow(directed_edge))
directed_edge$weight <<- directed_edge_weight
}

output <- function(outputprefix) {
con_file_name <- paste(outputprefix, "xgmml", sep=".") #output file
network_csv_file <- paste(outputprefix, "network", "csv", sep=".")

#### Writing Network file
kera_gingiva_net_file <- rbind(directed_edge,undirected_edge)
print("PEARSON...")
#### Writing xgmml file for cytosacpe-backend

data_file <- read.csv(data_file_name,row.names=1,header=TRUE)
#data_file <- as.numeric(data.frame(data_file))

#data_file <- apply(data_file, 1, as.numeric)
#data_file <- t(data_file)

#### Correlation Data
cor_value <- cor(data_file, method = "pearson")
cv_col <- colnames(cor_value)
cv_row <- rownames(cor_value)

correlation <- c()

print(nrow(kera_gingiva_net_file))
print(length(cv_col))
print(length(cv_row))

for(i in 1:nrow(kera_gingiva_net_file))
{
	print(i)
  for (j in 1: length(cv_col))
    
  {
    print(kera_gingiva_net_file[i,1])
    print(cv_col[j])
    if (kera_gingiva_net_file[i,1]==cv_col[j])
    {
      x = j
    }
    
  }  
  
  for (k in 1: length(cv_row))
    
  {
    if (kera_gingiva_net_file[i,2]==cv_row[k])
    {
      y = k
    }
    
  } 
  
  correlation <- append(correlation,cor_value[x,y])
  
}

comp_bn_cor <- kera_gingiva_net_file[,1:4]
comp_bn_cor$pearson <- correlation

###
print("SPEARMAN...")
spearman_cor <- cor(data_file, method = "spearman")
spearman <- c()

for(i in 1:nrow(kera_gingiva_net_file))
{
	print(i)
  for (j in 1: length(cv_col))
    
  {
    if (kera_gingiva_net_file[i,1]==cv_col[j])
    {
      x = j
    }
    
  }  
  
  for (k in 1: length(cv_row))
    
  {
    if (kera_gingiva_net_file[i,2]==cv_row[k])
    {
      y = k
    }
    
  } 
  
  spearman <- append(spearman,spearman_cor[x,y])
  
}
comp_bn_cor$spearman <- spearman

###
print("KENDALL...")
kendall_cor <- cor(data_file, method = "kendall")
kendall <- c()

for(i in 1:nrow(kera_gingiva_net_file))
{
	print(i)
  for (j in 1: length(cv_col))
    
  {
    if (kera_gingiva_net_file[i,1]==cv_col[j])
    {
      x = j
    }
    
  }  
  
  for (k in 1: length(cv_row))
    
  {
    if (kera_gingiva_net_file[i,2]==cv_row[k])
    {
      y = k
    }
    
  } 
  
  kendall <- append(kendall,kendall_cor[x,y])
  
}
comp_bn_cor$kendall <- kendall
write.csv(comp_bn_cor,paste(outputprefix,"correlation","csv",sep="."),row.names = F)

#library(bnlearn)
#data_file_name <- "asv.transpose.norm.csv" #input file

#data_file <- read.csv(data_file_name,row.names=1,header=TRUE)

data_file_colname <- colnames(data_file)
data_colname <- data.frame(data_file_colname)

node_id <- c()

for (i in 1:nrow(data_colname))
{
  node_id <- append(node_id,i+1000)
}

node_data <- data.frame(data_colname,node_id)

edge_id <- c()
source_id <- c()
target_id <- c()

for (i in 1:nrow(kera_gingiva_net_file))
{
  edge_id <- append(edge_id,i+100)
  for(j in 1:nrow(node_data))
  {
    if(kera_gingiva_net_file[i,1]==node_data[j,1])
    {
      source_id <- append(source_id,node_data[j,2])
      
    }
  }
}

for (i in 1:nrow(kera_gingiva_net_file))
{
  for(j in 1:nrow(node_data))
  {
    if(kera_gingiva_net_file[i,2]==node_data[j,1])
    {
      target_id <- append(target_id,node_data[j,2])
      
    }
  }
}

kera_gingiva_net_file$edgeid <- edge_id
kera_gingiva_net_file$sourceid<- source_id
kera_gingiva_net_file$targetid<- target_id

kera_gingiva_net_file <- na.omit(kera_gingiva_net_file)

##
## xml file writing

name <- network_csv_file

#### -------------------CoN--------------------------------------------
bn_df <- read.csv(data_file_name,row.names=1,header=TRUE)


sum_of_col <- data.frame(colnames(bn_df),colSums(bn_df))
col_name <- colnames(bn_df)
sink(con_file_name)
cat("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n")
cat(sprintf("<graph id=\"%d\" label=\"%s\" directed=\"1\" cy:documentVersion=\"3.0\" xmlns:dc=\"http://purl.org/dc/elements/1.1/\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns:cy=\"http://www.cytoscape.org\" xmlns=\"http://www.cs.rpi.edu/XGMML\">",graph_id,name))
cat("\n<att name=\"networkMetadata\">")
cat("\n<rdf:RDF>")
cat("\n<rdf:Description rdf:about=\"http://www.cytoscape.org/\">")
cat("\n<dc:type>Protein-Protein Interaction</dc:type>")
cat("\n<dc:description>N/A</dc:description>")
cat("\n<dc:identifier>N/A</dc:identifier>")
cat("\n<dc:date>2018-02-28 19:57:24</dc:date>\n")
cat(sprintf("<dc:title>%s</dc:title>",name))
cat("\n<dc:source>http://www.cytoscape.org/</dc:source>")
cat("\n<dc:format>Cytoscape-XGMML</dc:format>")
cat("\n</rdf:Description>")
cat("\n</rdf:RDF>")
cat("\n </att>")
cat(sprintf("\n<att name=\"shared name\" value=\"%s\" type=\"string\"/>",name))
cat(sprintf("\n<att name=\"name\" value=\"%s\" type=\"string\"/>",name))
cat("\n<att name=\"selected\" value=\"1\" type=\"boolean\"/>")
cat("\n<att name=\"__Annotations\" type=\"list\">")
cat("\n</att>")
cat("\n<att name=\"layoutAlgorithm\" value=\"Hierarchical Layout\" type=\"string\" cy:hidden=\"1\"/>")
cat("\n
    <graphics>
    \n\t<att name=\"NETWORK_CENTER_Y_LOCATION\" value=\"200.0\" type=\"string\"/>
    \n\t<att name=\"NETWORK_SCALE_FACTOR\" value=\"0.3942098543284541\" type=\"string\"/>
    \n\t<att name=\"NETWORK_DEPTH\" value=\"0.0\" type=\"string\"/>
    \n\t<att name=\"NETWORK_WIDTH\" value=\"833.0\" type=\"string\"/>
    \n\t<att name=\"NETWORK_EDGE_SELECTION\" value=\"true\" type=\"string\"/>
    \n\t<att name=\"NETWORK_NODE_SELECTION\" value=\"true\" type=\"string\"/>
    \n\t<att name=\"NETWORK_HEIGHT\" value=\"600.0\" type=\"string\"/>
    \n\t<att name=\"NETWORK_CENTER_X_LOCATION\" value=\"400\" type=\"string\"/>
    \n\t<att name=\"NETWORK_CENTER_Z_LOCATION\" value=\"0.0\" type=\"string\"/>
    \n\t<att name=\"NETWORK_BACKGROUND_PAINT\" value=\"#FFFFFF\" type=\"string\"/>
    \n\t<att name=\"NETWORK_TITLE\" value=\"\" type=\"string\"/>
    \n</graphics>
    ")

#### Writing node information
for(i in 1:nrow(node_data))
{
  cat(sprintf("\n<node id=\"%d\" label=\"%s\">",node_data[i,2],node_data[i,1]))
  cat(sprintf("\n<att name=\"shared name\" value=\"%s\" type=\"string\"/>",node_data[i,2]))
  cat(sprintf("\n<att name=\"name\" value=\"%s\" type=\"string\"/>",node_data[i,1]))
  cat("\n<att name=\"selected\" value=\"0\" type=\"boolean\"/>")
  
  ## group1
  for (j in 1:length(group1))
  {
    if (node_data[i,1] == group1[j])
    {
      for(k in 1:nrow(sum_of_col))
      {
        if (sum_of_col[k:k,1:1]==node_data[i,1])
        {
          abundance <- sum_of_col[k:k,2:2] 
          
        }
      }
      
      if (abundance >= 0)
        abun <- log(abundance,2)*12
      if(abun<=10)
        abun <- 10
      
      cat(sprintf("\n<graphics  outline=\"#9400D3\"  h=\"%f\" w=\"%f\" fill=\"#FFFFFF\" type=\"ELLIPSE\" width=\"5.0\">",abun,abun))
      
      cat(sprintf("
                  
                  \n\t<att name=\"NODE_NESTED_NETWORK_IMAGE_VISIBLE\" value=\"true\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL_TRANSPARENCY\" value=\"255\" type=\"string\"/>
                  \n\t<att name=\"NODE_VISIBLE\" value=\"true\" type=\"string\"/>
                  \n\t<att name=\"NODE_DEPTH\" value=\"0.0\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL_WIDTH\" value=\"200.0\" type=\"string\"/>
                  \n\t<att name=\"NODE_SELECTED_PAINT\" value=\"#FFFF00\" type=\"string\"/>
                  \n\t<att name=\"NODE_BORDER_TRANSPARENCY\" value=\"255\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL_COLOR\" value=\"#9400D3\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL_FONT_SIZE\" value=\"20\" type=\"string\"/>
                  \n\t<att name=\"NODE_BORDER_STROKE\" value=\"SOLID\" type=\"string\"/>
                  \n\t<att name=\"NODE_TRANSPARENCY\" value=\"0\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL\" value=\"%s\" type=\"string\"/>
                  
                  ", node_data[i,1]))
      
      cat("\n</graphics>")
      
    }
    
  }
  
  ## group2
  for (j in 1:length(group2))
  {
    if (node_data[i,1] == group2[j])
    {
      for(k in 1:nrow(sum_of_col))
      {
        if (sum_of_col[k:k,1:1]==node_data[i,1])
        {
          abundance <- sum_of_col[k:k,2:2] 
          
        }
      }
      
      if (abundance >= 0)
        abun <- log(abundance,2)*12
      if(abun<=10)
        abun <- 10
      
      cat(sprintf("\n<graphics  outline=\"#FF0033\"  h=\"%f\" w=\"%f\" fill=\"#FFFFFF\" type=\"ELLIPSE\" width=\"5.0\">",abun,abun))
      
      cat(sprintf("
                  
                  \n\t<att name=\"NODE_NESTED_NETWORK_IMAGE_VISIBLE\" value=\"true\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL_TRANSPARENCY\" value=\"255\" type=\"string\"/>
                  \n\t<att name=\"NODE_VISIBLE\" value=\"true\" type=\"string\"/>
                  \n\t<att name=\"NODE_DEPTH\" value=\"0.0\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL_WIDTH\" value=\"200.0\" type=\"string\"/>
                  \n\t<att name=\"NODE_SELECTED_PAINT\" value=\"#FFFF00\" type=\"string\"/>
                  \n\t<att name=\"NODE_BORDER_TRANSPARENCY\" value=\"255\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL_COLOR\" value=\"#FF0033\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL_FONT_SIZE\" value=\"20\" type=\"string\"/>
                  \n\t<att name=\"NODE_BORDER_STROKE\" value=\"SOLID\" type=\"string\"/>
                  \n\t<att name=\"NODE_TRANSPARENCY\" value=\"0\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL\" value=\"%s\" type=\"string\"/>
                  
                  ", node_data[i,1]))
      
      cat("\n</graphics>")
      
    }
    
  }
  
  
  # setdiff(x, y)
  
  rest_set <- setdiff(col_name,group1)
  rest_set <- setdiff(rest_set,group2)
  
  
  ## Rest
  for (j in 1:length(rest_set))
  {
    if (node_data[i,1] == rest_set[j])
    {
      for(k in 1:nrow(sum_of_col))
      {
        if (sum_of_col[k:k,1:1]==node_data[i,1])
        {
          abundance <- sum_of_col[k:k,2:2] 
          
        }
      }
      
      # if (abundance >= 0)
      #   abun <- 25
      # if(abun<=10)
      #   abun <- 25
      
      abun <- 25
      
      
      cat(sprintf("\n<graphics  outline=\"#000000\"  h=\"%f\" w=\"%f\" fill=\"#FFFFFF\" type=\"ELLIPSE\" width=\"5.0\">",abun,abun))
      
      cat(sprintf("
                  
                  \n\t<att name=\"NODE_NESTED_NETWORK_IMAGE_VISIBLE\" value=\"true\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL_TRANSPARENCY\" value=\"255\" type=\"string\"/>
                  \n\t<att name=\"NODE_VISIBLE\" value=\"true\" type=\"string\"/>
                  \n\t<att name=\"NODE_DEPTH\" value=\"0.0\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL_WIDTH\" value=\"200.0\" type=\"string\"/>
                  \n\t<att name=\"NODE_SELECTED_PAINT\" value=\"#FFFF00\" type=\"string\"/>
                  \n\t<att name=\"NODE_BORDER_TRANSPARENCY\" value=\"255\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL_COLOR\" value=\"#000000\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL_FONT_SIZE\" value=\"20\" type=\"string\"/>
                  \n\t<att name=\"NODE_BORDER_STROKE\" value=\"SOLID\" type=\"string\"/>
                  \n\t<att name=\"NODE_TRANSPARENCY\" value=\"0\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL\" value=\"%s\" type=\"string\"/>
                  
                  ", node_data[i,1]))
      
      cat("\n</graphics>")
      
    }
    
  }
  
  
  ### End of node graphics
  
  cat("\n</node>")
}

#### Writing edge information

cor_wt <- comp_bn_cor$pearson
kera_gingiva_net_file$pearson <- cor_wt
kera_gingiva_net_file <- na.omit(kera_gingiva_net_file)


for (i in 1:nrow(kera_gingiva_net_file))
{
  cat(sprintf("\n<edge id=\"%d\" label=\"%s (pp) %s\" source=\"%d\" target=\"%d\" cy:directed=\"1\">",kera_gingiva_net_file[i,5],kera_gingiva_net_file[i,1],kera_gingiva_net_file[i,2],kera_gingiva_net_file[i,6],kera_gingiva_net_file[i,7]))
  cat(sprintf("\n\t<att name=\"shared name\" value=\"%s (pp) %s\" type=\"string\"/>",kera_gingiva_net_file[i,1],kera_gingiva_net_file[i,2]))   
  cat("\n\t<att name=\"shared interaction\" value=\"pp\" type=\"string\"/>")
  cat(sprintf("\n\t<att name=\"name\" value=\"%s (pp) %s\" type=\"string\"/>",kera_gingiva_net_file[i,1],kera_gingiva_net_file[i,2]))
  cat("\n\t<att name=\"selected\" value=\"0\" type=\"boolean\"/>")
  cat("\n\t<att name=\"interaction\" value=\"pp\" type=\"string\"/>")
  
  pos_flag <- 0
  neg_flag <- 0
  
  if(kera_gingiva_net_file[i,3]==TRUE)
  {
    cat("\n\t<att name=\"Directed\" value=\"1\" type=\"boolean\"/>")
    
  }
  if(kera_gingiva_net_file[i,3]==FALSE)
  {
    cat("\n\t<att name=\"Directed\" value=\"0\" type=\"boolean\"/>")
  }
  
  
  if (kera_gingiva_net_file[i,8]>0)
  {
    edge_wt <- kera_gingiva_net_file[i,8]
    pos_flag <-1
    
  }
  
  else if (kera_gingiva_net_file[i,8]<0)
  {
    edge_wt <- kera_gingiva_net_file[i,8]*(-1)
    neg_flag <- 1
  }
  
  
  cat(sprintf("\n\t<att name=\"weight\" value=\"%f\" type=\"real\"/>",255))
  
  
  if (pos_flag == 1)
    #cat(sprintf("\n<graphics width=\"%f\" fill=\"#00CC00\">", edge_wt*12))
    cat(sprintf("\n<graphics width=\"%f\" fill=\"#00CC00\">", 2.75))
  
  if (neg_flag==1 )
    #cat(sprintf("\n<graphics width=\"%f\" fill=\"#FF0000\">", edge_wt*12))
    cat(sprintf("\n<graphics width=\"%f\" fill=\"#FF0000\">", 2.75))
  
  if(kera_gingiva_net_file[i,3]==TRUE)
  {
    
    cat("\n\t<att name=\"EDGE_TARGET_ARROW_SHAPE\" value=\"DELTA\" type=\"string\"/>")
    cat("\n\t<att name=\"EDGE_TARGET_ARROW_UNSELECTED_PAINT\" value=\"#000000\" type=\"string\"/>")
    
    
  }
  
  #<att name="EDGE_TRANSPARENCY" value="100" type="string"/>
  #cat(sprintf("\n\t<att name=\"EDGE_TRANSPARENCY\" value=\"%f\" type=\"string\"/>",kera_gingiva_net_file[i,4]*255))
  
  
  # cat(sprintf("\n\t<att Transparency=\"%f\"/>",kera_gingiva_net_file[i,4]*255))
  
  cat("\n</graphics>")
  
  cat("\n</edge>")
  
}

cat("\n</graph>")

sink()

}

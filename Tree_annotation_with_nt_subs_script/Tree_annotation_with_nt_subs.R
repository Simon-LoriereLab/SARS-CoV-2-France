library(tidyverse)
library(jsonlite)

jtree <- fromJSON(txt = "../auspice/ncov_2020-04-22.json", flatten = T)
jree_tree <- jtree$tree

#####  function for extracting nt subs info from nextstrain build in auspice directory
rec_tree <- function(tree, df){
  count_nodes = 1
  node_col <- ""
  nuc_col <- ""
  if (is.null(tree[["name"]])== TRUE){
    node_col <- ""
  }else{
    for (node in list(tree[["name"]])){
      #print(node)
        if (count_nodes  == 1){
          node_col <- node
          count_nodes = count_nodes + 1
        }else{
          node_col <- c(node_col, node)
          count_nodes = count_nodes + 1
        }
    }
  }
  count_nucs = 1
  if (is.null(tree[["branch_attrs.mutations.nuc"]]) == FALSE){
      for (nuc in tree[["branch_attrs.mutations.nuc"]]){
        if (is.null(nuc)== TRUE){
          nuc <- ""
        }
        if (count_nucs == 1){
          nuc_col <- c(paste(nuc, collapse = "_"))
          count_nucs = count_nucs + 1
        }else{
          nuc_col <- c(nuc_col, paste(nuc, collapse = "_"))
          count_nucs = count_nucs + 1
        }
      }
  }else{
    nuc_col <- rep(x = "", length(node_col))
  }
  df_temp <- tibble(node_col, nuc_col)
  df <- rbind(df, df_temp)
  if (class(tree) == "data.frame"){
    child_tree <- tree[["children"]]
    if (is.null(child_tree) == FALSE){
      ##  recursive
      df <- rec_tree(child_tree, df)
    }
  }else if (class(tree) == "list"){
    for (i in 1:length(tree)){
      child_tree <- tree[[i]]
      if (is.null(child_tree) == FALSE){
        ##  recursive
        df <- rec_tree(child_tree, df)
      }
    }
  }
  return(df)
}

# create df with first line
node_col <- "root" 
nuc_col <- ""
df <- tibble(node_col, nuc_col)

###### extract data
new_df <- rec_tree(jree_tree[["children"]], df)


##### remove empty rows and rows without nt subs
new_df1 <- new_df[which(!new_df$node_col == ""),]
new_df1[which(new_df1$nuc_col == ""),]$nuc_col <- "none"

##### keep only nt subs leading to internal nodes and not to tips
new_df2 <- new_df1[which(grepl("NODE_" , new_df1$node_col) == TRUE),] #### do not run this line if wish to leave these nt subs in
new_df3 <- new_df2


##### import the tree that has been already opened and saved with FigTree.
my_tree <- read_file("divergence_tree_saved_after_FigTree_modifications.tre")

##### annotate tree
for (i in 1:length(new_df3$node_col)){
  if (i ==1){
    if (grepl("NODE_" , new_df3[i,]$node_col) == TRUE){
      my_new_tree <- str_replace_all(string = my_tree, 
                                     pattern = paste0('&label=\"',new_df3[i,]$node_col,'\"'), 
                                     replacement = paste0('&label=\"',new_df3[i,]$nuc_col,'\"'))
    }else{
      my_new_tree <- str_replace_all(string = my_tree, 
                                     pattern = paste0("\'",new_df3[i,]$node_col,"\'\\[&!color"), 
                                     replacement = paste0("\'",new_df3[i,]$node_col
                                                          ,"\'\\[&label=",'\"',
                                                          new_df3[i,]$nuc_col,
                                                          '\"',",!color"))
    }
  }else{
    if (grepl("NODE_" , new_df3[i,]$node_col) == TRUE){
      my_new_tree <- str_replace_all(string = my_new_tree, 
                                     pattern = paste0('&label=\"',new_df3[i,]$node_col,'\"'), 
                                     replacement = paste0('&label=\"',new_df3[i,]$nuc_col,'\"'))
    }else{
      my_new_tree <- str_replace_all(string = my_new_tree, 
                                     pattern = paste0("\'",new_df3[i,]$node_col,"\'\\[&!color"), 
                                     replacement = paste0("\'",new_df3[i,]$node_col
                                                          ,"\'\\[&label=",'\"',
                                                          new_df3[i,]$nuc_col,
                                                          '\"',",!color"))
    }
    
  }
  
}

##### save it in the current directory
write_file(x = my_new_tree, path = "tree1.tre", append = FALSE )

########    Clean up from unneccessary annotations
##### import the new tree
my_new_tree1 <- read_file("tree1.tre")

##### remove the travel history annotations
tree_str <- str_split(string = my_new_tree1, pattern = '\\[&label=\"')
count_matches = 1
for (label in tree_str[[1]]){
  match  <- str_match(string = label, pattern = '^[[:print:]]+_travel_history')
  if (is.na(match) == FALSE){
    to_remove <- paste0(match[1,1])
    #print(to_remove)
    if (count_matches == 1){
      my_new_tree2 <- str_replace(my_new_tree1, pattern = to_remove, replacement = "")
      count_matches = count_matches  + 1
    }else{
      my_new_tree2 <- str_replace(my_new_tree2, pattern = to_remove, replacement = "")
      count_matches = count_matches  + 1
    }
    
  }
}
my_new_tree3 <- str_replace_all(my_new_tree2, pattern = "none", replacement = "")


##### save tree
write_file(x = my_new_tree3, path = "tree2.tre", append = FALSE )







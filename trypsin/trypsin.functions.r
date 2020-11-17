read.protein<- function(wild.path, mutated.path){
  
  wild.protein <- read.fasta(file = wild.path, as.string = T,
                             seqtype = "AA",whole.header = T,
                             strip.desc = T)
  
  mutated.protein <- read.fasta(file = mutated.path, as.string = T,
                                seqtype = "AA",whole.header = T)
  
  
  
  return(list("wild" = wild.protein, "mutated" = mutated.protein))
}

do.label <- function(protein.seq, mutation.pos){
  
  substr(protein.seq, mutation.pos, mutation.pos) <- "@"
  
  return(protein.seq)
}

replace.label <- function(protein.seq, AA){
  
  protein.seq <- str_replace_all(protein.seq,"@",AA)
  
  return(protein.seq)
}

do.Cleavages <- function(protein.seq, protein.header , 
                         enzyme="trypsin", missedCleavages = 0,
                         label = "@"){
  
  cleavages <- cleave(protein.seq, enzym=enzyme, 
                      missedCleavages = missedCleavages)
  
  
  cleavages <- data.frame(cleavages)
  
  colnames(cleavages) <- "peptides"
  
  cleavages <- cleavages[grepl(label, cleavages$peptides, fixed = T),]
  
  return(cleavages)
  
}
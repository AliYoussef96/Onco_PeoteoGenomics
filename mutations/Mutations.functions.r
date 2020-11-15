# Reading the mutations table
read.mutat.table <- function(path.mutation.file){
  mutations.table <- read.csv(path.mutation.file, sep = "\t")
  return(mutations.table)
}

# for bad request

stop_for_status.new <- function (x, task = NULL) 
{
  if (status_code(x) < 300) {
    return(TRUE)
  }
  call <- sys.call(-1)
  return(FALSE)
}

# get the protein seq. from ensemble

get.protein <- function(ENST.id){
  ENST.id <- strsplit(ENST.id, ".", fixed = T)[[1]][1]
  
  server <- "https://rest.ensembl.org"
  ext <- paste0("/sequence/id/",ENST.id,"?species=human;type=protein;db=core")
  
  get.fasta <- GET(paste(server, ext, sep = ""), content_type("x-fasta"),timeout(30))
  
  if (stop_for_status.new(get.fasta)){
    
    fasta <- content(get.fasta)
    
    return(fasta$seq)
    
  }else{
    return(ENST.id)
    
  }
  
}

# You got the true protein?

true.protein <- function(protein.seq, mutation.pos){
  mutation.pos <- strsplit(mutation.pos, "p.", fixed = T)[[1]][[2]]
  
  mutation.gusb <- gsub("[0-9]+", "" ,mutation.pos)
  
  normal <- substring(mutation.gusb , 1,1)
  mutated <- substring(mutation.gusb , 2,2)
  
  normal.loc <- str_locate(mutation.pos, normal)[[1]]
  mutated.loc <- str_locate(mutation.pos, mutated)[[1]]
  
  
  
  pos <- substring(mutation.pos , normal.loc+1, mutated.loc-1)
  
  if (substring(protein.seq, pos, pos) == normal){
    return(TRUE)
  }else{return(FALSE)}
}

# Do mutation

do.mutation <- function(protein.seq, mutation.pos){
  mutation.pos <- strsplit(mutation.pos, "p.", fixed = T)[[1]][[2]]
  
  mutation.gusb <- gsub("[0-9]+", "" ,mutation.pos)
  
  normal <- substring(mutation.gusb , 1,1)
  mutated <- substring(mutation.gusb , 2,2)
  
  normal.loc <- str_locate(mutation.pos, normal)[[1]]
  mutated.loc <- str_locate(mutation.pos, mutated)[[1]]
  
  
  
  pos <- substring(mutation.pos , normal.loc+1, mutated.loc-1)
  
  mutated.protein <- protein.seq
  substr(mutated.protein, pos, pos) <- mutated
  
  if (substr(protein.seq,pos,pos) == normal){
    if(substr(mutated.protein,pos,pos) == mutated){
      return(mutated.protein)
    }
  }
}

# Mutation done?

mutation.done <- function(protein.seq, mutated.protein, mutation.pos){
  mutation.pos <- strsplit(mutation.pos, "p.", fixed = T)[[1]][[2]]
  
  mutation.gusb <- gsub("[0-9]+", "" ,mutation.pos)
  
  normal <- substring(mutation.gusb , 1,1)
  mutated <- substring(mutation.gusb , 2,2)
  
  normal.loc <- str_locate(mutation.pos, normal)[[1]]
  mutated.loc <- str_locate(mutation.pos, mutated)[[1]]
  
  
  
  pos <- substring(mutation.pos , normal.loc+1, mutated.loc-1)
  
  if (substr(protein.seq,pos,pos) == normal){
    if(substr(mutated.protein,pos,pos) == mutated){
      return(TRUE)
    }
  }else{return(FALSE)}
}

# write fasta formate

write.protein <- function(protein.seq, gene.name, ensembl_id,
                          mutation.id, gene.cds.length,
                          file.out){
  
  mutation.id <- strsplit(mutation.id, "p.", fixed = T)[[1]][[2]]
  
  protein.length <- nchar(protein.seq)
  
  header.name = paste0("gene_name:",gene.name,"|",
                       "ensembl_id:",ensembl_id,"|",
                       "mutation_id:",mutation.id,"|",
                       "gene_cds_length:",gene.cds.length,"|",
                       "protein_length:",protein.length)
  
  write.fasta(protein.seq, as.string = T, names= header.name,
              file.out = file.out, open = "a")
}

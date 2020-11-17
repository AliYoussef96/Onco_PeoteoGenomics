library("cleaver")
library("seqinr")
library("stringr")
source("trypsin.functions.r")



fasta.files <- read.protein("Mutationsfiles1/normal.proteins1_106164.fasta", 
                            "Mutationsfiles1/mutated.proteins1_106164.fasta")

wild.fasta <- fasta.files$wild
mutated.fasta <- fasta.files$mutated

mCleavages <- 1

loading <- 0
total <- length(names(wild.fasta))
for(i_wild in names(wild.fasta)){
  loading <- loading + 1
  
  print(paste0("wild ", loading, "..." , " From...", total))
  
  protein.wild <- as.character(wild.fasta[[i_wild]])
  
  mutation <- str_split_fixed(i_wild, fixed("|"), 
                                  str_count(i_wild,fixed("|")) +1)[[3]]
  
  mutation <- str_replace_all(mutation, fixed("mutation_id:"), "")
  
  mutation.pos <- substr(mutation, 2 , nchar(mutation)-1)
  
  mutation.wild <- substr(mutation,1,1)
  
  protein.label <- do.label(protein.wild,mutation.pos =  mutation.pos)
  
  peptides <- do.Cleavages(protein.label, i_wild ,missedCleavages = mCleavages, label = "@",
                           enzym = "trypsin")
  
  if(length(peptides) != 0){
  for(i_peptide in peptides){
  
  final.wild <- replace.label(i_peptide, mutation.wild)
    #write fasta with same header
  
  header.name <- paste0(i_wild, "|", "mCleavages = ", mCleavages)
  
  write.fasta(final.wild, as.string = T, names= header.name,
              file.out = "Mutationsfiles1/Wild peptides mCleavages = 1.fasta", open = "a")
  

  }
  }
}

###############################
############
###############################

loading <- 0
total <- length(names(mutated.fasta))
for(i_wild in names(mutated.fasta)){
  loading <- loading + 1  
  print(paste0("mutated ", loading, "..." , " From...", total))
  
  protein.wild <- as.character(mutated.fasta[[i_wild]])
  
  mutation <- str_split_fixed(i_wild, fixed("|"), 
                              str_count(i_wild,fixed("|")) +1)[[3]]
  
  mutation <- str_replace_all(mutation, fixed("mutation_id:"), "")
  
  mutation.pos <- substr(mutation, 2 , nchar(mutation)-1)
  
  mutation.wild <- substr(mutation,1,1)
  
  mutation.mutation <- substr(mutation,nchar(mutation),nchar(mutation))
  
  protein.label <- do.label(protein.wild,mutation.pos =  mutation.pos)
  
  peptides <- do.Cleavages(protein.label, i_wild ,missedCleavages = mCleavages, label = "@",
                           enzym = "trypsin")
  
  if (length(peptides) != 0){
  for(i_peptide in peptides){
    

    fina.mutated <- replace.label(i_peptide, mutation.mutation)
    #write fasta with same header
    
    header.name <- paste0(i_wild, "|", "mCleavages = ", mCleavages)
    
    
    write.fasta(fina.mutated, as.string = T, names= header.name,
                file.out = "Mutationsfiles1/Mutated peptides mCleavages = 1.fasta", open = "a")
  }
  }
}

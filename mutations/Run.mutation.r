
Run.mutation <- function(path.mutation.table){
  
  mutation.table <- read.mutat.table(path.mutation.table)
  
  erros.df <- data.frame()
  for(i_mutation.table in 1:nrow(mutation.table)){
    
    print (paste0("Loading...", i_mutation.table, " From...", nrow(mutation.table)))
    mutation.table. <- mutation.table[i_mutation.table,]
    
    
    protein.seq <- get.protein(mutation.table.$ACCESSION_NUMBER)
    
    if (protein.seq != strsplit(mutation.table.$ACCESSION_NUMBER, ".", fixed = T)[[1]][1]){
      
      if (true.protein(protein.seq, mutation.table.$MUTATION_AA )){
        mutated.seq <- do.mutation(protein.seq, mutation.table.$MUTATION_AA)
        
      }else{
        error.ids <- paste0(mutation.table.$ACCESSION_NUMBER, "|" ,
                            mutation.table.$MUTATION_AA)
        erros.df <- rbind(erros.df , data.frame(error = error.ids))
        #write.csv(error.ids , "Protein.Errors.csv")
      }
      
      if(mutation.done(protein.seq, mutated.seq, mutation.table.$MUTATION_AA )){
        write.protein(protein.seq, mutation.table.$GENE_NAME,
                      mutation.table.$ACCESSION_NUMBER,
                      mutation.table.$MUTATION_AA,
                      mutation.table.$GENE_CDS_LENGTH,
                      "normal.proteins.fasta")
        
        write.protein(mutated.seq, mutation.table.$GENE_NAME,
                      mutation.table.$ACCESSION_NUMBER,
                      mutation.table.$MUTATION_AA,
                      mutation.table.$GENE_CDS_LENGTH,
                      "mutated.proteins.fasta")
      }else{
        error.ids <- paste0(mutation.table.$ACCESSION_NUMBER, "|" ,
                            mutation.table.$MUTATION_AA)
        erros.df <- rbind(erros.df , data.frame(error = error.ids))
      }
    }else{
      error.ids <- paste0(mutation.table.$ACCESSION_NUMBER, "|" ,
                          mutation.table.$MUTATION_AA)
      erros.df <- rbind(erros.df , data.frame(error = error.ids))
    }
    
  }
  write.csv(erros.df, "Protein.Errors.csv")
}
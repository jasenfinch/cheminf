
checkEntries <- function(entries){
  necessary_names <- c('ID','NAME','SMILES')
  
  entry_names <- colnames(entries)
  
  presence <- necessary_names %in% entry_names
  
  if (FALSE %in% presence) {
    stop(paste0('The table of metabolite entries should contain the following column names: ',
          paste(necessary_names,collapse = ', '),
          '.'))
  }
  
  invisible()
}

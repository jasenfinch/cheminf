#' SMARTS substructure search
#' @description SMARTS substructure searching for SMILES chemical structures.
#' @param SMILES a valid SMILES structure string
#' @param SMARTS a valid SMARTS symbol
#' @examples
#' smartsSearch("C[C@@H](C(=O)O)N","[OX2H]")
#' @importFrom ChemmineOB smartsSearch_OB
#' @export

smartsSearch <- function(SMILES,SMARTS){
  
  if (length(SMILES) > 1){
    stop('Argument `SMILES` should be of length 1.',
         call. = FALSE)
  }
  
  molRefs = forEachMol("SMILES",SMILES,identity)
  smartsSearch_OB(molRefs,SMARTS)
}

#' Chemical descriptors
#' @param SMILES a character vector of valid SMILES
#' @importFrom dplyr mutate relocate
#' @importFrom purrr map
#' @importFrom furrr future_map_dfr future_map_int future_map_chr furrr_options future_map_dbl
#' @export
#' @examples
#' data(amino_acids)
#' chemicalDescriptors(amino_acids$SMILES)

chemicalDescriptors <- function(SMILES){
  
  desc <- c('HBA1',
            'HBA2',
            'HBD',
            'logP',
            'TPSA'
  )
  
  descs <- SMILES %>% 
      future_map_dfr(~{
        molRefs = forEachMol("SMILES",.x,identity)
        
        prop_OB(molRefs) %>% 
          {.[,colnames(.) %in% desc]} %>% 
          as_tibble() %>% 
          mutate(MF = smilesToMF(.x))
      },
      .options = furrr_options(seed = TRUE)) %>% 
    mutate(SMILES = SMILES) %>% 
    relocate(SMILES:MF,.before = 'HBA1')
  
  Fgroups <- tibble(Name = c('Negative_Charge',
                             'Positive_Charge',
                             'NHH',
                             'OH',
                             'COOH',
                             'COO'),
                    String = c('[-]',
                               '[+]',
                               "[NX3;H2]",
                               "[OX2H]",
                               "[CX3](=O)[OX2H1]",
                               "[CX3](=O)[OX1H0-]")
  )
  
  
  groups <- Fgroups %>%
    split(1:nrow(Fgroups)) %>%
    map(~{
      string <- .$String
      g <- future_map_int(SMILES,~{
        s <- .
        s %>%
          smartsSearch(string)
      },
      .options = furrr_options(seed = TRUE)) %>%
        as_tibble()
      names(g) <- .$Name
      return(g)
    }) %>%
    bind_cols()
  
  desc <- bind_cols(descs,groups) %>%
    mutate(Total_Charge = -Negative_Charge + Positive_Charge,
           MF = future_map_chr(SMILES,smilesToMF,
                               .options = furrr_options(seed = TRUE)),
           `Accurate_Mass` = future_map_dbl(SMILES,smilesToAccurateMass,
                                            .options = furrr_options(seed = TRUE)) %>% round(5)
           ) %>% 
    relocate(Accurate_Mass,.after = 'MF')
  
  return(desc)
}

#' Calculate ionisation products for a SMILES structure
#' @description Calculate electrospray ionisation product *m/z* for a SMILES structure.
#' @param SMILES a valid SMILES string
#' @param adduct_rules_table table of adduct rules. Defaults to adducts()
#' @return A tibble containing the mass to charge ratios of electrospray ionisation products for the specified SMILES structure.
#' @examples 
#' ionisationProducts(amino_acids$SMILES[1])
#' @importFrom dplyr ungroup rowwise bind_cols
#' @importFrom tibble as_tibble
#' @importFrom mzAnnotation adduct_rules calcMZ adductTransformMF count.elements
#' @export

ionisationProducts <- function(SMILES,adduct_rules_table = adduct_rules()){
  
  desc <- chemicalDescriptors(SMILES)
  
  adduct_rules_table %>%
    select(Name,Rule) %>%
    bind_cols(desc) %>%
    rowwise() %>%
    mutate(Possible = eval(parse(text = Rule)),
           `m/z` = calcMZ(Accurate_Mass,
                          Name,
                          adduct_rules_table = adduct_rules_table),
           MF = adductTransformMF(MF,
                                  Name,
                                  adduct_rules_table = adduct_rules())) %>%
    select(Adduct = Name,`m/z`,MF,Possible) %>%
    ungroup()
}

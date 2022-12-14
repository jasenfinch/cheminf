#' Convert chemical notation
#' @description convert between SMILES and Inchi and to InchiKey
#' @param input a valid SMILE or Inchi
#' @param input_type either "smiles" or "inchi", denoting the input type
#' @param output_type either "smiles", "inchi" or "inchikey", denoting the output type
#' @examples
#' convert("C[C@@H](C(=O)O)N",'SMILES','INCHI')
#' @importFrom ChemmineOB convertFormat
#' @importFrom stringr str_remove_all
#' @export

convert <- function(input, 
                    input_type = c('SMILES','INCHI'), 
                    output_type = c('INCHI','INCHIKEY','SMILES')) {

  input_type <- match.arg(
    input_type,
    choices = c('SMILES','INCHI')
  )

  output_type <- match.arg(
    output_type,
    choices = c('INCHI','INCHIKEY','SMILES')
  )
  
  output <- convertFormat(input_type,output_type,input,
                          options = data.frame()) %>% 
    str_remove_all('\\n')
  
  return(output)
}

#' Convert a SMILES structure to a molecular formula
#' @description Convert a SMILES structure to a molecular formula.
#' @param SMILES a valid SMILES structure string
#' @examples
#' smilesToMF("C[C@@H](C(=O)O)N")
#' @importFrom ChemmineOB prop_OB
#' @export

smilesToMF <- function(SMILES){
  
  if (length(SMILES) > 1){
    stop('Argument `SMILES` should be of length 1.',
         call. = FALSE)
  }
  
  molRefs = forEachMol("SMILES",SMILES,identity)
  
  prop_OB(molRefs)$formula
}

#' Convert a SMILES structure to accurate mass
#' @description convert a smile to an accurate mass
#' @param SMILES a valid SMILES structure string
#' @examples
#' smilesToAccurateMass("C[C@@H](C(=O)O)N")
#' @importFrom ChemmineOB forEachMol exactMass_OB
#' @export

smilesToAccurateMass <- function(SMILES){
  
  if (length(SMILES) > 1){
    stop('Argument `SMILES` should be of length 1.',
         call. = FALSE)
  }
  
  molRefs = forEachMol("SMILES",SMILES,identity)
  exactMass_OB(molRefs)
}

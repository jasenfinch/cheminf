#' Convert chemical structure notation
#' @description Convert between SMILES and InChI and to InchiKey chemical structure notations.
#' @param input a valid SMILE or InChI
#' @param input_type either `"SMILES"` or `"INCHI"`, denoting the input type
#' @param output_type either `"SMILES"`, `"INCHI"` or `"INCHIKEY"`, denoting the output type
#' @details This functionality is not currently supported on Windows.
#' @return The converted chemical structure.
#' @examples
#' if (Sys.info()["sysname"] != 'Windows'){
#'   convert("C[C@@H](C(=O)O)N",'SMILES','INCHI')
#' }
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

#' Convert a SMILES structure to a molecular formula or accurate mass
#' @rdname smilesTo
#' @description Convert a SMILES structure to a molecular formula or accurate mass.
#' @param SMILES a valid SMILES structure string
#' @examples
#' ## Convert a SMILES structure to a molecular formula
#' smilesToMF("C[C@@H](C(=O)O)N")
#' 
#' ## Convert a SMILES structure to an accurate mass
#' smilesToAccurateMass("C[C@@H](C(=O)O)N")
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

#' @rdname smilesTo
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

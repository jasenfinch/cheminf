#' MetaboliteDatabase S4 class
#' @rdname MetaboliteDatabase-class
#' @description An S4 class for a metabolite electrospray ionisation database.
#' @slot entries database metabolite entries
#' @slot descriptors chemical descriptors of metabolite database entries
#' @importFrom tibble tibble
#' @export

setClass(
  'MetaboliteDatabase',
  slots = list(
    entries = 'tbl_df',
    descriptors = 'tbl_df'
  ),
  prototype = list(
    entries = tibble(),
    descriptors = tibble()
  ))

setMethod('show',signature = 'MetaboliteDatabase',
          function(object){
            entries <- object %>%
              entries() %>%
              nrow()
            cat('\nMetaboliteDatabase object containing ',entries,' entries\n\n',sep = '')
          }
)

#' Create a metabolite ionisation database
#' @description Construct a metabolite electrospray ionisation database.
#' @param entries tibble containing accession information. See details.
#' of the table containing the accession information within the SQL database.
#' @details 
#' The `entries` tibble should contain at least the following columns:
#' * ID - Metabolite identification numbers
#' * NAME - Metabolite names
#' * SMILES - The metabolite SMILES structures
#' @return An object of S4 class `MetaboliteDatabase`
#' @examples 
#' db <- metaboliteDB(amino_acids)
#' db
#' @importFrom dplyr tbl
#' @importFrom purrr map_chr
#' @importFrom methods new
#' @export

metaboliteDB <- function(entries){
  
  checkEntries(entries)
  
  metabolite_descriptors <- chemicalDescriptors(entries$SMILES)
  
  entries <- as_tibble(entries)
  
  db <- new(
    'MetaboliteDatabase',
    entries = entries,
    descriptors = bind_cols(
      select(entries,ID),
      as_tibble(metabolite_descriptors)
    )
  )
  
  return(db)
}

#' Metabolite database utility methods
#' @rdname utilities
#' @description Utilities for working with metabolite databases.
#' @param db S4 object of class `MetaboliteDatabase`
#' @param IDs a numeric vector of entry IDs
#' @param mf a molecular formula to filter
#' @param rule a filtering expression
#' @param lower lower mass boundary
#' @param upper upper mass boundary
#' @details 
#' * `entries` - return the metabolite entries
#' * `descriptors` - return the metabolite entry chemical descriptors
#' * `nEntries` - return the number of metabolite entries
#' * `filterEntries` - filter the metabolite entries based on a vector of metabolite IDs
#' * `filterMR` - filter the metabolite entries based on a mass range
#' * `filterER` - filter the metabolite entries based on an element frequency rule
#' * `filterIP` - filter the metabolite entries based on an ionisation product rule
#' * `filterMF` - filter the metabolite entries based on a molecular formula
#' @return A tibble containing metabolite entry information, the number of metabolite entires or an 
#' S4 object of class `MetaboliteDatabase` depending on the method used.
#' @examples 
#' ## Create a metablite ionisation database using the example amino acid data
#' metabolite_database <- metaboliteDB(amino_acids)
#' 
#' ## Return the entries
#' entries(metabolite_database)
#' 
#' ## Return the chemical descriptors
#' descriptors(metabolite_database)
#' 
#' ## Return the number of database entries
#' nEntries(metabolite_database)
#' 
#' ## Filter database entries
#' filterEntries(metabolite_database,c(1:5))
#' 
#' ## Filter database using a mass range
#' filterMR(metabolite_database,100,120)
#' 
#' ## Filter the database by an element frequency rule
#' filterER(metabolite_database,C > 2)
#' 
#' ## Filter the database by an ionisation product rule
#' filterIP(metabolite_database,HBA2>0 & Total_Charge==0)
#' 
#' ## Filter a database by a molecular formula
#' filterMF(metabolite_database,"C3H7NO2")
#' @export

setGeneric('entries',function(db)
  standardGeneric('entries'))

#' @rdname utilities

setMethod('entries',signature = 'MetaboliteDatabase',
          function(db){
            db@entries
          }
)

setGeneric('entries<-',function(db,value)
  standardGeneric('entries<-'))

setMethod('entries<-',signature = 'MetaboliteDatabase',
          function(db,value){
            db@entries <- value
            return(db)
          }
)

#' @rdname utilities
#' @export

setGeneric('descriptors',function(db) 
  standardGeneric('descriptors'))

#' @rdname utilities

setMethod('descriptors',signature = 'MetaboliteDatabase',
          function(db){
            db@descriptors
          }
)

setGeneric('descriptors<-',function(db,value)
  standardGeneric('descriptors<-'))

setMethod('descriptors<-',signature = 'MetaboliteDatabase',
          function(db,value){
            db@descriptors <- value
            return(db)
          }
)

#' @rdname utilities
#' @export

setGeneric("nEntries", function(db) {
  standardGeneric("nEntries")
})

#' @rdname utilities

setMethod('nEntries',signature = 'MetaboliteDatabase',
          function(db){
          nrow(entries(db))
          }
)

#' @rdname utilities
#' @export

setGeneric("filterMR", function(db,lower,upper) {
  standardGeneric("filterMR")
})

#' @rdname utilities

setMethod('filterMR',signature = 'MetaboliteDatabase',
          function(db,lower,upper){
            desc <- descriptors(db) %>%
              filter(Accurate_Mass > lower,
                     Accurate_Mass < upper)
            
            acc <- entries(db) %>%
              filter(SMILES %in% desc$SMILES)
            
            descriptors(db) <- desc
            entries(db) <- acc
            return(db)
          }
)

#' @rdname utilities
#' @export

setGeneric("filterER", function(db,rule) {
  standardGeneric("filterER")
})

#' @rdname utilities
#' @importFrom rlang enexpr expr_text
#' @importFrom stringr str_extract_all

setMethod('filterER',signature = 'MetaboliteDatabase',
          function(db,rule){
            
            rule <- enexpr(rule)
            
            ef <- elementFreq(db)
            
            if (all(str_extract_all(expr_text(rule),'[:alpha:]')[[1]] %in% 
                    colnames(ef))) {
              ef <- ef %>%
                filter(!!rule)   
            } else {
              ef[0,]
            }
           
            db <- filterEntries(db,
                                IDs = ef)
            
            return(db)
          }
)

#' @rdname utilities
#' @export

setGeneric("filterIP", function(db,rule) {
  standardGeneric("filterIP")
})

#' @rdname utilities

setMethod('filterIP',signature = 'MetaboliteDatabase',
          function(db,rule){
            
            rule <- enexpr(rule)
            
            desc <- descriptors(db) %>%
              filter(!!rule)
            
            db <- filterEntries(
              db,
              desc$ID)
            
            return(db)
          }
)

#' @rdname utilities
#' @export

setGeneric("filterEntries", function(db,IDs) {
  standardGeneric("filterEntries")
})

#' @rdname utilities

setMethod('filterEntries',signature = 'MetaboliteDatabase',
          function(db,IDs){
            entries(db) <- db %>%
              entries() %>%
              filter(ID %in% IDs)
            
            descriptors(db) <- db %>%
              descriptors() %>%
              filter(ID %in% IDs)
            
            return(db)
          }
)

#' @rdname utilities
#' @export

setGeneric('filterMF', function(db,mf){
  standardGeneric('filterMF')
})

#' @rdname utilities

setMethod('filterMF',signature = 'MetaboliteDatabase',
          function(db,mf){
            mfs <- descriptors(db) %>% 
              filter(MF == mf)
            
            db <- filterEntries(
              db,
              mfs$ID
            )
            
            return(db)
          }
)

#' Ionisation product searches of a metabolite database 
#' @rdname products
#' @description Methods for metabolite database ionisation product searches
#' @param db S4 object of class `MetaboliteDatabase`
#' @param id Database entry ID
#' @param mz mass to charge ratio to search
#' @param adduct the search adduct as available in `mzAnnotation::adduct_names()`
#' @param ppm the parts per million error search threshold
#' @param isotope the search isotope as available in `mzAnnotation::isotope_names()`
#' @param adduct_rules_table table containing adduct formation rules. Defaults to `adduct_rules()`.
#' @param isotope_rules_table table containing isotope rules. Defaults to `isotope_rules()`.
#' @details 
#' * `PIPsearch` - perform a putative ionisation product search on a metabolite database
#' * `calcAdducts` - calculate adduct mass to charge ratios for a metabolite database entry
#' @return 
#' A tibble containing the relevant product search results depending on the method used.
#' @examples 
#' metabolite_database <- metaboliteDB(amino_acids)
#' 
#' ## Perform a putative ionisation product search
#' PIPsearch(
#'   metabolite_database,
#'   mz = 133.03358,
#'   adduct = '[M-H]1-',
#'   isotope = '13C')
#' 
#' ## Calculate adduct m/z for a database entry
#' calcAdducts(metabolite_database,1)
#' 
#' @importFrom mzAnnotation calcM isotope_rules ppmRange ppmError
#' @importFrom tibble deframe
#' @importFrom rlang parse_expr
#' @importFrom dplyr left_join right_join
#' @export

setGeneric('calcAdducts',function(db,id,adduct_rules_table = adduct_rules())
  standardGeneric('calcAdducts'))

#' @rdname products

setMethod('calcAdducts',signature = 'MetaboliteDatabase',
          function(db,id,adduct_rules_table = adduct_rules()){
            
            smiles <- db %>%
              descriptors() %>%
              filter(ID == id) %>%
              select(SMILES) %>%
              deframe()
            
            smiles %>%
              ionisationProducts(adduct_rules_table = adduct_rules_table)
          })

#' @rdname products
#' @importFrom dplyr bind_rows select filter
#' @export

setGeneric('PIPsearch',function(db,
                                mz,
                                adduct,
                                ppm = 6,
                                isotope = NA, 
                                adduct_rules_table = adduct_rules(),
                                isotope_rules_table = isotope_rules()) 
  standardGeneric('PIPsearch')
)

#' @rdname products
#' @importFrom rlang expr

setMethod('PIPsearch',signature = 'MetaboliteDatabase',
          function(db,
                   mz,
                   adduct,
                   ppm = 6,
                   isotope = NA, 
                   adduct_rules_table = adduct_rules(),
                   isotope_rules_table = isotope_rules()){
            M <- calcM(mz,
                       adduct = adduct,
                       isotope = isotope,
                       adduct_rules_table = adduct_rules_table,
                       isotope_rules_table =  isotope_rules_table)
            mr <- ppmRange(M,ppm)
            
            res <- db %>%
              filterMR(mr$lower,mr$upper)
            
            if (!is.na(isotope) & nEntries(res) > 0) {
              isoRule <- isotope_rules_table$Rule[isotope_rules_table$Isotope == isotope] %>% 
                parse_expr()
              
              res <- res %>%
                filterER(!!isoRule)
            }
            
            addRule <- adduct_rules_table$Rule[adduct_rules_table$Name == adduct] %>% 
              parse_expr()
            
            res <- res %>%
              filterIP(!!addRule)
            
            res <- res %>%
              {left_join(entries(.),
                         descriptors(.),
                         by = c("ID", "SMILES"))} %>%
              select(ID:Accurate_Mass) %>%
              mutate(Isotope = isotope,
                     Adduct = adduct,
                     `Measured m/z` = mz,
                     `Theoretical m/z` = calcMZ(Accurate_Mass,adduct,isotope),
                     `PPM Error` = ppmError(`Measured m/z`,`Theoretical m/z`)
              ) 
            return(res)
          })


setGeneric("elementFreq", function(db) {
  standardGeneric("elementFreq")
})

#' @importFrom dplyr everything

setMethod('elementFreq',signature = 'MetaboliteDatabase',
          function(db){
            MFs <- db %>%
              descriptors() %>%
              .$MF %>%
              unique() %>%
              map(~{
                mf <- .
                mf %>%
                  count.elements() %>%
                  as.list() %>%
                  as_tibble()
              })
            names(MFs) <- db %>%
              descriptors() %>%
              .$MF %>% 
              unique()
            MFs <- MFs %>% 
              bind_rows(.id = 'MF') %>%
              right_join(db %>%
                           descriptors() %>%
                           select(ID,MF), 
                         by = "MF",
                         multiple = 'all') %>%
              select(ID,everything())
            return(MFs)
          })

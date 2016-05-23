
#' Find Duplicates Mix
#'
#' looks for duplicates in binary mixture codes
#'
#' @param data.object input data frame
#' @param case
#' @param molsingle
#' @param conc0
#'
#' @return
#' @export
#'
#' @examples
FindDuplicates_mix <- function(data.object,
                           case = "upper",    # transform to "upper", "lower" or "common" case, (common not implemented yet); or "test" only
                           molsingle = F,     # find misspelled MOL recordings (only for single recordings, DONT USE WITH TWINPAL), greps for "mol" and replaces all hits with "mol"
                           #usualsuspects = T, # correct for some ususal made misspellings like BJOT BIOT, ALOT AIOT
                           conc0 = T         # set NOConc for "AIR" and "MOL" to 0

                           )
{

  # remove spaces
  spaces <- levels(data.object$TOdour)[grep("[ ]",levels(data.object$TOdour))]
  if(length(spaces) > 0)
    {
    cat("found odor names containing spaces (");cat(spaces);cat(") now correcting...\n")
    data.object$TOdour <- factor(gsub("[ ]","",data.object$TOdour))
    }

  duplicates <- levels(data.object$TOdour)[duplicated(toupper(levels(data.object$TOdour)))]

  if (molsingle == T)
    {
      hits <- levels(data.object$TOdour)[grep("mol_mol",levels(data.object$TOdour), ignore.case = T)]
      cat(paste(length(hits), "hits for \"mol_mol\": (")); cat(hits); cat(") renaming all to \"mol_mol\"\n")
      levels(data.object$TOdour)[grep("mol_mol",levels(data.object$TOdour), ignore.case = T)] <- "mol_mol"
    }


  if (case == "upper")  levels(data.object$TOdour) <- toupper(levels(data.object$TOdour))
  if (case == "lower")  levels(data.object$TOdour) <- tolower(levels(data.object$TOdour))
  if (case == "common") print("not implemented yet, please choose either 'upper' or 'lower' as option")
  #if (case == "test")   {cat(paste("Found ", length(duplicates), " duplicates (", sep="")); cat(duplicates); cat(paste(") DID NOTHING only had a look.\n", sep=""))}
  if (case != "test" & case != "common") {cat(paste("Found ", length(duplicates), " duplicates (", sep="")); cat(duplicates); cat(paste(") transformed all to ", toupper(case), "-CASE\n", sep=""))}




  if (conc0 == T)
    {
    ref.pos <- which(data.object$TOdour %in% c("MOL_MOL","AIR_AIR","N2ONLY_N2ONLY"))
    wrong <- ref.pos[which(data.object$NOConc[ref.pos] != 0)]
    data.object$NOConc[wrong] <- 0
    cat("Found and corrected ");cat(length(wrong)); cat(" concentrations different from 0 for \"MOL_MOL\", \"AIR_AIR\" or \"N2ONLY_N2ONLY\" ")
    }

  return(data.object)
}




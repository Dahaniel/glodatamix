#' find duplicated odor names
#'
#' checks 4 letter odor codes for duplications
#'
#' @param data.object input data in gloDatamix format
#' @param case transform all to "upper", "lower" or "common" case, (common not implemented yet); or "test" only
#' @param molsingle find misspelled MOL recordings (only for single recordings, DONT USE WITH TWINPAL), greps for "mol" and replaces all hits with "mol"
#' @param usualsuspects correct for some ususal made misspellings like BJOT BIOT, ALOT AIOT
#' @param conc0 set NOConc for "AIR" and "MOL" to 0
#' @param cut.odors remove everything beyond 4 letters from the odor names
#' @param conc.from.bc read concentration from barcode -- !!handle with care!!
#' @param conc.from.odor read concentration from odorfield, assuming is is the last 2 digits -- !!!handle with even more care!!!
#'
#' @author Daniel Münch <daniel@@muench.bio>

findDuplicates <- function(data.object, case = "upper", molsingle = F, usualsuspects = T, conc0 = T, cut.odors = F, conc.from.bc = F, conc.from.odor = F) {

  # remove spaces
  spaces <- levels(data.object$TOdour)[grep("[ ]",levels(data.object$TOdour))]
  if(length(spaces) > 0) {
    cat("found odor names containing spaces ("); cat(spaces); cat(") now correcting...\n")
    data.object$TOdour <- factor(gsub("[ ]","",data.object$TOdour))
    }

  duplicates <- levels(data.object$TOdour)[duplicated(toupper(levels(data.object$TOdour)))]

  if (molsingle == T)  {
    hits <- levels(data.object$TOdour)[grep("mol",levels(data.object$TOdour), ignore.case = T)]
    cat(paste(length(hits), "hits for \"mol\": (")); cat(hits); cat(") renaming all to \"mol\"\n")
    levels(data.object$TOdour)[grep("mol",levels(data.object$TOdour), ignore.case = T)] <- "mol"
    }

  if (case == "upper")  levels(data.object$TOdour) <- toupper(levels(data.object$TOdour))
  if (case == "lower")  levels(data.object$TOdour) <- tolower(levels(data.object$TOdour))
  if (case == "common") print("not implemented yet, please choose either 'upper' or 'lower' as option")
  #if (case == "test")   {cat(paste("Found ", length(duplicates), " duplicates (", sep="")); cat(duplicates); cat(paste(") DID NOTHING only had a look.\n", sep=""))}
  if (case != "test" & case != "common") {
    cat(paste("Found ", length(duplicates), " duplicates (", sep=""))
    cat(duplicates)
    cat(paste(") transformed all to ", toupper(case), "-CASE\n", sep=""))
  }

  if (conc.from.bc ==T) {
    skip <- grep("mol|air|n2only|barcode|missing",data.object$TOdour,ignore.case=T)
    data.object$NOConc[-skip] <- - as.numeric(substr(data.object$TOdour,5,5)[-skip])
    data.object$NOConc[skip] <- 0
    cat("read concentration from barcode, set mol|air|n2only|barcode|missing to 0 - please recheck!\n")
  }

  if (conc.from.odor ==T) {
    refs <- levels(factor(data.object$TOdour[grep("R_",data.object$TOdour)]))
    skip.3 <- which(data.object$TOdour %in% refs[which(nchar(refs) < 8)]) # refs without concentration information
    refs <- refs[which(nchar(refs) == 8)]                               # refs with concentration information
    refs.pos <- which(data.object$TOdour %in% refs)

    skip.1   <- grep("mol|air|n2only|barcode|missing",data.object$TOdour,ignore.case=T)
    skip.2 <- which(nchar(as.character(data.object$TOdour)) == 4)
    skip   <- sort(c(skip.1,skip.2,skip.3,refs.pos))


    data.object$NOConc[-skip] <- as.numeric(substr(data.object$TOdour,5,6)[-skip])
    data.object$NOConc[refs.pos] <- as.numeric(substr(data.object$TOdour,7,8)[refs.pos])
    data.object$NOConc[skip.1] <- 0
    cat("read concentration from odor field, set mol|air|n2only|barcode|missing to 0 - please recheck!\n")
  }


  if (cut.odors == T) {
    refs <- levels(data.object$TOdour)[grep("R_",levels(data.object$TOdour))]
    skip <- c("MISSING", "BARCODE", "N2ONLY", refs)
    dontskip <- which(levels(data.object$TOdour) %ni% skip)

    levels(data.object$TOdour)[dontskip] <- strtrim(levels(data.object$TOdour)[dontskip],4)
    refs.pos <- which(levels(data.object$TOdour) %in% refs) # do only after renaming the first set of levels, as number of levels changes
    levels(data.object$TOdour)[refs.pos] <- strtrim(levels(data.object$TOdour)[refs.pos],6)

    cat("shortened all odorant names except ");cat(skip);cat(" to 4 letters\n")
    cat("shortened ");cat(refs);cat(" to 6 letters\n")
  }

  if (usualsuspects == T) {
    wrong <- c("BJOT", "BLOT", "ALOT", "AJOT", "LUFT", "ELVE", "OKTK")
    right <- c("BIOT", "BIOT", "AIOT", "AIOT", "AIR" , "EIVE", "OCTK")
    us.hits <- numeric(); for(i in 1:length(wrong)){ us.hit <- which(wrong[i] == levels(data.object$TOdour)); us.hits <- c(us.hits, us.hit)}
    us.hits <- levels(data.object$TOdour)[us.hits]
    cat("Found the following ");cat(length(us.hits)); cat(" usual suspects: ");cat(us.hits);cat(" ...now correcting...\n")
    for(i in 1:length(wrong)) levels(data.object$TOdour)[grep(wrong[i],levels(data.object$TOdour), ignore.case=T)] <- right[i]
  }

  if (conc0 == T) {
    ref.pos <- which(data.object$TOdour %in% c("MOL","AIR","N2ONLY"))
    wrong   <- ref.pos[which(data.object$NOConc[ref.pos] != 0)]
    data.object$NOConc[wrong] <- 0
    cat("Found and corrected "); cat(length(wrong)); cat(" concentrations different from 0 for \"MOL\", \"AIR\" or \"N2ONLY\" ")
  }

  return(data.object)
}

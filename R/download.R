##' Download predicted microRNA targets
##' original NAR paper resources:
##' http://genie.weizmann.ac.il/pubs/mir07/catalogs/PITA_targets_mm9_0_0_TOP.tab.gz
##' ftp://ftp.ebi.ac.uk/pub/databases/microcosm/v5/arch.v5.gff.mus_musculus.zip
##' http://www.targetscan.org//mmu_60//mmu_60_data_download/Conserved_Site_Context_Scores.txt.zip
##' updated resources:
##' http://www.targetscan.org/vert_61/vert_61_data_download/Conserved_Site_Context_Scores.txt.zip
##' @title Download predicted microRNA targets
##' @param destDir directory where you want to store the predicted microRNA target files
##' @return NULL
##' @author Maarten van Iterson
##' @export
##' @examples
##' \dontrun{
##' ##requires write access to the installed package directory
##' ##if you don't have define your own directory 
##' library(miRNAmRNA)
##' dir.create(file.path(path.package("miRNAmRNA"),"extdata")) 
##' downloadTargets(file.path(path.package("miRNAmRNA"),"extdata"))}
##' 
downloadTargets <- function(destDir){    
    ##original NAR paper resources
    pita <- "http://genie.weizmann.ac.il/pubs/mir07/catalogs/PITA_targets_mm9_0_0_TOP.tab.gz"
    microcosm <- "ftp://ftp.ebi.ac.uk/pub/databases/microcosm/v5/arch.v5.gff.mus_musculus.zip"
    targetscan <- "http://www.targetscan.org//mmu_60//mmu_60_data_download/Conserved_Site_Context_Scores.txt.zip"
    ##updated resources
    targetscan <- "http://www.targetscan.org/vert_61/vert_61_data_download/Conserved_Site_Context_Scores.txt.zip"

    download.file(pita, file.path(destDir, basename(pita)))
    download.file(microcosm, file.path(destDir, basename(microcosm)))
    unzip(file.path(destDir, basename(microcosm)), exdir=destDir)
    file.remove(file.path(destDir, basename(microcosm)))
    download.file(targetscan, file.path(destDir, basename(targetscan)))
    unzip(file.path(destDir, basename(targetscan)), exdir=destDir)
    invisible(file.remove(file.path(destDir, basename(targetscan))))
}





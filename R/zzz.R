.onAttach <- function(libname, pkgname){
    .BAMCOEFS <- new.env()
    ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), 
                    fields=c("Version", "Date"))
    packageStartupMessage(paste(pkgname, ver[1], "\t", ver[2],
        "\nEstimates based on Solymos et al. 2013",
        "\nUse 'load_BAM_QPAD()' to load estimates"))
    invisible(NULL)
}

.onUnload <- function(libpath){
    if (exists(".BAMCOEFS", envir=.GlobalEnv))
        rm(list=".BAMCOEFS", envir=.GlobalEnv)
    invisible(NULL)
}

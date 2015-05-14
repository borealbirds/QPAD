.onAttach <- function(libname, pkgname){
    ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), 
                    fields=c("Version", "Date"))

    source(system.file("v3/QPAD_v3.R", package="QPAD"))
    packageStartupMessage(paste(pkgname, ver[1], "\t", ver[2]))
    
    invisible(NULL)
}

.onUnload <- function(libpath){
    invisible(NULL)
}


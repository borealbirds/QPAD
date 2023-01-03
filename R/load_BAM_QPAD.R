unload_BAM_QPAD <- function() {
    if (!exists(".BAMCOEFS", envir=.GlobalEnv)) {
        cat("There was nothing to unload.\n",
            "Use 'load_BAM_QPAD()' to load estimates")
    } else {
        rm(list=".BAMCOEFS", envir=.GlobalEnv)
    }
    invisible(NULL)
}

checkBAMestimates <- function() {
    if (!exists(".BAMCOEFS", envir=.GlobalEnv))
        stop("Use 'load_BAM_QPAD()' to load estimates")
    invisible(TRUE)
}

load_BAM_QPAD <- function(version) {
    if (missing(version)) {
        cat("Select version of BAM QPAD estimates\n")
        cat("2 - version 2 (2013-02-26)\n")
        cat("3 - version 3 (2015-05-14)\n")
        cat("4 - version 4 (2023-01-06)\n")
        cat("Enter selection: ")
        version <- readline()
    } else {
        version <- as.character(version)
    }
    while (!(version %in% c("2","3", "4"))) { # available versions
        cat("Value out of range. Enter selection: ")
        version <- readline()
    }
    fn <- system.file(paste0("estimates/QPAD_v", version, ".R"), 
        package="QPAD")
    source(fn)
    cat("BAM QPAD parameter estimates loaded, version", version, "\n")
    invisible(NULL)
}

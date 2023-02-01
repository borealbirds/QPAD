plot_BAM_QPAD <-
function(spp, type=c("AIC", "BIC")){
  
  type <- match.arg(type)

    if (!exists(".BAMCOEFS"))
        stop("Use 'load_BAM_QPAD()' to load estimates")

    SPP <- getBAMspecieslist()
    if (!(spp %in% SPP))
        stop("species info not available")

    jd <- seq(0.35, 0.55, 0.01) # JDAY
    ts <- seq(-0.25, 0.5, 0.01) # TSSR
    ls <- seq(0, 0.25, len=length(jd)) # DSLS
    tm <- c("PC", "1SPT", "1SPM") # TM

    xp1 <- expand.grid(JDAY=jd,
        TSSR=ts)
    xp1$JDAY2 <- xp1$JDAY^2
    xp1$TSSR2 <- xp1$TSSR^2
    xp1$Jday <- xp1$JDAY * 365
    xp1$Tssr <- xp1$TSSR * 24

    xp2 <- expand.grid(DSLS=ls,
        TSSR=ts)
    xp2$DSLS2 <- xp2$DSLS^2
    xp2$TSSR2 <- xp2$TSSR^2
    xp2$Dsls <- xp2$DSLS * 365
    xp2$Tssr <- xp2$TSSR * 24

    Xp1 <- model.matrix(~., xp1)
    Xp2 <- model.matrix(~., xp2)

    if (getBAMversion() < 3) {
        lc <- as.factor(1:5)
        tr <- seq(0, 1, 0.01)
        xq <- expand.grid(LCC=as.factor(lc),
            TREE=tr)
    } else {
        lc <- factor(c("DecidMixed", "Conif", "Open", "Wet"),
            c("DecidMixed", "Conif", "Open", "Wet"))
        tr <- seq(0, 1, 0.1)
        xq <- expand.grid(LCC4=lc, TREE=tr)
        xq$LCC2 <- as.character(xq$LCC4)
        xq$LCC2[xq$LCC2 %in% c("DecidMixed", "Conif")] <- "Forest"
        xq$LCC2[xq$LCC2 %in% c("Open", "Wet")] <- "OpenWet"
        xq$LCC2 <- factor(xq$LCC2, c("Forest", "OpenWet"))
    }
    Xq0 <- model.matrix(~., xq)
    
    t <- seq(0, 10, 0.1)
    r <- seq(0, 4, 0.05)
    
    if(getBAMversion()>3){
      if(as.numeric(bestmodelBAMspecies(spp, type=type, TM=1)$sra > 14)){
        cf <- coefBAMspecies(spp, 15, 0)
        p <- data.frame(PC=sra_fun(t, exp(cf$sra[1])),
                        SPT = sra_fun(t, exp(sum(cf$sra[c(1,2)]))),
                        SPM = sra_fun(t, exp(sum(cf$sra[c(1,3)]))))
      } else {
        cf <- coefBAMspecies(spp, 0, 0)
        p <- data.frame(All = sra_fun(t, exp(cf$sra[1])))
      }
    }
    if(getBAMversion() < 4){
      cf <- coefBAMspecies(spp, 0, 0)
      p <- data.frame(All = sra_fun(t, exp(cf$sra[1])))
    }
    q <- edr_fun(r, exp(cf$edr))
     cfall <- exp(t(sapply(SPP, function(spp)
         unlist(coefBAMspecies(spp, 0, 0)))))

    ## model weights
     if(type=="AIC"){
       wp <- selectmodelBAMspecies(spp)$sra$wAIC
       wq <- selectmodelBAMspecies(spp)$edr$wAIC
     }
    if(type=="BIC"){
      wp <- selectmodelBAMspecies(spp)$sra$wBIC
      wq <- selectmodelBAMspecies(spp)$edr$wBIC
    }
    names(wp) <- rownames(selectmodelBAMspecies(spp)$sra)
    names(wq) <- rownames(selectmodelBAMspecies(spp)$edr)
    nsra <- selectmodelBAMspecies(spp)$sra$nobs[1]
    nedr <- selectmodelBAMspecies(spp)$edr$nobs[1]

    ## covariate effects
    if(getBAMversion() > 3){
      mi <- bestmodelBAMspecies(spp, type=type, TM=1)
    }
    if(getBAMversion() < 4){
      mi <- bestmodelBAMspecies(spp, type=type, TM=0)
    }
    if(as.numeric(mi$sra) > 14){
      mi$sra <- as.numeric(mi$sra)-15
    }
    
    cfi <- coefBAMspecies(spp, mi$sra, mi$edr)
    vci <- vcovBAMspecies(spp, mi$sra, mi$edr)

    Xp <- if (getBAMversion() >= 3 & mi$sra %in% c("9","10","11","12","13","14")) Xp2 else Xp1
    
    if (getBAMversion() < 3)
        colnames(Xp)[1] <- "INTERCEPT"
    Xp <- Xp[,names(cfi$sra),drop=FALSE]
    lphi1 <- drop(Xp %*% cfi$sra)
    pmat <- matrix(exp(lphi1), length(jd), length(ts))
    pmax <- sra_fun(10, max(exp(lphi1)))
    pmat <- sra_fun(3, pmat)
    pmax <- 1

    if (getBAMversion() < 3)
        colnames(Xq0)[1] <- "INTERCEPT"
    Xq <- Xq0[,names(cfi$edr),drop=FALSE]
    ltau1 <- drop(Xq %*% cfi$edr)
    qmat <- matrix(exp(ltau1), length(lc), length(tr))
    qmax <- edr_fun(0.5, max(exp(ltau1)))
    qmat <- edr_fun(1, qmat)
    qmax <- 1

    op <- par(las=1, mfrow=c(3,2))

    barplot(wp, space=0, col=grey(1-wp), border="grey", ylim=c(0,1),
        main=paste0(spp, " (n=", nsra, ") v", getBAMversion()),
        ylab="Model weight", xlab="Time-removal model ID")
    barplot(wq, space=0, col=grey(1-wq), border="grey", ylim=c(0,1),
        main=paste0(spp, " (n=", nedr, ") v", getBAMversion()),
        ylab="Model weight", xlab="Distance sampling model ID")

    plot(t, p[,1], type="n", ylim=c(0,1),
         xlab="Point count duration (min)",
         ylab="Probability of singing")
    if(getBAMversion() < 4){
      lines(t, p$All, col="black", lwd=2)
    }
    if(getBAMversion() > 3){
      lines(t, p$PC, col="black", lwd=2)
      lines(t, p$SPT, col="grey30", lwd=2)
      lines(t, p$SPM, col="grey70", lwd=2)
      legend("bottomright", legend=c("Human point count", "ARU - 1 tag/task", "ARU - 1 tag/minute"), col=c("black", "grey30", "grey70"), lwd=2, cex=0.8)
    }

    plot(r*100, q, type="n", ylim=c(0,1),
         xlab="Point count radius (m)",
         ylab="Probability of detection")
    lines(r*100, q, col=1, lwd=2)
    abline(v=cfall[spp,2]*100, lty=2)
    rug(cfall[,2]*100, side=1, col="grey")
    box()

    xval <- if (mi$sra %in% c("9","10","11","12","13","14"))
        ls*365 else jd*365
    image(xval, ts*24, pmat,
        col = rev(grey(seq(0, pmax, len=12))),
        xlab=ifelse(mi$sra %in% c("9","10","11","12","13","14"),
            "Days since local spring", "Julian days"),
        ylab="Hours since sunrise",
        main=paste("Best model:", mi$sra))
    box()
    if (getBAMversion() < 3)
        image(1:nlevels(xq$LCC), tr*100, qmat,
          col = rev(grey(seq(0, qmax, len=12))), axes=FALSE,
          xlab="Land cover types", ylab="Percent tree cover",
          main=paste("Best model:", mi$edr))
    if (getBAMversion() > 2)
        image(1:nlevels(xq$LCC4), tr*100, qmat,
              col = rev(grey(seq(0, qmax, len=12))), axes=FALSE,
              xlab="Land cover types", ylab="Percent tree cover",
              main=paste("Best model:", mi$edr))
    if (getBAMversion() < 3)
        axis(1, 1:5, c("DConif","DDecid","SConif","SDecid","Open"))
    if (getBAMversion() > 2)
        axis(1, 1:nlevels(xq$LCC4), levels(xq$LCC4))
    axis(2)
    box()

    par(op)
    invisible(NULL)
}


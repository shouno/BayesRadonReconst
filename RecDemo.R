source( 'MyHPEstimate.R' )
source( 'MyRadon.R' )

#library( Matrix )
#library( fields )
library( gplots )


ShowReconst <- function( obs, L=64 ){
    Nth <- nrow( obs )
    Ns <- ncol( obs )

    Lx <- L
    Ly <- L
    Rmin <- -0.5*((2*round(sqrt(2*L^2)/2)+1)-1)
    Slen <- (2*abs(Rmin)+1)
    ds <- Slen / Ns
    dth <- pi / Nth
    thetas <- ( seq(0,Nth-1) * pi / Nth )
    lftseq <- seq( 0, (Ns-1)/2 )
    rgtseq <- rev( seq( 1, Ns - Ns/2 ) )
    skabs <- c( lftseq, rgtseq ) / Slen
    eta <- 0.0001

    wseq <- seq(-Ns/2, Ns/2-1)
    hanw <- 0.5 - 0.5 * cos(2*pi*wseq/Ns)
    cutoff <- 0.2
    cutw <- abs(wseq) > floor( (1-cutoff) * Ns/2 )
    hanw <- cutw * hanw

    tau <- obs
    # avgtau <- mean(array(tau))
    # sdtau <- sd( array(tau) )
    # tau <- (tau-avgtau)/sdtau * 40.0  # 適当にスケーリング
    Ttau <- t(mvfft(t(tau)))

    lnBeta <- 0
    lnGamma <- 0
    lnH <- 0

    #logfname <- sprintf( "logs/ReconstHP.log" )
    #sink( file=logfname )

    ITER <- 500

    for ( t in seq(1, ITER) ){
        Beta <- exp(lnBeta)
        Gamma <- exp(lnGamma)  
        H <- exp(lnH)
      
        xx <- lnP( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma )
        cat( t, Beta, H, Gamma, xx, "\n" )    

        DBeta <- DlnPDBeta( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma )
        lnBeta <-  lnBeta + eta * DBeta * Beta
  
        DGamma <- DlnPDGamma( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma )
        lnGamma <- lnGamma + eta * DGamma * Gamma

        DH <- DlnPDH( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma )
        lnH <- lnH + eta * DH * H
    }

    FBPfilt <- skabs
    F <- (Beta*skabs^2+H)*skabs + Gamma
    Bayesfilt <- (Gamma/F) * skabs
    HANfilt <- hanw * skabs
    #sink()

    # LS filter
    FBPfltd <- t( FBPfilt * t(Ttau) ) 
    gFBP <- Re(1/Ns * t(mvfft(t(FBPfltd), inverse=TRUE)) )
    FBPReconst <- FBPT( Ly, Lx, gFBP )

    # Han filter
    HANfltd <- t( HANfilt * t(Ttau) ) 
    hFBP <- Re(1/Ns * t(mvfft(t(FBPfltd), inverse=TRUE)) )
    HANReconst <- FBPT( Ly, Lx, hFBP )

    # Bayes fiter
    Bayesfltd <- t( Bayesfilt * t(Ttau) ) 
    gBayes <- Re(1/Ns * t(mvfft(t(Bayesfltd), inverse=TRUE)) )
    BayesReconst <- FBPT( Ly, Lx, gBayes )


    par( mfrow=c(1,3) )
    image( t(FBPReconst)[,Ly:1], col=rich.colors(120), axes=FALSE )
    title( main='FBP Reconst' )

    image( t(HANReconst)[,Ly:1], col=rich.colors(120), axes=FALSE )
    title( main='HAN filtered Reconst' )

    image( t(BayesReconst)[,Ly:1], col=rich.colors(120), axes=FALSE )
    title( main='Bayes filtered Reconst' )

    par( mfrow=c(1,1) )

    return( list( FBPReconst=FBPReconst, HANReconst=HANReconst, BayesReconst=BayesReconst ) )
}


load( 'P1Radon.Rdata' )
ShowReconst( tauNormal )

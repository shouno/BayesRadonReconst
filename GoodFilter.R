#
# Initial Environment is moved to Environ.R
#

library( PET )
source( 'MyRadon.R' )
source( 'MyHPEstimate.R' )
source( 'Environ.R' )


Rmax <- 1
Rmin <- -1

Lx <- 256
Ly <- 256
dx <- (Rmax-Rmin)/Lx
dy <- (Rmax-Rmin)/Ly

Nth <- 256
Ns <- 256
Slen <- Rmax-Rmin
ds <- Slen/Ns
dth <- pi / Nth

P1 <- phantom()
P1 <- P1[1:Ly,1:Lx]

#tauTrue    <- RadonT( Nth, Ns, P1, pflg=0, periodic=TRUE )
#tauTrue    <- RadonT( Nth, Ns, P1, pflg=0 )
#TtauTrue    <- t( mvfft( t( tauTrue ) ) )


#tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=0.1, periodic=TRUE )
#tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=1.0, periodic=TRUE )
#tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=0.5, periodic=TRUE )
#tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=2.0, periodic=TRUE )
#tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=0.5 )
#tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=1.0 )
#tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=1.5 )
#tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=2.0 )
tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=2.5 )
#tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=3.0 )
#tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=3.5 )
#tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=4.0 )

TtauNormal  <- t( mvfft( t( tauNormal ) ) ) 

skabs <- c( seq(0,(Ns/2)-1), seq(Ns-Ns/2, 1) )/Slen
filt <- skabs



#fTtauTrue    <- t(filt * t(TtauTrue))
fTtauNormal  <- t(filt * t(TtauNormal))

#gTrue    <- Re( 1/Ns * t(mvfft(t(fTtauTrue), inverse=TRUE)) )
gNormal  <- Re( 1/Ns * t(mvfft(t(fTtauNormal), inverse=TRUE)) )

#recImageTrue    <- FBPT( Ly, Lx, gTrue )
recImageNormal  <- FBPT( Ly, Lx, gNormal )


##
## tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=2.0, periodic=TRUE )
## t: 1000 Energy(Beta <- 0.004792214 , Gamma <- 301.283 , H <- 0.04609314 ):  9703238
##

##
## tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=1.5, periodic=TRUE )
##t: 1000 Energy(Beta <- 0.004974799 , Gamma <- 558.3766 , H <- 0.04460882 ):  18114549 
##

##
## tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=1.0, periodic=TRUE )
## t: 1000 Energy(Beta <- 0.005209066 , Gamma <- 1346.014 , H <- 0.04279941 ):  43902298 
##

##
## tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=0.5, periodic=TRUE )
## t: 1000 Energy(Beta <- 0.00551309 , Gamma <- 5770.693 , H <- 0.04050712 ):  188848859 
##

##
## tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=0.1, periodic=TRUE )
## t: 1000 Energy(Beta <- 0.005849808 , Gamma <- 12530.78 , H <- 0.03823728 ):  410340577 
##

## tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=0.5 )
##t: 1000 Energy(Beta <- 0.008542396 , Gamma <- 16502.64 , H <- 0.03280429 ):  540493021 

## tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=1.0 )
## t: 1000 Energy(Beta <- 0.007799392 , Gamma <- 2757.482 , H <- 0.03616925 ):  90143460 

## tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=1.5 )
## t: 1000 Energy(Beta <- 0.007425619 , Gamma <- 970.1216 , H <- 0.03798115 ):  31601839 

## tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=2.0 )
## t: 1000 Energy(Beta <- 0.00714325 , Gamma <- 496.1433 , H <- 0.03943584 ):  16085182 

## tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=2.5 )
## t: 1000 Energy(Beta <- 0.006949347 , Gamma <- 300.199 , H <- 0.0404853 ):  9674219 

## tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=3.0 )
## t: 1000 Energy(Beta <- 0.006826676 , Gamma <- 200.6026 , H <- 0.04117587 ):  6417728 

## tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=3.5 )
## t: 1000 Energy(Beta <- 0.006758782 , Gamma <- 143.2556 , H <- 0.04157583 ):  4544023 

## tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=4.0 )
## t: 1000 Energy(Beta <- 0.006732937 , Gamma <- 107.3075 , H <- 0.04174843 ):  3370418 

##
## Reconstruction with Posterior 
##

Beta <- 0.006949347
Gamma <- 300.199
H <- 0.0404853

F <- (Beta * skabs^2 + H) * skabs + Gamma

plot( seq(0,Ns-1), filt, type="l", col="red" )
lines( seq(0,Ns-1), filt*(Gamma/F), col="blue" )

Ttau <- TtauNormal
filtered <- t(filt * Gamma/F * t(Ttau) )

g <- Re( 1/Ns * t(mvfft(t(filtered), inverse=TRUE)) )
recImageBayes <- FBPT( Ly, Lx, g )


#save( file="GoodFiltered.Rdata", recImageTrue, recImageNormal, recImageBayes )

colmap <- gray( 0:255/256 )

circ <- MakeCircMask( Lx, R=(Lx/2-2) )

#saveNormalizedTiff( img=recImageTrue*circ, fname="TrueP1.tiff" )
saveNormalizedTiff( img=recImageNormal*circ, fname="NormalP1.sd_2.5.tiff" )
saveNormalizedTiff( img=recImageBayes*circ, fname="BayesP1.sd_2.5.tiff" )

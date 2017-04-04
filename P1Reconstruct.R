source( 'MyRadon.R' )
source( 'MyGaussImage.R' )
library( adimpro )
library( PET )
library( biOps )

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

#
# set half contrast
#
# P1 <- P1 * .5

tauTrue    <- RadonT( Nth, Ns, P1, pflg=0, periodic=TRUE )
tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=1.0, periodic=TRUE )
tauPoisson <- RadonT( Nth, Ns, P1, pflg=2, periodic=TRUE )

TtauTrue    <- t( mvfft( t( tauTrue ) ) )
TtauNormal  <- t( mvfft( t( tauNormal ) ) ) 
TtauPoisson <- t( mvfft( t( tauPoisson ) ) ) 

filt <- c( seq(0,(Ns/2)-1), seq(Ns-Ns/2, 1) ) / Slen
fTtauTrue    <- t(filt * t(TtauTrue))
fTtauNormal  <- t(filt * t(TtauNormal))
fTtauPoisson <- t(filt * t(TtauPoisson))

gTrue    <- Re( 1/Ns * t(mvfft(t(fTtauTrue), inverse=TRUE)) )
gNormal  <- Re( 1/Ns * t(mvfft(t(fTtauNormal), inverse=TRUE)) )
gPoisson <- Re( 1/Ns * t(mvfft(t(fTtauPoisson), inverse=TRUE)) )

recImageTrue    <- FBPT( Ly, Lx, gTrue )
recImageNormal  <- FBPT( Ly, Lx, gNormal )
recImagePoisson <- FBPT( Ly, Lx, gPoisson )

save( tauTrue, tauNormal, tauPoisson,
     TtauTrue, TtauNormal, TtauPoisson,
     gTrue, gNormal, gPoisson,
     recImageTrue, recImageNormal, recImagePoisson, file="P1ReconstPeriodic.Rdata" )

colmap <- gray(0:255/256)

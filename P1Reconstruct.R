source( 'MyRadon.R' )

library( PET )
library( gplots )

Rmax <- 1
Rmin <- -1

Lx <- 64
Ly <- 64
dx <- (Rmax-Rmin)/Lx
dy <- (Rmax-Rmin)/Ly

Nth <- 128
Ns <- 256
Slen <- Rmax-Rmin
ds <- Slen/Ns
dth <- pi / Nth

P1 <- phantom(n=Ly)
P1 <- P1[1:Ly,1:Lx]

#
# set half contrast
#

tauTrue    <- RadonT( Nth, Ns, P1, pflg=0 )
tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=2.0 )
tauPoisson <- RadonT( Nth, Ns, P1, pflg=2 )

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

par(mfrow=c(1,3))
image( t(recImageTrue)[,Ly:1], col=rich.colors(120), axes=FALSE )
title( main="Noiseless Reconst." )
image( t(recImageNormal)[,Ly:1], col=rich.colors(120), axes=FALSE )
title( main="AWGN sd=2.0" )
image( t(recImagePoisson)[,Ly:1], col=rich.colors(120), axes=FALSE )
title( main="Poisson" )

par(mfrow=c(1,1))

source( 'MyRadon.R' )

library( PET )
library( gplots )

Rmax <- 1
Rmin <- -1

Lx <- 64
Ly <- 64
dx <- (Rmax-Rmin)/Lx
dy <- (Rmax-Rmin)/Ly

Nth <- 256
Ns <- 256
Slen <- Rmax-Rmin
ds <- Slen/Ns
dth <- pi / Nth

# using Phantom image
P1 <- phantom(n=Ly)
P1 <- P1[1:Ly,1:Lx]

tauTrue    <- RadonT( Nth, Ns, P1, pflg=0 )
tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=2.0 )
tauPoisson <- RadonT( Nth, Ns, P1, pflg=2 )

# Save the generated sinogram
save( tauTrue, tauNormal, tauPoisson, Lx, Ly, Nth, Ns, file="P1Radon.Rdata" )

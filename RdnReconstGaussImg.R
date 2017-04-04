source( 'MyRadon.R' )
source( 'MyGaussImage.R' )
source( 'Environ.R' )
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

tbeta <- 1.0^2
th <- 0.1
tgamma <- 10^5
P1 <- GaussImage( Ly, tbeta, th )
save( tbeta, th, P1, file="GaussImg.Rdata" )

tauTrue    <- RadonT( Nth, Ns, P1, pflg=0 )
tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=sqrt(1/(8*pi^2*tgamma)) )

TtauTrue    <- t( mvfft( t( tauTrue ) ) )
TtauNormal  <- t( mvfft( t( tauNormal ) ) ) 

filt <- c( seq(0,(Ns/2)-1), seq(Ns-Ns/2, 1) ) / Slen
fTtauTrue    <- t(filt * t(TtauTrue))
fTtauNormal  <- t(filt * t(TtauNormal))

gTrue    <- Re( 1/Ns * t(mvfft(t(fTtauTrue), inverse=TRUE)) )
gNormal  <- Re( 1/Ns * t(mvfft(t(fTtauNormal), inverse=TRUE)) )

recImageTrue    <- FBPT( Ly, Lx, gTrue )
recImageNormal  <- FBPT( Ly, Lx, gNormal )

save( tauTrue, tauNormal, 
     TtauTrue, TtauNormal, 
     gTrue, gNormal, 
     recImageTrue, recImageNormal, file="GaussReconstImg.Rdata" )

colmap <- gray(0:255/256)

ebeta <- tbeta
egamma <- tgamma
eh <- th
skabs <- c( seq(0,(Ns/2)-1), seq(Ns-Ns/2, 1) )/Slen
partF <- (ebeta * skabs^2 + eh) * skabs
F <- partF + egamma
#
myfilt <- filt*(egamma/F)
myfTtauTrue    <- t(myfilt * t(TtauTrue))
myfTtauNormal  <- t(myfilt * t(TtauNormal))
mygTrue    <- Re( 1/Ns * t(mvfft(t(myfTtauTrue), inverse=TRUE)) )
mygNormal  <- Re( 1/Ns * t(mvfft(t(myfTtauNormal), inverse=TRUE)) )
myrecImageTrue    <- FBPT( Ly, Lx, mygTrue )
myrecImageNormal  <- FBPT( Ly, Lx, mygNormal )

old.par <- par( no.readonly = TRUE )
circ <- MakeCircMask( L=Lx )
par( mfrow = c( 2, 2 ) )
plot( filt, type="l", col="red" )
lines( seq(1,Ns), myfilt, col="blue" )
image( circ*P1, col=colmap )
image( circ*recImageNormal, col=colmap )
image( circ*myrecImageNormal, col=colmap )

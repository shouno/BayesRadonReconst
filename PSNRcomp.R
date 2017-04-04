#
#
#

source( 'Environ.R' )

L <- 256
circ <- MakeCircMask( L )

truefname <- 'Image/TrueP1.tiff'
TP1 <- readTiff( truefname )
TP1 <- TP1[,,1]  ## since type is 'rgb'

Dirname <- 'Image'

seqnum <- seq(1:8)
sdseq <- seqnum/2
BayesSN <- rep( 0, length( sdseq ) )
NormalSN <- rep( 0, length( sdseq ) )

for( ss in seqnum ){
  sdval <- sdseq[ss]
  
  compfnameN <- sprintf( '%s/NormalP1.sd_%3.1f.tiff', Dirname, sdval )
  NP1 <- readTiff( compfnameN )
  NP1 <- NP1[,,1]  ## since type is 'rgb'

  compfnameB <- sprintf( '%s/BayesP1.sd_%3.1f.tiff', Dirname, sdval )
  BP1 <- readTiff( compfnameB )
  BP1 <- BP1[,,1]  ## since type is 'rgb'

  NormalSN[ss] <- MaskedPSNR( TP1, NP1, circ )
  BayesSN[ss]  <- MaskedPSNR( TP1, BP1, circ )
}

d <- data.frame( sdseq, NormalSN, BayesSN )

save( d, file='PSNRcomp.Rdata' )

plot( x=d$sdseq, y=d$NormalSN, t='l', col='blue' )
lines( x=d$sdseq, y=d$BayesSN, t='l', col='red' )

pdf( 'PSNRdiff.pdf', family='Japan1')
plot( x=d$sdseq, y=d$NormalSN, t='b', col='blue' )
lines( x=d$sdseq, y=d$BayesSN, t='b', pch=19, col='red' )
dev.off()

ShiftArray <- function( x ){
  len <- length(x)
  return( c( x[(len/2+1):len], x[1:(len/2)] ) );
  remove(len)
}


Len <- 32
RLH <- c( seq(0,Len/2-1), seq(Len/2,1) )
RLh <- fft( RLH, inverse=TRUE )

windx <- seq(-Len/2,Len/2-1)

pdf( "RLHw.pdf", family="Japan1")
plot( windx, ShiftArray(RLH), type='l', xlab="Freq. Index", ylab="" )
points( windx, ShiftArray(RLH) )
dev.off()

pdf( "RLht.pdf", family="Japan1")
plot( windx, ShiftArray(Re(RLh)), type='l', xlab="Pos. Index", ylab="" )
points( windx, ShiftArray(Re(RLh)) )
dev.off()


RhoH <- Len/2
LSH[1:(Len/2)] <- 2*RhoH/pi * abs(sin(pi/2 * seq(0,Len/2-1)/RhoH ) )
LSH[(Len/2+1):Len] <- 2*RhoH/pi * abs(sin(pi/2 * (seq(Len/2,1)-Len)/RhoH ))
LSh <- fft( LSH, inverse=TRUE )

pdf( "LSHw.pdf", family="Japan1")
plot( windx, ShiftArray(RLH), type='l', lty="dashed", xlab="Freq. Index", ylab="" )
lines( windx, ShiftArray(LSH) )
points( windx, ShiftArray(LSH) )
dev.off()

pdf( "LSht.pdf", family="Japan1")
plot( windx, ShiftArray(Re(LSh)), type='l', xlab="Pos. Index", ylab="" )
points( windx, ShiftArray(Re(LSh)) )
dev.off()

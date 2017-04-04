source( 'Environ.R' )

Rmax <- 1
Rmin <- -1
Ns <- 128
Slen <- (Rmax-Rmin)

skabs <- c( seq(0,(Ns/2)-1), seq(Ns-Ns/2, 1) )/Slen

filt <- skabs
  
##
## Result is described in GoodFilter.R as comment
##
sdseq <- c( 0.5, 1.0, 2.0, 4.0 )
betaseq <- c( 8.54e-3, 7.80e-3, 7.10e-3, 6.73e-3 )
hseq <- c( 3.28e-2, 3.62e-2, 3.94e-2, 4.17e-2 )
gammaseq <- c( 1.65e4, 2.76e3, 4.96e2, 1.07e2 )

pdf("MyFilters.pdf", family="Japan1")

oldpar <- par()
par( mfrow=c(2,2) )
for( i in seq( 1, length( sdseq ) ) ){
  Beta <- betaseq[i]
  H <- hseq[i]
  Gamma <- gammaseq[i]
  txt <- sprintf( "s.d.=%3.1f", sdseq[i] )
  
  F <- (Beta * skabs^2 + H) * skabs + Gamma
  MyFiltHw <- filt*(Gamma/F)
  windx <- seq(-Ns/2,Ns/2-1)

  plot( windx, ShiftArray(filt), lwd=2, type="l", col="red", lty="dashed", xlab="Frequency Index", ylab="" )
  lines( windx, ShiftArray(MyFiltHw), lwd=3, col="blue" )
  text( 48, 32, txt )
  ##points( windx, ShiftArray(MyFiltHw), col="blue" )
  ##dev.off()

  ##pdf( "MyFilterht.pdf", family="Japan1" )
  #MyFiltht <- Re(fft( MyFiltHw, inverse=TRUE))
  #plot( windx, ShiftArray(MyFiltht), type="l", col="blue" )
  #points(windx, ShiftArray(MyFiltht), col="blue" )
  ##dev.off()
}

dev.off()

library(biOps)

MakeCircMask <- function( L, R=L/2 ){
  R2 <- R^2
  ones <- matrix( rep(1,L), ncol=L )
  y <- (seq(1,L)-(L+1)/2)^2 %*% ones
  x <- t(y)
  circ <- as.numeric( x + y <= R2 )
  circ <- matrix( circ, nrow=L, ncol=L )
  return( circ )
}


saveNormalizedTiff <- function( img, fname, tmax=max(img), tmin=min(img) ){
  tiffmax <- 255
  tiffmin <- 0
  
  imgtmp <- imagedata( ((tiffmax-tiffmin)/(tmax-tmin))*(img-tmin) )
  writeTiff( fname, imgtmp )
}


MaskedMSE <- function( mat1, mat2, mask ){
  pixels <- sum(mask)
  
  return( sum( ( mask * ( mat1 - mat2 ) )^2 )/pixels )
}


MaskedPSNR <- function( mat1, mat2, mask, lev=255 ){
  ## lev <- 255 # 8bit max
  maxmat1 <- max( mat1*mask )
  minmat1 <- min( mat1*mask )
  scmat1 <- ( mat1*mask - minmat1 ) / (maxmat1 - minmat1) * lev

  maxmat2 <- max( mat2*mask )
  minmat2 <- min( mat2*mask )
  scmat2 <- ( mat2*mask - minmat2 ) / (maxmat2 - minmat2) * lev

  mse <- MaskedMSE( scmat1, scmat2, mask )

  return( 20 * log10( lev/ sqrt(mse) ) )
}


ShiftArray <- function( x ){
  len <- length(x)
  return( c( x[(len/2+1):len], x[1:(len/2)] ) );
  remove(len)
}


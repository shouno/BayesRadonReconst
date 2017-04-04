GaussImage <- function( L, beta, h )
{
  eps <- 1e-12
  d <- 2 * pi / L

  coord <- c( seq(0,L/2-1), seq(-L/2,-1) )
  ones  <- rep(1,L)
  k1cos <- 2 - 2 * cos(coord*d) %*% t(ones)
  k2cos <- t(k1cos)
  Gk <- k1cos + k2cos

  Gk[1, 1] <- eps
  acc <- beta * Gk + h
  xr <- 1.0/sqrt( 2.0 * acc ) * matrix(rnorm(L*L),nrow=L)
  xi <- 1.0/sqrt( 2.0 * acc ) * matrix(rnorm(L*L),nrow=L)
  
  x <- matrix( complex( real=xr, imaginary=xi), nrow=L )

  x1 <- x[ seq(2, L/2), seq(2, L/2) ]
  x2 <- x[ seq(2, L/2), seq(L/2+2, L) ];
  x3 <- Conj(x2[ seq(L/2-1,1), seq(L/2-1,1) ]);
  x4 <- Conj(x1[ seq(L/2-1,1), seq(L/2-1,1) ]);
 
  xp <- matrix( x[ 1, seq(2,L/2) ], nrow=1 )
  xm <- matrix( Conj( xp[ 1, seq(L/2-1, 1) ] ), nrow=1 )
  yp <- matrix( x[ seq(2,L/2), 1 ], ncol=1 )
  ym <- matrix( Conj( yp[ seq(L/2-1,1), 1 ] ), ncol=1 )
  xA <- matrix( x[ L/2+1, seq(2,L/2) ], nrow=1 )
  xB <- matrix( Conj( xA[1, seq(L/2-1, 1)] ), nrow=1 )
  yA <- matrix( x[ seq(2,L/2), L/2+1 ], ncol=1 )
  yB <- matrix( Conj( yA[ seq(L/2-1, 1), 1] ), ncol=1 )
  z0 <- matrix( 0.0+0.0i, nrow=1, ncol=1 )
  z1 <- matrix( complex(real=xr[1, L/2+1], imaginary=0), nrow=1, ncol=1 )
  z2 <- matrix( complex(real=xr[L/2+1, 1], imaginary=0), nrow=1, ncol=1 )
  z3 <- matrix( complex(real=xr[L/2+1, L/2+1], imaginary=0), nrow=1, ncol=1 )

  X <- rbind(
             cbind( z0, xp, z1, xm ),
             cbind( yp, x1, yA, x2 ),
             cbind( z2, xA, z3, xB ),
             cbind( ym, x3, yB, x4 )
             )

  x <- Re( FFT2( X, inverse=TRUE) ) /L/L
  return( x )
}


FFT2 <- function( x, inverse=FALSE )
{
  return( t( mvfft( t( mvfft( x, inverse ) ), inverse ) ) )
}

ShiftFFT <- function( x )
{
  L <- ncol( x )
  if( L != nrow( x ) ){
    stop( "nrow(x) != ncol(x)" )
  }

  P1 <- x[ seq(1, L/2), seq(1, L/2) ]
  P2 <- x[ seq(1, L/2), seq(L/2+1, L) ]
  P3 <- x[ seq(L/2+1, L), seq(1, L/2) ]
  P4 <- x[ seq(L/2+1, L), seq(L/2+1, L) ]

  x <- rbind(
             cbind( P4, P3 ),
             cbind( P2, P1 )
             )
}

dyn.load( './Cobj/RadonTrans.so' )

RadonT <- function( Nth, Ns, inimg,
                   left = -1.0, right = 1.0,
                   top = 1.0, bottom = -1.0,
                   pflg = 0, psd = 1.0,
                   nh = 4096,
                   periodic = FALSE
                   )
{
  Nx <- ncol(inimg)
  Ny <- nrow(inimg)
  pmu <- 0.0
  ds <- (right-left)/Ns;
  if( periodic ){
    pp <- 1
  }
  else{
    pp <- 0
  }
  ans <- .C( "RadonInterface",
            arg1 = double(Nth*Ns),
            arg2 = as.integer( Nth ),
            arg3 = as.integer( Ns ),
            arg4 = as.double( ds ),
            arg5 = as.double(as.vector(t(inimg))),
            arg6 = as.integer( nrow( inimg ) ),
            arg7 = as.integer( ncol( inimg ) ),
            arg8 = as.double(left),
            arg9 = as.double(right),
            arga = as.double(top),
            argb = as.double(bottom),
            argc = as.double(pmu),
            argd = as.double(psd),
            arge = as.integer(nh),
            argf = as.integer(pflg),
            argg = as.integer(pp)
            )$arg1
  return( t(matrix(ans, nrow = Ns, ncol = Nth) ) )
}


FBPT <- function( Ny, Nx, g,
                 left = -1.0, right = 1.0,
                 top = 1.0, bottom = -1.0 )
{
  Nth <- nrow( g )
  Ns <- ncol( g )
  ds <- (right - left) / Ns;

  ans <- .C( "ReconstructInterface",
            arg1 = double(Nx*Ny),
            arg2 = as.integer(Ny),
            arg3 = as.integer(Nx),
            arg4 = as.double(left),
            arg5 = as.double(right),
            arg6 = as.double(top),
            arg7 = as.double(bottom),
            arg8 = as.double( as.vector(t(g)) ),
            arg9 = as.integer(Nth),
            arga = as.integer(Ns),
            argb = as.double(ds)
            )$arg1
  return( t(matrix(ans, nrow = Ny, ncol = Nx ) ) )
}

NormalizeImg <- function( mat ) {
  tmax <- max( mat )
  tmin <- min( mat )
  return( 255/(tmax-tmin)*(mat-tmin) )
}


lnP <- function( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma ){

#  lnZnoise <- -0.5 * Nth * Ns * log( Gamma )

#  lnTmp <- log(Beta * (skabs^2) + H)
#  lnTmp[1] <- 0 # 直流成分は常に 0 にしとく
#  lnZpri <- -0.5 * Nth * sum( lnTmp )

#  TtauAbsSums <- colSums( abs(Ttau)^2 )  # Theta 方向に和をとっても良い
#  lnFsk <- log(Fsk)
#  lnZpost <- - 0.5 * Nth * sum( lnFsk )
#             - (4*pi^2)*ds*dth/Ns * sum( Gamma*( 1 - Gamma/Fsk ) * TtauAbsSums )

  a <- lnZpost( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma )
  b <- lnZpri( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma )
  c <- lnZnoise( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma )
#  return( lnZpost - lnZpri - lnZnoise )
  return( a - b - c )
}

lnZnoise <- function( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma ){
  return( -0.5 * Nth * Ns * log( Gamma ) )
}

lnZpri <- function( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma ){
  lnTmp <- log(Beta * (skabs^2) + H)
  lnTmp[1] <- 0 # 直流成分は常に 0 にしとく
  return( -0.5 * Nth * sum( lnTmp ) )
}

lnZpost <- function( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma ){
  TtauAbsSums <- colSums( abs(Ttau)^2 )  # Theta 方向に和をとっても良い
  Fsk <- (Beta * skabs^2 + H)*skabs + Gamma
  return(  - 0.5 * Nth * sum( log(Fsk) ) 
         - (4*pi^2)*dth/Ns * sum( Gamma*( 1 - Gamma/Fsk ) * TtauAbsSums ) )
#         - (4*pi^2)*ds*dth/Ns * sum( Gamma*( 1 - Gamma/Fsk ) * TtauAbsSums ) )
}


DlnZnoiseDGamma <- function( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma ){
  return( -0.5 * Nth * Ns  / Gamma )
}


DlnZpriorDBeta <- function( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma ){
  skabs2 <- skabs^2
  tmp <- skabs2 / (Beta * skabs2 + H)
  tmp[1] <- 0 # 直流成分は0にしておく
  return( -0.5 * Nth * sum( tmp ) )
}

DlnZpriorDH <- function( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma ){
  skabs2 <- skabs^2
  tmp <- 1 / (Beta * skabs2 + H)
  tmp[1] <- 0 # 直流成分は0にしておく
  return( -0.5 * Nth * sum( tmp ) )
}


DlnZpostDBeta <- function( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma ){
  skabs3 <- skabs^3
  TtauAbsSums <- colSums( abs(Ttau)^2 )  # Theta 方向に和をとっても良い
  Fsk <- (Beta * skabs^2 + H)*skabs + Gamma
  return( -0.5 * Nth * sum(skabs3/Fsk) 
         - (4*pi^2)*dth/Ns * Gamma^2 * sum( skabs3/(Fsk^2) * TtauAbsSums )  )
#         - (4*pi^2)*ds*dth/Ns * Gamma^2 * sum( skabs3/(Fsk^2) * TtauAbsSums )  )
}

DlnZpostDH <- function( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma ){
  TtauAbsSums <- colSums( abs(Ttau)^2 )  # Theta 方向に和をとっても良い
  Fsk <- (Beta * skabs^2 + H)*skabs + Gamma
  return( -0.5 * Nth * sum(skabs/Fsk)
         - (4*pi^2)*dth/Ns * (Gamma^2) * sum( skabs/(Fsk^2) * TtauAbsSums )  )
#        - (4*pi^2)*ds*dth/Ns * (Gamma^2) * sum( skabs/(Fsk^2) * TtauAbsSums )  )
}

DlnZpostDGamma <- function( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma ){
  TtauAbsSums <- colSums( abs(Ttau)^2 )  # Theta 方向に和をとっても良い
  Fsk <- (Beta * skabs^2 + H)*skabs + Gamma
  return( -0.5 * Nth * sum(1/Fsk)
         - (4*pi^2)*dth/Ns * sum( (1-Gamma/Fsk)^2 * TtauAbsSums )  )
#         - (4*pi^2)*ds*dth/Ns * sum( (1-Gamma/Fsk)^2 * TtauAbsSums )  )
}

DlnPDBeta <- function( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma ){
  post <- DlnZpostDBeta( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma )
  prio <- DlnZpriorDBeta( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma )
#  cat( "# prior, post:", prio, post, "\n" )
  return( post - prio )
}

DlnPDH <- function( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma ){
  post <- DlnZpostDH( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma )
  prio <- DlnZpriorDH( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma )
#  cat( "# prior, post:", prio, post, "\n" )
  return( post - prio )
}

DlnPDGamma <- function( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma ){
  post <- DlnZpostDGamma( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma )
  nois <- DlnZnoiseDGamma( Ttau, Nth, Ns, ds, dth, skabs, Beta, H, Gamma )
  return( post - nois )
}

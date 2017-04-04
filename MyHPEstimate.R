#
# Initial Environment is moved to Environ.R
#

lnP <- function( Ttau, Nth, Ns, dth, skabs, Fsk, Beta, H, Gamma ){

  lnZnoise <- -0.5 * Gamma * Nth * Ns

  lnTmp <- log(Beta * skabs^2 + H)
  lnTmp[1] <- 0 # 直流成分は常に 0 にしとく
  lnZpri <- -0.5 * Nth * sum( lnTmp )

  TtauAbsSums <- colSums( abs(Ttau)^2 )  # Theta 方向に和をとっても良い
  lnFsk <- log(Fsk)
  lnZpost <- -0.5 * Nth * sum( lnFsk ) - dth/Ns * sum( Gamma*( 1 - Gamma/Fsk ) * TtauAbsSums )
  
  return( lnZpost - lnZpri - lnZnoise )
}


DlnZnoiseDGamma <- function( Ttau, Nth, Ns, dth, skabs, Fsk, Beta, H, Gamma ){
  return( -0.5 * Nth * Ns  / Gamma )
}

DlnZpriorDBeta <- function( Ttau, Nth, Ns, dth, skabs, Fsk, Beta, H, Gamma ){
  skabs2 <- skabs^2
  tmp <- skabs2 / (Beta * skabs2 + H)
  tmp[1] <- 0 # 直流成分は0にしておく
  return( -0.5 * Nth * sum( tmp ) )
}

DlnZpriorDH <- function( Ttau, Nth, Ns, dth, skabs, Fsk, Beta, H, Gamma ){
  skabs2 <- skabs^2
  tmp <- 1 / (Beta * skabs2 + H)
  tmp[1] <- 0 # 直流成分は0にしておく
  return( -0.5 * Nth * sum( tmp ) )
}


DlnZpostDBeta <- function( Ttau, Nth, Ns, dth, skabs, Fsk, Beta, H, Gamma ){
  skabs3 <- skabs^3
  TtauAbsSums <- colSums( abs(Ttau)^2 )  # Theta 方向に和をとっても良い
  return( -0.5 * Nth * sum(skabs3/Fsk) - (2*pi)^2*dth/Ns * Gamma^2 * sum( skabs3/(Fsk^2) * TtauAbsSums )  )
}

DlnZpostDH <- function( Ttau, Nth, Ns, dth, skabs, Fsk, Beta, H, Gamma ){
  TtauAbsSums <- colSums( abs(Ttau)^2 )  # Theta 方向に和をとっても良い
  return( -0.5 * Nth * sum(skabs/Fsk) - (2*pi)^2*dth/Ns * Gamma^2 * sum( skabs/(Fsk^2) * TtauAbsSums )  )
}

DlnZpostDGamma <- function( Ttau, Nth, Ns, dth, skabs, Fsk, Beta, H, Gamma ){
  TtauAbsSums <- colSums( abs(Ttau)^2 )  # Theta 方向に和をとっても良い
  return( -0.5 * Nth * sum(1/Fsk) -  (2*pi)^2*dth/Ns * sum( (1-Gamma/Fsk)^2 * TtauAbsSums )  )
}

DlnPDBeta <- function( Ttau, Nth, Ns, dth, skabs, Fsk, Beta, H, Gamma ){
  post <- DlnZpostDBeta( Ttau, Nth, Ns, dth, skabs, Fsk, Beta, H, Gamma )
  prio <- DlnZpriorDBeta( Ttau, Nth, Ns, dth, skabs, Fsk, Beta, H, Gamma )
#  cat( "# prior, post:", prio, post, "\n" )
  return( post - prio )
}

DlnPDGamma <- function( Ttau, Nth, Ns, dth, skabs, Fsk, Beta, H, Gamma ){
  post <- DlnZpostDGamma( Ttau, Nth, Ns, dth, skabs, Fsk, Beta, H, Gamma )
  nois <- DlnZnoiseDGamma( Ttau, Nth, Ns, dth, skabs, Fsk, Beta, H, Gamma )
  return( post - nois )
}

DlnPDH <- function( Ttau, Nth, Ns, dth, skabs, Fsk, Beta, H, Gamma ){
  post <- DlnZpostDH( Ttau, Nth, Ns, dth, skabs, Fsk, Beta, H, Gamma )
  prio <- DlnZpriorDH( Ttau, Nth, Ns, dth, skabs, Fsk, Beta, H, Gamma )
#  cat( "# prior, post:", prio, post, "\n" )
  return( post - prio )
}

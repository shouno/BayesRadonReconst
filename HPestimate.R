#
# Initial Environment is moved to MyHPEstimate.R
#
library( PET )
library( biOps )
source( 'MyRadon.R' )
source( 'MyHPEstimate.R' )


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

P1 <- phantom()
P1 <- P1[1:Ly,1:Lx]
#tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=0.1, periodic=TRUE )
#tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=1.0, periodic=TRUE )
#tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=0.5, periodic=TRUE )
#tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=2.0, periodic=TRUE )
#tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=0.5 )
#tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=1.0 )
#tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=1.5 )
#tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=2.0 )
#tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=2.5 )
#tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=3.0 )
#tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=3.5 )
tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=4.0 )
TtauNormal  <- t( mvfft( t( tauNormal ) ) ) 



#tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=1.0, periodic=TRUE )
#load( 'P1ReconstPeriodic.Rdata' )
# loaded variable of file="P1ReconstPeriodic.Rdata" )
# tauTrue, tauNormal, tauPoisson                 =====> Sinograms
# TtauTrue, TtauNormal, TtauPoisson              =====> 1D Fourier Transformed
# gTrue, gNormal, gPoisson                       =====> Filtered Images
# recImageTrue, recImageNormal, recImagePoisson  =====> Reconstructed by FBP


thetas <- ( seq(0,Nth-1) * pi / Nth )
skabs <- c( seq(0,(Ns/2)-1), seq(Ns-Ns/2, 1) ) / Slen
filt <- skabs

eta <- 0.00001

Ttau <- TtauNormal

# initial value for Hyper parameters.
lnBeta <- 0
lnGamma <- 0
lnH <- 0

for ( trial in seq(1,1000) ){
  Beta <- exp(lnBeta)
  Gamma <- exp(lnGamma)  
  H <- exp(lnH)
    
  Fk <- (Beta * skabs^2 + H) * skabs + Gamma
  xx <- lnP( Ttau, Nth, Ns, dth, skabs, Fk, Beta, H, Gamma )
  cat( "t:", trial, "Energy(Beta <-", Beta, ", Gamma <-", Gamma, ", H <-", H, "): ", xx, "\n" )
    
  DBeta <- DlnPDBeta( Ttau, Nth, Ns, dth, skabs, Fk, Beta, H, Gamma )
  lnBeta <-  lnBeta + eta * DBeta * Beta
    
  DGamma <- DlnPDGamma( Ttau, Nth, Ns, dth, skabs, Fk, Beta, H, Gamma )
  lnGamma <- lnGamma + eta * DGamma * Gamma

  DH <- DlnPDH( Ttau, Nth, Ns, dth, skabs, Fk, Beta, H, Gamma )
#    cat( "# DH =", DH, "\n" )
  lnH <- lnH + eta * DH * H
}
cat( "t:", trial, "Energy(Beta <-", Beta, ", Gamma <-", Gamma, ", H <-", H, "): ", xx, "\n" )

#
# tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=2.0, periodic=TRUE )
# t: 1    Energy(Beta <- 1 , Gamma <- 1 , H <- 1 ):  -72042.51 
# t: 100  Energy(Beta <- 0.004874861 , Gamma <- 282.2839 , H <- 0.05857081 ):  9082896 
# t: 500  Energy(Beta <- 0.004790266 , Gamma <- 301.3651 , H <- 0.04625332 ):  9705923 
# t: 1000 Energy(Beta <- 0.004792214 , Gamma <- 301.283 , H <- 0.04609314 ):  9703238
#

#
# tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=1.5, periodic=TRUE )
#t: 1    Energy(Beta <- 1 , Gamma <- 1 , H <- 1 ):  -72041.16 
#t: 100  Energy(Beta <- 0.005161583 , Gamma <- 475.1064 , H <- 0.05778482 ):  15391076 
#t: 500 Energy(Beta <- 0.004972841 , Gamma <- 558.6185 , H <- 0.04478154 ):  18122464 
#t: 1000 Energy(Beta <- 0.004974799 , Gamma <- 558.3766 , H <- 0.04460882 ):  18114549 
#

#
# tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=1.0, periodic=TRUE )
# t: 1    Energy(Beta <- 1 , Gamma <- 1 , H <- 1 ):  -72040.19 
# t: 100  Energy(Beta <- 0.005634413 , Gamma <- 825.5476 , H <- 0.05653882 ):  26862958 
# t: 500  Energy(Beta <- 0.005208921 , Gamma <- 1345.231 , H <- 0.04290959 ):  43876657 
# t: 1000 Energy(Beta <- 0.005209066 , Gamma <- 1346.014 , H <- 0.04279941 ):  43902298 
#

#
# tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=0.5, periodic=TRUE )
# t: 1    Energy(Beta <- 1 , Gamma <- 1 , H <- 1 ):  -72039.62 
# t: 100  Energy(Beta <- 0.006207916 , Gamma <- 1305.476 , H <- 0.05525508 ):  42579685 
# t: 500  Energy(Beta <- 0.005602374 , Gamma <- 4519.47 , H <- 0.03988514 ):  147856700 
# t: 1000 Energy(Beta <- 0.00551309 , Gamma <- 5770.693 , H <- 0.04050712 ):  188848859 
#

#
# tauNormal  <- RadonT( Nth, Ns, P1, pflg=1, psd=0.1, periodic=TRUE )
# t: 1 Energy(Beta <- 1 , Gamma <- 1 , H <- 1 ):  -72039.44 
# t: 100 Energy(Beta <- 0.006487406 , Gamma <- 1562.047 , H <- 0.05470817 ):  50983244 
# t: 500 Energy(Beta <- 0.00594528 , Gamma <- 7754.351 , H <- 0.03775313 ):  253842212 
# t: 1000 Energy(Beta <- 0.005849808 , Gamma <- 12530.78 , H <- 0.03823728 ):  410340577 
#

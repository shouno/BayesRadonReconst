library( adimpro )
library( PET )

Lx <- 256
Ly <- 256
P1 <- phantom()
P1 <- P1[1:Ly,1:Lx]

pdf( "P1hist.pdf", family="Japan1")
hist( P1 )
dev.off()


rm ( list=ls() )
source( './efa.R') 

S     <- 20
f_epi <- .05 ## hivert 2021 AJHG f_epi ~~ 20% (h2=.208, eta2=.055 [albet w huge sd])
Ns    <- c( 1e2, 1e3, 1e4 )
h2s   <- c( .1, .5, .8 )

r2fxn <- function(U,Uhat)
  max( cor( c(Uhat), c(U) )^2, cor( c(Uhat[,2:1]), c(U) )^2 )

# expectation using purely random guesses
baseline <- mean( replicate( 1e4, r2fxn(matrix(rnorm(S*2),S,2),matrix(rnorm(S*2),S,2)) ) )

nit <- 20
out <- array( NA, dim=c( 4, length(Ns), length(h2s), nit ) )
for( it in 1:nit )
{
  for( j in seq(along=Ns) )
    for( k in seq(along=h2s) )
  {
    print( j )
    N <- Ns[j]
    h2<- h2s[k]

    G     <- scale( matrix( rnorm(N*S), N, S ) )
    betas <- sqrt(h2/S) * rnorm(S)

    Z <- sample( 0:1, S, rep=T )
    U <- cbind( betas * Z, betas*(1-Z) ) ### optimistic because U1 \perp U2
    P <- G %*% U

    h2x <- h2 * ( f_epi / (1-f_epi) ) 
    lambda <-  sqrt( h2x/mean( as.numeric( P[,1]*P[,2] )^2 ) )

    y <- P[,1] + P[,2] + lambda * P[,1] * P[,2] + rnorm(N) * sqrt(1-h2-h2x)

    efaout  <-  efa( y, G, Utrue=U, lam=lambda, U_init=NA, U=U ) 
    out[1,j,k,it] <- r2fxn( efaout$U, U )

    efaout  <-  efa( y, G, Utrue=U, U_init='ols+mask' ) 
    out[2,j,k,it] <- r2fxn( efaout$U, U )

    efaout  <-  efa( y, G, Utrue=U, U_init='ols+mask', nit=30 ) 
    out[3,j,k,it] <- r2fxn( efaout$U, U )

    efaout  <-  efa( y, G, Utrue=U, U_init='rand' )
    out[4,j,k,it] <- r2fxn( efaout$U, U )
  }

  pdf( 'r2s.pdf', width=3*4, height=4 )
  par( mfrow=c(1,3) )
  for( k in seq(along=h2s) ){
    xs  <- log10(Ns)
    ys  <- apply( out[,,k,], c(1,2), mean, na.rm=T )
    plot( range(xs), 0:1, type='n', axes=F, xlab='N (log scale)', ylab='R2(Uhat,U)' )
    axis( 1, at=xs, lab=Ns ); axis( 2 )
    for( i in 1:4 ){
      points( xs, ys[i,], col=i+1, pch=16, cex=1.5 )
      lines(  xs, ys[i,], col=i+1, lty=1 )
    }
    abline( h=1 )
    abline( h=baseline, col='grey', lwd=2 )
    if(k==1)
    legend( 'topleft'     , cex=1.2, bg='white', fill=c( 1:4+1, 'grey' ), leg=c( 'Init @ Oracle', 'Init @ OLS', 'Init @ OLS,30 iters', 'Init @ 0', 'Pure noise' ) )
    legend( 'bottomright' , cex=1.5, bg='white', leg=paste0( 'h2 = ', h2s[k] ) )
  }
  dev.off()
}

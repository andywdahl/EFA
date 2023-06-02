efa <- function(y,G,U,lam,nit=300,Utrue,U_init){
  if( missing(U) )
    U <- init_U_fxn( y, G, U_init )
  P <- G %*% U
  if( missing( lam ) )
    lam <- coef( lm( ( y - rowSums(P) ) ~ -1 + P[,1]:P[,2] ) )[1]

  out <- list()
  if(!missing( Utrue )){
    out$errs[1] <- max( cor( c(U), c(Utrue) )^2, cor( c(U[,2:1]), c(Utrue) )^2 )
  } else {
    out$errs <- NA
  }
  out$objs <- objfxn(y,G,U,lam,P=P)

  for( it in 1:nit ){
       
    U[,1] <- ols( y-P[,2], G + apply( G, 2, function(g)  lam * P[,2] * g ) )
    P[,1] <- as.numeric( G %*% U[,1] )
      
    U[,2] <- ols( y-P[,1], G + apply( G, 2, function(g)  lam * P[,1] * g ) )
    P[,2] <- as.numeric( G %*% U[,2] )
      
    lam <- ols( y - rowSums(P), P[,1] * P[,2])

    if(!missing( Utrue ))
      out$errs  <- c( out$errs    , max( cor( c(U), c(Utrue) )^2, cor( c(U[,2:1]), c(Utrue) )^2 ) )
    out$objs    <- c( out$objs    , objfxn(y,G,U,lam,P=P) )
    out$change  <- c( out$change  , ( rev(out$objs)[2] - rev(out$objs)[1] )/rev(out$objs)[2] )

    if( rev(out$change)[1] < 1e-10 ) break
  }
  list( P=P, U=U, lam=lam, obj=rev(out$objs)[1], out=out )
}

objfxn <- function(y,G,U,lam,P=G %*% U )
  mean(( y - ( rowSums(P) + lam*P[,1]*P[,2] ) )^2)

ols <- function(y,X)
  as.numeric( solve( t(X) %*% X ) %*% ( t(X) %*% y ) ) 

init_U_fxn <- function( y, G, U_init=c( 'ols+mask', 'rand' ) ){
  M     <- ncol(G)
  if(  U_init == 'ols+mask' ){
    b_ols <- ols( y, G )
    Z     <- sample( 0:1, M, rep=T )
    U     <- cbind( b_ols * Z, b_ols*(1-Z) )
  } else if( U_init == 'rand' ){
    U   <- cbind( rnorm(M)/1e3, rnorm(M)/1e3 )
  }
  U
}

#ecuacion de von bertalanffy


  
VB<-function(Li,k,t,t0){
  VB<-Li*(1-exp(-k*(t-t0)))
  return(VB)
}
curve(expr = VB(Li=183.2,k=0.3706,t0=-0.2,t=x),from=0, to=15)


clasificar <- function(n){
  x1 <- runif(n = 10000,min =  0,max =  15)
  t1 <- VB(Li = 183.2,
           k = 0.3706,
           t0 = -0.2,
           t = x1)
  x2 <- x1 + 1
  t2 <- VB(Li = 183.2,
           k = 0.3706,
           t0 = -0.2,
           t = x2)
  c1 <- sapply(X=t1, FUN = function(X){
    val <- numeric()
    if (X <= 150){
      val <- "juv"
    } else if (X > 150 & X <= 170){
      val <- "sub" #subadulto
    } else if ( X > 170){
      val <- "adu"
    }
    return(val)
  })
  
  c2 <- sapply(X=t2, FUN = function(X){
    val <- numeric()
    if (X <= 150){
      val <- "juv"
    } else if (X > 150 & X <= 170){
      val <- "sub" #subadulto
    } else if ( X > 170){
      val <- "adu"
    }
    return(val)
  })
  
  S00 <- S01 <- S02 <- S11 <- S12 <- S22 <- numeric()
  
  S00 <- sum(ifelse(c1 == "juv" & c2 == "juv", 1, 0))
  S01 <- sum(ifelse(c1 == "juv" & c2 == "sub", 1, 0))
  S02 <- sum(ifelse(c1 == "juv" & c2 == "adu", 1, 0))
  S11 <- sum(ifelse(c1 == "sub" & c2 == "sub", 1, 0))
  S12 <- sum(ifelse(c1 == "sub" & c2 == "adu", 1, 0))
  S22 <- sum(ifelse(c1 == "adu" & c2 == "adu", 1, 0))
  
  ctab1 <- table(c1)
  
  lst <- list()
  
  for( i in 1:n){
    S00. <- rbeta(1, shape1 = S00, ctab1["juv"] - S00)
    S01. <- rbeta(1, shape1 = S01, ctab1["juv"] - S01)
    S02. <- rbeta(1, shape1 = S02, ctab1["juv"] - S02)
    S11. <- rbeta(1, shape1 = S11, ctab1["sub"] - S11)
    S12. <- rbeta(1, shape1 = S12, ctab1["sub"] - S12)
    S22. <- rbeta(1, shape1 = S22, ctab1["adu"] - S22)
  
  lst[[i]] <- matrix(c(S00., NA, NA,
                       S01., S11., 0,
                       S02., S12., S22.), byrow = TRUE, ncol=3, nrow=3)
  }
  
  
  return(list(
    stochMatriz = lst,
    M = matrix(
      c(S00, NA, NA,
        S01, S11, 0,
        S02, S12, S22),
      byrow = TRUE,
      ncol = 3,
      nrow = 3
    )
  ))
}
(res <- clasificar(3))


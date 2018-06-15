library(popbio)

#-----abundancia por remocion-----

##Cargar paquetes
library(unmarked)


## Leer los datos
juvenil <- read.csv("juvenil.csv", header = FALSE, row.names = NULL)
juvenil <- as.matrix(juvenil)


umf1 <- unmarkedFrameGMM(y= juvenil, numPrimary=1, type="removal")# remocion solo se realizo una vez

m1 <- gmultmix(~1, ~1, ~1, data=umf1, K=30) #no hay covariables 
backTransform(m1, type="lambda") # Individuals per plot
#backTransform(m1, type="phi") # Probability of being avilable
#(p <- backTransform(m1, type="det")) # Probability of detection
#p <- coef(p)
# Multinomial cell probabilities under removal design
#c(p, (1-p) * p, (1-p)^2 * p)
# Or more generally:
#head(getP(m1))
# Empirical Bayes estimates of super-population size
re <- ranef(m1)
#plot(re, layout=c(5,5), xlim=c(-1,20), subset=site%in%1:25)



## Leer los datos Subadultos
Subadulto <- read.csv("C:/localforks/grupal2018/informe/Subadulto.csv", header = FALSE, row.names = NULL)
Subadulto <- as.matrix(Subadulto)


umf2 <- unmarkedFrameGMM(y= Subadulto, numPrimary=1, type="removal")

(m2 <- gmultmix(~1, ~1,~1, data=umf2, K=30)) #no hay covariables 
backTransform(m2, type="lambda") # Individuals per plot
#backTransform(m2, type="phi") # Probability of being avilable
#(p <- backTransform(m1, type="det")) # Probability of detection
#p <- coef(p)
# Multinomial cell probabilities under removal design
#c(p, (1-p) * p, (1-p)^2 * p)
# Or more generally:
#head(getP(m1))
# Empirical Bayes estimates of super-population size
#re <- ranef(m2)
#plot(re, layout=c(5,5), xlim=c(-1,20), subset=site%in%1:25)


#-----ecuacion de von bertalanffy----

VB<-function(Li,k,t,t0){
  VB<-Li*(1-exp(-k*(t-t0)))
  return(VB)
}
curve(expr = VB(Li=183.2,k=0.3706,t0=-0.2,t=x),from=0, to=15)

#-----Matrices que ocupa popbio-----

Matrices <- function(n){
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
  
  lst[[i]] <- matrix(c(S00., 1.007, 1.03,
                       S01., S11., 0,
                       S02., S12., S22.), byrow = TRUE, ncol=3, nrow=3)
  }
  #buscar las fecundidades de sub (los inventamos) adultos y adultos que corresponden a las NA anteriores
  
  return(stochMatriz = lst)
}

(res <- Matrices(3))

#-----Funcion viabilidad poblacional-----

# Funcion para generar viabilidad poblacional

transStochMat <- setRefClass("transStochMat",
   fields = list(
     matriz = "list",
     n0 = "numeric",
     t = "numeric",
     p = "matrix",
     umbral = "numeric",
     probExt = "matrix",
     incluirEtapa = "numeric"
   ),
   methods = list(
     initialize = function(matriz,
                           n0,
                           t,
                           p,
                           umbral,
                           incluirEtapa,
                           probExt
     ) {
       .self$matriz <- matriz
       .self$n0 <- n0
       .self$t <- t
       setanddone <- FALSE
       if (missing(incluirEtapa)){
         .self$incluirEtapa <- rep(1, length(n0))
         setanddone <- TRUE
       } else {
         .self$incluirEtapa <- incluirEtapa
         setanddone <- TRUE
       }
       if (missing(umbral)){
         .self$umbral <- 0.05*(sum(n0))
       } else {
         .self$umbral <- umbral
       }
       if (missing(p)){
         require(popbio)
         .self$p <- stoch.projection(matrices = .self$matriz,
                                     nreps = 5000,
                                     n0 = n0,
                                     tmax = t)
       }
       if (setanddone & missing(probExt)){
         require(popbio)
         .self$probExt <-
           stoch.quasi.ext(
             matrices = matriz,
             n0 = n0,
             Nx = .self$umbral,
             nreps = 5000,
             maxruns = 10,
             tmax = t,
             sumweight = .self$incluirEtapa,
             verbose = FALSE
           )
         
       }
     },
     
     plotN = function(){
       assign("op", par())
       original <- sum(n0*incluirEtapa)
       iters <- rowSums(p%*%incluirEtapa)
       xLabThis <- paste0("Número de individuos en t = ",t)
       
       hist(
         iters,
         main = "Tamaño de población",
         xlab = xLabThis,
         ylab = "Frecuencia",
         lwd = 2, las = 1,
         xlim = c(min(c(original,iters)), max(c(original,iters)))
       )
       abline(v=original,col="gray")
       abline(v=umbral,col="red",lwd=2)
       suppressWarnings(suppressMessages(par(op)))
     },
     
     darR0 = function(){
       require(popbio)
       valR0 <- stoch.growth.rate(matriz)
       return(list(
         approx = exp(valR0$approx),
         sim = exp(valR0$sim),
         simCI = exp(valR0$sim.CI)
       ))
     },
     
     plotExtProb = function(){
       matplot(
         probExt,
         xlab = "Tiempo",
         ylab = "Probabilidad de quasi-extinción",
         type = "l",
         lty = 1,
         col = rainbow(10)
       )
     },
     extProb = function(){
       pop <- p%*%incluirEtapa#solo etapas incluidas
       prob <- mean(ifelse(rowSums(pop) < umbral, 1, 0))
       return(prob)
     }
   )#methods
)#class

#-----Viabilidad poblacional-----
#funcion para calcular abundancia con modelos de rmocion

Area <- 1161252
(n_0 <- floor((2.25 / 100) * Area))
(n_1 <- floor((1.56 / 100) * Area))
n_2 <- 0
Nx <- floor(Area * (0.005))#bichos por metro
NT <- 5
n <- c(n_0, n_1, n_2)

Viabilidadcambute <- transStochMat$new(
  matriz = Matrices(500), n0 = n, 
  t=NT,
  umbral=Nx,# numero de plantas con potencial reproduvtivo
  incluirEtapa=c(0,1,1)
)




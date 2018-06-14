#ecuacion de von bertalanffy

  
VB<-function(Li,k,t,t0){
  VB<-Li*(1-exp(-k*(t-t0)))
  return(VB)
}
VB(Li=183.2, k=0.3706, t=2, t0=0.2)

curve(expr =VB(
      Li=183.2,
      k=0.3706,
      t0=-0.2,
      t=x),from = 0, to=15)

Transform2ab <- function (d_a, c){
  b <- -(c+1)/(c-1)
  a <- (d_a/sqrt(2))*sqrt(1+b^2)
  return( list(
    a = a,
    b = b
  ) )
}

GetLimits <- function (d_a, c)
{
  if (c < 0.0) {
    LL <-  d_a/4/c*sqrt(1+c^2)
    UL <-  10.0
  }
  
  if (c == 0.0) 
  {
    LL  <-  -10.0
    UL <-  10.0
  }
  
  if (c > 0.0) 
  {
    LL <-  -10.0
    UL <-  d_a/4/c*sqrt(1+c^2)
  }
  return( list(
    LL = LL,
    UL = UL
  ) )
  
}



Heaviside <- function (c)
{
  if (c < 0.0) rval <-  0.0 else rval  <-  1.0
  return (rval)
}

FalsePositiveFraction <- function (vc, d_a, c)
{
  
  arg1 <-  -(1-c)*vc-d_a/2*sqrt(1+c^2)
  arg2 <-  -(1-c)*vc+d_a/2/c*sqrt(1+c^2)
  
  rval = pnorm (arg1) + pnorm (arg2) - Heaviside (c)
  
  return (rval)
}



TruePositiveFraction <- function (vc, d_a, c)
{
  arg1 = -(1+c)*vc+d_a/2*sqrt(1+c^2)
  arg2 = -(1+c)*vc+d_a/2/c*sqrt(1+c^2)
  
  rval = pnorm (arg1) + pnorm (arg2) - Heaviside (c)
  
  return (rval)
}




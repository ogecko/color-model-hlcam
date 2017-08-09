## convert XYZ to CCT using McCamy's Approximation https://en.wikipedia.org/wiki/Color_temperature#Approximation
XYZtoCCT <- function(Xw, Yw, Zw) {
  xw = Xw / (Xw + Yw + Zw)
  yw = Yw / (Xw + Yw + Zw)
  n = (xw - 0.3320) / (yw - 0.1858)
  CCT = -449*n^3 + 3525*n^2 - 6823.3*n + 5520.33
  return (CCT)
}

## conversion matrices
## MCAT02 - MCAT02 transformation matrix (calculated using the CIE 1931 2° standard colorimetric observer)
## MH - Hunt–Pointer–Estévez space transformation matrix from XYZc to LMSp
## MHMCAT02 - transformation matrix from LMSc to LMSp
MCAT02 <- matrix(c(0.7328, 0.4296, -0.1624, -0.7036, 1.6975, 0.0061, 0.0030, 0.0136, 0.9834), ncol=3, byrow=TRUE)
MH <- matrix(c(0.38971, 0.68898, -0.07868, -0.22981, 1.18340, 0.04641, 0.0, 0.0, 1), ncol=3, byrow=TRUE)
MHMCAT02 <- MH %*% solve(MCAT02)
MCAT02INV <- solve(MCAT02)
MCAT02INVY <- MCAT02INV[2,]
Heh <- matrix(c(0.0, 100.0, 200.0, 300.0, 400.0, 0.8, 0.7, 1.0, 1.2, 0.8, 20.14, 90.0, 164.25, 237.53, 380.14), ncol=3)

## convert XYZ to LMS space using MCAT02 transform
## returns an (n x 3) matrix
XYZtoLMS <- function(X,Y,Z) {
  XYZ <- matrix(c(X,Y,Z), ncol=3)
  LMS <- t(MCAT02 %*% t(XYZ))
  return (LMS)
}

## convert LMS to XYZ space using inverse MCAT02 transform
## returns an (n x 3) matrix
LMStoXYZ <- function(LMS) {
  XYZ <- t(MCAT02INV %*% t(LMS))
  return (XYZ)
}

## convert LMS to XYZ space using inverse MCAT02 transform
## returns only the Y component
LMStoY <- function(LMS) {
  Y <- t(MCAT02INVY %*% t(LMS))
  return (Y)
}

## calculates the degree of adaption based on Ambient surround condition and La
## A = Ambient surround condition (Average, Dim, Dark)
## D = degree of adaptation (0 for no adaptation, 1 for full adaption)
calculateD <- function(A, La) {
  F[A=='dark'] <- 0.8
  F[A=='dim'] <- 0.9
  F[A=='average'] <- 1.0
  D <- F * (1 - (1 / 3.6) * exp((- La - 42) / 92))
  return (D)
}

## use CAT02 to color adapt a sample from the illumant adopted white to a reference white used in CIECAM02
## D = degree of adaptation (0 for no adaptation, 1 for full adaption)
## LMSs = LMS of the color sample to adapt
## LMSw = LMS of the white point (ie illuminant adopted white)
## LMSc = LMS of the adapted color sample
LMStoLMSc <- function(D, LMSs,LMSw) {
  Yw <- LMStoY(LMSw)
  Lc <- LMSs[,1] * (D * Yw / LMSw[,1] + 1 - D)
  Mc <- LMSs[,2] * (D * Yw / LMSw[,2] + 1 - D)
  Sc <- LMSs[,3] * (D * Yw / LMSw[,3] + 1 - D)
  LMSc <- matrix(c(Lc, Mc, Sc), ncol=3)
  return (LMSc)  
}

## perform post adaption conversion to cone responses using Hunt-Pointer-Esteves
## LMSc = LMS of the adapted color sample
## LMSp = Post Adaption Hunt-Point_estevez space
LMSctoLMSp <- function(LMSc) {
  LMSp = t(MHMCAT02 %*% t(LMSc))
  return (LMSp)
}

## convert to absolute cone responses
## LMSp = Post Adaption Hunt-Point_estevez space
## LMSr = Absolute cone response
LMSptoLMSr <- function(LMSp, La) {
  nc = 0.57
  Lp <- LMSp[,1]^nc / (LMSp[,1]^nc + La^nc)
  Mp <- LMSp[,2]^nc / (LMSp[,2]^nc + La^nc)
  Sp <- LMSp[,3]^nc / (LMSp[,3]^nc + La^nc)
  LMSr <- matrix(c(Lp, Mp, Sp), ncol=3)
  return(LMSr)
}

## Calculate achromatic signal by averaging across cone responses in 40:20:1 ratio
calculateA <- function(LMSr) {
  A <- (40 * LMSr[,1] + 20 * LMSr[,2] + LMSr[,3]) / 61
  return (A)
}


## Calculate chromatic signals based on psychophysical results from Vos and Walraven
calculateab <- function(LMSr) {
  a <- 1/11 * (11 * LMSr[,1] - 12 * LMSr[,2] + LMSr[,3])
  b <- 1/9 * (LMSr[,1] + LMSr[,2] - 2 * LMSr[,3])
  ab  <- matrix(c(a, b), ncol=2)
  return (ab)
}

## Calculate perceived Lightness based on achromatic signals and physical stimulus
## Inverse sigmoidal function to convert As (achromtic signal) to J (perceived lightness)
## As, Aw achromatic signal for color sample As and for adopted white Aw
## E physical stimulus 1.0=LCD, 1.2175=transparent media, 1.4572=CRT, 1.7526=reflective paper
calcLightness <- function(As, Aw, E) {
  X <- As / Aw
  gX <- ((-(X-0.24)*0.65^3.65)/(X-0.24-0.89))^(1/3.65)
  gX[gX<0] <- 0
  J <- 100 * (E * (gX - 1) + 1)
  return (J)
}

# calculate Chroma based on chromatic signals
calcChroma <- function(a, b) {
  C <- 456.5 * (sqrt(a*a + b*b))^0.62
  return (C)
}

# calculate Hue Angle based on chromatic signals
calcHue <- function(a, b) {
  H <- 180 / pi * atan2(b, a)
  H <- H + 360 * (H < 0)
  return (H)
}

## Calculate Hue Quadrature (Moroney et al 2002) based on chromatic signals (to factor in ecentricy)
## Uses a linear interpolation based on the matrix Heh lookup table
## i is index in Heh[i,] matrixi which contains three values Hi, ei, hi
## See https://en.wikipedia.org/wiki/CIECAM02#Appearance_correlates
calcHueQ <- function(a, b) {
  h <- calcHue(a, b)
  h <- h + 360 * (h < 20.14)
  i = 1 + (h>=90) + (h>=164.25) + (h>=237.53) 
  ## H = Hi + 100 * ((h - hi) / ei)/((h - hi)/ei + (hi+1 - h)/ei+1)
  H = Heh[i,1] + 100 * ((h - Heh[i,3]) / Heh[i,2])/((h - Heh[i,3])/Heh[i,2] + (Heh[i+1,3] - h)/Heh[i+1,2])
  return (H)  
}

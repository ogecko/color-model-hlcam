## conversion matrices
## MCAT02 - MCAT02 transformation matrix (calculated using the CIE 1931 2° standard colorimetric observer)
## MH - Hunt–Pointer–Estévez space transformation matrix from XYZc to LMSp
## MHMCAT02 - transformation matrix from LMSc to LMSp
MCAT02 <- matrix(c(0.7328, 0.4296, -0.1624, -0.7036, 1.6975, 0.0061, 0.0030, 0.0136, 0.9834), ncol=3, byrow=TRUE)
MH <- matrix(c(0.38971, 0.68898, -0.07868, -0.22981, 1.18340, 0.04641, 0.0, 0.0, 1), ncol=3, byrow=TRUE)
MCAT02INV <- solve(MCAT02)
MCAT02INVY <- MCAT02INV[2,]
MHMCAT02 <- MH %*% MCAT02INV
MHMCAT02INV <- solve(MHMCAT02)
MSIG02 <- matrix(c(40.0/61.0, 20.0/61.0, 1.0/61.0, 1.0, -12.0/11.0, 1.0/11.0, 1.0/9.0, 1.0/9.0, -2.0/9.0), ncol=3, byrow=TRUE)
MSIG02INV <- solve(MSIG02)
Heh <- matrix(c(0.0, 100.0, 200.0, 300.0, 400.0, 0.8, 0.7, 1.0, 1.2, 0.8, 20.14, 90.0, 164.25, 237.53, 380.14), ncol=3)

## convert XYZ to CCT using McCamy's Approximation https://en.wikipedia.org/wiki/Color_temperature#Approximation

# forward convert CCT to XYZ: CCT -> XYZ
CCTtoXYZ <- function(CCT) {

}

# reverse convert XYZ to CCT: XYZ -> CCT
XYZtoCCT <- function(Xw, Yw, Zw) {
  xw = Xw / (Xw + Yw + Zw)
  yw = Yw / (Xw + Yw + Zw)
  n = (xw - 0.3320) / (yw - 0.1858)
  CCT = -449*n^3 + 3525*n^2 - 6823.3*n + 5520.33
  return (CCT)
}


## convert between XYZ and LMS spaces using the MCAT02 transform
## returns an (n x 3) matrix, or just the Y component

## forward convert XYZ to LMS space using MCAT02 transform: XYZ -> LMS
XYZtoLMS <- function(X,Y,Z) {
  XYZ <- matrix(c(X,Y,Z), ncol=3)
  LMS <- t(MCAT02 %*% t(XYZ))
  return (LMS)
}

## reverse convert LMS to XYZ space using inverse MCAT02 transform: LMS -> XYZ
LMStoXYZ <- function(LMS) {
  XYZ <- t(MCAT02INV %*% t(LMS))
  return (XYZ)
}

## reverse convert LMS to XYZ space using inverse MCAT02 transform: LMS -> Y
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

## use CAT02 to color adapt a sample based on the illumant adopted white and the degree of adaption
## D = degree of adaptation (0 for no adaptation, 1 for full adaption)
## LMSs = LMS of the color sample to adapt
## LMSw = LMS of the white point (ie illuminant adopted white)
## LMSc = LMS of the adapted color sample (based on a reference white used in CIECAM02)

# forward raw LMS to color adapted LMS: LMS, D, LMSw -> LMSc
LMStoLMSc <- function(D, LMSs, LMSw) {
  Yw <- LMStoY(LMSw)
  Lc <- LMSs[,1] * (D * Yw / LMSw[,1] + 1 - D)
  Mc <- LMSs[,2] * (D * Yw / LMSw[,2] + 1 - D)
  Sc <- LMSs[,3] * (D * Yw / LMSw[,3] + 1 - D)
  LMSc <- matrix(c(Lc, Mc, Sc), ncol=3)
  return (LMSc)  
}

# reverse color adapted LMS to raw LMS: LMSc, D, LMSw -> LMS
LMSctoLMS <- function(LMSc, D, LMSw) {
  Yw <- LMStoY(LMSw)
  L <- LMSc[,1] / (D * Yw / LMSw[,1] + 1 - D)  
  M <- LMSc[,2] / (D * Yw / LMSw[,2] + 1 - D)  
  S <- LMSc[,3] / (D * Yw / LMSw[,3] + 1 - D)  
  LMS <- matrix(c(L, M, S), ncol=3)
  return (LMS)  
}

## perform post adaption conversion to cone responses using Hunt-Pointer-Esteves
## LMSc = LMS of the adapted color sample
## LMSp = Post Adaption Hunt-Point_estevez space

# forward adapted LMS to Hunt-Point_estevez space: LMSc -> LMSp
LMSctoLMSp <- function(LMSc) {
  LMSp = t(MHMCAT02 %*% t(LMSc))
  return (LMSp)
}

# reverse Hunt-Point_estevez space to LMS adapted: LMSp -> LMSc
LMSptoLMSc <- function(LMSp) {
  LMSc = t(MHMCAT02INV %*% t(LMSp))
  return (LMSc)
}

## convert to absolute cone responses based on physiological experiments by Valenton and Norren
## they showed cone response of primates follows a sigmoidal curve influenced by level of adaption
## LMSp = Post Adaption Hunt-Point_estevez space
## LMSr = Absolute cone response

# forward convert Post Adaption space to absolute cone response: LMSp, La -> LMSr
LMSptoLMSr <- function(LMSp, La) {
  nc = 0.57
  Lr <- LMSp[,1]^nc / (LMSp[,1]^nc + La^nc)
  Mr <- LMSp[,2]^nc / (LMSp[,2]^nc + La^nc)
  Sr <- LMSp[,3]^nc / (LMSp[,3]^nc + La^nc)
  LMSr <- matrix(c(Lr, Mr, Sr), ncol=3)
  return(LMSr)
}

# reverse convert absolute cone response to Post Adaption space: LMSr, La -> LMSp
LMSrtoLMSp <- function(LMSr, La) {
  nc = 0.57
  Lp <- (-La^nc * LMSr[,1] / (LMSr[,1] - 1))^(1/nc)
  Mp <- (-La^nc * LMSr[,2] / (LMSr[,2] - 1))^(1/nc)
  Sp <- (-La^nc * LMSr[,3] / (LMSr[,3] - 1))^(1/nc)
  LMSp <- matrix(c(Lp, Mp, Sp), ncol=3)
  return(LMSp)
}



## Calculate achromatic signal by averaging across cone responses in 40:20:1 ratio
## Calculate chromatic signals based on psychophysical results from Vos and Walraven
## Linear transform is stored in MSIG02

# forward absolute cone response to achromatic signals: LMSr -> A 
LMSrtoA <- function(LMSr) {
  A <- t(MSIG02[1,] %*% t(LMSr))
  return (A)
}

# forward absolute cone response to chromatic signals: LMSr -> ab 
LMSrtoAB <- function(LMSr) {
  ab <- t(MSIG02[2-3,] %*% t(LMSr))
  return (ab)
}

# forward absolute cone response to both chromatic/achromatic signals: LMSr -> Aab 
LMSrtoAab <- function(LMSr) {
  Aab <- t(MSIG02 %*% t(LMSr))
  return (Aab)
}

# reverse chromatic/achromatic signals to absolute cone response: Aab -> LMSr 
AabtoLMSr <- function(As, ab) {
  Aab <- matrix(c(As, ab[,1], ab[,2]), ncol=3)
  LMSr <- t(MSIG02INV %*% t(Aab))
  return (LMSr)
}

## Calculate perceived Lightness based on achromatic signals and physical stimulus
## Inverse sigmoidal function to convert As (achromtic signal) to J (perceived lightness)
## As, Aw achromatic signal for color sample As and for adopted white Aw
## E physical stimulus 1.0=LCD, 1.2175=transparent media, 1.4572=CRT, 1.7526=reflective paper

# forward achromatic signal to lightness: As, Aw, E -> J
AtoLightness <- function(As, Aw, E) {
  X <- As / Aw
  gX <- ((-(X-0.24)*0.65^3.65)/(X-0.24-0.89))^(1/3.65)
  gX[gX<0] <- 0
  J <- 100 * (E * (gX - 1) + 1)
  return (J)
}

# reverse lightness to achromatic signal: J, Aw, E -> As
LightnesstoA <- function(J, Aw, E) {
  gX <- (J/100 - 1)/E + 1
  As <- Aw * (0.89 * gX^3.65 / (gX^3.65 + 0.65^3.65) + 0.24)
  return (As)
}

## Calculate perceived Chroma and Hue based on chromatic signals ab

# forward chromatic signals to Chroma: ab -> C
ABtoChroma <- function(ab) {
  C <- 456.5 * (sqrt(ab[,1]*ab[,1] + ab[,2]*ab[,2]))^0.62
  return (C)
}

# forward chromatic signals to Hue Quadrature: ab -> H
ABtoHueQ <- function(ab) {
  h <- 180 / pi * atan2(ab[,2], ab[,1])
  h <- h + 360 * (h < 0)
  H <- htoHueQ(h)
  return (H)
}

# reverse Hue Quadrature and Chroma to chromatic signals: H,C -> ab
HueQChromatoAB <- function(H, C) {
  h <- HueQtoh(H)
  a <- cos(pi*h/180)*(C/456.5)^(1/0.62) 
  b <- sin(pi*h/180)*(C/456.5)^(1/0.62) 
  ab  <- matrix(c(a, b), ncol=2)
  return(ab)  
}


## Calculate Hue Quadrature (Moroney et al 2002) based on Hue Angle (to factor in ecentricy)
## Uses a linear interpolation based on the matrix Heh lookup table
## i is index in Heh[i,] matrixi which contains three values Hi, ei, hi
## See https://en.wikipedia.org/wiki/CIECAM02#Appearance_correlates

## forward Hue Angle to Hue Quadrature: h -> H
htoHueQ <- function(h) {
  h <- h + 360 * (h < 20.14)
  i = 1 + (h>=90) + (h>=164.25) + (h>=237.53) 
  ## H = Hi + 100 * ((h - hi) / ei)/((h - hi)/ei + (hi+1 - h)/ei+1)
  H <- Heh[i,1] + 100 * ((h - Heh[i,3]) / Heh[i,2])/((h - Heh[i,3])/Heh[i,2] + (Heh[i+1,3] - h)/Heh[i+1,2])
  return (H)  
}

## reverse Hue Quadrature to Hue Angle: H -> h
HueQtoh <- function(H) {
  i <- 1 + (H>=100.0) + (H>=200.0) + (H>=300.0) 
  h <- ((H - Heh[i,1]) * (Heh[i+1,2]*Heh[i,3] - Heh[i,2]*Heh[i+1,3]) - 100 * Heh[i,3] * Heh[i+1,2]) / ((H - Heh[i,1])*(Heh[i+1,2] - Heh[i,2]) - 100 * Heh[i+1,2])
  h <- h - 360 * (h > 360)
  return (h)
}

# forward convert XYZ to JCH: XYZs, XYZw, La, A, E -> JHC
XYZtoJHC <- function(Xs, Ys, Zs, Xw, Yw, Zw, La, Ambient, E) {
  ## convert (color sample) and (adopted white) to LMS space
  LMSs <- XYZtoLMS(Xs,Ys,Zs)
  LMSw <- XYZtoLMS(Xw,Yw,Zw)
  
  ## use CAT02 to adapt the sample from the adopted white illumant to a reference white used in CIECAM02
  D <- calculateD(Ambient, La)
  LMSsc <- LMStoLMSc(D, LMSs, LMSw)
  LMSwc <- LMStoLMSc(D, LMSw, LMSw)
  
  ## perform post adaption conversion to cone responses using Hunt-Pointer-Esteves
  LMSsp <- LMSctoLMSp(LMSsc)
  LMSwp <- LMSctoLMSp(LMSwc)
  
  ## convert to absolute cone responses
  LMSsr <- LMSptoLMSr(LMSsp, La)
  LMSwr <- LMSptoLMSr(LMSwp, La)
  
  ## Calculate achromatic and chromatic signals 
  As <- LMSrtoA(LMSsr)
  Aw <- LMSrtoA(LMSwr)
  ab <- LMSrtoAB(LMSsr)
  
  ## calculate the perceived lightness (0 to 100), Hue Quadrature (0 to 400) and Chroma (0 to 100)
  J <- AtoLightness(As, Aw, E)
  H <- ABtoHueQ(ab)
  C <- ABtoChroma(ab)

  JHC <- matrix(c(J, H, C), ncol=3)
  return(JHC)
}

# reverse convert JCH to XYZ: JHC, XYZw, La, A, E -> XYZs
JHCtoXYZ <- function(Js, Hs, Cs, Xw, Yw, Zw, La, Ambient, E) {
  ## calculate the achromatic white point for the target device
  RLMSw <- XYZtoLMS(all$Xw,all$Yw,all$Zw)
  RD <- calculateD(all$Ambient, all$La)
  RLMSwc <- LMStoLMSc(RD, RLMSw, RLMSw)
  RLMSwp <- LMSctoLMSp(RLMSwc)
  RLMSwr <- LMSptoLMSr(RLMSwp, all$La)
  RAw <- LMSrtoA(RLMSwr)
  
  # calculate the acrhomatic signal A from the sample lightness J
  RAs <- LightnesstoA(all$Js1, RAw, 1.0)
  
  # calculate the chromatic ab signal from the Hue Quadrant and Chroma 
  Rab <- HueQChromatoAB(all$Hs1, all$Cs1)
  
  # compute cone response LMS from the achromatic signal A and opponents a & b
  RLMSsr <- AabtoLMSr(RAs, Rab)
  
  # compute cone signal LMS from the LMS cone response and La
  RLMSsp <- LMSrtoLMSp(RLMSsr, all$La)
  
  # compute the color adapted LMS
  RLMSsc <- LMSptoLMSc(RLMSsp)
  
  # compute the raw LMS
  RLMS <- LMSctoLMS(RLMSsc, RD, RLMSw)
  
  # compute the XYZ
  XYZ <- LMStoXYZ(RLMS)

  return(XYZ)
}  
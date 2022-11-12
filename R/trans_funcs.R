#' @title Transfer models from underwater remote-sensing ratio (rrs) to
#'   above-water remote-sensing reflectance (Rrs)
#' @name trans_function
#' @param rrs remote sensing reflectance just beneath the water surface
#' @param wavelen wavelength in nm
#' @param Theta Viewing zenith angle in deg
#' @param Phi Relative azimuth angle in deg
#' @param suntheta Solar zenith angle in deg
#' @param windspd Wind speed in meter/s default as 5 m/s
#' @param cloud Cloud fraction default as 0 (a clear sky) and 1 is a overcast
#'   sky
#' @param at total absorption in m-1
#' @param bbt total backward scattering in m-1
#' @param bbw water backward scattering in m-1
#' @param nw Refractive index of water
#' @export
#'
#' @references
#'
#' - Lee, ZhongPing, Kendall L. Carder, and Robert A. Arnone. “Deriving Inherent
#'   Optical Properties from Water Color: A Multiband Quasi-Analytical Algorithm
#'   for Optically Deep Waters.” Applied Optics 41, no. 27 (September 20, 2002):
#'   5755. https://doi.org/10.1364/AO.41.005755.
#'
#' - Lee, Zhongping, Keping Du, Kenneth J. Voss, Giuseppe Zibordi, Bertrand
#'   Lubac, Robert Arnone, and Alan Weidemann. “An
#'   Inherent-Optical-Property-Centered Approach to Correct the Angular Effects in
#'   Water-Leaving Radiance.” Applied Optics 50, no. 19 (July 1, 2011): 3155.
#'   https://doi.org/10.1364/AO.50.003155.
#'
#' - Loisel, Hubert. “Rrs(0+) -> Rrs(0-) & Water Coefficients.” Ocean Color
#'   Forum, Inherent OPtical Properties Workshop, 2008.
#'   http://oceancolor.gsfc.nasa.gov/forum/oceancolor/topic_show.pl?tid=2657.
#'
#' - Mobley, Curtis D. “Estimation of the Remote-Sensing Reflectance from
#'   above-Surface Measurements.” Applied Optics 38, no. 36 (December 20, 1999):
#'   7442. https://doi.org/10.1364/AO.38.007442.
#'
#' - Morel, André, and Bernard Gentili. “Diffuse Reflectance of Oceanic Waters
#'   III Implication of Bidirectionality for the Remote-Sensing Problem.” Applied
#'   Optics 35, no. 24 (August 20, 1996): 4850. https://doi.org/10.1364/AO.35.004850.
#'
#' @examples
#'
#' # parameters preparation
#' wavelen  = c(410, 440, 560, 620, 670, 800)
#' rrs      = c(0.00557, 0.00834, 0.0208, 0.01035, 0.00624, 0.00132)
#' Rrs      = c(0.00297, 0.0045, 0.01158, 0.00566, 0.00339, 0.00072)
#' at       = c(0.7576, 0.5058, 0.1964, 0.3547, 0.5321, 2.2494)
#' bbt      = c(0.04685, 0.04536, 0.04089, 0.03851, 0.0364, 0.03578)
#' bbw      = c(0.00306, 0.00226, 0.00082, 0.00054, 0.00039, 0.00019)
#' suntheta = 59.90437
#' Theta    = 40
#' Phi      = 135
#' windspd  = 6
#' cloud    = 0
#' Temp     = 1.1
#' Sal      = 31.1
#'
#' # run models
#' nw = WOPP(Temp, Sal, wavelen)$nw
#' Rrs_Lee02 = trans_Lee02(rrs)
#' Rrs_Loisel08 = trans_Loisel08(rrs, suntheta)
#' Rrs_Lee11 = trans_Lee11(at, bbt, bbw, Theta, Phi, suntheta)
#' Rrs_Bi22 = trans_Bi22(rrs, Theta, Phi, suntheta, windspd, cloud, nw = nw)
#' Rrs_Mobley99 = trans_Mobley99(rrs)
#' Rrs_Morel96 = trans_Morel96(rrs, wavelen, Theta, suntheta)
#'
#' # plot
#' pal_col = rainbow(n=6)
#' plot(wavelen, Rrs, col = "black", type = "l",
#'      xlim = c(400, 810), ylim = c(0, 0.012),
#'      xlab = "Wavelength [nm]", ylab = "Rrs [1/sr]")
#' lines(wavelen, Rrs_Bi22, col = pal_col[1], type = "l")
#' lines(wavelen, Rrs_Lee02, col = pal_col[2], type = "l")
#' lines(wavelen, Rrs_Loisel08, col = pal_col[3], type = "l")
#' lines(wavelen, Rrs_Lee11, col = pal_col[4], type = "l")
#' lines(wavelen, Rrs_Morel96, col = pal_col[5], type = "l")
#' lines(wavelen, Rrs_Mobley99, col = pal_col[6], type = "l")
#' lines(wavelen, Rrs, col = "black", type = "l", lty = 2)
#' legend(x = 810, y = 0.012, xjust = 1, yjust = 1,
#'        c("HydroLight", "Bi22", "Lee02", "Loisel08",
#'          "Lee11", "Morel96", "Mobley99"),
#'        col = c("black", pal_col), lty = c(1, 2, rep(1, 5)),
#'        x.intersp = 0.5, y.intersp = 0.8)
#'

trans_Lee02 <- function(rrs) {

  a = 0.52
  b = 1.7
  trans_abfun_rrs(rrs, a, b)

}

#' @rdname trans_function
#' @export
trans_Mobley99 <- function(rrs) {

  a = 0.54
  b = 0
  trans_abfun_rrs(rrs, a, b)

}

#' @rdname trans_function
#' @export
trans_Loisel08 <- function(rrs, suntheta) {

  if(suntheta < 0 | suntheta > 90) {
    stop("suntheta should be between 0 and 90")
  }
  if(suntheta > 60) suntheta = 60

  s_seq = c(0, 30, 60)
  a_seq = c(0.5236, 0.5169, 0.4933)
  b_seq = c(2.1941, 2.3001, 2.6796)
  names(a_seq) <- names(b_seq) <- s_seq

  bd_s <- get_boundary(suntheta, s_seq)
  mat_A <- matrix(c(1, 1, bd_s[1], bd_s[2]), ncol = 2, nrow = 2)
  coef_a <- solve(mat_A, c(a_seq[as.character(bd_s)]))
  coef_b <- solve(mat_A, c(b_seq[as.character(bd_s)]))
  a <- as.numeric(t(c(1, suntheta)) %*% coef_a)
  b <- as.numeric(t(c(1, suntheta)) %*% coef_b)

  trans_abfun_rrs(rrs, a, b)

}

#' @rdname trans_function
#' @export
trans_Lee11 <- function(at, bbt, bbw, Theta, Phi, suntheta) {

  g <- lut_interp_Lee11(Theta, Phi, suntheta)
  gw0 = g["gw0"]
  gw1 = g["gw1"]
  gp0 = g["gp0"]
  gp1 = g["gp1"]
  kai = at + bbt
  bbp = bbt - bbw
  as.numeric((gw0 + gw1 * bbw/kai) * (bbw/kai) + (gp0 + gp1 * bbp/kai) * (bbp/kai))

}

#' @rdname trans_function
#' @export
trans_Bi22 <- function(rrs, Theta, Phi, suntheta, windspd, cloud = NULL, nw) {

  ab <- lut_interp_Bi22(Theta, Phi, suntheta, windspd, cloud)

  a = as.numeric(ab["ap"] / nw^2)
  b = as.numeric(ab["b"])

  trans_abfun_rrs(rrs, a, b)

}

#' @rdname trans_function
#' @export
trans_Morel96 <- function(rrs, wavelen, Theta, suntheta) {

  if(length(rrs) != length(wavelen)) {
    stop("Error - length of rrs and wavelen should match!")
  }

  rho_bar = 0.043
  nw <- Quan_Fry(wavelen)

  xa = suntheta / 180 * pi
  na = 1
  xw = asin(sin(xa)/nw)
  Rs = ((na*cos(xw)-nw*cos(xa))/(na*cos(xw)+nw*cos(xa)))^2
  Rp = ((na*cos(xa)-nw*cos(xw))/(na*cos(xa)+nw*cos(xw)))^2
  rho_wa = (Rs + Rp)/2

  a0 <- with(LUT_Morel96, approx(wv, a0, xout = wavelen))$y
  a1 <- with(LUT_Morel96, approx(wv, a1, xout = wavelen))$y
  a2 <- with(LUT_Morel96, approx(wv, a2, xout = wavelen))$y
  a3 <- with(LUT_Morel96, approx(wv, a3, xout = wavelen))$y

  Qvalue <- a0 + a1 * suntheta + a2 * suntheta^2 + a3 * suntheta^3

  (1 - rho_bar) * (1 - rho_wa) / ((1-rho_bar * Qvalue * rrs) * nw^2) * rrs

}


trans_abfun_rrs <- function(x, a, b) {
  a * x / (1 - b * x)
}

trans_abfun_big_Rrs <- function(x, a, b) {
  x / (a + b * x)
}






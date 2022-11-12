lut_interp_Lee11 <- function(Theta, Phi, suntheta) {

  d <- LUT_Lee11

  if(Theta < 0 | Theta > 90) stop("Theta should be between 0 and 90")
  if(Theta > max(d$Theta)) Theta = max(d$Theta)
  if(Theta < 10) {
    Theta = 0
    Phi = 0
  }

  if(Phi < 0 | Phi > 180) stop("Phi should be between than 0 and 180")
  if(Phi > max(d$Phi)) Phi = max(d$Phi)

  if(suntheta < 0 | suntheta > 90) stop("suntheta should be between 0 and 90")
  if(suntheta > max(d$suntheta)) suntheta = max(d$suntheta)

  bd_T <- get_boundary(Theta, sort(unique(d$Theta)))
  bd_P <- get_boundary(Phi, sort(unique(d$Phi)))
  bd_s <- get_boundary(suntheta, sort(unique(d$suntheta)))

  mat_parm <- expand.grid(Theta = bd_T, Phi = bd_P, suntheta = bd_s)
  mat <- data.table::merge.data.table(mat_parm, d, by = c("Theta", "Phi", "suntheta"),
                                      all.x = TRUE, all.y = FALSE)

  if(0 %in% (mat$Theta + mat$Phi)) {
    mat[mat$Theta == 0 & mat$Phi != 0, -c(1:3)] <- mat[mat$Theta == 0 & mat$Phi == 0, -c(1:3)]
  } else {
    mat[mat$Theta == 0 & mat$Phi != 0, -c(1:3)] <- mat[!(mat$Theta == 0 & mat$Phi != 0), -c(1:3)]
  }

  mat <- data.table::as.data.table(mat)

  g_values <- numeric(4) * NA
  names(g_values) <- c("gp0", "gp1", "gw0", "gw1")

  for(gvar in c("gp0", "gp1", "gw0", "gw1")) {

    mat_g <- array(dim = rep(2, 3), dimnames = list(bd_T, bd_P, bd_s))
    for(iT in as.character(bd_T)) {
      for(iP in as.character(bd_P)) {
        for(is in as.character(bd_s)) {
          mat_g[iT, iP, is] <-
            mat[Theta==iT & Phi==iP & suntheta==is, gvar, with = FALSE][[1]]
        }
      }
    }

    g_values[gvar] <- approx3d(mat_g, bd_T, bd_P, bd_s, Theta, Phi, suntheta)

  }

  g_values

}



lut_interp_Bi22 <- function(Theta, Phi, suntheta, windspd, cloud) {

  d <- LUT_Bi22

  if(Theta < 0 | Theta > 90) stop("Theta should be between 0 and 90")
  if(Theta > max(d$Theta)) Theta = max(d$Theta)

  if(Phi < 0 | Phi > 180) stop("Phi should be between than 0 and 180")
  if(Phi > max(d$Phi)) Phi = max(d$Phi)

  if(suntheta < 0 | suntheta > 90) stop("suntheta should be between 0 and 90")
  if(suntheta > max(d$suntheta)) suntheta = max(d$suntheta)

  if(windspd < 0) windspd = 0
  if(windspd > 10) windspd = 10

  if(is.null(cloud)) {
    cloud = 0
  } else {
    if(cloud < 0) cloud = 0
    if(cloud > 1) cloud = 1
  }

  bd_T <- get_boundary(Theta, sort(unique(d$Theta)))
  bd_P <- get_boundary(Phi, sort(unique(d$Phi)))
  bd_s <- get_boundary(suntheta, sort(unique(d$suntheta)))
  bd_w <- get_boundary(windspd, sort(unique(d$windspd)))
  bd_c <- get_boundary(cloud, sort(unique(d$cloud)))

  mat_parm <- expand.grid(Theta = bd_T, Phi = bd_P, suntheta = bd_s, windspd = bd_w, cloud = bd_c)
  mat <- data.table::merge.data.table(mat_parm, d, by = c("Theta", "Phi", "suntheta", "windspd", "cloud"),
                                      all.x = TRUE, all.y = FALSE)

  if(0 %in% (mat$Theta + mat$Phi)) {
    mat[mat$Theta == 0 & mat$Phi != 0, -c(1:3)] <- mat[mat$Theta == 0 & mat$Phi == 0, -c(1:3)]
  } else {
    mat[mat$Theta == 0 & mat$Phi != 0, -c(1:3)] <- mat[!(mat$Theta == 0 & mat$Phi != 0), -c(1:3)]
  }

  mat <- data.table::as.data.table(mat)

  ab_values <- numeric(2) * NA
  names(ab_values) <- c("ap", "b")

  for(ab_var in c("ap", "b")) {

    mat_ab <- array(dim = rep(2, 5), dimnames = list(bd_T, bd_P, bd_s, bd_w, bd_c))
    for(iT in as.character(bd_T)) {
      for(iP in as.character(bd_P)) {
        for(is in as.character(bd_s)) {
          for(iw in as.character(bd_w)) {
            for(ic in as.character(bd_c)) {
              mat_ab[iT, iP, is, iw, ic] <-
                mat[Theta==iT & Phi==iP & suntheta==is & windspd == iw & cloud == ic, ab_var, with = FALSE][[1]]
            }
          }
        }
      }
    }

    ab_values[ab_var] <- approx5d(mat_ab, bd_T, bd_P, bd_s, bd_w, bd_c,
                                  Theta, Phi, suntheta, windspd, cloud)

  }

  ab_values

}





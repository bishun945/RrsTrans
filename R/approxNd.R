#' get_boundary
#'
#' @param x Points to find boundary
#' @param vec Unique sorted index values
#'
#' @return Numeric two values. Left and right boundary of \code{x} in \code{vec}
#' @noRd
#' @author Shun Bi
get_boundary <- function(x, vec) {

  sort(vec[sort.int(abs(as.numeric(vec) - x), index.return = TRUE)$ix[1:2]])

}



approx2d <- function(mat, x_seq, y_seq, x, y) {

  Ndim = 2
  Q = expand.grid(1:2, 1:2)
  colnames(Q) = c("x", "y")

  mat_C = numeric(2^Ndim)
  for(i in 1:length(mat_C)) {
    mat_C[i] = mat[Q[i,1],Q[i,2]]
  }

  # combn(c("x", "y"), 2)

  mat_A = diag(2^Ndim) * NA
  for(i in 1:nrow(mat_A)) {
    xx = x_seq[Q[i,"x"]]
    yy = y_seq[Q[i,"y"]]
    mat_A[i,] = c(1, xx, yy, xx*yy)
  }

  input <- c(1, x, y, x*y)

  a <- solve(mat_A, mat_C)
  r  <- as.numeric(t(input) %*% a)

  return(r)

}



approx3d <- function(mat, x_seq, y_seq, z_seq, x, y, z) {

  Ndim = 3
  Q = expand.grid(1:2, 1:2, 1:2)
  colnames(Q) = c("x", "y", "z")

  mat_C = numeric(2^Ndim)
  for(i in 1:length(mat_C)) {
    mat_C[i] = mat[Q[i,1],Q[i,2],Q[i,3]]
  }

  # combn(c("x", "y", "z"), 2)

  mat_A = diag(2^Ndim) * NA
  for(i in 1:nrow(mat_A)) {
    xx = x_seq[Q[i,"x"]]
    yy = y_seq[Q[i,"y"]]
    zz = z_seq[Q[i,"z"]]
    mat_A[i,] = c(1, xx, yy, zz, xx*yy, xx*zz, yy*zz, xx*yy*zz)
  }

  input <- c(1, x, y, z, x*y, x*z, y*z, x*y*z)

  a <- solve(mat_A, mat_C)
  r  <- as.numeric(t(input) %*% a)

  return(r)

}


#' approx5d
#'
#' @param mat mat
#' @param x_seq x_seq
#' @param y_seq y_seq
#' @param z_seq z_seq
#' @param u_seq u_seq
#' @param v_seq v_seq
#' @param x x
#' @param y y
#' @param z z
#' @param u u
#' @param v v
#'
#' @note
#'
#' this one is more general than 2d... others should be revised...
#' f(x,y,...,n) = ao + a1x + a2y + ... + a[]xy + ... + a[]xyz + ... + # a\[n-1\]xy..n
#' 2d approx needs 4 points surrounded -> a0, a1, ..., a3
#' 3d needs 8 points                   -> a0, a1, ..., a7
#' 5d need 2^5=32 points               -> a0, a1, ..., a31
#'
#' @return Interpolated value
#' @author Shun Bi Nov-12-22
#' @noRd
#'
#' @importFrom stats approx
#'
approx5d <- function(mat, x_seq, y_seq, z_seq, u_seq, v_seq, x, y, z, u, v) {

  Q = expand.grid(1:2, 1:2, 1:2, 1:2, 1:2)
  colnames(Q) = c("x", "y", "z", "u", "v")

  mat_C = numeric(2^5)
  for(i in 1:length(mat_C)) {
    mat_C[i] = mat[Q[i,1],Q[i,2],Q[i,3],Q[i,4],Q[i,5]]
  }

  mat_A = diag(2^5) * NA
  for(i in 1:nrow(mat_A)) {
    xx = x_seq[Q[i,"x"]]
    yy = y_seq[Q[i,"y"]]
    zz = z_seq[Q[i,"z"]]
    uu = u_seq[Q[i,"u"]]
    vv = v_seq[Q[i,"v"]]
    mat_A[i,] = c(
      1, xx, yy, zz, uu, vv,                                                    # a5
      xx*yy, xx*zz, xx*uu, xx*vv, yy*zz, yy*uu, yy*vv, zz*uu, zz*vv, uu*vv,     # a15
      xx*yy*zz, xx*yy*uu, xx*yy*vv, xx*zz*uu, xx*zz*vv, xx*uu*vv, yy*zz*uu,     # a22
      yy*zz*vv, yy*uu*vv, zz*uu*vv,                                             # a25
      xx*yy*zz*uu, xx*yy*zz*vv, xx*yy*uu*vv, xx*zz*uu*vv, yy*zz*uu*vv,          # a30
      xx*yy*zz*uu*vv                                                            # a31
    )
  }

  input <- c(
    1, x, y, z, u, v,                                                           # a5
    x*y, x*z, x*u, x*v, y*z, y*u, y*v, z*u, z*v, u*v,                           # a15
    x*y*z, x*y*u, x*y*v, x*z*u, x*z*v, x*u*v, y*z*u,                            # a22
    y*z*v, y*u*v, z*u*v,                                                        # a25
    x*y*z*u, x*y*z*v, x*y*u*v, x*z*u*v, y*z*u*v,                                # a30
    x*y*z*u*v                                                                   # a31
  )

  a <- solve(mat_A, mat_C)
  r  <- as.numeric(t(input) %*% a)

  return(r)

}



#' approxNd
#'
#' @param mat Array on N dimensions
#' @param list_seq Boundary seq on N dimensions
#' @param list_q Query values on N dimensions
#'
#' @return A numeric value
#' @export
#'
#' @importFrom utils combn
#'
#' @examples
#' list_seq = list(1:2, 3:4, 5:6)
#' list_q   = list(1.5, 3.5, 5.5)
#' mat = array(runif(8),
#'             dim = lapply(list_seq, length),
#'             dimnames = list_seq)
#' approxNd(mat, list_seq, list_q)
#' # the interpolated value should be identical to the mean values
#' mean(mat)
#'
#' list_seq = list(1:2, 3:4, 5:6, 7:8, 9:10, 11:12)
#' list_q   = list(1.5, 3.5, 5.5, 7.5, 9.5, 11.5)
#' mat = array(runif(2^length(list_seq)),
#'             dim = lapply(list_seq, length),
#'             dimnames = list_seq)
#' approxNd(mat, list_seq, list_q)
#' mean(mat)
#'
approxNd <- function(mat, list_seq, list_q) {

  N <- length(dim(mat))
  VN <- LETTERS[1:N]
  eval(parse(text = sprintf("Q = expand.grid(%s)", paste(rep("1:2", N), collapse = ", "))))
  names(Q) <- VN

  mat_C = numeric(2^N)
  for(i in 1:length(mat_C)) {
    eval(parse(text = sprintf("mat_C[i]=mat[%s]", paste(sprintf("Q[i,%s]", 1:N), collapse = ","))))
  }

  multi_vec <- lapply(1:N, function(iN) {
    apply(combn(paste0(VN, VN), iN), 2, paste, collapse = "*")
  })
  multi_vec <- unlist(multi_vec)
  multi_str <- paste(c(1, multi_vec), collapse = ", ")

  mat_A = diag(2^N) * NA
  for(i in 1:nrow(mat_A)) {

    for(j in 1:N) {
      eval(parse(text = sprintf("%s%s=list_seq[[j]][Q[i,'%s']]",VN[j],VN[j],VN[j])))
    }

    eval(parse(text = sprintf("mat_A[i,]=c(%s)", multi_str)))

  }

  # generate input vec
  multi_vec <- lapply(1:N, function(iN) {
    apply(combn(paste0(VN), iN), 2, paste, collapse = "*")
  })
  multi_vec <- unlist(multi_vec)
  multi_str <- paste(c(1, multi_vec), collapse = ", ")

  for(j in 1:N) {
    eval(parse(text = sprintf("%s=list_q[[j]]",VN[j])))
  }

  input <- NULL
  eval(parse(text = sprintf("input=c(%s)", multi_str)))

  mat_B <- solve(mat_A, mat_C)

  r <- as.numeric(t(input) %*% mat_B)

  return(r)

}






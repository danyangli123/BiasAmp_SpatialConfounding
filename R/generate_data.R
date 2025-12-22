#covariance-function
mat_corfun_matern <- function(phi, D) {
  Sigma <- (1+sqrt(3)*D/phi)*exp(-sqrt(3)*D/phi) #p=1
  Sigma
}
mat_corfun_exp <- function(phi, D) {
  Sigma <- exp(-D/phi)
  Sigma
}

gen_data_comps <- function(size, n_obs = 1000, phi_seq, sigsq_y_true_seq, rand_points = TRUE, mat_corfun = mat_corfun_exp, scale_proc = T) {

  if (rand_points) {
    #n_obs=200
    locs <- matrix(runif(2*n_obs, 0, sqrt(size)),ncol=2)
  } else {
    ## regular grid of points
    npt <- ceiling(sqrt(n_obs))
    sq <- (1:npt)/(npt + 1)
    locs <- crossing(p0 = sq, p1 = sq)
  }
  df <- locs %>%
    set_colnames(c("p0", "p1")) %>%
    as_tibble()

  D <- as.matrix(fields::rdist(as.matrix(locs)))
  proc_c <- proc_u <- matrix(NA, n_obs, length(phi_seq))
  for (j in seq_along(phi_seq)) {
    Sigma <- mat_corfun(phi = phi_seq[j], D = D)
    if (n_obs <= 10000) {
      proc <- t(mvtnorm::rmvnorm(2,sigma=Sigma,method="svd"))
    } else {
      print("larger")
      cholS <- chol(Sigma)
      z <- matrix(rnorm(2*n_obs),ncol = 2)
      proc <- t(cholS)%*%z
    }
    if (scale_proc) proc %<>% scale()
    proc_c[, j] <- proc[, 1]
    proc_u[, j] <- proc[, 2]
  }
  colnames(proc_c) <- paste0("proc_c_", c(0.04,0.15,0.6))
  colnames(proc_u) <- paste0("proc_u_", c(0.04,0.15,0.6))
  proc_c %<>% as_tibble()
  proc_u %<>% as_tibble()

  #non-spatial component (contructed from a std. normal  variable) in exposure and error term in outcome
  sigsq_seq <- c(1, sigsq_y_true_seq)
  mat <- matrix(NA, nrow = n_obs, ncol = length(sigsq_seq))
  colnames(mat) <- c("proc_ns_std",paste0("err_y_", sigsq_y_true_seq))
  for (j in seq_along(sigsq_seq)) {
    mat[, j] <- rnorm(n = n_obs, 0, sqrt(sigsq_seq[j]))
  }
  mat %<>% as_tibble()

  df %<>%
    bind_cols(proc_c) %>%
    bind_cols(proc_u) %>%
    bind_cols(mat)

  df_comps <- df
  df_comps
}

gen_data_comps_phiY <- function(size, n_obs = 1000, phi_u, phi_c, phi_y, sigsq_y_true_seq, rand_points = TRUE, mat_corfun = mat_corfun_exp, scale_proc = T) {

  if (rand_points) {
    #n_obs=200
    locs <- matrix(runif(2*n_obs, 0, sqrt(size)),ncol=2)
  } else {
    ## regular grid of points
    npt <- ceiling(sqrt(n_obs))
    sq <- (1:npt)/(npt + 1)
    locs <- crossing(p0 = sq, p1 = sq)
  }
  df <- locs %>%
    set_colnames(c("p0", "p1")) %>%
    as_tibble()

  D <- as.matrix(fields::rdist(as.matrix(locs)))
  Sigma_c <- mat_corfun(phi = phi_c, D = D)
  proc_c <- t(mvtnorm::rmvnorm(1,sigma=Sigma_c,method="svd"))
  Sigma_u <- mat_corfun(phi = phi_u, D = D)
  proc_u <- t(mvtnorm::rmvnorm(1,sigma=Sigma_u,method="svd"))
  Sigma_y <- mat_corfun(phi = phi_y, D = D)
  proc_y <- t(mvtnorm::rmvnorm(1,sigma=Sigma_y,method="svd"))
  if(scale_proc){
    proc_c <- matrix(scale(proc_c[,1]), ncol=1); proc_u <- matrix(scale(proc_u[,1]), ncol=1); proc_y <- matrix(scale(proc_y[,1]), ncol=1)
  }
  colnames(proc_c) <- paste0("proc_c_", phi_c)
  colnames(proc_u) <- paste0("proc_u_", phi_u)
  colnames(proc_y) <- paste0("proc_y_", phi_y)

  proc_c %<>% as_tibble()
  proc_u %<>% as_tibble()
  proc_y %<>% as_tibble()

  #non-spatial component (contructed from a std. normal  variable) in exposure and error term in outcome
  sigsq_seq <- c(1, sigsq_y_true_seq)
  mat <- matrix(NA, nrow = n_obs, ncol = length(sigsq_seq))
  colnames(mat) <- c("proc_ns_std",paste0("err_y_", sigsq_y_true_seq))
  for (j in seq_along(sigsq_seq)) {
    mat[, j] <- rnorm(n = n_obs, 0, sqrt(sigsq_seq[j]))
  }
  mat %<>% as_tibble()

  df %<>%
    bind_cols(proc_c) %>%
    bind_cols(proc_u) %>%
    bind_cols(proc_y) %>%
    bind_cols(mat)

  df_comps <- df
  df_comps
}

#function to independently generate unmeasured spatial effect modifier "proc_e"
gen_data_comps_ze <- function(size, n_obs = 1000, phi_seq, sigsq_y_true_seq, rand_points = TRUE, mat_corfun = mat_corfun_exp, scale_proc = T) {

  if (rand_points) {
    #n_obs=200
    locs <- matrix(runif(2*n_obs, 0, sqrt(size)),ncol=2)
  } else {
    ## regular grid of points
    npt <- ceiling(sqrt(n_obs))
    sq <- (1:npt)/(npt + 1)
    locs <- crossing(p0 = sq, p1 = sq)
  }
  df <- locs %>%
    set_colnames(c("p0", "p1")) %>%
    as_tibble()

  D <- as.matrix(fields::rdist(as.matrix(locs)))
  proc_c <- proc_u <- proc_e <- matrix(NA, n_obs, length(phi_seq))
  for (j in seq_along(phi_seq)) {
    Sigma <- mat_corfun(phi = phi_seq[j], D = D)
    if (n_obs <= 10000) {
      proc <- t(mvtnorm::rmvnorm(3,sigma=Sigma,method="svd"))
    } else {
      print("larger")
      cholS <- chol(Sigma)
      z <- matrix(rnorm(2*n_obs),ncol = 2)
      proc <- t(cholS)%*%z
    }
    if (scale_proc) proc %<>% scale()
    proc_c[, j] <- proc[, 1]
    proc_u[, j] <- proc[, 2]
    proc_e[, j] <- proc[, 3]
  }
  colnames(proc_c) <- paste0("proc_c_", c(0.04,0.15,0.6))
  colnames(proc_u) <- paste0("proc_u_", c(0.04,0.15,0.6))
  colnames(proc_e) <- paste0("proc_e_", c(0.04,0.15,0.6))
  proc_c %<>% as_tibble()
  proc_u %<>% as_tibble()
  proc_e %<>% as_tibble()

  #non-spatial component (contructed from a std. normal  variable) in exposure and error term in outcome
  sigsq_seq <- c(1, sigsq_y_true_seq)
  mat <- matrix(NA, nrow = n_obs, ncol = length(sigsq_seq))
  colnames(mat) <- c("proc_ns_std",paste0("err_y_", sigsq_y_true_seq))
  for (j in seq_along(sigsq_seq)) {
    mat[, j] <- rnorm(n = n_obs, 0, sqrt(sigsq_seq[j]))
  }
  mat %<>% as_tibble()

  df %<>%
    bind_cols(proc_c) %>%
    bind_cols(proc_u) %>%
    bind_cols(proc_e) %>%
    bind_cols(mat)

  df_comps <- df
  df_comps
}

## deltaE is a 9X3 df
gen_data_from_comps <- function(df_comps, beta, deltaE, coef_gen_y, coef_gen_c, coef_gen_u, coef_inter, power_effect=1, thresh, Pns, phi_c, phi_u, sigsq_y_true) {

  if (phi_c == 0) {
    z_c <- 0
  } else {
    z_c <- with(df_comps, get(paste0("proc_c_", phi_c)))
  }

  if (phi_u == 0) {
    z_u <- 0
  } else {
    z_u <- with(df_comps, get(paste0("proc_u_", phi_u)))
  }

  z_n <- with(df_comps, get("proc_ns_std"))

  #no confounding process or no independent process, the coef is 3, otherwise it is sqrt(9/2)
  #coef_gen_x = ifelse(phi_c == 0 | phi_u == 0, 3, sqrt(9/2))
  theta = sqrt(Pns*(coef_gen_c^2+coef_gen_u^2)/(1-Pns))
  x_lat = coef_gen_c*z_c + coef_gen_u*z_u + theta*z_n
  X_indv <- ifelse(x_lat > thresh, 1, 0)

  eps_y <- with(df_comps, get(paste0("err_y_", sigsq_y_true)))
  #deltae <- deltaE %>% filter(phic == phi_c, phiu == phi_u)
  #deltae <- as.matrix(deltae)[which(c(0,0.5,0.9)==Pns)+2]
  #beta_lat <- beta/deltae

  df <- df_comps %>%
    mutate(z_c = z_c,
           z_u = z_u,
           z_n_mod = theta*z_n,
           x_latent = x_lat,
           X_indv = X_indv,
           mu_y_indv = beta*X_indv + coef_gen_y*z_c + coef_inter*X_indv*(z_c^power_effect),
           effect = beta + coef_inter*z_c^power_effect,
           y_indv = mu_y_indv + eps_y
           #y_lat = beta_lat*x_lat + coef_gen_y*z_c + eps_y
           ) %>% dplyr::select(-contains("proc_"))

  df
}

gen_data_from_comps_Zy <- function(df_comps, beta, coef_y_c, coef_x_c, coef_x_u, coef_y_u, coef_inter, power_effect=1, thresh, Pns, phi_c, phi_u, phi_y, sigsq_y_true) {

  if (phi_c == 0) {
    z_c <- 0
  } else {
    z_c <- with(df_comps, get(paste0("proc_c_", phi_c)))
  }

  if (phi_u == 0) {
    z_u <- 0
  } else {
    z_u <- with(df_comps, get(paste0("proc_u_", phi_u)))
  }

  if (phi_y == 0) {
    z_y <- 0
  } else {
    z_y <- with(df_comps, get(paste0("proc_y_", phi_y)))
  }

  z_n <- with(df_comps, get("proc_ns_std"))

  #no confounding process or no independent process, the coef is 3, otherwise it is sqrt(9/2)
  #coef_gen_x = ifelse(phi_c == 0 | phi_u == 0, 3, sqrt(9/2))
  theta = sqrt(Pns*(coef_x_c^2+coef_x_u^2)/(1-Pns))
  x_lat = coef_x_c*z_c + coef_x_u*z_u + theta*z_n
  X_indv <- ifelse(x_lat > thresh, 1, 0)

  eps_y <- with(df_comps, get(paste0("err_y_", sigsq_y_true)))

  df <- df_comps %>%
    mutate(z_c = z_c,
           z_u = z_u,
           z_y = z_y,
           z_n_mod = theta*z_n,
           x_latent = x_lat,
           X_indv = X_indv,
           mu_y_indv = beta*X_indv + coef_y_c*z_c + coef_y_u*z_y + coef_inter*X_indv*(z_c^power_effect),
           effect = beta + coef_inter*z_c^power_effect,
           y_indv = mu_y_indv + eps_y
           #y_lat = beta_lat*x_lat + coef_gen_y*z_c + eps_y
    ) %>% dplyr::select(-contains("proc_"))

  df
}

gen_data_from_comps_ze <- function(df_comps, beta, coef_y_c, coef_y_e, thresh, theta, phi_c, phi_u, phi_e, sigsq_y_true) {

  if (phi_c == 0) {
    z_c <- 0
  } else {
    z_c <- with(df_comps, get(paste0("proc_c_", phi_c)))
  }

  if (phi_u == 0) {
    z_u <- 0
  } else {
    z_u <- with(df_comps, get(paste0("proc_u_", phi_u)))
  }

  if (phi_e == 0) {
    z_e <- 0
  } else {
    z_e <- with(df_comps, get(paste0("proc_e_", phi_e)))
  }

  z_n <- with(df_comps, get("proc_ns_std"))

  #no confounding process or no independent process, the coef is 3, otherwise it is sqrt(9/2)
  coef_gen_x = ifelse(phi_c == 0 | phi_u == 0, 3, sqrt(9/3))
  x_lat = coef_gen_x*(z_c + z_u) - coef_gen_x*z_e + theta*z_n
  X_indv <- ifelse(x_lat > thresh, 1, 0)

  eps_y <- with(df_comps, get(paste0("err_y_", sigsq_y_true)))

  df <- df_comps %>%
    mutate(z_c = z_c,
           z_u = z_u,
           z_e = z_e,
           z_ctot = coef_y_c*z_c + coef_y_e*z_e,
           z_n_mod = theta*z_n,
           x_latent = x_lat,
           X_indv = X_indv,
           mu_y_indv = beta*X_indv + coef_y_c*z_c + coef_y_e*z_e,
           effect = beta + coef_y_c*z_c + coef_y_e*z_e,
           y_indv = mu_y_indv + eps_y) %>% dplyr::select(-contains("proc_"))

  df
}


gen_data_from_comps_binary <- function(df_comps, beta, coef_gen_y, coef_inter, power_effect=1, thresh_x, thresh_y, theta, phi_c, phi_u, sigsq_y_true) {

  if (phi_c == 0) {
    z_c <- 0
  } else {
    z_c <- with(df_comps, get(paste0("proc_c_", phi_c)))
  }

  if (phi_u == 0) {
    z_u <- 0
  } else {
    z_u <- with(df_comps, get(paste0("proc_u_", phi_u)))
  }

  z_n <- with(df_comps, get("proc_ns_std"))

  #no confounding process or no independent process, the coef is 3, otherwise it is sqrt(9/2)
  coef_gen_x = ifelse(phi_c == 0 | phi_u == 0, 3, sqrt(9/2))
  x_lat = coef_gen_x*(z_c + z_u) + theta*z_n
  X_indv <- ifelse(x_lat > thresh_x, 1, 0)

  eps_y <- with(df_comps, get(paste0("err_y_", sigsq_y_true)))
  y_lat = beta*X_indv + coef_gen_y*z_c + coef_inter*X_indv*(z_c^power_effect) + eps_y
  y_indv = ifelse(y_lat > thresh_y, 1, 0)

  df <- df_comps %>%
    mutate(z_c = z_c,
           z_u = z_u,
           z_n_mod = theta*z_n,
           x_latent = x_lat,
           X_indv = X_indv,
           y_latent = y_lat,
           y_indv = y_indv) %>% dplyr::select(-contains("proc_"))

  df
}



#Direct Match using Macth Package
# 1-to-Match, covariates is list argument consisting of covariates names in df
#clp_trype: "covariates", "distance", "None"
#match_obj = m$bin_att; caliper = 0.25
###Match_distance <- function(df, match_obj, caliper, case){
  #trtp <- df[match_obj$index.treated, c("p0","p1")]
  #crlp <- df[match_obj$index.control, c("p0","p1")]
  #dist <- fields::rdist.vec(trtp, crlp)
  #set the threshold r=0.25, pair_id stores the indices of pairs within the threshold
  #pair_id <- which(dist <= caliper)
  #matched_trt <- match_obj$index.treated[pair_id]
  #matched_crl <- match_obj$index.control[pair_id]
  #est
  #if (case == "bin"){
    #est = mean(df$y_binom[matched_trt]-df$y_binom[matched_crl])
    #n.treated = sum(df$X_binom)
  #} else if (case == "indv") {
    #est = mean(df$y_indv[matched_trt]- df$y_indv[matched_crl])
    #n.treated = sum(df$X_indv)
  #} else {
    #cat("not available")
  #}
  #se = 0
  #return(list(est = est, se = se, index.treated = matched_trt,
              #index.control = matched_crl,
              #ndrops.matches = sum(dist > caliper), orig.treated.nobs = n.treated))
#}
#sum(matched_crl != match_obj$index.control)
#trt <- m$bin_att$index.treated
#crl <- m$bin_att$index.control
#mean(df$y_binom[trt] - df$y_binom[crl])
#clp_type: distance, covariates or None
#est = "ATE" or "ATC"
#Match_list <- function(df, M, covariates, clp_typ, clp=NULL, weight=NULL, est, ties=TRUE){
  #if (clp_typ == "covariates") {
    #caliper = clp
    #if (length(covariates) == 1) {
      ##ties = FALSE
      #weight = NULL
    #}
  #} else {caliper = NULL}
  ##standard match
  ##m_bin_att <- Match(Y=df$y_binom, Tr=df$X_binom, X=df[,covariates], M=M,
  ##caliper = caliper, Weight.matrix = weight, ties = ties)
  ##m_bin_ate <- Match(Y=df$y_binom, Tr=df$X_binom, X=df[,covariates], M=M,
  ##caliper = caliper, Weight.matrix = weight, estimand = est, ties = ties)
  #m_indv_att <- Match(Y=df$y_indv, Tr=df$X_indv, X=df[,covariates], M=M,
                      #caliper = caliper, Weight.matrix = weight, ties = ties)
  #m_indv_ate <- Match(Y=df$y_indv, Tr=df$X_indv, X=df[,covariates], M=M,
                      #caliper = caliper, Weight.matrix = weight, estimand = est, ties = ties)
  #if (clp_typ == "distance") {
    #caliper_dist = sqrt(4*(clp^2)*sd(df$p0)*sd(df$p1)/pi)
    ##caliper_dist = clp
    ##m_bin_att <- Match_distance(df, m_bin_att, caliper=caliper_dist, "bin")
    ##m_bin_ate <- Match_distance(df, m_bin_ate, caliper=caliper_dist, "bin")
    #m_indv_att <- Match_distance(df, m_indv_att, caliper=caliper_dist, "indv")
    #m_indv_ate <- Match_distance(df, m_indv_ate, caliper=caliper_dist, "indv")
    ##return(list(bin_att = m_bin_att, bin_ate = m_bin_ate,
    #return(list(indv_att = m_indv_att , indv_ate = m_indv_ate))
  #}
  #return(list(bin_att = m_bin_att, bin_ate = m_bin_ate,
  #return(list(indv_att = m_indv_att , indv_ate = m_indv_ate))
#}

Match_prop <- function(df, M, clp_typ, clp=NULL, est, ties=TRUE){
  if(clp_typ == "distance"){
    caliper = clp
  }else if (clp_typ == "None") {
    caliper = NULL
  }else {
    cat("not available")
  }
  ##standard match
  glm <- glm(X_indv~p0 + p1, family = binomial, data = df)
  m <- Match(Y=df$y_indv, Tr=df$X_indv, X=glm$fitted, M=M,
             caliper = caliper, estimand = est, ties = ties)
  return(list(indv_ate = m))
}

# matching performance
#distance BETWEEN MATCHED PAIRS
#dist_MB <- function(df, m){
  ##trtp_bin_unm <- df %>% filter(X_binom == 1) %>% dplyr::select(p0, p1)
  #trtp_indv_unm <- df %>% filter(X_indv == 1) %>% dplyr::select(p0, p1)
  ##crlp_bin_unm <- df %>% filter(X_binom == 0) %>% dplyr::select(p0,p1)
  #crlp_indv_unm <- df %>% filter(X_indv == 0) %>% dplyr::select(p0,p1)
  ##dbin_unm_mat <- as.matrix(fields::rdist(x1 = trtp_bin_unm, x2 = crlp_bin_unm))
  #dindv_unm_mat <- as.matrix(fields::rdist(x1 = trtp_indv_unm, x2 = crlp_indv_unm))
  ##att
  ##trtp_bin_att <- df[m$bin_att$index.treated, c("p0","p1")]
  #trtp_indv_att <- df[m$indv_att$index.treated, c("p0","p1")]
  ##crlp_bin_att <- df[m$bin_att$index.control, c("p0","p1")]
  #crlp_indv_att <- df[m$indv_att$index.control, c("p0","p1")]
  ##dbin_att <- fields::rdist.vec(x1 = trtp_bin_att, x2 = crlp_bin_att)
  #dindv_att <- fields::rdist.vec(x1 = trtp_indv_att, x2 = crlp_indv_att)
  ##ate
  ##trtp_bin_ate <- df[m$bin_ate$index.treated, c("p0","p1")]
  #trtp_indv_ate <- df[m$indv_ate$index.treated, c("p0","p1")]
  ##crlp_bin_ate <- df[m$bin_ate$index.control, c("p0","p1")]
  #crlp_indv_ate <- df[m$indv_ate$index.control, c("p0","p1")]
  ##dbin_ate <- fields::rdist.vec(x1 = trtp_bin_ate, x2 = crlp_bin_ate)
  #dindv_ate <- fields::rdist.vec(x1 = trtp_indv_ate, x2 = crlp_indv_ate)
  #return(list(#dbin_before = dbin_unm_mat,
    #dindv_before = dindv_unm_mat, #dbin_att = dbin_att,
    #dindv_att = dindv_att, #dbin_ate = dbin_ate,
    #dindv_ate = dindv_ate))
#}

dist_MB_inter <- function(df, m){
  trtp_indv_unm <- df %>% filter(X_indv == 1) %>% dplyr::select(p0, p1)
  #crlp_bin_unm <- df %>% filter(X_binom == 0) %>% dplyr::select(p0,p1)
  crlp_indv_unm <- df %>% filter(X_indv == 0) %>% dplyr::select(p0,p1)
  #dbin_unm_mat <- as.matrix(fields::rdist(x1 = trtp_bin_unm, x2 = crlp_bin_unm))
  dindv_unm_mat <- as.matrix(fields::rdist(x1 = trtp_indv_unm, x2 = crlp_indv_unm))

  trtp_indv_att <- df[m$index.treated, c("p0","p1")]
  #crlp_bin_att <- df[m$bin_att$index.control, c("p0","p1")]
  crlp_indv_att <- df[m$index.control, c("p0","p1")]
  #dbin_att <- fields::rdist.vec(x1 = trtp_bin_att, x2 = crlp_bin_att)
  dindv_att <- fields::rdist.vec(x1 = trtp_indv_att, x2 = crlp_indv_att)
  return(list(#dbin_before = dbin_unm_mat,
    dindv_before = dindv_unm_mat, #dbin_att = dbin_att,
    dindv_att = dindv_att)) #dbin_ate = dbin_ate,
    #dindv_ate = dindv_ate))
}

#dist_mom <- function(dist_mat){
  ##dbin_mean_before <- mean(colMeans(dist_mat$dbin_before))
  #dindv_mean_before <- mean(colMeans(dist_mat$dindv_before))

  ##dbin_mean_att <- mean(dist_mat$dbin_att)
  ##dbin_sd_att <- sd(dist_mat$dbin_att)
  ##dbin_mean_ate <- mean(dist_mat$dbin_ate)
  ##dbin_sd_ate <- sd(dist_mat$dbin_ate)

  #dindv_mean_att <- mean(dist_mat$dindv_att)
  #dindv_sd_att <- sd(dist_mat$dindv_att)
  #dindv_mean_ate <- mean(dist_mat$dindv_ate)
  #dindv_sd_ate <- sd(dist_mat$dindv_ate)

  #return(list(#dbin_before = dbin_mean_before,dbin_att = dbin_mean_att,dbin_att_sd = dbin_sd_att,
    ##dbin_ate = dbin_mean_ate,dbin_ate_sd = dbin_sd_ate,
    #dindv_before = dindv_mean_before,dindv_att = dindv_mean_att,dindv_att_sd = dindv_sd_att,
    #dindv_ate = dindv_mean_ate, dindv_ate_sd = dindv_sd_ate))
#}

dist_mom_inter <- function(dist_mat) {
  dindv_mean_before <- mean(colMeans(dist_mat$dindv_before))

  #dbin_mean_att <- mean(dist_mat$dbin_att)
  #dbin_sd_att <- sd(dist_mat$dbin_att)
  #dbin_mean_ate <- mean(dist_mat$dbin_ate)
  #dbin_sd_ate <- sd(dist_mat$dbin_ate)

  dindv_mean_att <- mean(dist_mat$dindv_att)
  dindv_sd_att <- sd(dist_mat$dindv_att)
  #dindv_mean_ate <- mean(dist_mat$dindv_ate)
  #dindv_sd_ate <- sd(dist_mat$dindv_ate)

  return(list(#dbin_before = dbin_mean_before,dbin_att = dbin_mean_att,dbin_att_sd = dbin_sd_att,
    #dbin_ate = dbin_mean_ate,dbin_ate_sd = dbin_sd_ate,
    dindv_before = dindv_mean_before,
    dindv_att = dindv_mean_att,
    dindv_att_sd = dindv_sd_att))
    #dindv_ate = dindv_mean_ate, dindv_ate_sd = dindv_sd_ate))
}

#p0, p1
MB_p0p1 <- function(df_cpy, m) {

  #bin att
  bin_att_ksbefore0 <- MatchBalance(X_binom~p0_scale+p1_scale, data = df_cpy,
                                  match.out = m[[1]])$BeforeMatching[[1]]$ks$ks.boot.pvalue
  bin_att_ksbefore1 <- MatchBalance(X_binom~p0_scale+p1_scale, data = df_cpy,
                                   match.out = m[[1]])$BeforeMatching[[2]]$ks$ks.boot.pvalue

  bin_att_ksafter0 <- MatchBalance(X_binom~p0_scale+p1_scale, data = df_cpy,
                                 match.out = m[[1]])$AfterMatching[[1]]$ks$ks.boot.pvalue
  bin_att_ksafter1 <- MatchBalance(X_binom~p0_scale+p1_scale, data = df_cpy,
                                   match.out = m[[1]])$AfterMatching[[2]]$ks$ks.boot.pvalue
  bin_att <- data.frame(p0 = c(bin_att_ksbefore0, bin_att_ksafter0), p1=c(bin_att_ksbefore1, bin_att_ksafter1),
                        row.names=c("before","after"))
  #bin ate
  bin_ate_ksbefore0 <- MatchBalance(X_binom~p0_scale+p1_scale, data = df_cpy,
                                   match.out = m[[2]])$BeforeMatching[[1]]$ks$ks.boot.pvalue
  bin_ate_ksbefore1 <- MatchBalance(X_binom~p0_scale+p1_scale, data = df_cpy,
                                    match.out = m[[2]])$BeforeMatching[[2]]$ks$ks.boot.pvalue

  bin_ate_ksafter0 <- MatchBalance(X_binom~p0_scale+p1_scale, data = df_cpy,
                                  match.out = m[[2]])$AfterMatching[[1]]$ks$ks.boot.pvalue
  bin_ate_ksafter1 <- MatchBalance(X_binom~p0_scale+p1_scale, data = df_cpy,
                                  match.out = m[[2]])$AfterMatching[[2]]$ks$ks.boot.pvalue
  bin_ate <- data.frame(p0 = c(bin_ate_ksbefore0, bin_ate_ksafter0), p1=c(bin_ate_ksbefore1, bin_ate_ksafter1),
                        row.names=c("before","after"))
  #indv att
  indv_att_ksbefore0 <- MatchBalance(X_indv~p0_scale+p1_scale, data = df_cpy,
                                    match.out = m[[3]])$BeforeMatching[[1]]$ks$ks.boot.pvalue
  indv_att_ksbefore1 <- MatchBalance(X_indv~p0_scale+p1_scale, data = df_cpy,
                                    match.out = m[[3]])$BeforeMatching[[2]]$ks$ks.boot.pvalue

  indv_att_ksafter0 <- MatchBalance(X_indv~p0_scale+p1_scale, data = df_cpy,
                                   match.out = m[[3]])$AfterMatching[[1]]$ks$ks.boot.pvalue
  indv_att_ksafter1 <- MatchBalance(X_indv~p0_scale+p1_scale, data = df_cpy,
                                   match.out = m[[3]])$AfterMatching[[2]]$ks$ks.boot.pvalue
  indv_att <- data.frame(p0 = c(indv_att_ksbefore0, indv_att_ksafter0), p1=c(indv_att_ksbefore1, indv_att_ksafter1),
                        row.names=c("before","after"))
  #indv ate
  indv_ate_ksbefore0 <- MatchBalance(X_indv~p0_scale+p1_scale, data = df_cpy,
                                    match.out = m[[4]])$BeforeMatching[[1]]$ks$ks.boot.pvalue
  indv_ate_ksbefore1 <- MatchBalance(X_indv~p0_scale+p1_scale, data = df_cpy,
                                    match.out = m[[4]])$BeforeMatching[[2]]$ks$ks.boot.pvalue

  indv_ate_ksafter0 <- MatchBalance(X_indv~p0_scale+p1_scale, data = df_cpy,
                                   match.out = m[[4]])$AfterMatching[[1]]$ks$ks.boot.pvalue
  indv_ate_ksafter1 <- MatchBalance(X_indv~p0_scale+p1_scale, data = df_cpy,
                                   match.out = m[[4]])$AfterMatching[[2]]$ks$ks.boot.pvalue
  indv_ate <- data.frame(p0 = c(indv_ate_ksbefore0, indv_ate_ksafter0), p1=c(indv_ate_ksbefore1, indv_ate_ksafter1),
                         row.names=c("before","after"))

  return(list(bin_att, bin_ate, indv_att, indv_ate))

}

#estimate
#extract_ests_match <- function(m, modname){
  #if (modname == "bin") {
    #summary <- c(m$bin_ate$est, m$bin_ate$se, m$bin_ate$orig.treated.nobs, m$bin_ate$ndrops.matches,
                 #m$bin_att$est, m$bin_att$se, m$bin_att$orig.treated.nobs, m$bin_att$ndrops.matches) %>%
      #set_names(c("bin_ate_est","bin_ate_se","bin_ate_n.treat","bin_ate_n.drop",
       #         "bin_att_est","bin_att_se","bin_att_n.treat","bin_att_n.drop"))%>%
      #as.list() %>% as_tibble()
  #} else if (modname == "indv") {
    #summary <- c(m$indv_ate$est, m$indv_ate$se, m$indv_ate$orig.treated.nobs, m$indv_ate$ndrops.matches,
     #            m$indv_att$est, m$indv_att$se, m$indv_att$orig.treated.nobs, m$indv_att$ndrops.matches) %>%
    #  set_names(c("indv_ate_est","indv_ate_se","indv_ate_n.treat","indv_ate_n.drop",
      #          "indv_att_est","indv_att_se","indv_att_n.treat","indv_att_n.drop"))%>%
      #as.list() %>% as_tibble()
  #} else {
   # cat("not available")
  #}
  #summary
#}
#bin <- extract_ests_match(m,"bin")
#m$bin_att

balance_zc_wait <- function(m, df, modname){
  if (modname == "bin") {
    bin_att <- MatchBalance(X_binom~z_c, match.out = m$bin_att, data = df)
    bin_att_sdiff <- bin_att$AfterMatching[[1]]$sdiff/100
    bin_att_diff <- abs(bin_att$AfterMatching[[1]]$mean.Tr - bin_att$AfterMatching[[1]]$mean.Co)
    bin_att_ks <- bin_att$AfterMatching[[1]]$ks$ks.boot.pvalue

    bin_ate <- MatchBalance(X_binom~z_c, match.out = m$bin_ate, data = df)
    bin_ate_sdiff <- bin_ate$AfterMatching[[1]]$sdiff/100
    bin_ate_diff <- abs(bin_ate$AfterMatching[[1]]$mean.Tr - bin_ate$AfterMatching[[1]]$mean.Co)
    bin_ate_ks <- bin_ate$AfterMatching[[1]]$ks$ks.boot.pvalue
    balance <- c(bin_ate_sdiff, bin_ate_diff, bin_ate_ks,
                 bin_att_sdiff, bin_att_diff, bin_att_ks) %>%
      set_names(c("bin_ate_sdiff", "bin_ate_diff","bin_ate_ks",
                  "bin_att_sdiff", "bin_att_diff", "bin_att_ks"))%>%
      as.list() %>% as_tibble()
  } else {
    indv_att <- MatchBalance(X_indv~z_c, match.out = m$indv_att, data = df)
    indv_att_sdiff <- indv_att$AfterMatching[[1]]$sdiff/100
    indv_att_diff <- abs(indv_att$AfterMatching[[1]]$mean.Tr - indv_att$AfterMatching[[1]]$mean.Co)
    indv_att_ks <- indv_att$AfterMatching[[1]]$ks$ks.boot.pvalue

    indv_ate <- MatchBalance(X_indv~z_c, match.out = m$indv_ate, data = df)
    indv_ate_sdiff <- indv_ate$AfterMatching[[1]]$sdiff/100
    indv_ate_diff <- abs(indv_ate$AfterMatching[[1]]$mean.Tr - indv_ate$AfterMatching[[1]]$mean.Co)
    indv_ate_ks <- indv_ate$AfterMatching[[1]]$ks$ks.boot.pvalue
    balance <- c(indv_ate_sdiff, indv_ate_diff, indv_ate_ks,
                 indv_att_sdiff, indv_att_diff, indv_att_ks) %>%
      set_names(c("indv_ate_sdiff", "indv_ate_diff", "indv_ate_ks",
                  "indv_att_sdiff", "indv_att_diff", "indv_att_ks"))%>%
      as.list() %>% as_tibble()
  }

  return(balance)
}


#balance_zc <- function(m, df, case) {
  #if (case == "bin") {
    #std_bin <- sd(df$z_c[df$X_binom == 1])
    #trt_zc_att <- df$z_c[m$bin_att$index.treated]
    #crl_zc_att <- df$z_c[m$bin_att$index.control]
    #trt_zc_ate <- df$z_c[m$bin_ate$index.treated]
    #crl_zc_ate <- df$z_c[m$bin_ate$index.control]
    ##att
    #sdiff_att <- abs(mean(trt_zc_att) - mean(crl_zc_att))/std_bin
    #vr_att <- var(trt_zc_att)/var(crl_zc_att)
    #ks_att <- ks.boot(trt_zc_att, crl_zc_att)
    #ks_att_stat <- ks_att$ks$statistic; ks_att_p <- ks_att$ks$p.value
    ##ate
    #sdiff_ate <- abs(mean(trt_zc_ate) - mean(crl_zc_ate))/std_bin
    #vr_ate <- var(trt_zc_ate)/var(crl_zc_ate)
    #ks_ate <- ks.boot(trt_zc_ate, crl_zc_ate)
    #ks_ate_stat <- ks_ate$ks$statistic; ks_ate_p <- ks_ate$ks$p.value
  #} else if (case == "indv") {
    #s_trt_att <- var(df$z_c[m$indv_att$index.treated]); s_crl_att <- var(df$z_c[m$indv_att$index.control])
    #s_trt_ate <- var(df$z_c[m$indv_ate$index.treated]); s_crl_ate <- var(df$z_c[m$indv_ate$index.control])
    ##the number of treated is approximately equaling to that of untreated
    #std_att <- sqrt((s_trt_att + s_crl_att)/2)
    #std_ate <- sqrt((s_trt_ate + s_crl_ate)/2)
    #trt_zc_att <- df$z_c[m$indv_att$index.treated]
    #crl_zc_att <- df$z_c[m$indv_att$index.control]
    #trt_zc_ate <- df$z_c[m$indv_ate$index.treated]
    #crl_zc_ate <- df$z_c[m$indv_ate$index.control]
    ##att
    #sdiff_att <- abs(mean(trt_zc_att) - mean(crl_zc_att))/std_att
    #diff_att <- abs(mean(trt_zc_att) - mean(crl_zc_att))
    #vr_att <- var(trt_zc_att)/var(crl_zc_att)
    #ks_att <- ks.boot(trt_zc_att, crl_zc_att)
    #ks_att_stat <- ks_att$ks$statistic; ks_att_p <- ks_att$ks$p.value
    ##ate
    #sdiff_ate <- abs(mean(trt_zc_ate) - mean(crl_zc_ate))/std_ate
    #diff_ate <- abs(mean(trt_zc_ate) - mean(crl_zc_ate))
    #vr_ate <- var(trt_zc_ate)/var(crl_zc_ate)
    #ks_ate <- ks.boot(trt_zc_ate, crl_zc_ate)
    #ks_ate_stat <- ks_ate$ks$statistic; ks_ate_p <- ks_ate$ks$p.value
  #} else {
   # cat("not available")
  #}
  #return(list(sdiff_att = sdiff_att, diff_att = diff_att, vr_att = vr_att, ks_att_stat = ks_att_stat, ks_att_p = ks_att_p,
  #            sdiff_ate = sdiff_ate, diff_ate = diff_ate, vr_ate = vr_ate, ks_ate_stat = ks_ate_stat, ks_ate_p = ks_ate_p))
#}

balance_zc_inter <- function(m, df) {
  s_trt_att <- var(df$z_c[m$index.treated]); s_crl_att <- var(df$z_c[m$index.control])
  #s_trt_ate <- var(df$z_c[m$indv_ate$index.treated]); s_crl_ate <- var(df$z_c[m$indv_ate$index.control])
  #the number of treated is approximately equaling to that of untreated
  std_att <- sqrt((s_trt_att + s_crl_att)/2)
  #std_ate <- sqrt((s_trt_ate + s_crl_ate)/2)
  trt_zc_att <- df$z_c[m$index.treated]
  crl_zc_att <- df$z_c[m$index.control]
  #trt_zc_ate <- df$z_c[m$indv_ate$index.treated]
  #crl_zc_ate <- df$z_c[m$indv_ate$index.control]
  #att
  sdiff_att <- abs(mean(trt_zc_att) - mean(crl_zc_att))/std_att
  diff_att <- abs(mean(trt_zc_att) - mean(crl_zc_att))
  vr_att <- var(trt_zc_att)/var(crl_zc_att)
  ks_att <- ks.boot(trt_zc_att, crl_zc_att)
  ks_att_stat <- ks_att$ks$statistic; ks_att_p <- ks_att$ks$p.value
  #ate
  #sdiff_ate <- abs(mean(trt_zc_ate) - mean(crl_zc_ate))/std_ate
  #diff_ate <- abs(mean(trt_zc_ate) - mean(crl_zc_ate))
  #vr_ate <- var(trt_zc_ate)/var(crl_zc_ate)
  #ks_ate <- ks.boot(trt_zc_ate, crl_zc_ate)
  #ks_ate_stat <- ks_ate$ks$statistic; ks_ate_p <- ks_ate$ks$p.value

  return(list(sdiff_att = sdiff_att, diff_att = diff_att, vr_att = vr_att, ks_att_stat = ks_att_stat, ks_att_p = ks_att_p))
              #sdiff_ate = sdiff_ate, diff_ate = diff_ate, vr_ate = vr_ate, ks_ate_stat = ks_ate_stat, ks_ate_p = ks_ate_p))
}


balance_ndrop <- function(res,phi_seq,te,case){
  ndrop <- c(0,0,0,0)
  for (c in seq_along(phi_seq)) {
    dat_phic <- res %>% filter(phi_c == phi_seq[c])
    if (te == 0 & phi_seq[c]==0) {
      #remove 0 from list for phi_i
      phi_seq_i = phi_seq[2:4]
    } else {
      phi_seq_i = phi_seq
    }
    #n.match = 1:row
    for (i in seq_along(phi_seq_i)) {
      dat_phici <- dat_phic %>% filter(phi_i == phi_seq_i[i], theta == te)
      if (case == "bin") {
        #cat("phi_c=",phi_seq[c], "phi_i=", phi_seq_i[i],
        #"\nunique bin_ate.n.match=",unique(dat_phici$bin_ate_n.match),
        #"\nunique bin_att.n.match=",unique(dat_phici$bin_att_n.match),"\n")
        Cov_att <- as.data.frame(dat_phici %>% filter(caliper == "covariates")
                                 %>% dplyr::select(bin_att_n.drop, bin.att_dm, bin.att_dse, bin.att_sdiff,
                                                   bin.att_vr, bin.att_ks_stat, bin.att_ks_p))
        Cov_ate <- as.data.frame(dat_phici %>% filter(caliper == "covariates")
                                 %>% dplyr::select(bin_ate_n.drop, bin.ate_dm, bin.ate_dse, bin.ate_sdiff,
                                                   bin.ate_vr, bin.ate_ks_stat, bin.ate_ks_p))
        dist_att <- as.data.frame(dat_phici %>% filter(caliper == "distance")
                                  %>% dplyr::select(bin_att_n.drop, bin.att_dm, bin.att_dse, bin.att_sdiff,
                                                    bin.att_vr, bin.att_ks_stat, bin.att_ks_p))
        dist_ate <-  as.data.frame(dat_phici %>% filter(caliper == "distance")
                                   %>% dplyr::select(bin_ate_n.drop, bin.ate_dm, bin.ate_dse, bin.ate_sdiff,
                                                     bin.ate_vr, bin.ate_ks_stat, bin.ate_ks_p))
        ndropi <- data.frame(Cov_att = mean(Cov_att$bin_att_n.drop), Cov_ate = mean(Cov_ate$bin_ate_n.drop),
                             dist_att = mean(dist_att$bin_att_n.drop), dist_ate = mean(dist_ate$bin_ate_n.drop))
      } else if (case == "indv") {
        #cat("phi_c=",phi_seq[c], "phi_i=", phi_seq_i[i],
        #"\nunique indv_ate.n.match=",unique(dat_phici$indv_ate_n.match),
        #"\nunique indv_att.n.match=",unique(dat_phici$indv_att_n.match),"\n")
        Cov_att <- as.data.frame(dat_phici %>% filter(caliper == "covariates")
                                 %>%dplyr::select(indv_att_n.drop, indv.att_dm, indv.att_dse, indv.att_sdiff,
                                                  indv.att_vr, indv.att_ks_stat, indv.att_ks_p))
        Cov_ate <- as.data.frame(dat_phici %>% filter(caliper == "covariates")
                                 %>% dplyr::select(indv_ate_n.drop, indv.ate_dm, indv.ate_dse, indv.ate_sdiff,
                                                   indv.ate_vr, indv.ate_ks_stat, indv.ate_ks_p))
        dist_att <- as.data.frame(dat_phici %>% filter(caliper == "distance")
                                  %>%dplyr::select(indv_att_n.drop, indv.att_dm, indv.att_dse, indv.att_sdiff,
                                                   indv.att_vr, indv.att_ks_stat, indv.att_ks_p))
        dist_ate <- as.data.frame(dat_phici %>% filter(caliper == "distance")
                                  %>% dplyr::select(indv_ate_n.drop, indv.ate_dm, indv.ate_dse, indv.ate_sdiff,
                                                    indv.ate_vr, indv.ate_ks_stat, indv.ate_ks_p))
        ndropi <- data.frame(Cov_att = mean(Cov_att$indv_att_n.drop), Cov_ate = mean(Cov_ate$indv_ate_n.drop),
                             dist_att = mean(dist_att$indv_att_n.drop), dist_ate = mean(dist_ate$indv_ate_n.drop))
      } else {
        cat("not available")
      }
      #combine the row for each scenario
      ndrop <- rbind(ndrop, ndropi)
    }
  }
  ndrop[2:nrow(ndrop),]
}

match_bias <- function(method, case, phic, phii, te, typ, caliper, est_method, k=NULL) {
  tabble <- NULL
  if(method == "M1") {
    for (i in 1:100) {
      df_comps <- readRDS(paste0(output_root, "/datasets_demo/",i,".rds"))
      df <- gen_data_from_comps(df_comps, seed_bin = unique(df_comps$seed_dataset),
                                intercept_x = 0.5, coef_gen_y = 1,
                                theta=te, phi_c=phic, phi_i=phii, sigsq_y_true = 4)
      m1 <- Match_list(df, M=1, c("p0", "p1"), clp_typ=typ, clp=caliper, weight=diag(2), est="ATE")
      m2 <- Match_list(df, M=1, c("p0", "p1"), clp_typ=typ, clp=caliper, weight=diag(2), est="ATC")
      #naive
      if (case == "bin") {
        df_att <- rbind(df[m1$bin_att$index.treated,], df[m1$bin_att$index.control,])
        df_atc <- rbind(df[m2$bin_ate$index.treated,], df[m2$bin_ate$index.control,])
        df_ate <- rbind(df[m1$bin_ate$index.treated,], df[m1$bin_ate$index.control,])
        #m2$bin_ate is actually atc
        p <- length(m1$bin_att$index.treated)/(length(m1$bin_att$index.treated) + length(m2$bin_ate$index.control))
        #att
        d_sigy_att <- mean(df$err_y_4[m1$bin_att$index.treated]) - mean(df$err_y_4[m1$bin_att$index.control])
        d_zc_att <- mean(df$z_c[m1$bin_att$index.treated]) - mean(df$z_c[m1$bin_att$index.control])
        #atc
        d_sigy_atc <- mean(df$err_y_4[m2$bin_ate$index.treated])- mean(df$err_y_4[m2$bin_ate$index.control])
        d_zc_atc <- mean(df$z_c[m2$bin_ate$index.treated]) - mean(df$z_c[m2$bin_ate$index.control])

        #analysis stage
        if (est_method == "naive") {
          fit_att <- lm(y_binom~X_binom,data=df_att)
          fit_atc <- lm(y_binom~X_binom,data=df_atc)
          fit_ate <- lm(y_binom~X_binom,data=df_ate)
        } else if(est_method == "PS") {
          fit_att <- gam(y_binom ~ X_binom + s(p0, p1, k = k), data = df_att)
          fit_atc <- gam(y_binom ~ X_binom + s(p0, p1, k = k), data = df_atc)
          fit_ate <- gam(y_binom ~ X_binom + s(p0, p1, k = k), data = df_ate)
        } else if (est_method == "fix_df") {
          fit_att <- gam(y_binom ~ X_binom + s(p0, p1, k = k, fx=TRUE), data = df_att)
          fit_atc <- gam(y_binom ~ X_binom + s(p0, p1, k = k, fx=TRUE), data = df_atc)
          fit_ate <- gam(y_binom ~ X_binom + s(p0, p1, k = k, fx=TRUE), data = df_ate)
        } else {cat("not available")}
        tabble <- rbind(tabble, c(phic,phii,i,p, m1$bin_att$est, d_zc_att, d_sigy_att, fit_att$coefficients["X_binom"],
                                  m2$bin_ate$est, d_zc_atc, d_sigy_atc, fit_atc$coefficients["X_binom"],
                                  m1$bin_ate$est, fit_ate$coefficients["X_binom"]))
      } else if (case == "indv"){
        df_att <- rbind(df[m1$indv_att$index.treated,], df[m1$indv_att$index.control,])
        df_atc <- rbind(df[m2$indv_ate$index.treated,], df[m2$indv_ate$index.control,])
        df_ate <- rbind(df[m1$indv_ate$index.treated,], df[m1$indv_ate$index.control,])
        p <- length(m1$indv_att$index.treated)/(length(m1$indv_att$index.treated) + length(m2$indv_ate$index.control))
        #att
        d_sigy_att <- mean(df$err_y_4[m1$indv_att$index.treated]) - mean(df$err_y_4[m1$indv_att$index.control])
        d_zc_att <- mean(df$z_c[m1$indv_att$index.treated]) - mean(df$z_c[m1$indv_att$index.control])
        #atc
        d_sigy_atc <- mean(df$err_y_4[m2$indv_ate$index.treated])- mean(df$err_y_4[m2$indv_ate$index.control])
        d_zc_atc <- mean(df$z_c[m2$indv_ate$index.treated]) - mean(df$z_c[m2$indv_ate$index.control])

        if (est_method == "naive") {
          fit_att <- lm(y_indv~X_indv,data=df_att)
          fit_atc <- lm(y_indv~X_indv,data=df_atc)
          fit_ate <- lm(y_indv~X_indv,data=df_ate)
        } else if(est_method == "PS") {
          #fit_att <- gam(y_indv ~ X_indv + s(p0, p1, k = k), data = df_att)
          #fit_atc <- gam(y_indv ~ X_indv + s(p0, p1, k = k), data = df_atc)
          fit_ate <- gam(y_indv ~ X_indv + s(p0, p1, k = k), data = df_ate)
        } else if (est_method == "fix_df") {
          #fit_att <- gam(y_indv ~ X_indv + s(p0, p1, k = k, fx=TRUE), data = df_att)
          #fit_atc <- gam(y_indv ~ X_indv + s(p0, p1, k = k, fx=TRUE), data = df_atc)
          fit_ate <- gam(y_indv ~ X_indv + s(p0, p1, k = k, fx=TRUE), data = df_ate)
        } else {cat("not available")}
        tabble <- rbind(tabble, c(phic,phii,i,p, m1$indv_att$est, d_zc_att, d_sigy_att,#fit_att$coefficients["X_indv"],
                                  m2$indv_ate$est, d_zc_atc, d_sigy_atc,#fit_atc$coefficients["X_indv"],
                                  m1$indv_ate$est, fit_ate$coefficients["X_indv"]))
      } else {
        cat("not available")
      }
    }
    colnames(tabble) <- c("phi_c", "phi_i", "dataset","proportion","ATT","diff_zc_att","diff_sig_att",#"est_att",
                          "ATc","diff_zc_atc","diff_sig_atc",#"est_atc",
                          "ATE", "est_ate")
    df <- as.data.frame(tabble)
  }
  return(df)
}

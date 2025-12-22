target_match <- function(df_comb, a, b, caliper, clp=NULL, match_method){
  target_id <- which(df_comb$X_indv == a)
  comp_id <- which(df_comb$X_indv == b)

  match.matrix<-matrix(NA, nrow = length(target_id), ncol = 1)
  for (i in seq_along(target_id)) {
    d2 <- (df_comb$p0[target_id[i]]-df_comb$p0[comp_id])^2 + (df_comb$p1[target_id[i]]-df_comb$p1[comp_id])^2
    if(caliper != "None") {
      match_id <- ifelse(sqrt(min(d2))<clp,comp_id[which.min(d2)],NA)
      target_id[i] <- ifelse(sqrt(min(d2))<clp,target_id[i],NA)
    } else {
      match_id <- comp_id[which.min(d2)]
    }
    match.matrix[i,] <- match_id
  }
  rownames(match.matrix) <- target_id
  return(list(match.matrix = match.matrix))
}
#match_method = "nearest
#m <- Matchit_list(df, caliper= res0$caliper, clp=res0$clp,
 #            rep = res0$match == "nearest",
  #           estimand=estm, match_method = res0$match)
cluster_nearest <- function(m) {
  df_pair <- data.frame(trt_id = as.character(m$index.treated),
                        crl_id = as.character(m$index.control),
                        cluster_id = 0)
  g = graph_from_data_frame(df_pair, directed = F,
                            vertices = as.character(1:1000))
  trt_left <- as.character(m$index.treated)
  i=1
  while (length(trt_left)>0) {
    trt <- trt_left[1]
    connect <- ego(g, 1000, trt)[[1]]
    df_pair$cluster_id[which(df_pair$trt_id %in% connect)] = i
    i = i + 1
    # remove connected trts
    trt_left <- setdiff(trt_left, connect)
  }
  # duplicated data
  m$df_dup$cluster_id <- rep(df_pair$cluster_id,2)
  # weighted data
  unique_cluster <- m$df_dup %>% group_by(obs.id) %>% summarise(cluster_id = mean(cluster_id))
  m$df_weight %<>% left_join(unique_cluster, by="obs.id")
  return(m)
}

Match_large <- function(df, trt_value, crl_value, caliper, clp, rep = TRUE, est="ATE", match_method, formula){
  if(match_method == "genetic"){
    if(caliper == "None"){clp <- NULL}
    stdclp <- ifelse(caliper == "None", TRUE, FALSE)
    m.out_ATT <-  matchit(as.formula(formula) , data = df, method = "genetic",
                          replace = T, distance = "mahalanobis", caliper=clp, std.caliper=stdclp,
                          pop.size = length(unique(df$p0))/10)
    m.out_ATC <- NA
    if(est == "ATE"){
      m.out_ATC <- matchit(as.formula(formula) , data = df, method = "genetic",
                           replace = T, distance = "mahalanobis", estimand = "ATC", caliper=clp, std.caliper=stdclp,
                           pop.size = length(unique(df$p0))/10)
    }
  } else {
    df$obs.id <- 1:dim(df)[1]
    trt_id <- which(df$X_indv==trt_value); crl_id <- which(df$X_indv==crl_value)
    if(match_method != "full"){
      m.out_ATT <- target_match(df_comb = df, a=trt_value, b=crl_value, caliper, clp, match_method)

      m.out_ATC <- NA
      if(est == "ATE"){
        m.out_ATC <- target_match(df_comb = df, a=crl_value, b=trt_value, caliper, clp, match_method)
      }
    } else {
      cat("not applicable")
    }
  }
  trt_mATT <- as.numeric(rownames(m.out_ATT$match.matrix))
  crl_mATC <- as.numeric(rownames(m.out_ATC$match.matrix))
  rownames(m.out_ATT$match.matrix) <- NULL
  rownames(m.out_ATC$match.matrix) <- NULL
  id.treated <- c(trt_mATT,m.out_ATC$match.matrix[,1])
  id.control <- c(m.out_ATT$match.matrix[,1],crl_mATC)
  id.mtrt <- id.treated[id.treated != 0]
  id.mcrl <- id.control[id.control != 0]

  df_dup <- df[c(id.mtrt,id.mcrl),]
  wt <- table(id.mtrt)/mean(table(id.mtrt))
  id.trt <- as.numeric(names(wt))
  wc <- table(id.mcrl)/mean(table(id.mcrl))
  id.crl <- as.numeric(names(wc))
  df_weight <- df[c(id.trt, id.crl),]
  df_weight$weights <- c(wt, wc)
  return(list(df_dup =df_dup,
              df_weight = df_weight,
              index.treated = id.mtrt,
              index.control = id.mcrl,
              #for full match and est = ATE, ATT_obj is ATE_obj
              ATT_obj = m.out_ATT, ATC_obj = m.out_ATC))
}
#caliper = "distance"
#match_method = "nearest"
#clp=0.08
Matchit_list <- function(df, caliper, clp, rep = TRUE, est="ATE", match_method) {
  df$obs.id <- 1:dim(df)[1]
  trt_id <- which(df$X_indv==1); crl_id <- which(df$X_indv==0)
  dist_tc <-  optmatch::match_on(X_indv~p0+p1,method = "euclidean", data=df)
  if(est == "ATE") {
    df1 <- df
    #change treated/control labels
    df1[trt_id,]$X_indv = 0
    df1[crl_id,]$X_indv = 1
    dist_ct <- optmatch::match_on(X_indv~p0+p1,method = "euclidean", data=df1)
  }
  if(caliper == "None"){clp <- NULL}
  stdclp <- ifelse(caliper == "None", TRUE, FALSE)
  #ATT & ATC
  if(match_method == "full"){
    m.out_ATT <- matchit(X_indv~p0+p1, data = df, method = match_method,
                         estimand = est, distance = dist_tc, std.caliper=stdclp, caliper = clp)
    #For full matching, the estimand can only be ATE
    m.out_ATC <- NA
  } else {
    m.out_ATT <- matchit(X_indv~p0+p1, data = df, method = match_method,
                         replace = rep, distance = dist_tc, std.caliper=stdclp, caliper = clp)
    if(est == "ATE"){
      m.out_ATC <- matchit(X_indv~p0+p1, data = df1, method = match_method,
                           replace = rep, distance = dist_ct, std.caliper=stdclp, caliper = clp)
    } else {
      m.out_ATC <- NA
    }
  }
  df_dup <- NA; id.treated <- id.control <- 0
  md_ATT <- match.data(m.out_ATT)
  df_weight <- md_ATT; trt_ATC <- NA
  if(match_method != "full"){
    crl_mATT <- as.numeric(m.out_ATT$match.matrix[,1])
    trt_mATT <- trt_id[!is.na(crl_mATT)]
    id.treated <- trt_mATT
    id.control <- na.omit(crl_mATT)
    df_dup <- df[c(trt_mATT,crl_mATT),] %>% drop_na()
    df_dup$pair_id <- rep(seq_along(id.treated),times=2)
    if(est == "ATE"){
      trt_mATC <- as.numeric(m.out_ATC$match.matrix[,1])
      crl_mATC <- crl_id[!is.na(trt_mATC)]
      id.treated <- c(trt_mATT, na.omit(trt_mATC))
      id.control <- c(na.omit(crl_mATT),crl_mATC)
      df_dup <- df[c(id.treated,id.control),]
      df_dup$pair_id <- rep(seq_along(id.treated),times=2)
      trt_ATC <- trt_mATC
      # df weighted
      weights_df <-  data.frame(obs.id = c(names(table(id.treated)), names(table(id.control))),
                                weights = c(table(id.treated)/mean(table(id.treated)),
                                            table(id.control)/mean(table(id.control))))
      weights_df$obs.id = as.numeric(weights_df$obs.id)
      df_weight <- df %>% right_join(weights_df, by="obs.id")
    }
  }

  #fit_ATE <- lm(y_indv~X_indv, data = df_new)
  return(list(df_dup =df_dup,
              df_weight = df_weight,
              index.treated = id.treated,
              index.control = id.control,
              #for full match and est = ATE, ATT_obj is ATE_obj
              ATT_obj = m.out_ATT, ATC_obj = m.out_ATC))
}

#trt_mATC
dist_MB <- function(df, m){
  #trtp_bin_unm <- df %>% filter(X_binom == 1) %>% dplyr::select(p0, p1)
  trtp_indv_unm <- df %>% filter(X_indv == 1) %>% dplyr::select(p0, p1)
  #crlp_bin_unm <- df %>% filter(X_binom == 0) %>% dplyr::select(p0,p1)
  crlp_indv_unm <- df %>% filter(X_indv == 0) %>% dplyr::select(p0,p1)
  #dbin_unm_mat <- as.matrix(fields::rdist(x1 = trtp_bin_unm, x2 = crlp_bin_unm))
  dindv_unm_mat <- as.matrix(fields::rdist(x1 = trtp_indv_unm, x2 = crlp_indv_unm))
  #ate
  #trtp_bin_ate <- df[m$bin_ate$index.treated, c("p0","p1")]
  trtp_indv_ate <- df[m$index.treated, c("p0","p1")]
  #crlp_bin_ate <- df[m$bin_ate$index.control, c("p0","p1")]
  crlp_indv_ate <- df[m$index.control, c("p0","p1")]
  #dbin_ate <- fields::rdist.vec(x1 = trtp_bin_ate, x2 = crlp_bin_ate)
  dindv_ate <- fields::rdist.vec(x1 = trtp_indv_ate, x2 = crlp_indv_ate)
  return(list(dindv_before = dindv_unm_mat, dindv_ate = dindv_ate))
}

dist_mom <- function(dist_mat){
  #dbin_mean_before <- mean(colMeans(dist_mat$dbin_before))
  dindv_mean_before <- mean(colMeans(dist_mat$dindv_before))

  #dbin_mean_att <- mean(dist_mat$dbin_att)
  #dbin_sd_att <- sd(dist_mat$dbin_att)
  #dbin_mean_ate <- mean(dist_mat$dbin_ate)
  #dbin_sd_ate <- sd(dist_mat$dbin_ate)
  dindv_mean_ate <- mean(dist_mat$dindv_ate)
  dindv_sd_ate <- sd(dist_mat$dindv_ate)

  return(list(dindv_before = dindv_mean_before,
              dindv_ate = dindv_mean_ate, dindv_ate_sd = dindv_sd_ate))
}

#estimate
extract_ests_match <- function(m){
  #original treated/ control units number
  crl_num <- dim(m$ATC_obj$match.matrix)[1]
  trt_num <- dim(m$ATT_obj$match.matrix)[1]
  summary <- c(m$ATE_est, crl_num, m$ATT_crl_num, trt_num, m$ATC_trt_num) %>%
    set_names(c("ATE_est", "crl_num_orig", "crl_num_matched", "trt_num_orig", "trt_num_matched"))%>%
    as.list() %>% as_tibble()
  summary
}

balance_zc <- function(m, df) {
  s_trt_ate <- var(df$z_c[m$index.treated]); s_crl_ate <- var(df$z_c[m$index.control])
  #the number of treated is approximately equaling to that of untreated
  std_ate <- sqrt((s_trt_ate + s_crl_ate)/2)
  trt_zc_ate <- df$z_c[m$index.treated]
  crl_zc_ate <- df$z_c[m$index.control]
  #ate
  sdiff_ate <- abs(mean(trt_zc_ate) - mean(crl_zc_ate))/std_ate
  diff_ate <- abs(mean(trt_zc_ate) - mean(crl_zc_ate))
  vr_ate <- s_trt_ate/s_crl_ate
  ks_ate <- ks.boot(trt_zc_ate, crl_zc_ate)
  ks_ate_stat <- ks_ate$ks$statistic; ks_ate_p <- ks_ate$ks$p.value
  return(list(sdiff_ate = sdiff_ate, diff_ate = diff_ate, vr_ate = vr_ate, ks_ate_stat = ks_ate_stat, ks_ate_p = ks_ate_p))
}

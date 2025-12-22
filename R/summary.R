#function that seperate the result and summerises mean, quantitles and empirical standard error of estimated coefficient
#summary_regmodel(NSlatcomb,phi_seq, phi_seq, 9, coef_gen_x=1)
summary_regmodel <- function(res,phi_seq_c, phi_seq_u, te, coef_gen_x) {
  stat_constant <- NULL
  for(c in seq_along(phi_seq_c)){
    dat_phic <- res %>% filter(phi_c == phi_seq_c[c])
    if (te == 0 & phi_seq_c[c] == 0) {
      phi_i_seq = setdiff(phi_seq_u,0)
    } else {
      phi_i_seq = phi_seq_u
    }
    for (i in seq_along(phi_i_seq)) {
      dat <- dat_phic %>% filter(phi_u == phi_i_seq[i], theta == te)
      mean_indv_est <- mean(dat$ATE_est); med_indv_est <- median(dat$ATE_est)
      indv_emp_se <- sd(dat$ATE_est);
      Sk2_indv_est <- 3*abs(mean_indv_est - med_indv_est)/indv_emp_se
      if (Sk2_indv_est > 0.3669) {
        cat(unique(dat$model),": the distribution of coef. estimates under phi_c=",phi_seq_c[c],
            "phi_i=", phi_i_seq[i], "theta=",te,"is significantly skewed with Sk2=", Sk2_indv_est, "\n")
      }
      indv_est25 <- quantile(dat$ATE_est, 0.25); indv_est75 <- quantile(dat$ATE_est, 0.75)

      mean_indv_se <- sqrt(mean((dat$ATE_se)^2)); med_indv_se <- sqrt(median((dat$ATE_se)^2))
      mean1 <- mean(dat$ATE_se); med1 <- median(dat$ATE_se)
      Sk2_indv <- 3*abs(mean1 - med1)/sd(dat$ATE_se)
      if (Sk2_indv > 0.3669) {
        cat(unique(dat$model),": the distribution of estimated se under phi_c=",phi_seq_c[c],"phi_i=",
            phi_i_seq[i], "theta=",te,"is significantly skewed with Sk2=", Sk2_indv, "\n")
        indv_est_se <- med_indv_se
      } else { indv_est_se <- mean_indv_se }

      #get rid of se_rob
      if(FALSE){
        mean_indv_se_rob <- mean(dat$coef_se_rob); med_indv_se_rob <- median(dat$coef_se_rob)
        Sk2_indv_rob <- 3*abs(mean_indv_se_rob - med_indv_se_rob)/sd(dat$coef_se_rob)
        if (Sk2_indv_rob > 0.3669) {
          cat(unique(dat$model),": the distribution of estimated robust se under phi_c=",phi_seq_c[c],"phi_i=",
              phi_i_seq[i], "theta=",te,"is significantly skewed with Sk2=", Sk2_indv_rob, "\n")
          indv_est_se_rob <- med_indv_se_rob
        } else { indv_est_se_rob <- mean_indv_se_rob }
      }

      indv_rMSE <- sqrt(mean((dat$ATE_est - coef_gen_x)^2))

      stat_constant <- rbind(stat_constant, c(phi_seq_c[c], phi_i_seq[i], te, mean_indv_est, indv_est25, med_indv_est, indv_est75,
                                              indv_est_se, mean_indv_se, med_indv_se, indv_emp_se, indv_rMSE))
    }
  }
  colnames(stat_constant) <- c("phi_c","phi_u","theta","ATE_est_mean","25%ATE_est","50%ATE_est","75%ATE_est",
                               "ATE_se_est","ATE_se_est_mean", "ATE_se_est_med","ATE_se_emp","rMSE")
  return(stat_constant)
}
#res <- nearest0.08PS_
#phi_seq_c = phi_seq_u = phi_seq; te =0; coef_gen_x = 1
#c=i=1
summary_combMR <- function(res, phi_seq_c, phi_seq_u, te, coef_gen_x, rob = FALSE) {
  stat_constant <-  NULL
  for(c in seq_along(phi_seq_c)){
    dat_phic <- res %>% filter(phi_c == phi_seq_c[c])
    if (te == 0 & phi_seq_c[c] == 0) {
      phi_i_seq = setdiff(phi_seq_u,0)
    } else {
      phi_i_seq = phi_seq_u
    }
    for (i in seq_along(phi_i_seq)) {
      dat <- dat_phic %>% filter(phi_u == phi_i_seq[i], theta == te)
      mean_est <- mean(dat$ATE_est); med_est <- median(dat$ATE_est)
      emp_se <- sd(dat$ATE_est)
      Sk2_est <- 3*abs(mean_est - med_est)/emp_se

      if (Sk2_est > 0.3669) {
        cat(unique(dat$save_file),"the distribution of estimated ATE under phi_c=",phi_seq_c[c],
            "phi_u=", phi_i_seq[i], "is significantly skewed with Sk2=", Sk2_est, "\n")
      }

      est_ate25 <- quantile(dat$ATE_est, 0.25); est_ate75 <- quantile(dat$ATE_est, 0.75)

      mean1 <- mean(dat$ATE_se); med1 <- median(dat$ATE_se)
      mean_indv_se <- sqrt(mean((dat$ATE_se)^2)); med_indv_se <- sqrt(median((dat$ATE_se)^2))
      Sk2_indv <- 3*abs(mean1 - med1)/sd(dat$ATE_se)
      if (Sk2_indv > 0.3669) {
        cat(unique(dat$save_file),": the distribution of estimated se under phi_c=",phi_seq_c[c],"phi_u=",
            phi_i_seq[i], "theta=",te,"is significantly skewed with Sk2=", Sk2_indv, "\n")
        indv_est_se <- med_indv_se
      } else { indv_est_se <- mean_indv_se }

      #WLS
      if(FALSE){
        mean_WLS_est <- mean(dat$WLS_ATE_est); med_WLS_est <- median(dat$WLS_ATE_est)
        emp_WLS_se <- sd(dat$WLS_ATE_est)
        Sk2_WLS_est <- 3*abs(mean_WLS_est - med_WLS_est)/emp_WLS_se

        if (Sk2_WLS_est > 0.3669) {
          cat(unique(dat$save_file), "the distribution of WLS ate estimates for M1comb under phi_c=",phi_seq_c[c],
              "phi_U=", phi_i_seq[i], "is significantly skewed with Sk2=", Sk2_WLS_est, "\n")
        }
        est_WLS_ate25 <- quantile(dat$WLS_ATE_est, 0.25); est_WLS_ate75 <- quantile(dat$WLS_ATE_est, 0.75)

        mean_WLS1 <- mean(dat$WLS_ATE_se); med_WLS1 <- median(dat$WLS_ATE_se)
        mean_WLS_se <- sqrt(mean((dat$WLS_ATE_se)^2)); med_WLS_se <- sqrt(median((dat$WLS_ATE_se)^2))
        Sk2_WLS_se <- 3*abs(mean_WLS1 - med_WLS1)/sd(dat$WLS_ATE_se)
        if (Sk2_WLS_se > 0.3669) {
          cat(unique(dat$save_file),": the distribution of WLS estimated se under phi_c=",phi_seq_c[c],"phi_U=",
              phi_i_seq[i], "theta=",te,"is significantly skewed with Sk2=", Sk2_WLS_se, "\n")
          WLS_est_se <- med_WLS_se
        } else { WLS_est_se <- mean_WLS_se }
      }

      rMSE <- sqrt(mean((dat$ATE_est - coef_gen_x)^2))
      #rMSE_WLS <- sqrt(mean((dat$WLS_ATE_est - coef_gen_x)^2))
      time <- difftime(dat$time1, dat$time0)
      #time <- time_WLSonly <- 0
      #time_WLSonly <- difftime(dat$time3, dat$time2)

      mean_cluster1 = mean(dat$coef_se_CL); med_cluster1 <- median(dat$coef_se_CL)
      mean_cluster_se <- sqrt(mean((dat$coef_se_CL)^2)); med_cluster_se <- sqrt(median((dat$coef_se_CL)^2))
      Sk2_cl_se <- 3*abs(mean_cluster1 - med_cluster1)/sd(dat$coef_se_CL)
      if (Sk2_cl_se > 0.3669) {
        cat(unique(dat$save_file),": the distribution of clustered se under phi_c=",phi_seq_c[c],"phi_U=",
            phi_i_seq[i], "theta=",te,"is significantly skewed with Sk2=", Sk2_cl_se, "\n")
        cluster_est_se <- med_cluster_se
      } else { cluster_est_se <- mean_cluster_se }

      #mean_cluster_se_pair = mean(dat$coef_se_CL_pair); med_cluster_se_pair <- median(dat$coef_se_CL_pair)
      #Sk2_cl_se_pair <- 3*abs(mean_cluster_se_pair - med_cluster_se_pair)/sd(dat$coef_se_CL_pair)
      #if (Sk2_cl_se_pair > 0.3669) {
       # cat(unique(dat$save_file),": the distribution of paired clustered se under phi_c=",phi_seq_c[c],"phi_U=",
        #    phi_i_seq[i], "theta=",te,"is significantly skewed with Sk2=", Sk2_cl_se_pair, "\n")
        #cluster_est_se_pair <- med_cluster_se_pair
      #} else { cluster_est_se_pair <- mean_cluster_se_pair }

      #mean_cluster_WLS1 = mean(dat$WLS_coef_se_CL); med_cluster_WLS1 <- median(dat$WLS_coef_se_CL)
      #mean_cluster_se_WLS <- sqrt(mean((dat$WLS_coef_se_CL)^2)); med_cluster_se_WLS <- sqrt(median((dat$WLS_coef_se_CL)^2))
      #Sk2_cl_se_WLS <- 3*abs(mean_cluster_WLS1 - med_cluster_WLS1)/sd(dat$WLS_coef_se_CL)
      #if (Sk2_cl_se_WLS > 0.3669) {
       # cat(unique(dat$save_file),": the distribution of clustered se of weighted version under phi_c=",phi_seq_c[c],"phi_U=",
        #    phi_i_seq[i], "theta=",te,"is significantly skewed with Sk2=", Sk2_cl_se_WLS, "\n")
        #cluster_est_se_WLS <- med_cluster_se_WLS
      #} else { cluster_est_se_WLS <- mean_cluster_se_WLS }

      res_stat <- c(phi_seq_c[c], phi_i_seq[i], te, mean_est, est_ate25, med_est, est_ate75, indv_est_se,
                    mean_indv_se, med_indv_se, emp_se, rMSE, median(dat$edf),  median(dat$sp), #mean_WLS_est,
                    #est_WLS_ate25, med_WLS_est, est_WLS_ate75, WLS_est_se, mean_WLS_se, med_WLS_se, emp_WLS_se,
                    #rMSE_WLS, median(dat$WLS_edf), median(dat$WLS_sp),
                    median(time), quantile(time, c(0.25,0.75)), #median(time_WLSonly), quantile(time_WLSonly, c(0.25,0.75)),
                    cluster_est_se, mean_cluster_se, med_cluster_se)
                    #cluster_est_se_pair, mean_cluster_se_pair, med_cluster_se_pair,
                    #cluster_est_se_WLS, mean_cluster_se_WLS, med_cluster_se_WLS)

      if (rob == TRUE) {
        mean_indv_se_rob <- mean(dat$coef_se_rob); med_indv_se_rob <- median(dat$coef_se_rob)
        Sk2_indv_rob <- 3*abs(mean_indv_se_rob - med_indv_se_rob)/sd(dat$coef_se_rob)
        if (Sk2_indv_rob > 0.3669) {
          cat(unique(dat$save_file),": the distribution of estimated robust se under phi_c=",phi_seq_c[c],"phi_u=",
              phi_i_seq[i], "theta=",te,"is significantly skewed with Sk2=", Sk2_indv_rob, "\n")
          indv_est_se_rob <- med_indv_se_rob
        } else { indv_est_se_rob <- mean_indv_se_rob }

        #mean_WLS_se_rob <- mean(dat$WLS_coef_se_rob); med_WLS_se_rob <- median(dat$WLS_coef_se_rob)
        #Sk2_WLS_rob <- 3*abs(mean_WLS_se_rob - med_WLS_se_rob)/sd(dat$WLS_coef_se_rob)
        #if (Sk2_WLS_rob > 0.3669) {
         # cat(unique(dat$save_file),": the distribution of WLS estimated robust se under phi_c=",phi_seq_c[c],"phi_U=",
          #    phi_i_seq[i], "theta=",te,"is significantly skewed with Sk2=", Sk2_WLS_rob, "\n")
          #WLS_est_se_rob <- med_WLS_se_rob
        #} else { WLS_est_se_rob <- mean_WLS_se_rob}
        res_stat <- c(res_stat, indv_est_se_rob) #WLS_est_se_rob)
      }


      stat_constant <- rbind(stat_constant, res_stat)

    }
  }
  stat_name <- c("phi_c","phi_u","theta",
                 "ATE_est_mean", "25%ATE_est", "50%ATE_est", "75%ATE_est", "ATE_se_est", "ATE_se_est_mean",
                 "ATE_se_est_med", "ATE_se_emp","rMSE", "edf_med", "sp_med", #"ATE_est_WLS_mean", "25%ATE_est_WLS",
                 #"50%ATE_est_WLS", "75%ATE_est_WLS", "ATE_se_WLS_est", "ATE_se_WLS_est_mean", "ATE_se_WLS_est_med",
                 #"ATE_se_WLS_emp", "rMSE_WLS", "edf_WLS_med", "sp_WLS_med",
                 "medtime","25%time","75%time",#"medtimeWLSonly","25%timeWLSonly","75%timeWLSonly",
                 "ATE_cluster_se","ATE_cluster_se_mean","ATE_cluster_se_med")
                 #"ATE_cluster_se_pair","ATE_cluster_se_pair_mean","ATE_cluster_se_pair_med",
                 #"ATE_cluster_WLS_se","ATE_cluster_WLS_se_mean","ATE_cluster_WLS_se_med")
  if(rob == TRUE) { stat_name <- c(stat_name,"ATE_se_rob")}#"ATE_se_WLS_rob") }

  colnames(stat_constant) <- stat_name

  return(stat_constant)
}

#summary_fixed(FDF5latcomb, phi_seq, phi_seq, 0, coef_gen_x=1)
summary_fixed <- function(res, phi_seq_c, phi_seq_u, te, coef_gen_x) {
  stat_constant <- NULL
  for(c in seq_along(phi_seq_c)){
    dat_phic <- res %>% filter(phi_c == phi_seq_c[c])
    if (te == 0 & phi_seq_c[c] == 0) {
      phi_i_seq = setdiff(phi_seq_u,0)
    } else {
      phi_i_seq = phi_seq_u
    }
    for (i in seq_along(phi_i_seq)) {
      dat <- dat_phic %>% filter(phi_u == phi_i_seq[i], theta == te)

      #indv
      mean_indv_est <- mean(dat$ATE_est); med_indv_est <- median(dat$ATE_est)
      indv_emp_se <- sd(dat$ATE_est)
      Sk2_indv_est <- 3*abs(mean_indv_est - med_indv_est)/indv_emp_se
      if (Sk2_indv_est > 0.3669) {
        cat(unique(dat$model), ": the distribution of coef. estimates under phi_c=",phi_seq_c[c],
            "phi_i=", phi_i_seq[i], "\ntheta=", te, "k_nots=",unique(dat$n_knots), "is significantly skewed with Sk2=", Sk2_indv_est, "\n")
      }
      indv_est25 <- quantile(dat$ATE_est, 0.25); indv_est75 <- quantile(dat$ATE_est, 0.75)

      mean_indv_se <- sqrt(mean((dat$ATE_se)^2)); med_indv_se <- sqrt(median((dat$ATE_se)^2))
      mean1 <- mean(dat$ATE_se); med1 <- median(dat$ATE_se)
      Sk2_indv <- 3*abs(mean1 - med1)/sd(dat$ATE_se)
      if (Sk2_indv > 0.3669) {
        cat(unique(dat$model), ": the distribution of estimated se under phi_c=",phi_seq_c[c],"phi_i=",
            phi_i_seq[i], "\ntheta=", te, "k_nots=", unique(dat$n_knots), "is significantly skewed with Sk2=", Sk2_indv, "\n")
        indv_est_se <- med_indv_se
      } else { indv_est_se <- mean_indv_se }

      #get rid of se_rob
      if(FALSE){
        mean_indv_se_rob <- mean(dat$coef_se_rob); med_indv_se_rob <- median(dat$coef_se_rob)
        Sk2_indv_rob <- 3*abs(mean_indv_se_rob - med_indv_se_rob)/sd(dat$coef_se_rob)
        if (Sk2_indv_rob > 0.3669) {
          cat(unique(dat$model),": the distribution of estimated robust se under phi_c=",phi_seq_c[c],"phi_i=",
              phi_i_seq[i], "theta=",te,"k_nots=", unique(dat$n_knots), "is significantly skewed with Sk2=", Sk2_indv_rob, "\n")
          indv_est_se_rob <- med_indv_se_rob
        } else { indv_est_se_rob <- mean_indv_se_rob }
      }

      indv_rMSE <- sqrt(mean((dat$ATE_est - coef_gen_x)^2))
      time = difftime(dat$time1, dat$time0)
      stat_constant <- rbind(stat_constant, c(phi_seq_c[c], phi_i_seq[i], te, unique(dat$n_knots), mean_indv_est, indv_est25, med_indv_est,
                                              indv_est75, indv_est_se, mean_indv_se, med_indv_se,indv_emp_se, indv_rMSE,
                                              median(time), quantile(time, c(0.25,0.75))))

    }


  }
  colnames(stat_constant) <- c("phi_c","phi_u", "theta", "k", "ATE_est_mean","25%ATE_est", "50%ATE_est", "75%ATE_est",
                           "ATE_se_est","ATE_se_est_mean","ATE_se_est_med","ATE_se_emp","rMSE", "medtime","25%time","75%time")
  return(stat_constant)
}

#res = EPS100; phi_seq_c = 0.15; phi_seq_u=0.6; te=3; modname = "est_df_x"
#coef_gen_x = 1; c=i=1
#summary_est(PS500lat,phi_seq, phi_seq, te=0, modname="est_df_y",1)
summary_est <- function(res,phi_seq_c,phi_seq_u, te, modname, coef_gen_x) {
  stat_constant <- NULL
  for(c in seq_along(phi_seq_c)){
    dat_phic <- res %>% filter(phi_c == phi_seq_c[c])
    if (te == 0 & phi_seq_c[c] == 0) {
      phi_i_seq = setdiff(phi_seq_u,0)
    } else {
      phi_i_seq = phi_seq_u
    }
    for (i in seq_along(phi_i_seq)) {
      dat <- dat_phic %>% filter(phi_u == phi_i_seq[i], theta == te)
      mean_indv_est <- mean(dat$ATE_est); med_indv_est <- median(dat$ATE_est)
      indv_emp_se <- sd(dat$ATE_est)
      Sk2_indv_est <- 3*abs(mean_indv_est - med_indv_est)/indv_emp_se
      if (Sk2_indv_est > 0.3669) {
        cat(unique(dat$model), ": the distribution of coef. estimates under phi_c=",phi_seq_c[c],
            "phi_i=", phi_i_seq[i], "\ntheta=", te, "k_nots=", unique(dat$n_knots), "is significantly skewed with Sk2=", Sk2_indv_est, "\n")
      }
      indv_est25 <- quantile(dat$ATE_est, 0.25); indv_est75 <- quantile(dat$ATE_est, 0.75)

      mean1 <- mean(dat$ATE_se); med1 <- median(dat$ATE_se)
      mean_indv_se <- sqrt(mean((dat$ATE_se)^2)); med_indv_se <- sqrt(median((dat$ATE_se)^2))
      Sk2_indv <- 3*abs(mean1 - med1)/sd(dat$ATE_se)
      if (Sk2_indv > 0.3669) {
        cat(unique(dat$model), ": the distribution of estimated se under phi_c=",phi_seq_c[c],"phi_i=",
            phi_i_seq[i], "\ntheta=", te, "k_nots=", unique(dat$n_knots), "is significantly skewed with Sk2=", Sk2_indv, "\n")
        indv_est_se <- med_indv_se
      } else { indv_est_se <- mean_indv_se }

      # now get rid of se_rob
      if(FALSE){
        mean_indv_se_rob <- mean(dat$coef_se_rob); med_indv_se_rob <- median(dat$coef_se_rob)
        Sk2_indv_rob <- 3*abs(mean_indv_se_rob - med_indv_se_rob)/sd(dat$coef_se_rob)
        if (Sk2_indv_rob > 0.3669) {
          cat(unique(dat$model),": the distribution of estimated robust se under phi_c=",phi_seq_c[c],"phi_i=",
              phi_i_seq[i], "theta=",te, "k_nots=", unique(dat$n_knots), "is significantly skewed with Sk2=", Sk2_indv_rob, "\n")
          indv_est_se_rob <- med_indv_se_rob
        } else { indv_est_se_rob <- mean_indv_se_rob }
      }

      indv_rMSE <- sqrt(mean((dat$ATE_est - coef_gen_x)^2))
      time = difftime(dat$time1, dat$time0)

      stat_constant <- rbind(stat_constant, c(phi_seq_c[c], phi_i_seq[i], te, unique(dat$n_knots), mean_indv_est, indv_est25, med_indv_est, indv_est75,
                                              indv_est_se, mean_indv_se, med_indv_se, indv_emp_se, indv_rMSE, median(dat$edf), median(dat$sp),
                                              median(time), quantile(time, c(0.25,0.75))))
    }

  }
  colnames(stat_constant) <- c("phi_c", "phi_u", "theta", "k", "ATE_est_mean","25%ATE_est", "50%ATE_est", "75%ATE_est",
                               "ATE_se_est","ATE_se_est_mean","ATE_se_est_med","ATE_se_emp","rMSE", "edf_med", "sp_med", "medtime","25%time","75%time")

  return(stat_constant)
}
#summary_Mmatch(res = res_list_constant[[1]],
  #             phi_seq_c = phi_seq, phi_seq_u = phi_seq,
 #              te=theta_seq[1],coef_gen_x=1)
#c=i=1
summary_Mmatch <- function(res, phi_seq_c, phi_seq_u, te, coef_gen_x){
  stat_constant <- NULL
  for (c in seq_along(phi_seq_c)) {
    dat_phic <- res %>% filter(phi_c == phi_seq_c[c])
    if (te == 0 & phi_seq_c[c]==0) {
      phi_i_seq = setdiff(phi_seq_u,0)
    } else {
      phi_i_seq = phi_seq_u
    }
    for (i in seq_along(phi_i_seq)) {
      dat <- dat_phic %>% filter(phi_u == phi_i_seq[i], theta == te)
      mean_indv_ate <- mean(dat$ATE_est);med_indv_ate <- median(dat$ATE_est)
      indv_ate_emp_se <- sd(dat$ATE_est)
      Sk2_indv_ate <- 3*abs(mean_indv_ate - med_indv_ate)/indv_ate_emp_se
      if (Sk2_indv_ate > 0.3669) {
        cat(unique(dat$save_file), ": the distribution of ate estimates under phi_c=",phi_seq_c[c],"phi_i=",phi_i_seq[i],
            "\ntheta=", te, "is significantly skewed with Sk2=",Sk2_indv_ate,"\n")
      }
      indv_ate25 <- quantile(dat$ATE_est, 0.25); indv_ate75 <- quantile(dat$ATE_est, 0.75)

      mean1 <- mean(dat$ATE_se); med1 <- median(dat$ATE_se)
      mean_indv_se <- sqrt(mean((dat$ATE_se)^2)); med_indv_se <- sqrt(median((dat$ATE_se)^2))
      Sk2_indv <- 3*abs(mean1 - med1)/sd(dat$ATE_se)
      if (Sk2_indv > 0.3669) {
        cat(unique(dat$save_file),": the distribution of estimated se under phi_c=",phi_seq_c[c],"phi_i=",
            phi_i_seq[i], "theta=",te,"is significantly skewed with Sk2=", Sk2_indv, "\n")
        indv_est_se <- med_indv_se
      } else { indv_est_se <- mean_indv_se }

      #mean_indv_se_rob <- mean(dat$coef_se_rob); med_indv_se_rob <- median(dat$coef_se_rob)
      #Sk2_indv_rob <- 3*abs(mean_indv_se_rob - med_indv_se_rob)/sd(dat$coef_se_rob)
      #if (Sk2_indv_rob > 0.3669) {
       # cat(unique(dat$save_file),": the distribution of robust estimator of se under phi_c=",phi_seq_c[c],"phi_i=",
        #    phi_i_seq[i], "theta=",te,"is significantly skewed with Sk2=", Sk2_indv_rob, "\n")
        #indv_est_se_rob <- med_indv_se_rob
      #} else { indv_est_se_rob <- mean_indv_se_rob }

      #WLS
      mean1 <- mean(dat$WLS_ATE_se); med1 <- median(dat$WLS_ATE_se)
      mean_WLS_se <- sqrt(mean((dat$WLS_ATE_se)^2)); med_WLS_se <- sqrt(median((dat$WLS_ATE_se)^2))
      Sk2_WLS_se <- 3*abs(mean1 - med1)/sd(dat$WLS_ATE_se)
      if (Sk2_WLS_se > 0.3669) {
        cat(unique(dat$save_file),": the distribution of WLS estimated se under phi_c=",phi_seq_c[c],"phi_i=",
            phi_i_seq[i], "theta=",te,"is significantly skewed with Sk2=", Sk2_WLS_se, "\n")
        WLS_est_se <- med_WLS_se
      } else { WLS_est_se <- mean_WLS_se }

      #mean_WLS_se_rob <- mean(dat$coef_se_rob_wls); med_WLS_se_rob <- median(dat$coef_se_rob_wls)
      #Sk2_WLS_rob <- 3*abs(mean_WLS_se_rob - med_WLS_se_rob)/sd(dat$coef_se_rob_wls)
      #if (Sk2_WLS_rob > 0.3669) {
       # cat(unique(dat$save_file),": the distribution of WLS estimated robust se under phi_c=",phi_seq_c[c],"phi_i=",
        #    phi_i_seq[i], "theta=",te,"is significantly skewed with Sk2=", Sk2_WLS_rob, "\n")
        #WLS_est_se_rob <- med_WLS_se_rob
      #} else { WLS_est_se_rob <- mean_WLS_se_rob }

      #MP est se
      MP_est_se <- 0
      if(mean(dat$MP_se) != 0){
        mean1 <- mean(dat$MP_se); med1 <- median(dat$MP_se)
        mean_MP_se <- sqrt(mean((dat$MP_se)^2)); med_MP_se <- sqrt(median((dat$MP_se)^2))
        Sk2_MP_se <- 3*abs(mean1 - med1)/sd(dat$MP_se)
        if (Sk2_MP_se > 0.3669) {
          cat(unique(dat$save_file),": the distribution of MP estimated se under phi_c=",phi_seq_c[c],"phi_i=",
              phi_i_seq[i], "theta=",te,"is significantly skewed with Sk2=", Sk2_MP_se, "\n")
          MP_est_se <- med_MP_se
        } else { MP_est_se <- mean_MP_se }
      } else {
        mean_MP_se <- sqrt(mean((dat$MP_se)^2))
        med_MP_se <- sqrt(median((dat$MP_se)^2))
      }

      #MP_new est se
      MP_est_se_new <- 0
      if(mean(dat$MP_se_new) != 0){
        mean1 <- mean(dat$MP_se_new); med1 <- median(dat$MP_se_new)
        mean_MP_se_new <- sqrt(mean((dat$MP_se_new)^2)); med_MP_se_new <- sqrt(median((dat$MP_se_new)^2))
        Sk2_MP_se <- 3*abs(mean1 - med1)/sd(dat$MP_se_new)
        if (Sk2_MP_se > 0.3669) {
          cat(unique(dat$save_file),": the distribution of MP new estimated se under phi_c=",phi_seq_c[c],"phi_i=",
              phi_i_seq[i], "theta=",te,"is significantly skewed with Sk2=", Sk2_MP_se, "\n")
          MP_est_se_new <- med_MP_se_new
        } else { MP_est_se_new <- mean_MP_se_new }
      } else {
        mean_MP_se_new <- sqrt(mean((dat$MP_se_new)^2))
        med_MP_se_new <- sqrt(median((dat$MP_se_new)^2))
      }

      #Two-sample est se
      TS_est_se <- 0
      if(mean(dat$TS_se) != 0){
        mean1 <- mean(dat$TS_se); med1 <- median(dat$TS_se)
        mean_TS_se <- sqrt(mean((dat$TS_se)^2)); med_TS_se <- sqrt(median((dat$TS_se)^2))
        Sk2_TS_se <- 3*abs(mean1-med1)/sd(dat$TS_se)
        if (Sk2_TS_se > 0.3669) {
          cat(unique(dat$model),": the distribution of TS estimated se under phi_c=",phi_seq_c[c],"phi_i=",
              phi_i_seq[i], "theta=",te,"is significantly skewed with Sk2=", Sk2_TS_se, "\n")
          TS_est_se <- med_TS_se
        } else { TS_est_se <- mean_TS_se }
      } else {
        mean_TS_se <- sqrt(mean((dat$TS_se)^2)); med_TS_se <- sqrt(median((dat$TS_se)^2))
      }

      indv_ate_rMSE <- sqrt(mean((dat$ATE_est - coef_gen_x)^2))
      time = difftime(dat$time1, dat$time0)
      #indv_ate_dm <- mean(dat$indv.ate_dm); indv_ate_dse <- mean(dat$indv.ate_dse)
      #indv_ate_sdiff <- mean(dat$sdiff_ate); indv_ate_diff <- mean(dat$diff_ate); indv_ate_vr <- mean(dat$vr_ate)
      #indv_ate_ks_stat <- mean(dat$ks_ate_stat); indv_ate_ks_p <- mean(dat$ks_ate_p)

      #summary
      #stat_bin_ate <- rbind(stat_bin_ate, c(phi_seq_c[c], phi_i_seq[i], te, M[j], mean_bin_ate, bin_ate25, med_bin_ate, bin_ate75,
      #bin_ate_est_se, bin_ate_emp_se, bin_ate_rMSE, bin_ate_dm, bin_ate_dse, bin_ate_sdiff,
      #bin_ate_vr, bin_ate_ks_stat, bin_ate_ks_p,
      #mean(dat$bin_ate_n.treat), mean(dat$bin_ate_n.drop)))
      #stat_bin_att <- rbind(stat_bin_att, c(phi_seq_c[c], phi_i_seq[i], te, M[j], mean_bin_att, bin_att25, med_bin_att, bin_att75,
      #bin_att_est_se, bin_att_emp_se, bin_att_rMSE, bin_att_dm, bin_att_dse, bin_att_sdiff,
      #bin_att_vr, bin_att_ks_stat, bin_att_ks_p,
      #mean(dat$bin_att_n.treat), mean(dat$bin_att_n.drop)))

      stat_constant <- rbind(stat_constant, c(phi_seq_c[c], phi_i_seq[i], te, mean_indv_ate, indv_ate25, med_indv_ate, indv_ate75,
                                              indv_est_se, mean_indv_se, med_indv_se, indv_ate_emp_se, indv_ate_rMSE, WLS_est_se, mean_WLS_se, med_WLS_se,
                                              MP_est_se, mean_MP_se, med_MP_se, MP_est_se_new, mean_MP_se_new, med_MP_se_new,
                                              TS_est_se, mean_TS_se, med_TS_se, median(time), quantile(time, c(0.25,0.75))))
                                              #median(time), quantile(time, c(0.25,0.75)))) #indv_ate_dm, indv_ate_dse, indv_ate_sdiff,
                                              #indv_ate_diff, indv_ate_vr, indv_ate_ks_stat, indv_ate_ks_p,
                                              #mean(dat$crl_num_orig),mean(dat$crl_num_matched),
                                              #mean(dat$trt_num_orig),mean(dat$trt_num_matched)))
    }

  }

  colnames(stat_constant) <- c("phi_c", "phi_u", "theta", "ATE_est_mean","25%ATE_est", "50%ATE_est", "75%ATE_est",
                               "ATE_se_est","ATE_se_est_mean","ATE_se_est_med","ATE_se_emp","rMSE", "ATE_se_WLS_est","ATE_se_WLS_est_mean","ATE_se_WLS_est_med",
                              "ATE_se_MP","ATE_se_MP_mean","ATE_se_MP_med","ATE_se_MP_new","ATE_se_MP_new_mean","ATE_se_MP_new_med",
                              "ATE_se_TS","ATE_se_TS_mean","ATE_se_TS_med","medtime","25%time","75%time")

  return(stat_constant)
}

#phi_seq_c=phi_seq_i=c(0.04,0.6);te=1;n_datasets=100;coef_gen_x=1
#c=i=1
summary_IPW <- function(res, phi_seq_c, phi_seq_i, te, coef_gen_x){
  stat_ate_IPW <- NULL
  stat_ate_IPWRA <- NULL
  for (c in seq_along(phi_seq_c)) {
    dat_phic <- res %>% filter(phi_c == phi_seq_c[c])
    if (te == 0 & phi_seq_c[c]==0) {
      #remove 0 from list for phi_i
      phi_i_seq = setdiff(phi_seq_i,0)
    } else {
      phi_i_seq = phi_seq_i
    }
    for (i in seq_along(phi_i_seq)) {
      dat <- dat_phic %>% filter(phi_i == phi_i_seq[i], theta == te)
      #indv_ipw
      mean_est_IPW <- mean(dat$ATE_IPW);med_est_IPW <- median(dat$ATE_IPW)
      emp_se_IPW <- sd(dat$ATE_IPW)
      Sk2_indv_ate <- 3*abs(mean_est_IPW - med_est_IPW)/emp_se_IPW
      if (Sk2_indv_ate > 0.3669) {
        cat(unique(dat$model), ": the distribution of IPW ate estimates under phi_c=",phi_seq_c[c],"phi_i=",phi_i_seq[i],
            "\ntheta=", te, "is significantly skewed with Sk2=",Sk2_indv_ate,"\n")
      }
      est_se_IPW <- 0

      indv_ate_IPW25 <- quantile(dat$ATE_IPW, 0.25); indv_ate_IPW75 <- quantile(dat$ATE_IPW, 0.75)
      indv_ate_IPW_rMSE <- sqrt(mean((dat$ATE_IPW - coef_gen_x)^2))

      #indv_ipwra
      mean_est_IPWRA <- mean(dat$ATE_IPWRA);med_est_IPWRA <- median(dat$ATE_IPWRA)
      emp_se_IPWRA <- sd(dat$ATE_IPWRA)
      Sk2_indv_ate <- 3*abs(mean_est_IPWRA - med_est_IPWRA)/emp_se_IPWRA
      if (Sk2_indv_ate > 0.3669) {
        cat(unique(dat$model), ": the distribution of IPWRA ate estimates under phi_c=",phi_seq_c[c],"phi_i=",phi_i_seq[i],
            "\ntheta=", te, "is significantly skewed with Sk2=",Sk2_indv_ate,"\n")
      }
      est_se_IPWRA <- 0

      indv_ate_IPWRA25 <- quantile(dat$ATE_IPWRA, 0.25); indv_ate_IPWRA75 <- quantile(dat$ATE_IPWRA, 0.75)
      indv_ate_IPWRA_rMSE <- sqrt(mean((dat$ATE_IPWRA - coef_gen_x)^2))
      time_IPW <- difftime(dat$time2, dat$time0)
      time_IPWRA <- difftime(dat$time1, dat$time0) + difftime(dat$time3, dat$time2)

      stat_ate_IPW <- rbind(stat_ate_IPW, c(phi_seq_c[c], phi_i_seq[i], te, mean_est_IPW, indv_ate_IPW25, med_est_IPW, indv_ate_IPW75,
                                            est_se_IPW, emp_se_IPW, indv_ate_IPW_rMSE,
                                            mean(dat$`5%weights`), mean(dat$`95%weights`),mean(dat$mean_weight), mean(dat$sd_weight),
                                            median(time_IPW), quantile(time_IPW, c(0.25,0.75))))
      stat_ate_IPWRA <- rbind(stat_ate_IPWRA, c(phi_seq_c[c], phi_i_seq[i], te, mean_est_IPWRA, indv_ate_IPWRA25, med_est_IPWRA, indv_ate_IPWRA75,
                                                est_se_IPWRA, emp_se_IPWRA, indv_ate_IPWRA_rMSE,
                                                mean(dat$`5%weights`), mean(dat$`95%weights`),mean(dat$mean_weight), mean(dat$sd_weight),
                                                median(time_IPWRA), quantile(time_IPWRA, c(0.25,0.75))))
    }

  }

  colnames(stat_ate_IPW) <- c("phi_c", "phi_i", "theta", "mean_est_IPW", "indv_ate_IPW25", "med_est_IPW", "indv_ate_IPW75",
                               "est_se_IPW", "emp_se_IPW", "indv_ate_IPW_rMSE",
                               "5%weights", "95%weights","mean_weight", "sd_weight","medtime_IPW","25time_IPW","75time_IPW")
  colnames(stat_ate_IPWRA) <- c("phi_c", "phi_i", "theta", "mean_est_IPWRA", "indv_ate_IPWRA25", "med_est_IPWRA", "indv_ate_IPWRA75",
                               "est_se_IPWRA", "emp_se_IPWRA", "indv_ate_IPWRA_rMSE",
                               "5%weights", "95%weights","mean_weight", "sd_weight","medtime_IPWRA","25time_IPWRA","75time_IPWRA")

  return(list(ATE_IPW = stat_ate_IPW, ATE_IPWRA = stat_ate_IPWRA))

}

summary_Mmatch_zc <- function(res, phi_seq_c, phi_seq_i, te, M, n_datasets=100, coef_gen_x=1){
  stat_bin_ate <- NULL
  stat_bin_att <- NULL
  stat_indv_ate <- NULL
  stat_indv_att <- NULL
  for (c in seq_along(phi_seq_c)) {
    dat_phic <- res %>% filter(phi_c == phi_seq_c[c])
    if (te == 0 & phi_seq_c[c]==0) {
      #remove 0 from list for phi_i
      phi_i_seq = setdiff(phi_seq_i,0)
    } else {
      phi_i_seq = phi_seq_i
    }
    for (i in seq_along(phi_i_seq)) {
      dat_phici <- dat_phic %>% filter(phi_i == phi_i_seq[i], theta == te)
      for (j in seq_along(M)) {
        dat <- dat_phici %>% filter(M == M[j])
        #bin_ate
        mean_bin_ate <- mean(dat$bin_ate_est);med_bin_ate <- median(dat$bin_ate_est)
        bin_ate_emp_se <- sd(dat$bin_ate_est)
        Sk2_bin_ate <- 3*abs(mean_bin_ate - med_bin_ate)/bin_ate_emp_se
        if (Sk2_bin_ate > 0.3669) {
          cat(unique(dat$model), ": the distribution of bin ate estimates under phi_c=",phi_seq_c[c],"phi_i=",phi_i_seq[i],
              "\ntheta=", te, "M=", M[j], "is significantly skewed with Sk2=",Sk2_bin_ate,"\n")
        } else {cat(c(c,i))}
        bin_ate25 <- quantile(dat$bin_ate_est, 0.25); bin_ate75 <- quantile(dat$bin_ate_est, 0.75)
        bin_ate_se_mean <- mean(dat$bin_ate_se); bin_ate_se_med <- median(dat$bin_ate_se)
        Sk2_bin_ate_se <- 3*abs(bin_ate_se_mean - bin_ate_se_med)/sd(dat$bin_ate_se)
        bin_ate_est_se <- ifelse(Sk2_bin_ate_se <= 0.3669, bin_ate_se_mean, bin_ate_se_med)
        bin_ate_rMSE <- sqrt(sum((dat$bin_ate_est - coef_gen_x)^2)/n_datasets)
        bin_ate_dm <- mean(dat$bin.ate_dm); bin_ate_dse <- mean(dat$bin.ate_dse)

        #bin_att
        mean_bin_att <- mean(dat$bin_att_est);med_bin_att <- median(dat$bin_att_est)
        bin_att_emp_se <- sd(dat$bin_att_est)
        Sk2_bin_att <- 3*abs(mean_bin_att - med_bin_att)/bin_att_emp_se
        if (Sk2_bin_att > 0.3669) {
          cat(unique(dat$model), ": the distribution of bin att estimates under phi_c=",phi_seq_c[c],"phi_i=",phi_i_seq[i],
              "\ntheta=", te, "M=", M[j], "is significantly skewed with Sk2=",Sk2_bin_att,"\n")
        }
        bin_att25 <- quantile(dat$bin_att_est, 0.25); bin_att75 <- quantile(dat$bin_att_est, 0.75)
        bin_att_se_mean <- mean(dat$bin_att_se); bin_att_se_med <- median(dat$bin_att_se)
        Sk2_bin_att_se <- 3*abs(bin_att_se_mean - bin_att_se_med)/sd(dat$bin_att_se)
        bin_att_est_se <- ifelse(Sk2_bin_att_se <= 0.3669, bin_att_se_mean, bin_att_se_med)
        bin_att_rMSE <- sqrt(sum((dat$bin_att_est - coef_gen_x)^2)/n_datasets)
        bin_att_dm <- mean(dat$bin.att_dm); bin_att_dse <- mean(dat$bin.att_dse)

        #indv_ate
        mean_indv_ate <- mean(dat$indv_ate_est);med_indv_ate <- median(dat$indv_ate_est)
        indv_ate_emp_se <- sd(dat$indv_ate_est)
        Sk2_indv_ate <- 3*abs(mean_indv_ate - med_indv_ate)/indv_ate_emp_se
        if (Sk2_indv_ate > 0.3669) {
          cat(unique(dat$model), ": the distribution of indv ate estimates under phi_c=",phi_seq_c[c],"phi_i=",phi_i_seq[i],
              "\ntheta=", te, "M=", M[j], "is significantly skewed with Sk2=",Sk2_indv_ate,"\n")
        }
        indv_ate25 <- quantile(dat$indv_ate_est, 0.25); indv_ate75 <- quantile(dat$indv_ate_est, 0.75)
        indv_ate_se_mean <- mean(dat$indv_ate_se); indv_ate_se_med <- median(dat$indv_ate_se)
        Sk2_indv_ate_se <- 3*abs(indv_ate_se_mean - indv_ate_se_med)/sd(dat$indv_ate_se)
        indv_ate_est_se <- ifelse(Sk2_indv_ate_se <= 0.3669, indv_ate_se_mean, indv_ate_se_med)
        indv_ate_rMSE <- sqrt(sum((dat$indv_ate_est - coef_gen_x)^2)/n_datasets)
        indv_ate_dm <- mean(dat$indv.ate_dm); indv_ate_dse <- mean(dat$indv.ate_dse)

        #indv_att
        mean_indv_att <- mean(dat$indv_att_est);med_indv_att <- median(dat$indv_att_est)
        indv_att_emp_se <- sd(dat$indv_att_est)
        Sk2_indv_att <- 3*abs(mean_indv_att - med_indv_att)/indv_att_emp_se
        if (Sk2_indv_att > 0.3669) {
          cat(unique(dat$model), ": the distribution of indv att estimates under phi_c=",phi_seq_c[c],"phi_i=",phi_i_seq[i],
              "\ntheta=", te, "M=", M[j], "is significantly skewed with Sk2=",Sk2_indv_att,"\n")
        }
        indv_att25 <- quantile(dat$indv_att_est, 0.25); indv_att75 <- quantile(dat$indv_att_est, 0.75)
        indv_att_se_mean <- mean(dat$indv_att_se); indv_att_se_med <- median(dat$indv_att_se)
        Sk2_indv_att_se <- 3*abs(indv_att_se_mean - indv_att_se_med)/sd(dat$indv_att_se)
        indv_att_est_se <- ifelse(Sk2_indv_att_se <= 0.3669, indv_att_se_mean, indv_att_se_med)
        indv_att_rMSE <- sqrt(sum((dat$indv_att_est - coef_gen_x)^2)/n_datasets)
        indv_att_dm <- mean(dat$indv.att_dm); indv_att_dse <- mean(dat$indv.att_dse)

        #summary
        stat_bin_ate <- rbind(stat_bin_ate, c(phi_seq_c[c], phi_i_seq[i], te, M[j], mean_bin_ate, bin_ate25, med_bin_ate, bin_ate75,
                                              bin_ate_est_se, bin_ate_emp_se, bin_ate_rMSE, bin_ate_dm, bin_ate_dse,
                                              mean(dat$bin_ate_n.treat), mean(dat$bin_ate_n.drop)))
        stat_bin_att <- rbind(stat_bin_att, c(phi_seq_c[c], phi_i_seq[i], te, M[j], mean_bin_att, bin_att25, med_bin_att, bin_att75,
                                              bin_att_est_se, bin_att_emp_se, bin_att_rMSE, bin_att_dm, bin_att_dse,
                                              mean(dat$bin_att_n.treat), mean(dat$bin_att_n.drop)))

        stat_indv_ate <- rbind(stat_indv_ate, c(phi_seq_c[c], phi_i_seq[i], te, M[j], mean_indv_ate, indv_ate25, med_indv_ate, indv_ate75,
                                                indv_ate_est_se, indv_ate_emp_se, indv_ate_rMSE, indv_ate_dm, indv_ate_dse,
                                                mean(dat$indv_ate_n.treat), mean(dat$indv_ate_n.drop)))
        stat_indv_att <- rbind(stat_indv_att, c(phi_seq_c[c], phi_i_seq[i], te, M[j], mean_indv_att, indv_att25, med_indv_att, indv_att75,
                                                indv_att_est_se, indv_att_emp_se, indv_att_rMSE, indv_att_dm, indv_att_dse,
                                                mean(dat$indv_att_n.treat), mean(dat$indv_att_n.drop)))
      }
    }

  }
  #set col and row names
  colnames(stat_bin_ate) <- c("phi_c", "phi_i", "theta", "M", "bin_ate_mean","25%bin_ate", "50%bin_ate", "75%bin_ate",
                              "bin_ate_est_se","bin_ate_emp_se","bin_ate_rMSE", "bin_ate_dm", "bin_ate_dse",
                              "bin_ate_n.treat", "bin_ate_n.drop")
  colnames(stat_bin_att) <- c("phi_c", "phi_i", "theta", "M", "bin_att_mean","25%bin_att", "50%bin_att", "75%bin_att",
                              "bin_att_est_se","bin_att_emp_se","bin_att_rMSE", "bin_att_dm", "bin_att_dse",
                              "bin_att_n.treat", "bin_att_n.drop")

  colnames(stat_indv_ate) <- c("phi_c", "phi_i", "theta", "M", "indv_ate_mean","25%indv_ate", "50%indv_ate", "75%indv_ate",
                               "indv_ate_est_se","indv_ate_emp_se","indv_ate_rMSE", "indv_ate_dm", "indv_ate_dse",
                               "indv_ate_n.treat", "indv_ate_n.drop")
  colnames(stat_indv_att) <- c("phi_c", "phi_i", "theta", "M", "indv_att_mean","25%indv_att", "50%indv_att", "75%indv_att",
                               "indv_att_est_se","indv_att_emp_se","indv_att_rMSE", "indv_att_dm", "indv_att_dse",
                               "indv_att_n.treat", "indv_att_n.drop")

  return(list(bin_ate = stat_bin_ate, bin_att = stat_bin_att,
              indv_ate = stat_indv_ate, indv_att = stat_indv_att))

}


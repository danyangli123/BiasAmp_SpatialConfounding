#function that organize stat of all
#require("rlist")
stat_all <- function(res_list, list, phi_seq_c, phi_seq_u, theta_seq, ATE, rob_need=FALSE) {
  #cat("go stat_all")
  stat <- list();
  for (j in seq_along(res_list)) {
    cat("j=",j,"\n")
    statj <- list()
    for (i in seq_along(theta_seq)) {
      if (list[j]=="true_Model"|list[j]=="Naive"|list[j]=="partial_true") {
        statij <- summary_regmodel(res_list[[j]], phi_seq_c, phi_seq_u, te=theta_seq[i], ATE)
        cat("true/Naive_done","\n")
      } else if (list[j]=="fdf") {
        statij <- summary_fixed(res_list[[j]], phi_seq_c, phi_seq_u, te=theta_seq[i], ATE)
        cat("fdf_done","\n")
      } else if (list[j]=="PS") {
        statij <- summary_est(res_list[[j]], phi_seq_c, phi_seq_u, te = theta_seq[i], modname = "est_df_y", ATE)
        cat("PS_done","\n")
      } else if (list[j]=="EPS") {
        statij <- summary_est(res_list[[j]], phi_seq_c, phi_seq_u, te=theta_seq[i], modname = "est_df_x", ATE)
        cat("EPS_done","\n")
      } else if (list[j]=="M1") {
        statij <- summary_Mmatch(res_list[[j]], phi_seq_c, phi_seq_u, te=theta_seq[i],ATE)
        statij <- as.data.frame(statij) %>% dplyr::select(-contains("WLS"))
        statij <- as.matrix(statij)
        cat("M1_done","\n")
      } else if (list[j]=="M1_wls") {
        statij <- summary_Mmatch(res_list[[j]], phi_seq_c, phi_seq_u, te=theta_seq[i], ATE)
        statij <- cbind(statij[,1:3],as.data.frame(statij) %>% dplyr::select(contains("WLS"), "ATE_se_emp"))
        colnames(statij) <- c("phi_c","phi_u","theta",
                              "ATE_se_est", "ATE_se_rob", "ATE_se_emp")
        statij <- as.matrix(statij)
        cat("M1_wls_done","\n")
      } else if (list[j]=="M1_weight") {
        statij <- summary_combMR(res_list[[j]], phi_seq_c, phi_seq_u, te=theta_seq[i], ATE, rob_need)
        statij <- cbind(statij[,1],statij[,2],statij[,3],
                        as.data.frame(statij) %>% dplyr::select(contains("WLS")))
        colnames(statij) <- c("phi_c","phi_u","theta",
                              "ATE_est_mean", "25%ATE_est", "50%ATE_est", "75%ATE_est", "ATE_se_est",
                              "ATE_se_est_mean","ATE_se_est_med","ATE_se_emp",
                              "rMSE", "edf_med", "sp_med", "medtime", "25%time", "75%time",
                              "ATE_cluster_se","ATE_cluster_se_mean","ATE_cluster_se_med")
        if(rob_need){
          colnames(statij) <- c(colnames(statij),"ATE_se_rob")
        }
        statij <- as.matrix(statij)
        cat("M1PS_weight_done with",theta_seq[i],"\n")
      } else if (list[j]=="M1_comb") {
        statij <- summary_combMR(res_list[[j]], phi_seq_c, phi_seq_u, te=theta_seq[i], ATE, rob_need)
        statij <- as.data.frame(statij) %>% dplyr::select(-contains("WLS"))
        statij <- as.matrix(statij)
        cat("M1PS_done with",theta_seq[i],"\n")
      } else if (list[j]=="IPW") {
        statij <- summary_IPW(res_list[[j]],phi_seq_c, phi_seq_i,te=theta_seq[i], coef_gen_x=ATE)
      } else {
        print("not available")
      }
      statj <- list.append(statj, statij)
    }
    stat <- list.append(stat, statj)
  }
  return(stat)
}

stat_all_ze <- function(res_list, list, phi_seq_c, phi_seq_e, phi_seq_u, theta_seq, ATE, rob_need=FALSE) {
  #cat("go stat_all")
  stat <- list();
  for (j in seq_along(res_list)) {
    cat("j=",j,"\n")
    statj <- list()
    for (i in seq_along(theta_seq)) {
      if (list[j]=="true_Model"|list[j]=="Naive") {
        statij <- summary_regmodel_ze(res_list[[j]], phi_seq_c, phi_seq_e, phi_seq_u, te=theta_seq[i], ATE)
        cat("true/Naive_done","te=",theta_seq[i],"\n")
      } else if (list[j]=="fdf") {
        statij <- summary_fixed_ze(res_list[[j]], phi_seq_c, phi_seq_e, phi_seq_u, te=theta_seq[i], ATE)
        cat("fdf_done","\n")
      } else {
        print("not available")
      }
      statj <- list.append(statj, statij)
    }
    stat <- list.append(stat, statj)
  }
  return(stat)
}

stat_all_betac <- function(res_list, list, phi_seq_c, phi_seq_u, theta_seq, betac_list, ATE, rob_need=FALSE) {
  #cat("go stat_all")
  stat <- list();
  for (j in seq_along(res_list)) {
    cat("j=",j,"\n")
    statj <- list()
    for (i in seq_along(theta_seq)) {
      if (list[j]=="true_Model"|list[j]=="Naive") {
        statij <- summary_regmodel_betac(res_list[[j]], phi_seq_c, phi_seq_u, te=theta_seq[i], betac_list, ATE)
        cat("true/Naive_done","\n")
      } else if (list[j]=="fdf") {
        statij <- summary_fixed_betac(res_list[[j]], phi_seq_c, phi_seq_u, te=theta_seq[i], betac_list, ATE)
        cat("fdf_done","\n")
      } else {
        print("not available")
      }
      statj <- list.append(statj, statij)
    }
    stat <- list.append(stat, statj)
  }
  return(stat)
}

biasRMSE <- function(stat, model_name, phi_c, phi_u, theta){
  res <- array(0, dim=c(length(theta), length(phi_c), length(phi_u), length(model_name), 4))
  for (t in seq_along(theta)) {
    for (c in seq_along(phi_c)) {
      for (u in seq_along(phi_u)) {
        for (m in seq_along(model_name)) {
          #med 1
          res[t,c,u,m,1] <- stat[[m]][[t]][(c-1)*length(phi_u)+u,"ATE_est_mean"]
          #lb 2
          res[t,c,u,m,2] <- stat[[m]][[t]][(c-1)*length(phi_u)+u,"25%ATE_est"]
          #ub 3
          res[t,c,u,m,3] <- stat[[m]][[t]][(c-1)*length(phi_u)+u,"75%ATE_est"]
          #rMSE 4
          res[t,c,u,m,4] <- stat[[m]][[t]][(c-1)*length(phi_u)+u,"rMSE"]
        }
      }
    }
  }

  size = prod(dim(res)[1:4])
  tesq <- phi_c_seq <- phi_u_seq <- model_seq <- meanest <- lb <- ub <- rMSE <- rep(0,size)
  l1 <- dim(res)[4]; l2 <- l1*dim(res)[3]; l3 <- l2*dim(res)[2]
  for (t in seq_along(theta)) {
    for (c in seq_along(phi_c)) {
      for (u in seq_along(phi_u)) {
        for (m in seq_along(model_name)) {
          tesq[(t-1)*l3+(c-1)*l2+(u-1)*l1+m] = theta[t]
          phi_c_seq[(t-1)*l3+(c-1)*l2+(u-1)*l1+m] = phi_c[c]
          phi_u_seq[(t-1)*l3+(c-1)*l2+(u-1)*l1+m] = phi_u[u]
          model_seq[(t-1)*l3+(c-1)*l2+(u-1)*l1+m] = model_name[m]
          meanest[(t-1)*l3+(c-1)*l2+(u-1)*l1+m] = res[t,c,u,m,1]
          lb[(t-1)*l3+(c-1)*l2+(u-1)*l1+m] = res[t,c,u,m,2]
          ub[(t-1)*l3+(c-1)*l2+(u-1)*l1+m] = res[t,c,u,m,3]
          rMSE[(t-1)*l3+(c-1)*l2+(u-1)*l1+m] = res[t,c,u,m,4]
        }
      }
    }
  }
  result <- data.frame(theta = factor(tesq), phic = factor(phi_c_seq,levels=phi_c),
                       phiu = factor(phi_u_seq,levels=phi_u),
                       model = factor(model_seq, levels = model_name),
                       mean = meanest, lb = lb, ub = ub, rMSE = rMSE)
  if("NS" %in% model_name){
    result$NS = rep(result[seq(which(model_name == "NS"),nrow(result),length(model_name)),"mean"],
                       each=length(model_name))
  }

  return(result)
}

biasRMSE_ze <- function(stat, model_name, phi_c, phi_u, phi_e, theta){
  res <- array(0, dim=c(length(theta), length(phi_c), length(phi_u), length(phi_e), length(model_name), 4))
  for (t in seq_along(theta)) {
    for (c in seq_along(phi_c)) {
      for (e in seq_along(phi_e)) {
        for (u in seq_along(phi_u)) {
          for (m in seq_along(model_name)) {
            #med 1
            res[t,c,u,e,m,1] <- stat[[m]][[t]][(c-1)*length(phi_e)*length(phi_u)+(u-1)*length(phi_e)+e,"ATE_est_mean"]
            #lb 2
            res[t,c,u,e,m,2] <- stat[[m]][[t]][(c-1)*length(phi_e)*length(phi_u)+(u-1)*length(phi_e)+e,"25%ATE_est"]
            #ub 3
            res[t,c,u,e,m,3] <- stat[[m]][[t]][(c-1)*length(phi_e)*length(phi_u)+(u-1)*length(phi_e)+e,"75%ATE_est"]
            #rMSE 4
            res[t,c,u,e,m,4] <- stat[[m]][[t]][(c-1)*length(phi_e)*length(phi_u)+(u-1)*length(phi_e)+e,"rMSE"]
          }
        }
      }
    }
  }

  size = prod(dim(res)[1:5])
  tesq <- phi_c_seq <- phi_e_seq <- phi_u_seq <- model_seq <- meanest <- lb <- ub <- rMSE <- rep(0,size)
  l1 <- dim(res)[5]; l2 <- l1*dim(res)[4]; l3 <- l2*dim(res)[3]; l4 <- l3*dim(res)[2]
  for (t in seq_along(theta)) {
    for (c in seq_along(phi_c)) {
      for (e in seq_along(phi_e)) {
        for (u in seq_along(phi_u)) {
          for (m in seq_along(model_name)) {
            tesq[(t-1)*l4+(c-1)*l3+(u-1)*l2+(e-1)*l1+m] = theta[t]
            phi_c_seq[(t-1)*l4+(c-1)*l3+(u-1)*l2+(e-1)*l1+m] = phi_c[c]
            phi_e_seq[(t-1)*l4+(c-1)*l3+(u-1)*l2+(e-1)*l1+m] = phi_e[e]
            phi_u_seq[(t-1)*l4+(c-1)*l3+(u-1)*l2+(e-1)*l1+m] = phi_u[u]
            model_seq[(t-1)*l4+(c-1)*l3+(u-1)*l2+(e-1)*l1+m] = model_name[m]
            meanest[(t-1)*l4+(c-1)*l3+(u-1)*l2+(e-1)*l1+m] = res[t,c,u,e,m,1]
            lb[(t-1)*l4+(c-1)*l3+(u-1)*l2+(e-1)*l1+m] = res[t,c,u,e,m,2]
            ub[(t-1)*l4+(c-1)*l3+(u-1)*l2+(e-1)*l1+m] = res[t,c,u,e,m,3]
            rMSE[(t-1)*l4+(c-1)*l3+(u-1)*l2+(e-1)*l1+m] = res[t,c,u,e,m,4]
          }
        }
      }
    }
  }
  result <- data.frame(theta = factor(tesq), phic = factor(phi_c_seq,levels=phi_c),
                       phiu = factor(phi_u_seq,levels=phi_u),
                       phie = factor(phi_e_seq,levels=phi_e),
                       model = factor(model_seq, levels = model_name),
                       mean = meanest, lb = lb, ub = ub, rMSE = rMSE)
  if("Naive" %in% model_name){
    result$Naive = rep(result[seq(which(model_name == "Naive"),nrow(result),length(model_name)),"mean"],
                       each=length(model_name))
  }

  return(result)
}


#i=t=1
#est_theta_se(stat, se_name = names(se_included), model_name,phi_c = phi_seq, phi_u = phi_seq, theta = theta_seq)
est_theta_se <- function(stat, se_name, model_name, phi_c, phi_u, theta){
  lc<- length(phi_c); lu <- length(phi_u); lt <- length(theta); lm <- length(model_name)
  data_se <- data.frame(matrix(NA, ncol=3+length(se_name),
                               nrow=lu*lc*lt*lm)) %>%
    set_names(c("theta","phi_c","phi_u",se_name))
  for (i in seq_along(stat)) {
    for(t in seq_along(theta)){
      statit <- as.data.frame(stat[[i]][[t]])
      se_included <- statit %>% dplyr::select(contains("_se"))
      se_notincluded <- setdiff(colnames(se_included),se_name)
      se_included %<>% dplyr::select(-se_notincluded)
      data_se[lu*lc*lt*(i-1)+(lu*lc*(t-1)+1):(lu*lc*t),c("theta",names(se_included))] <-
        cbind(theta[t],se_included)
    }
  }
  data_se$phi_c <- rep(stat[[1]][[1]][,"phi_c"], lt*lm)
  data_se$phi_u <- rep(stat[[1]][[1]][,"phi_u"], lt*lm)
  data_se$model <- rep(model_name, each = lu*lc*lt)
  return(data_se)
}

biasRMSE_betac <- function(stat, model_name, phi_c, betac, phi_u, theta){
  res <- array(0, dim=c(length(theta), length(phi_c), length(phi_u), length(betac), length(model_name), 4))
  for (t in seq_along(theta)) {
    for (c in seq_along(phi_c)) {
      for (b in seq_along(betac)) {
        for (u in seq_along(phi_u)) {
          for (m in seq_along(model_name)) {
            #med 1
            res[t,c,u,b,m,1] <- stat[[m]][[t]][(c-1)*length(betac)*length(phi_u)+(u-1)*length(betac)+b,"ATE_est_mean"]
            #lb 2
            res[t,c,u,b,m,2] <- stat[[m]][[t]][(c-1)*length(betac)*length(phi_u)+(u-1)*length(betac)+b,"25%ATE_est"]
            #ub 3
            res[t,c,u,b,m,3] <- stat[[m]][[t]][(c-1)*length(betac)*length(phi_u)+(u-1)*length(betac)+b,"75%ATE_est"]
            #rMSE 4
            res[t,c,u,b,m,4] <- stat[[m]][[t]][(c-1)*length(betac)*length(phi_u)+(u-1)*length(betac)+b,"rMSE"]
          }
        }
      }
    }
  }

  size = prod(dim(res)[1:5])
  tesq <- phi_c_seq <- Beta_c <- phi_u_seq <- model_seq <- meanest <- lb <- ub <- rMSE <- rep(0,size)
  l1 <- dim(res)[5]; l2 <- l1*dim(res)[4]; l3 <- l2*dim(res)[3]; l4 <- l3*dim(res)[2]
  for (t in seq_along(theta)) {
    for (c in seq_along(phi_c)) {
      for (b in seq_along(betac)) {
        for (u in seq_along(phi_u)) {
          for (m in seq_along(model_name)) {
            tesq[(t-1)*l4+(c-1)*l3+(u-1)*l2+(b-1)*l1+m] = theta[t]
            phi_c_seq[(t-1)*l4+(c-1)*l3+(u-1)*l2+(b-1)*l1+m] = phi_c[c]
            phi_u_seq[(t-1)*l4+(c-1)*l3+(u-1)*l2+(b-1)*l1+m] = phi_u[u]
            Beta_c[(t-1)*l4+(c-1)*l3+(u-1)*l2+(b-1)*l1+m] = betac[b]
            model_seq[(t-1)*l4+(c-1)*l3+(u-1)*l2+(b-1)*l1+m] = model_name[m]
            meanest[(t-1)*l4+(c-1)*l3+(u-1)*l2+(b-1)*l1+m] = res[t,c,u,b,m,1]
            lb[(t-1)*l4+(c-1)*l3+(u-1)*l2+(b-1)*l1+m] = res[t,c,u,b,m,2]
            ub[(t-1)*l4+(c-1)*l3+(u-1)*l2+(b-1)*l1+m] = res[t,c,u,b,m,3]
            rMSE[(t-1)*l4+(c-1)*l3+(u-1)*l2+(b-1)*l1+m] = res[t,c,u,b,m,4]
          }
        }
      }
    }
  }
  result <- data.frame(theta = factor(tesq), phic = factor(phi_c_seq,levels=phi_c),
                       phiu = factor(phi_u_seq,levels=phi_u),
                       beta_c = factor(Beta_c,levels=betac),
                       model = factor(model_seq, levels = model_name),
                       mean = meanest, lb = lb, ub = ub, rMSE = rMSE)
  if("Naive" %in% model_name){
    result$Naive = rep(result[seq(which(model_name == "Naive"),nrow(result),length(model_name)),"mean"],
                       each=length(model_name))
  }

  return(result)
}


edf_time <- function(stat,model_name, phi_c, phi_u, theta){
  lc<- length(phi_c); lu <- length(phi_u); lt <- length(theta); lm <- length(model_name)
  data <- data.frame(matrix(NA, ncol=8,
                            nrow=lu*lc*lt*lm)) %>%
    set_names(c("theta","phi_c","phi_u","edf","sp",
                "medtime","25%time","75%time"))
  for (i in seq_along(stat)) {
    for(t in seq_along(theta)){
      statit <- as.data.frame(stat[[i]][[t]])
      if(grepl("FDF", model_name[i], fixed = TRUE)){
        edf <- statit %>% dplyr::select(k)
        sp <- 0
      }else if(model_name[i] == "1NN"){
        edf <- sp <- 0
      }else{
        edf <- statit %>% dplyr::select(contains("edf"))
        sp <- statit %>% dplyr::select(contains("sp"))
      }
      time <- statit %>% dplyr::select(contains("time"))
      data[lu*lc*lt*(i-1)+(lu*lc*(t-1)+1):(lu*lc*t),
           c("theta","edf","sp","medtime","25%time","75%time")] <-
        cbind(theta[t],edf,sp,time)
    }
  }
  data$phi_c <- rep(stat[[1]][[1]][,"phi_c"], lt*lm)
  data$phi_u <- rep(stat[[1]][[1]][,"phi_u"], lt*lm)
  data$model <- rep(model_name, each = lu*lc*lt)
  return(data)
}

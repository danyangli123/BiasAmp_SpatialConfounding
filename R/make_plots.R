#df=df1; proc="z_c"; ini
surface_kc <- function(df, proc, ini, pred.grid, expression, file.path=NA){
  if(!is.na(file.path)){
    if(file.exists(file.path)){
      kc <- readRDS(file.path)
    }else{
      process <- as.geodata(df,coords.col=c(3,4),data.col=which(colnames(df)==proc))
      fit <- likfit(process,ini=ini)
      kc <- krige.conv(process,loc=pred.grid, krige=krige.control(obj.m=fit))
      saveRDS(kc, file.path)
    }
  }else{
    process <- as.geodata(df,coords.col=c(3,4),data.col=which(colnames(df)==proc))
    fit <- likfit(process,ini=ini)
    kc <- krige.conv(process,loc=pred.grid, krige=krige.control(obj.m=fit))
  }

  n <- sqrt(dim(pred.grid)[1])
  image.plot(x=pred.grid[["Var1"]][1:n], y=unique(pred.grid[["Var2"]]),
             z=matrix(kc$predict,nrow=n, ncol=n),col=terrain.colors(100),
             nlevel = 20, xlab=" ",ylab=" ",
             main=expression, cex.axis=1.5, cex.lab=1.8, cex.main=1.8, yaxt = "n", xaxt = "n")
  points(df$p0, df$p1, col = as.factor(df$X_indv), pch = 16, cex=1.5)
}

#constant effect surface default
positivity_visual <- function(output_dir, beta=1, coef_gen_y=1, coef_inter=0, thresh=0,
                              theta, phi_c, phi_u, sigsq_y_true=4){
  zc_seq <- coef_gen_y*seq(-3,3,0.1)
  emp_prop <- NULL
  for (i in 1:100) {
    df_comps <- readRDS(paste0(output_dir,i,".rds"))
    df <- gen_data_from_comps(df_comps, beta=beta, coef_gen_y=coef_gen_y, coef_inter=coef_inter,
                              thresh=thresh, theta=theta, phi_c=phi_c, phi_u=phi_u, sigsq_y_true=sigsq_y_true)
    emp_propi <- matrix(NA, nrow = length(zc_seq)-1, ncol = 4)
    for (v in 1:(length(zc_seq)-1)) {
      id_within <- (zc_seq[v] <= df$z_c) & (df$z_c < zc_seq[v+1])
      num_expv <- sum(df[id_within,"X_indv"])
      if(sum(id_within) == 0) {
        prop_expv <- NA
      } else {
        prop_expv <- num_expv/sum(id_within)
      }
      emp_propi[v,] <- c((zc_seq[v]+zc_seq[v+1])/2, prop_expv, sum(id_within), num_expv)
      #cat("v=",v,"\n")
    }
    emp_propi <- data.frame(zc = emp_propi[,1], emp_prob = emp_propi[,2], num_obs = emp_propi[,3], num_exp = emp_propi[,4])
    emp_prop <- rbind(emp_prop,emp_propi)
    cat(dim(emp_prop), "\n")
  }
  emp_prop %<>% filter(num_obs != 0) %>% group_by(zc) %>% summarise(prob_mean = mean(emp_prob), prob_sd = sd(emp_prob),
                                                                    valid_num = n(), obs_num_mean = mean(num_obs),
                                                                    exp_num_mean = mean(num_exp))
  g <- ggplot(emp_prop,aes(x=zc, y=prob_mean)) +
    geom_line()+geom_linerange(aes(ymin=prob_mean-1.96*prob_sd, ymax=prob_mean+1.96*prob_sd),color="red") +
    labs(title = paste0("phic=",phi_c,"phiu=",phi_u))

  return(list(emp_prop_df = emp_prop, gplot = g))
}

match_pairs <- function(df, typ, clp, rep_need, proc, ini, pred.grid) {
  #M1, ATT
  m <- Matchit_list(df, typ, clp, rep = rep_need)
  df_trt <- df[m$index.treated,]
  df_crl <- df[m$index.control,]
  surface_kc(df, proc, ini, pred.grid)
  points(df$p0, df$p1, col=as.factor(df$X_indv),xlab="p0",ylab='p1',main='indv')
  for (i in 1:length(m$index.treated)) {
    arrows(x0 = df_trt$p0[i], y0=df_trt$p1[i], x1=df_crl$p0[i], y1=df_crl$p1[i],length=0)
  }
}

combine_zc_matchparis <- function(df, cov, ini, pred.grid, typ, clp, weight, est, case) {
  m <- Match_list(df, M=1, covariates=cov, clp_typ=typ, clp=clp, weight=weight, est=est, ties=TRUE)
  if(case=="bin"){
    m <- m$bin_ate
  } else if (case == "indv") {
    m <- m$indv_ate
  } else {
    cat("not available")
  }
  df_trt <- df[m$index.treated,]
  df_crl <- df[m$index.control,]
  confounder_kc(df, "z_c", ini, pred.grid)
  if(case=="bin"){
    points(df$p0, df$p1, col=as.factor(df$X_binom),xlab="p0",ylab='p1',main='bin')
  } else if (case == "indv") {
    points(df$p0, df$p1, col=as.factor(df$X_indv),xlab="p0",ylab='p1',main='indv')
  } else {
    cat("not available")
  }
  #for (i in 1:length(m$index.treated)) {
    #arrows(x0 = df_trt$p0[i], y0=df_trt$p1[i], x1=df_crl$p0[i], y1=df_crl$p1[i],length=0)
  #}

}

#m is a Match_list obj.
AB_zc <- function(df, m, case) {
  c1 <- rgb(0,0,0, max = 255, alpha = 80, names = "gray")
  c2 <- rgb(255,192,203, max = 255, alpha = 255, names = "lt.pink")
  if (case == "bin") {
    crl0 <- hist(df$z_c[df$X_binom==0], freq = TRUE, xlab = "z_c")
    crl1 <- hist(df$z_c[m$bin_att$index.control], freq = TRUE, xlab = "z_c")
    crl2 <- hist(df$z_c[m$bin_ate$index.control], freq = TRUE, xlab = "z_c")
    par(mfrow=c(3,1))
    hist(df$z_c[df$X_binom==1], freq = TRUE, ylim = c(0,80), xlab = "z_c", col=c2, main = "Before match")
    plot(crl0, add = TRUE, col=c1)
    hist(df$z_c[m$bin_att$index.treated], freq = TRUE, ylim = c(0,100), xlab = "z_c", col=c2, main = "After match ATT")
    plot(crl1, add = TRUE, col=c1)
    hist(df$z_c[m$bin_ate$index.treated], freq = TRUE, ylim = c(0,120), xlab = "z_c", col=c2, main = "After match ATE")
    plot(crl2, add = TRUE, col=c1)
    cat("bin_att estimate=",m$bin_att$est, "\nbin_ate estimate=", m$bin_ate$est)
  } else if (case == "indv") {
    crl0 <- hist(df$z_c[df$X_indv==0], freq = TRUE, xlab = "z_c")
    crl1 <- hist(df$z_c[m$indv_att$index.control], freq = TRUE, xlab = "z_c")
    crl2 <- hist(df$z_c[m$indv_ate$index.control], freq = TRUE, xlab = "z_c")
    par(mfrow=c(3,1))
    hist(df$z_c[df$X_indv==1], freq = TRUE, ylim = c(0,80), xlab = "z_c", col=c2, main = "Before match")
    plot(crl0, add = TRUE, col=c1)
    hist(df$z_c[m$indv_att$index.treated], freq = TRUE, ylim = c(0,100), xlab = "z_c", col=c2, main = "After match ATT")
    plot(crl1, add = TRUE, col=c1)
    hist(df$z_c[m$indv_ate$index.treated], freq = TRUE, ylim = c(0,120), xlab = "z_c", col=c2, main = "After match ATE")
    plot(crl2, add = TRUE, col=c1)
    cat("indv_att estimate=",m$indv_att$est, "\nindv_ate estimate=", m$indv_ate$est)
  } else {
    cat("not available")
  }
}
#AB_zc(df, m, "indv")


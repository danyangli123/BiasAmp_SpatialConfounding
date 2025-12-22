#extract summary from gam object
extract_ests_gam <- function(mod, inter, coef, s = NULL, cluster = FALSE, match_data=NULL) {
  if (is.null(s)) {
    s <- summary(mod)
  }

  sp <- ifelse(is.null(mod$full.sp), ifelse(is.null(mod$sp), NA, mod$sp), mod$full.sp)
  if(inter){
    edf = c("X1_edf","X0_edf")
    colname12 <- c("coef_est", "coef_se")
  }else{
    edf =  c("edf")
    colname12 <- c("ATE_est", "ATE_se")
  }
  if(cluster){
    cl_stde <- sqrt(sandwich::vcovCL(mod, cluster = match_data$cluster_id)[coef,coef])
  }else{
    cl_stde <- 0
  }

  summ0 <- c(s$p.table[coef, ], s$edf, #sqrt(sandwich::sandwich(mod)[coef, coef]),
             cl_stde,
             BIC(mod), mod$gcv.ubre, sp) %>% set_names(c(colname12,"indv_t", "coef_p","edf", #"coef_se_rob",
                                                         "coef_se_CL","bic","gcv","sp"))%>%
    as.list() %>% as_tibble()

}


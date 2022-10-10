library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)

setwd('/home/huangqiangru/pcr')
TimeAndStrategy <- read.csv('TimeAndStrategy_pcr-2days.csv')
source("./initialize_pcr.R")
source("./functions_pcr.R")

set.seed(10031)  

stochastic_NH <- function(parms, Ns, delta_t, t, num_staff){

  beta = parms[["beta"]] 
  beta.s = parms[["beta.s"]] 
  beta.rm = parms[["beta.rm"]] 
  sigma1 = parms[["sigma1"]] 
  sigma2 = parms[["sigma2"]] 
  gamma = parms[["gamma1"]]
  gamma = parms[["gamma2"]]
  I.C_time=parms[["I.C_time"]] 
  I.C = ifelse(t>I.C_time,parms[["I.C"]],0) 
  k.HH = parms[["k.HH"]] 
  k.HR=parms[["k.HR"]]
  k.RH=parms[["k.RH"]]
  k.RR=parms[["k.RR"]] 
  alpha.r= parms[["alpha.r"]]
  alpha.hcw= parms[["alpha.hcw"]]
  ID= parms[["id.I"]] 
  int.r = parms[["int.r"]] 
  int.hcw = parms[["int.hcw"]] 
  int.rooms = parms[["int.rooms"]] 
  mu.C = parms[["mu.C"]] 
  mu.NC = parms[["mu.NC"]] 
  ppe = parms[["ppe"]]
  ppe_r = parms[["ppe_r"]]
  prop_rhcwR = parms[["prop_rhcwR"]] 
  VL.PCR_threshold = parms[["VL.PCR_threshold"]] 
  VL.Antigen_threshold = parms[["VL.Antigen_threshold"]]
  int_time = parms[["int_time"]] 
  VE_s1 = parms[["VE_s1"]] 
  VE_p1 = parms[["VE_p1"]]
  VE_i1 = parms[["VE_i1"]]
  VE_s2 = parms[["VE_s2"]] 
  VE_p2 = parms[["VE_p2"]]
  VE_i2 = parms[["VE_i2"]]
  VE_s3 = parms[["VE_s3"]] 
  VE_p3 = parms[["VE_p3"]]
  VE_i3 = parms[["VE_i3"]]
  VE_s4 = parms[["VE_s4"]] 
  VE_p4 = parms[["VE_p4"]]
  VE_i4 = parms[["VE_i4"]]
  V_doses = parms[["V_doses"]]
  boost_vax_coverage_r = parms[["boost_vax_coverage_r"]] 
  full_vax_coverage_r = parms[["full_vax_coverage_r"]] 
  refusal_r = parms[["refusal_r"]]
  boost_vax_coverage_hcw = parms[["boost_vax_coverage_hcw"]]
  full_vax_coverage_hcw = parms[["full_vax_coverage_hcw"]]
  refusal_hcw= parms[["refusal_hcw"]] 
  comm.vax = parms[["comm.vax"]] 


  S.rNC=Ns[["S.rNC"]] 
  E.rNC=Ns[["E.rNC"]]
  A.rNC=Ns[["A.rNC"]]
  I.rNC=Ns[["I.rNC"]]
  R.rNC=Ns[["R.rNC"]]
  I.rC=Ns[["I.rC"]]
  R.rC=Ns[["R.rC"]]
  S.hcwNC=Ns[["S.hcwNC"]]
  E.hcwNC=Ns[["E.hcwNC"]]
  A.hcwNC=Ns[["A.hcwNC"]]
  I.hcwNC=Ns[["I.hcwNC"]]
  R.hcwNC=Ns[["R.hcwNC"]]
  S.hcwC=Ns[["S.hcwC"]]
  E.hcwC=Ns[["E.hcwC"]]
  A.hcwC=Ns[["A.hcwC"]]
  I.hcwC=Ns[["I.hcwC"]]
  R.hcwC=Ns[["R.hcwC"]]
  I.hcwH=Ns[["I.hcwH"]]
  S.rNC3=Ns[["S.rNC3"]]
  S.hcwNC3=Ns[["S.hcwNC3"]]
  S.hcwC3=Ns[["S.hcwC3"]]
  death.r=Ns[["death.r"]]
  death.rC=Ns[["death.rC"]]
  cum_inc_r=Ns[["cum_inc_r"]]
  cum_inc_hcw=Ns[["cum_inc_hcw"]]
  inc_vax1_r=Ns[["inc_vax1_r"]] 
  inc_vax2_r=Ns[["inc_vax2_r"]]
  inc_vax1_hcw=Ns[["inc_vax1_hcw"]] 
  inc_vax2_hcw=Ns[["inc_vax2_hcw"]]
  cum_inc_community=Ns[["cum_inc_community"]]
  total=Ns[["total"]]
  total_hcw=Ns[["total_hcw"]]
  family=Ns[["family"]] 
  vax_count1_r=Ns[["vax_count1_r"]]
  vax_count2_r=Ns[["vax_count2_r"]]
  vax_count1_hcw=Ns[["vax_count1_hcw"]]
  vax_count2_hcw=Ns[["vax_count2_hcw"]]

  N.rNC <- nrow(S.rNC) + nrow(E.rNC) + nrow(A.rNC) + nrow(I.rNC) + nrow(R.rNC)
  N.rC <- ifelse(nrow(I.rC) + nrow(R.rC) > 0, nrow(I.rC) + nrow(R.rC), 1) 
  N.rNC + N.rC


  N.hcwNC <-  nrow(S.hcwNC) + nrow(E.hcwNC) + nrow(A.hcwNC) + nrow(I.hcwNC) + nrow(R.hcwNC)
  N.hcwC <- nrow(S.hcwC) + nrow(E.hcwC) + nrow(A.hcwC) + nrow(I.hcwC) + nrow(R.hcwC)
  N.hcwNC + N.hcwC


  recover.A.rNC <- recover_or_test_r(A.rNC,parms, "A",t,TimeAndStrategy)
  R.rNC <- rbind(R.rNC,recover.A.rNC[[1]])    
  I.rC <- rbind(I.rC,recover.A.rNC[[2]])     
  A.rNC <- recover.A.rNC[[3]]               
  
  family[(family$ResID %in% recover.A.rNC[[2]]$ID),"ResidentFlag"] =0
  

  recover.I.rNC <- recover_or_test_r(I.rNC,parms, "I",t,TimeAndStrategy)
  R.rNC <- rbind(R.rNC,recover.I.rNC[[1]])     # move I.rNC to R.rNC if recovered
  I.rC <- rbind(I.rC,recover.I.rNC[[2]])      # move I.rNC to I.rC if identified as posiitive
  I.rNC <- recover.I.rNC[[3]]                
  
  family[(family$ResID %in% recover.I.rNC[[2]]$ID),"ResidentFlag"] =0

  recover.I.rC <- recover_I.rC(I.rC, parms)

  if (int.r==1 & t>int_time){   
    R.rNC <- rbind(R.rNC,recover.I.rC[[1]]) 
    R.room <- recover.I.rC[[1]]$ID          # need room assignment
  } else {
    R.rC <- rbind(R.rC,recover.I.rC[[1]])   
    R.room <- NULL
  }
  I.rC <- recover.I.rC[[2]] # keep rest in I.rC


  ratio <- (N.hcwNC + N.hcwC)/(N.rNC + N.rC)

  recover.A.hcwNC <- recover_or_test_hcw(A.hcwNC,parms,total_hcw, "A", ratio,t,TimeAndStrategy)
  R.hcwNC <- rbind(R.hcwNC,recover.A.hcwNC[[1]]) # move A.hcwNC to R.hcwNC if recover
  I.hcwH <- rbind(I.hcwH,recover.A.hcwNC[[2]]) # move A.hcwNC to I.hcwH if test positive
  A.hcwNC <- recover.A.hcwNC[[3]] # keep rest of A.hcwNC in A.hcwNC
  E.hcwNC <- rbind(E.hcwNC,recover.A.hcwNC[[4]])
  S.hcwNC <- rbind(S.hcwNC,recover.A.hcwNC[[5]])
  R.hcwNC <- rbind(R.hcwNC,recover.A.hcwNC[[6]])
  total_hcw <- total_hcw - nrow(recover.A.hcwNC[[2]])+ nrow(recover.A.hcwNC[[4]]) + nrow(recover.A.hcwNC[[5]]) + nrow(recover.A.hcwNC[[6]]) #新进的工作人员
  
  vax_count1_hcw <- vax_count1_hcw + nrow(recover.A.hcwNC[[4]] %>% subset(V_doses==2))+ nrow(recover.A.hcwNC[[5]] %>% subset(V_doses==2))+ nrow(recover.A.hcwNC[[6]] %>% subset(V_doses==2)) #全程接种
  vax_count2_hcw <- vax_count2_hcw + nrow(recover.A.hcwNC[[4]] %>% subset(V_doses==3))+ nrow(recover.A.hcwNC[[5]] %>% subset(V_doses==3))+ nrow(recover.A.hcwNC[[6]] %>% subset(V_doses==3)) #加强针接种

  recover.I.hcwNC <- recover_or_test_hcw(I.hcwNC,parms,total_hcw, "I", ratio,t,TimeAndStrategy)
  R.hcwNC <- rbind(R.hcwNC,recover.I.hcwNC[[1]]) # move I.hcwNC to R.hcwNC if recover
  I.hcwH <- rbind(I.hcwH,recover.I.hcwNC[[2]]) # move I.hcwNC to I.hcwH if test positive
  I.hcwNC <- recover.I.hcwNC[[3]] # keep rest of I.hcwNC in I.hcwNC
  E.hcwNC <- rbind(E.hcwNC,recover.I.hcwNC[[4]])
  S.hcwNC <- rbind(S.hcwNC,recover.I.hcwNC[[5]])
  R.hcwNC <- rbind(R.hcwNC,recover.I.hcwNC[[6]])
  total_hcw <- total_hcw - nrow(recover.I.hcwNC[[2]])+ nrow(recover.I.hcwNC[[4]]) + nrow(recover.I.hcwNC[[5]]) + nrow(recover.I.hcwNC[[6]])
  vax_count1_hcw <- vax_count1_hcw + nrow(recover.I.hcwNC[[4]] %>% subset(V_doses==2))+ nrow(recover.I.hcwNC[[5]] %>% subset(V_doses==2))+ nrow(recover.I.hcwNC[[6]] %>% subset(V_doses==2)) #全程接种
  vax_count2_hcw <- vax_count2_hcw + nrow(recover.I.hcwNC[[4]] %>% subset(V_doses==3))+ nrow(recover.I.hcwNC[[5]] %>% subset(V_doses==3))+ nrow(recover.I.hcwNC[[6]] %>% subset(V_doses==3)) #加强针接种
  
  recover.A.hcwC <- recover_or_test_hcw(A.hcwC,parms,total_hcw, "A", ratio,t,TimeAndStrategy)
  R.hcwC <- rbind(R.hcwC,recover.A.hcwC[[1]]) # move A.hcwC to R.hcwC if recover
  I.hcwH <- rbind(I.hcwH,recover.A.hcwC[[2]]) # move A.hcwC to I.hcwH if test positive
  A.hcwC <- recover.A.hcwC[[3]] # keep rest of A.hcwC in A.hcwC
  E.hcwC <- rbind(E.hcwC,recover.A.hcwC[[4]])
  S.hcwC <- rbind(S.hcwC,recover.A.hcwC[[5]])
  R.hcwC <- rbind(R.hcwC,recover.A.hcwC[[6]])
  total_hcw <- total_hcw - nrow(recover.A.hcwC[[2]]) + nrow(recover.A.hcwC[[4]]) + nrow(recover.A.hcwC[[5]]) + nrow(recover.A.hcwC[[6]])
  vax_count1_hcw <- vax_count1_hcw + nrow(recover.A.hcwC[[4]] %>% subset(V_doses==2))+ nrow(recover.A.hcwC[[5]] %>% subset(V_doses==2))+ nrow(recover.A.hcwC[[6]] %>% subset(V_doses==2)) #全程接种
  vax_count2_hcw <- vax_count2_hcw + nrow(recover.A.hcwC[[4]] %>% subset(V_doses==3))+ nrow(recover.A.hcwC[[5]] %>% subset(V_doses==3))+ nrow(recover.A.hcwC[[6]] %>% subset(V_doses==3)) #加强针接种
  
  recover.I.hcwC <- recover_or_test_hcw(I.hcwC,parms,total_hcw, "I", ratio,t,TimeAndStrategy)
  R.hcwC <- rbind(R.hcwC,recover.I.hcwC[[1]]) # move I.hcwC to R.hcwC if recover
  I.hcwH <- rbind(I.hcwH,recover.I.hcwC[[2]]) # move I.hcwC to I.hcwH if test positive
  I.hcwC <- recover.I.hcwC[[3]] # keep rest of I.hcwC in I.hcwC
  E.hcwC <- rbind(E.hcwC,recover.I.hcwC[[4]])
  S.hcwC <- rbind(S.hcwC,recover.I.hcwC[[5]])
  R.hcwC <- rbind(R.hcwC,recover.I.hcwC[[6]])
  total_hcw <- total_hcw -nrow(recover.I.hcwC[[2]])+ nrow(recover.I.hcwC[[4]]) + nrow(recover.I.hcwC[[5]]) + nrow(recover.I.hcwC[[6]])
  vax_count1_hcw <- vax_count1_hcw + nrow(recover.I.hcwC[[4]] %>% subset(V_doses==2))+ nrow(recover.I.hcwC[[5]] %>% subset(V_doses==2))+ nrow(recover.I.hcwC[[6]] %>% subset(V_doses==2)) #全程接种
  vax_count2_hcw <- vax_count2_hcw + nrow(recover.I.hcwC[[4]] %>% subset(V_doses==3))+ nrow(recover.I.hcwC[[5]] %>% subset(V_doses==3))+ nrow(recover.I.hcwC[[6]] %>% subset(V_doses==3)) #加强针接种
  
  recover.home <- recover_home_hcw(I.hcwH)
  if (int.hcw==1 & t>int_time){
    R.hcwNC <- rbind(R.hcwNC,recover.home[[1]]) 
    if (nrow(recover.home[[1]])>0){
      num_remove <- nrow(recover.home[[1]])
      temp <- unlist(c(S.hcwNC %>% subset(ID>100000) %>% dplyr::select(ID),
                       S.hcwC %>% subset(ID>100000) %>% dplyr::select(ID),
                       E.hcwNC %>% subset(ID>100000) %>% dplyr::select(ID),
                       E.hcwC %>% subset(ID>100000) %>% dplyr::select(ID),
                       A.hcwNC %>% subset(ID>100000) %>% dplyr::select(ID),
                       A.hcwC %>% subset(ID>100000) %>% dplyr::select(ID),
                       I.hcwNC %>% subset(ID>100000) %>% dplyr::select(ID),
                       I.hcwC %>% subset(ID>100000) %>% dplyr::select(ID),
                       R.hcwNC %>% subset(ID>100000) %>% dplyr::select(ID),
                       R.hcwC %>% subset(ID>100000) %>% dplyr::select(ID)), use.names=FALSE)
      num_remove <- min(num_remove,length(temp))
      if (length(temp)>0){
        if (length(temp)>1){
          temp_remove <- cbind("ID"=sample(temp,num_remove,replace=FALSE))
        } else{
          temp_remove <- temp
        }
        S.hcwNC %>%
          subset(!(ID %in% temp_remove)) -> S.hcwNC
        S.hcwC %>%
          subset(!(ID %in% temp_remove)) -> S.hcwC
        E.hcwNC %>%
          subset(!(ID %in% temp_remove)) -> E.hcwNC
        E.hcwC %>%
          subset(!(ID %in% temp_remove)) -> E.hcwC
        A.hcwNC %>%
          subset(!(ID %in% temp_remove)) -> A.hcwNC
        A.hcwC %>%
          subset(!(ID %in% temp_remove)) -> A.hcwC
        I.hcwNC %>%
          subset(!(ID %in% temp_remove)) -> I.hcwNC
        I.hcwC %>%
          subset(!(ID %in% temp_remove)) -> I.hcwC
        R.hcwNC %>%
          subset(!(ID %in% temp_remove)) -> R.hcwNC
        R.hcwC %>%
          subset(!(ID %in% temp_remove)) -> R.hcwC
      }
    }
  } else{
    R.hcwC <- rbind(R.hcwC,recover.home[[1]])  
    if (nrow(recover.home[[1]])>0){
      num_remove <- nrow(recover.home[[1]])
      temp <- unlist(c(S.hcwNC %>% subset(ID>100000) %>% dplyr::select(ID),
                       S.hcwC %>% subset(ID>100000) %>% dplyr::select(ID),
                       E.hcwNC %>% subset(ID>100000) %>% dplyr::select(ID),
                       E.hcwC %>% subset(ID>100000) %>% dplyr::select(ID),
                       A.hcwNC %>% subset(ID>100000) %>% dplyr::select(ID),
                       A.hcwC %>% subset(ID>100000) %>% dplyr::select(ID),
                       I.hcwNC %>% subset(ID>100000) %>% dplyr::select(ID),
                       I.hcwC %>% subset(ID>100000) %>% dplyr::select(ID),
                       R.hcwNC %>% subset(ID>100000) %>% dplyr::select(ID),
                       R.hcwC %>% subset(ID>100000) %>% dplyr::select(ID)), use.names=FALSE)
      num_remove <- min(num_remove,length(temp))
      if (length(temp)>0){
        if (length(temp)>1){
          temp_remove <- cbind("ID"=sample(temp,num_remove,replace=FALSE))
        } else{
          temp_remove <- temp
        }
        S.hcwNC %>%
          subset(!(ID %in% temp_remove)) -> S.hcwNC
        S.hcwC %>%
          subset(!(ID %in% temp_remove)) -> S.hcwC
        E.hcwNC %>%
          subset(!(ID %in% temp_remove)) -> E.hcwNC
        E.hcwC %>%
          subset(!(ID %in% temp_remove)) -> E.hcwC
        A.hcwNC %>%
          subset(!(ID %in% temp_remove)) -> A.hcwNC
        A.hcwC %>%
          subset(!(ID %in% temp_remove)) -> A.hcwC
        I.hcwNC %>%
          subset(!(ID %in% temp_remove)) -> I.hcwNC
        I.hcwC %>%
          subset(!(ID %in% temp_remove)) -> I.hcwC
        R.hcwNC %>%
          subset(!(ID %in% temp_remove)) -> R.hcwNC
        R.hcwC %>%
          subset(!(ID %in% temp_remove)) -> R.hcwC
      }
    }
  }

  I.hcwH <- recover.home[[2]]  
  
  total_hcw <- total_hcw+ nrow(recover.home[[1]])

  # moving hcw
  prop.rNC <- (nrow(S.rNC) + nrow(E.rNC) + nrow(A.rNC) + nrow(I.rNC) + nrow(R.rNC))/
    (nrow(S.rNC) + nrow(E.rNC) + nrow(A.rNC) + nrow(I.rNC) + nrow(R.rNC) + nrow(I.rC) + nrow(R.rC))   ##NC队列中居民数占所有居民数的比例
  prop.hcwNC <- (nrow(S.hcwNC) + nrow(E.hcwNC) + nrow(A.hcwNC) +  nrow(I.hcwNC) + nrow(R.hcwNC))/
    (nrow(S.hcwNC) + nrow(E.hcwNC) + nrow(A.hcwNC) +  nrow(I.hcwNC) + nrow(R.hcwNC) +
       nrow(S.hcwC) + nrow(E.hcwC) + nrow(A.hcwC) + nrow(I.hcwC) + nrow(R.hcwC))  ##NC队列中hcw数占所有hcw数的比例

  moved_hcw <- move_hcw(S.hcwNC,
                        E.hcwNC,
                        A.hcwNC,
                        I.hcwNC,
                        R.hcwNC,
                        S.hcwC,
                        E.hcwC,
                        A.hcwC,
                        I.hcwC,
                        R.hcwC,
                        prop.rNC,
                        prop.hcwNC)
  S.hcwNC <- moved_hcw[[1]]
  E.hcwNC <- moved_hcw[[2]]
  A.hcwNC <- moved_hcw[[3]]
  I.hcwNC <- moved_hcw[[4]]
  R.hcwNC <- moved_hcw[[5]]
  S.hcwC <- moved_hcw[[6]]
  E.hcwC <- moved_hcw[[7]]
  A.hcwC <- moved_hcw[[8]]
  I.hcwC <- moved_hcw[[9]]
  R.hcwC <- moved_hcw[[10]]

  # E -> A/I
  infect_E.rNC <- E_to_I_r(E.rNC,parms)
  A.rNC <- rbind(A.rNC,infect_E.rNC[[1]])
  I.rNC <- rbind(I.rNC,infect_E.rNC[[2]])
  E.rNC <- infect_E.rNC[[3]]
  
  asymptomatic_inc_r <- nrow(infect_E.rNC[[1]]) 
  symptomatic_inc_r <- nrow(infect_E.rNC[[2]]) 


  infect_E.hcwNC <- E_to_I_hcw(E.hcwNC,parms)
  A.hcwNC <- rbind(A.hcwNC,infect_E.hcwNC[[1]])
  I.hcwNC <- rbind(I.hcwNC,infect_E.hcwNC[[2]])
  E.hcwNC <- infect_E.hcwNC[[3]]

  infect_E.hcwC <- E_to_I_hcw(E.hcwC,parms)
  A.hcwC <- rbind(A.hcwC,infect_E.hcwC[[1]])
  I.hcwC <- rbind(I.hcwC,infect_E.hcwC[[2]])
  E.hcwC <- infect_E.hcwC[[3]]
  
  asymptomatic_inc_hcw <- nrow(infect_E.hcwNC[[1]]) + nrow(infect_E.hcwC[[1]])  #新增无症状感染hcw
  symptomatic_inc_hcw <- nrow(infect_E.hcwNC[[2]]) + nrow(infect_E.hcwC[[2]])  #新增有症状感染hcw
  
  # update contact rates based on ratio
  k.HR <- k.HR*num_staff/(N.hcwNC + N.hcwC)
  #print(k.HR)

  # S -> E
  S.rNC2 <- S.rNC
  inc_vax1_r <- 0
  inc_vax2_r <- 0
  
  if(nrow(S.rNC)>=1){
    if(TimeAndStrategy[t+1,2]==1){k.RR = 0} else {k.RR = 1}
    for (i in 1:nrow(S.rNC)){
      roommate.exp <- expose_roommate(family,S.rNC$ID[i],A.rNC,I.rNC) #该函数输出结果为1/0
      S.rNC.exp <- rbinom(1,1, I.C*(1-S.rNC$VE_s[i]))
     
      # adjust betas based on vaccination status
      beta.rm.v <- beta.rm*(1-S.rNC$VE_s[i])
      beta.s.v <- beta.s*(1-S.rNC$VE_s[i])
      beta.v <- beta*(1-S.rNC$VE_s[i]) # reduce susceptibility by VEs (assuming leaky here)

      exp.prob <- beta.rm.v*roommate.exp +beta.v*ppe_r*S.rNC.exp+
        beta.v*ppe_r*(ifelse(N.rNC>0,k.RR*(sum(I.rNC$Infectiousness*(1-I.rNC$VE_i))+sum(A.rNC$Infectiousness*(1-A.rNC$VE_i)))/N.rNC,0)) +
        beta.s.v*ppe*(ifelse(N.hcwNC>0,k.RH*(sum(I.hcwNC$Infectiousness*(1-I.hcwNC$VE_i))+sum(A.hcwNC$Infectiousness*(1-A.hcwNC$VE_i)))/N.hcwNC,0))
      if (exp.prob >= 1){
        exp.prob <- 1
      } else if (exp.prob <0){
        exp.prob <- 0
      }

      S.rNC.exposed <- rbinom(1,1,exp.prob)  
      if (S.rNC.exposed==1){  
        S.rNC3 <- rbind(S.rNC3,cbind(ID=S.rNC$ID[i],VL=0, V_doses = S.rNC$V_doses[i], V_t=S.rNC$V_t[i],
                                   VE_i = S.rNC$VE_i[i],VE_s = S.rNC$VE_s[i],VE_p = S.rNC$VE_p[i],
                                   Inc.pd= t+round(runif(1, parms[["sigma1"]], parms[["sigma2"]])), Days=t))#改$$$$$$$
        E.rNC <- rbind(E.rNC,cbind(ID=S.rNC$ID[i],VL=0, V_doses = S.rNC$V_doses[i], V_t=S.rNC$V_t[i],
                                    VE_i = S.rNC$VE_i[i],VE_s = S.rNC$VE_s[i],VE_p = S.rNC$VE_p[i],
                                    Inc.pd= t+round(runif(1, parms[["sigma1"]], parms[["sigma2"]])), Days=t))#改$$$$$$$
        S.rNC2 %>%
          subset(ID != S.rNC$ID[i]) -> S.rNC2
      }
    }
  }

  expose.S.rNC <- nrow(S.rNC) - nrow(S.rNC2) 
  S.rNC <- S.rNC2 
  inc_vax1_r <- nrow(S.rNC3 %>% subset(V_doses == 2))
  inc_vax2_r <- nrow(S.rNC3 %>% subset(V_doses == 3))

  inc_community <- 0
  inc_vax1_hcw <- 0
  inc_vax2_hcw <- 0

  S.hcwNC2 <- S.hcwNC

  if(nrow(S.hcwNC2)>=1){
    for (i in 1:nrow(S.hcwNC)){

      S.hcwNC.exposed <- rbinom(1,1, I.C*(1-S.hcwNC$VE_s[i])) 

      if (S.hcwNC.exposed != 1){ 
        beta.s.v <- beta.s*(1-S.hcwNC$VE_s[i])
        beta.v <- beta*(1-S.hcwNC$VE_s[i]) # reduce susceptibility by VEs (assuming leaky here)
        exp.prob <- beta.v*ppe*(ifelse(N.hcwNC>0,k.HH*(sum(I.hcwNC$Infectiousness*(1-I.hcwNC$VE_i))+sum(A.hcwNC$Infectiousness*(1-A.hcwNC$VE_i)))/N.hcwNC,0)) +
        beta.s.v*ppe*(ifelse(N.rNC>0,k.HR*(sum(I.rNC$Infectiousness*(1-I.rNC$VE_i))+sum(A.rNC$Infectiousness*(1-A.rNC$VE_i)))/N.rNC,0))

        if (exp.prob >= 1){
          exp.prob <- 1} else if (exp.prob <0){
          exp.prob <- 0}

        S.hcwNC.exposed <- rbinom(1,1,exp.prob)

      } else{ inc_community <- inc_community + 1} # add to comm intro tracker

      if (S.hcwNC.exposed==1){
        S.hcwNC3 <- rbind(S.hcwNC3,cbind(ID=S.hcwNC$ID[i],VL=0,
                                       V_doses = S.hcwNC$V_doses[i],V_t = S.hcwNC$V_t[i],VE_i = S.hcwNC$VE_i[i],VE_s = S.hcwNC$VE_s[i],VE_p = S.hcwNC$VE_p[i],
                                       Inc.pd=t+round(runif(1, parms[["sigma1"]], parms[["sigma2"]])), Days=t))#改$$$$$$$$
        E.hcwNC <- S.hcwNC3
        S.hcwNC2 %>%
          subset(ID != S.hcwNC$ID[i]) -> S.hcwNC2 
      }
    }
  }

  expose.S.hcwNC <- nrow(S.hcwNC) - nrow(S.hcwNC2)
  S.hcwNC <- S.hcwNC2

  S.hcwC2 <- S.hcwC
  if(nrow(S.hcwC2)>=1){
    for (i in 1:nrow(S.hcwC)){
      S.hcwC.exposed <- rbinom(1,1, I.C*(1-S.hcwC$VE_s[i]))
      if (S.hcwC.exposed != 1){ 
        beta.s.v <- beta.s*(1-S.hcwC$VE_s[i])
        beta.v <- beta*(1-S.hcwC$VE_s[i]) 
        exp.prob <- beta.v*ppe*(ifelse(N.hcwC>0,k.HH*(sum(I.hcwC$Infectiousness*(1-I.hcwC$VE_i))+sum(A.hcwC$Infectiousness*(1-A.hcwC$VE_i)))/N.hcwC,0)) +
          beta.s.v*ppe*(ifelse(N.rC>0,k.HR*(sum(I.rC$Infectiousness*(1-I.rC$VE_i)))/N.rNC,0))
        if (exp.prob >= 1){
          exp.prob <- 1} else if (exp.prob <0){
          exp.prob <- 0}

        S.hcwC.exposed <- rbinom(1,1,exp.prob)

      } else{
        inc_community <- inc_community + 1}

      if (S.hcwC.exposed==1){
        S.hcwC3 <- rbind(S.hcwC3,cbind(ID=S.hcwC$ID[i],VL=0,
                                     V_doses = S.hcwC$V_doses[i],V_t = S.hcwC$V_t[i],VE_i = S.hcwC$VE_i[i],VE_s = S.hcwC$VE_s[i],VE_p = S.hcwC$VE_p[i],
                                     Inc.pd=t+round(runif(1, parms[["sigma1"]], parms[["sigma2"]])), Days=t))
        E.hcwC <- S.hcwC3
        S.hcwC2 %>%
          subset(ID != S.hcwC$ID[i]) -> S.hcwC2
      }
    }
  }

  expose.S.hcwC <- nrow(S.hcwC) - nrow(S.hcwC2)
  S.hcwC <- S.hcwC2

  inc_vax1_hcw <- nrow(S.hcwNC3 %>% subset(V_doses == 2))+nrow(S.hcwC3 %>% subset(V_doses == 2))
  inc_vax2_hcw <- nrow(S.hcwNC3 %>% subset(V_doses == 3))+nrow(S.hcwC3 %>% subset(V_doses == 3))
  

  # death
  death.S.rNC <- death(S.rNC,mu.NC,parms,total,t,dt)
  S.rNC <- death.S.rNC[[1]]
  death.r <- rbind(death.r,death.S.rNC[[3]]) 
  
  death.E.rNC <- death(E.rNC,mu.NC,parms,total,t,dt)
  E.rNC <- death.E.rNC[[1]]
  death.r <- rbind(death.r,death.E.rNC[[3]])

  death.A.rNC <- covid_death(A.rNC,mu.C,parms,total,t,dt)
  A.rNC <- death.A.rNC[[1]]
  death.rC <- rbind(death.rC,death.A.rNC[[3]])

  death.I.rNC <- covid_death(I.rNC,mu.C,parms,total,t,dt)
  I.rNC <- death.I.rNC[[1]]
  death.rC <- rbind(death.rC,death.I.rNC[[3]])

  death.R.rNC <- death(R.rNC,mu.NC,parms,total,t,dt)
  R.rNC <- death.R.rNC[[1]]
  death.r <- rbind(death.r,death.R.rNC[[3]]) 

  death.I.rC <- covid_death(I.rC,mu.C,parms,total,t,dt)  
  I.rC <- death.I.rC[[1]]
  death.rC <- rbind(death.rC,death.I.rC[[3]])

  death.R.rC <- death(R.rC,mu.NC,parms,total,t,dt)
  R.rC <- death.R.rC[[1]]
  death.r <- rbind(death.r,death.R.rC[[3]]) 


  new_dead <- c(death.S.rNC[[3]]$ID,death.E.rNC[[3]]$ID,death.A.rNC[[3]]$ID,death.I.rNC[[3]]$ID,death.R.rNC[[3]]$ID,
                    death.I.rC[[3]]$ID,death.R.rC[[3]]$ID) 
  
  
  family[(family$ResID %in% new_dead),"ResidentFlag"] =0

  if (length(R.room)>0){
    family <- assign_rooms(family,R.room,int.rooms,t, parms)
  }

  mortality <- death.S.rNC[[2]] + death.E.rNC[[2]] + death.A.rNC[[2]] + death.I.rNC[[2]] + death.R.rNC[[2]] +
    death.I.rC[[2]] + death.R.rC[[2]]
  total <- total-mortality


  inc_r <- expose.S.rNC
  inc_hcw <- expose.S.hcwNC + expose.S.hcwC

  cum_inc_r = (cum_inc_r + inc_r)
  cum_inc_hcw = (cum_inc_hcw + inc_hcw)
  cum_inc_community = (cum_inc_community + inc_community)


  final <- list("S.rNC"=S.rNC,
                "E.rNC"=E.rNC,
                "A.rNC"=A.rNC,
                "I.rNC"=I.rNC,
                "R.rNC"=R.rNC,
                "I.rC"=I.rC,
                "R.rC"=R.rC,
                "S.hcwNC"=S.hcwNC,
                "E.hcwNC"=E.hcwNC,
                "A.hcwNC"=A.hcwNC,
                "I.hcwNC"=I.hcwNC,
                "R.hcwNC"=R.hcwNC,
                "S.hcwC"=S.hcwC,
                "E.hcwC"=E.hcwC,
                "A.hcwC"=A.hcwC,
                "I.hcwC"=I.hcwC,
                "R.hcwC"=R.hcwC,
                "I.hcwH"=I.hcwH,
                "S.rNC3"=S.rNC3,
                "S.hcwNC3"=S.hcwNC3,
                "S.hcwC3"=S.hcwC3,
                "death.r"=death.r,
                "death.rC"=death.rC,
                "inc_r"=inc_r,
                "inc_hcw"=inc_hcw,
                "cum_inc_r"=cum_inc_r,
                "cum_inc_hcw"=cum_inc_hcw,
                "inc_vax1_r"=inc_vax1_r,
                "inc_vax2_r"=inc_vax2_r,
                "inc_vax1_hcw"=inc_vax1_hcw,
                "inc_vax2_hcw"=inc_vax2_hcw,
                "cum_inc_community"=cum_inc_community,
                "mortality"=mortality,
                "total"=total,
                "total_hcw"=total_hcw,
                "family"=family,  
                "vax_count1_r"=vax_count1_r,
                "vax_count2_r"=vax_count2_r, 
                "vax_count1_hcw"=vax_count1_hcw,
                "vax_count2_hcw"=vax_count2_hcw, 
                "asymptomatic_inc_r"=asymptomatic_inc_r,
                "symptomatic_inc_r"=symptomatic_inc_r,
                "asymptomatic_inc_hcw"=asymptomatic_inc_hcw,
                "symptomatic_inc_hcw"=symptomatic_inc_hcw)

  return(final)

}
inits=c(families=2000,
        S.rNC.init = 0,
        E.rNC.init = 33,
        A.rNC.init = 15,
        I.rNC.init = 2,
        R.rNC.init = 0,
        I.rC.init  = 0,
        R.rC.init = 0,
        S.hcwNC.init = 199,
        E.hcwNC.init = 0,
        A.hcwNC.init = 0,
        I.hcwNC.init = 0,
        R.hcwNC.init = 0,
        S.hcwC.init = 1,
        E.hcwC.init = 0,
        A.hcwC.init = 0,
        I.hcwC.init = 0,
        R.hcwC.init = 0,
        I.hcwH.init = 0, 
        group_1p=0.2539,
        group_2p=0.2968,
        group_3p=0.2099,
        group_4p=0.1317,
        group_5p=0.1076
        )

num_staff <- sum(inits[8:18])
ratio <- sum(inits[8:18])/sum(inits[1:7]) #hcw/res人数的比例


nsim <-50
j <-1

### run simulations ----
parms <- list(beta=0.0466,    ## beta
           beta.s=0.0466, ## beta for interactions between staff and residents
           beta.rm=(1-(1-0.0466)**10),  ## beta for roommates
           sigma1=2,      ## incubation period shortest
           sigma2=4,      ## incubation period longest
           gamma1=6,     ## infectious period
           gamma2=10,
           I.C=0,    ## probability of infection from community
           k.HH=2,       ## n contacts between staff  2
           k.RH=6,       ## n staff contacted by each resident, per day 6
           k.HR=2,       ## n residents contacted by each staff member, per day 2
           k.RR=0,       ## n daily contacts between residents besides roommates 0
           alpha.hcw = 0.907, # proportion hcw asymptomatic
           alpha.r = 0.907,   # proportion residents asymptomatic
           id.I= 2,         # duration of pre-symptomatic transmission
           ppe=0.05,        # reduction in beta with ppe
           ppe_r=0.82,
           prop_rhcwR=0.9,    # probability replacement hcw recovered
           mu.C=0.000008,       # COVID mortality
           mu.NC=0.000019,   # regular mortality
           int.r=1,         # 1 if go to new_recovered.rNC, 0 if go to new_recovered.rC
           int.rooms=1,     # if 1 put susceptible and recovered together in NC
           int.hcw=1,       # 1 if recovered work with NC and 0 if work with C
           VL.PCR_threshold=3,  # threshold for detectable VL (put arbitrary # in for now)
           VL.Antigen_threshold=5,
           test_PCR_delay=0, # staff test delay (days)
           test_Antigen_delay=0,   # resident test delay (days)
           int_time=1,      
           I.C_time=0,
           VE_s1 = 0.091,
           VE_p1 = 0.269,
           VE_i1 = 0.0,
           VE_s2 = 0.059,
           VE_p2 = 0.173,
           VE_i2 = 0.0,
           VE_s3 = 0.17,
           VE_p3 = 0.465,
           VE_i3 = 0.106,
           VE_s4 = 0.138,
           VE_p4 = 0.378,
           VE_i4 = 0.0,
           full_vax_coverage_r = 0.35,
           boost_vax_coverage_r = 0.55,
           refusal_r = 0.1, 
           full_vax_coverage_hcw = 0.35,
           boost_vax_coverage_hcw = 0.55,
           refusal_hcw = 0.1, 
           shortage = FALSE,
           shortage_threshold = 0.80*ratio,
           comm.vax = 0.9)

t_step <- 1

dt <- seq(0,180,t_step)

## loop over intervention strategies
intervention_scenarios <- as.data.frame(matrix(nrow=4, ncol=3))
colnames(intervention_scenarios) <- c("int.r", "int.rooms", "int.hcw")
rownames(intervention_scenarios) <- c("A) No intervention", "B) Resident intervention", "C) HCW intervention", "D) Both interventions")

intervention_scenarios[,1] <- c(0,1,0,1)
intervention_scenarios[,2] <- c(0,1,0,1)
intervention_scenarios[,3] <- c(0,0,1,1)


res_master <- NULL
VL_master <- NULL

I.C_list <- c(0.005, 0.002, 0.001, 0.0002, 0.00002,0.00001, 0)

ic=7
v=4      
s=1      
      parms[["I.C"]] <- I.C_list[ic]
      

      Intervention <- rownames(intervention_scenarios)[v]
      
      parms[["int.r"]] <- intervention_scenarios[Intervention, "int.r"]
      parms[["int.rooms"]] <- intervention_scenarios[Intervention, "int.rooms"]
      parms[["int.hcw"]] <- intervention_scenarios[Intervention, "int.hcw"]
      

      for (sim in 1:nsim){
        cat(sim,s,v,"\n")
        Ns <- initialize(inits,parms,dt)
        res <- as.data.frame(matrix(nrow=length(dt),ncol=length(Ns)))
        VLs <- NULL

        for(i in 1:length(dt)){

           cat(i,"\n")
           debug <- Ns
          final <- stochastic_NH(parms,Ns,t_step,(i-1)*t_step,num_staff)
          Ns <- final
          res[i,] <- c(i, nrow(Ns[["S.rNC"]]), nrow(Ns[["E.rNC"]]), nrow(Ns[["A.rNC"]]), nrow(Ns[["I.rNC"]]), nrow(Ns[["R.rNC"]]),
                       nrow(Ns[["I.rC"]]), nrow(Ns[["R.rC"]]),
                       nrow(Ns[["S.hcwNC"]]), nrow(Ns[["E.hcwNC"]]), nrow(Ns[["A.hcwNC"]]), nrow(Ns[["I.hcwNC"]]), nrow(Ns[["R.hcwNC"]]),
                       nrow(Ns[["S.hcwC"]]), nrow(Ns[["E.hcwC"]]), nrow(Ns[["A.hcwC"]]), nrow(Ns[["I.hcwC"]]), nrow(Ns[["R.hcwC"]]),
                       nrow(Ns[["I.hcwH"]]), nrow(Ns[["S.rNC3"]]),nrow(Ns[["S.hcwNC3"]]),nrow(Ns[["S.hcwC3"]]),
                       nrow(Ns[["death.r"]]),nrow(Ns[["death.rC"]]),Ns[["inc_r"]],  Ns[["inc_hcw"]],
                       Ns[["cum_inc_r"]], Ns[["cum_inc_hcw"]], 
                       Ns[["inc_vax1_r"]], Ns[["inc_vax2_r"]],Ns[["inc_vax1_hcw"]], Ns[["inc_vax2_hcw"]],
                       Ns[["cum_inc_community"]], Ns[["mortality"]],
                       Ns[["total"]],Ns[["total_hcw"]],Ns[["vax_count1_r"]],Ns[["vax_count2_r"]], Ns[["vax_count1_hcw"]],Ns[["vax_count2_hcw"]],
                       Ns[["asymptomatic_inc_r"]], Ns[["symptomatic_inc_r"]], Ns[["asymptomatic_inc_hcw"]], Ns[["symptomatic_inc_hcw"]])

          staff <- sum(res[i,c(9:18)])
          residents <- sum(res[i,2:8])
          
        }

        res %>%
          add_column("Sim"=sim) %>%
          add_column("Intervention"=Intervention) %>%
          add_column("I.C" = parms[["I.C"]]) %>%
          bind_rows(res_master) -> res_master

      }


colnames(res_master) <- c("time", "S.rNC", "E.rNC", "A.rNC", "I.rNC", "R.rNC",
                          "I.rC", "R.rC",
                          "S.hcwNC", "E.hcwNC", "A.hcwNC", "I.hcwNC", "R.hcwNC",
                          "S.hcwC", "E.hcwC", "A.hcwC", "I.hcwC", "R.hcwC", "I.hcwH","S.rNC3","S.hcwNC3","S.hcwC3", 
                          "death.r","death.rC",
                          "inc_r", "inc_hcw", "cum_inc_r", "cum_inc_hcw", 
                          "inc_vax1_r", "inc_vax2_r","inc_vax1_hcw", "inc_vax2_hcw",
                          "cum_inc_community", "mortality",
                          "total","total_hcw",
                          "vax_count1_r", "vax_count2_r", "vax_count1_hcw", "vax_count2_hcw",
                          "asymptomatic_inc_r", "symptomatic_inc_r",
                          "asymptomatic_inc_hcw", "symptomatic_inc_hcw",
                          "Sim","Intervention", 
                          "I.C")

write.csv(res_master, paste0(j,"pcr-2days-优化50-seed10031.csv"))
#


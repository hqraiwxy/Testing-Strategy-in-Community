initialize <- function(inits,params,dt){
  families = inits["families"]
  S.rNC.init = inits["S.rNC.init"]
  E.rNC.init = inits["E.rNC.init"]
  A.rNC.init = inits["A.rNC.init"]
  I.rNC.init = inits["I.rNC.init"]
  R.rNC.init = inits["R.rNC.init"]
  I.rC.init  = inits["I.rC.init"]
  R.rC.init = inits["R.rC.init"]
  S.hcwNC.init = inits["S.hcwNC.init"]
  E.hcwNC.init = inits["E.hcwNC.init"]
  A.hcwNC.init = inits["A.hcwNC.init"]
  I.hcwNC.init = inits["I.hcwNC.init"]
  R.hcwNC.init = inits["R.hcwNC.init"]
  S.hcwC.init = inits["S.hcwC.init"]
  E.hcwC.init = inits["E.hcwC.init"]
  A.hcwC.init = inits["A.hcwC.init"]
  I.hcwC.init = inits["I.hcwC.init"]
  R.hcwC.init = inits["R.hcwC.init"]
  I.hcwH.init = inits["I.hcwH.init"]
  I.hcwH.init = inits["I.hcwH.init"]
  group_1p=round(inits["group_1p"]*families)
  group_2p=round(inits["group_2p"]*families)
  group_3p=round(inits["group_3p"]*families)
  group_4p=round(inits["group_4p"]*families)
  group_5p=round(inits["group_5p"]*families)
  
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
  
  population=group_1p*1+group_2p*2+group_3p*3+group_4p*4+group_5p*5
  S.rNC.init=population-( E.rNC.init +  A.rNC.init +  I.rNC.init +  R.rNC.init +  I.rC.init +  R.rC.init)
  
  total <- 0
  if (S.rNC.init > 0){

    S.rNC=as.data.frame(cbind("ID"=1:S.rNC.init,"VL"=rep(NA,S.rNC.init)  ) )
    
    S.rNC %>%
      mutate(V_doses = sample(c(3,2,0),
                             size=S.rNC.init,
                             prob=c(parms[["boost_vax_coverage_r"]],parms[["full_vax_coverage_r"]],parms[["refusal_r"]]),
                             replace=TRUE),
             V_t = case_when(V_doses==2 ~ round(runif(1,1,15)),
                             V_doses==3 ~ round(runif(1,1,9)),
                             V_doses==0 ~ 1000,
                             TRUE ~ V_doses),
             VE_i = case_when(V_doses==2 & V_t <=6 ~ VE_i1, 
                              V_doses==2 & V_t > 6 ~ VE_i2,
                              V_doses==3 & V_t <=6 ~ VE_i3,
                              V_doses==3 & V_t > 6 ~ VE_i4,
                              V_doses==0 ~ 0
                              ),
             VE_s = case_when(V_doses==2 & V_t <=6 ~ VE_s1, 
                              V_doses==2 & V_t > 6 ~ VE_s2,
                              V_doses==3 & V_t <=6 ~ VE_s3,
                              V_doses==3 & V_t > 6 ~ VE_s4,
                              V_doses==0 ~ 0
                              ),
             VE_p = case_when(V_doses==2 & V_t <=6 ~ VE_p1, # Progression to symptoms
                              V_doses==2 & V_t > 6 ~ VE_p2,
                              V_doses==3 & V_t <=6 ~ VE_p3,
                              V_doses==3 & V_t > 6 ~ VE_p4,
                              V_doses==0 ~ 0
                              )
      ) -> S.rNC
    
    total <- total + S.rNC.init
  } else{
    S.rNC=as.data.frame(cbind("ID"=as.numeric(),"VL"=as.numeric(), 
                              "V_doses" = as.numeric(), #"V_refused" = as.numeric(),  
                              "V_t"=as.numeric(),                        
                              "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric())) 
  }
  
  if (E.rNC.init > 0){
    E.rNC=as.data.frame(cbind("ID"=(total+1):(total + E.rNC.init),"VL"=rep(0,E.rNC.init), 
                              #"V_doses" = rep(0,E.rNC.init), "V_refused" = rep(NA,E.rNC.init),
                              "Inc.pd"=round(runif(E.rNC.init, parms[["sigma1"]], parms[["sigma2"]])),
                              "Days"=rep(0,E.rNC.init)))
    E.rNC %>%
      mutate(
        V_doses = sample(c(3,2,0),
                         size=E.rNC.init,
                         prob=c(parms[["boost_vax_coverage_r"]],parms[["full_vax_coverage_r"]],parms[["refusal_r"]]),
                         replace=TRUE),
        V_t = case_when(V_doses==2 ~ round(runif(1,1,15)),
                        V_doses==3 ~ round(runif(1,1,9)),
                        V_doses==0 ~ 1000,
                        TRUE ~ V_doses),
        VE_i = case_when(V_doses==2 & V_t <=6 ~ VE_i1, 
                         V_doses==2 & V_t > 6 ~ VE_i2,
                         V_doses==3 & V_t <=6 ~ VE_i3,
                         V_doses==3 & V_t > 6 ~ VE_i4,
                         V_doses==0 ~ 0
        ),
        VE_s = case_when(V_doses==2 & V_t <=6 ~ VE_s1, # susceptibility (i.e. becoming infected)
                         V_doses==2 & V_t > 6 ~ VE_s2,
                         V_doses==3 & V_t <=6 ~ VE_s3,
                         V_doses==3 & V_t > 6 ~ VE_s4,
                         V_doses==0 ~ 0
        ),
        VE_p = case_when(V_doses==2 & V_t <=6 ~ VE_p1, # Progression to symptoms
                         V_doses==2 & V_t > 6 ~ VE_p2,
                         V_doses==3 & V_t <=6 ~ VE_p3,
                         V_doses==3 & V_t > 6 ~ VE_p4,
                         V_doses==0 ~ 0
        )
      ) -> E.rNC

    total <- total + E.rNC.init
  } else{
    E.rNC=as.data.frame(cbind("ID"=as.numeric(),"VL"=as.numeric(),
                              "V_doses" = as.numeric(), #"V_refused" = as.numeric(),  
                              "V_t"=as.numeric(), #"V2_t"=as.numeric(),                        
                              "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric(),
                              "Inc.pd"=as.numeric(),"Days"=as.numeric())) 
  }
  
  if (A.rNC.init > 0){
    A.rNC=as.data.frame(cbind("ID"=(total+1):(total+ A.rNC.init), "VL_rise" = round(runif(A.rNC.init, 1, runif(1,4.7,6.2))),
                              "removal.pd"=rep(NA, A.rNC.init),
                              #"VL"=abs(rnorm(A.rNC.init, 8, 1)*0.5), "VL_waning"=rep(8/runif(1, 15, 21),A.rNC.init), 
                              #"VL"=runif(A.rNC.init,1, 8), "VL_waning"=rep(0.35,A.rNC.init),
                              "VL"=abs(rnorm(A.rNC.init, 8, 1)*0.5),"VL_waning"=rep(0.35,A.rNC.init),
                              "Inc.pd"=rep(0,A.rNC.init), #新增
                              "Inf.pd"=round(runif(A.rNC.init, parms[["gamma1"]], parms[["gamma2"]])),
                              "Inf.days"=rep(0,A.rNC.init), # can change if we want to assume they are partway into inf period
                              "ID.days"=rep(0,A.rNC.init),
                              "Infectiousness"=rep(1,A.rNC.init)))
    
    A.rNC %>%
      mutate(
        V_doses = sample(c(3,2,0),
                         size=A.rNC.init,
                         prob=c(parms[["boost_vax_coverage_r"]],parms[["full_vax_coverage_r"]],parms[["refusal_r"]]),
                         replace=TRUE),
        V_t = case_when(V_doses==2 ~ round(runif(1,1,15)),
                        V_doses==3 ~ round(runif(1,1,9)),
                        V_doses==0 ~ 1000,
                        TRUE ~ V_doses),
        VE_i = case_when(V_doses==2 & V_t <=6 ~ VE_i1, 
                         V_doses==2 & V_t > 6 ~ VE_i2,
                         V_doses==3 & V_t <=6 ~ VE_i3,
                         V_doses==3 & V_t > 6 ~ VE_i4,
                         V_doses==0 ~ 0
        ),
        VE_s = case_when(V_doses==2 & V_t <=6 ~ VE_s1, # susceptibility (i.e. becoming infected)
                         V_doses==2 & V_t > 6 ~ VE_s2,
                         V_doses==3 & V_t <=6 ~ VE_s3,
                         V_doses==3 & V_t > 6 ~ VE_s4,
                         V_doses==0 ~ 0
        ),
        VE_p = case_when(V_doses==2 & V_t <=6 ~ VE_p1, # Progression to symptoms
                         V_doses==2 & V_t > 6 ~ VE_p2,
                         V_doses==3 & V_t <=6 ~ VE_p3,
                         V_doses==3 & V_t > 6 ~ VE_p4,
                         V_doses==0 ~ 0
        )
      ) -> A.rNC
    
    total <- total + A.rNC.init
  } else{
    A.rNC=as.data.frame(cbind("ID"=as.numeric(),"VL_rise"=as.numeric(), 
                              "removal.pd"=as.numeric(), "VL"=as.numeric(), "VL_waning"=as.numeric(),
                              "V_doses" = as.numeric(),
                              "V_t"=as.numeric(),                         
                              "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric(),
                              "Inc.pd"=as.numeric(),
                              "Inf.pd"=as.numeric(),"Inf.days"=as.numeric(),
                              "ID.days"=as.numeric(),
                              "Infectiousness"=as.numeric())) 
  }
  
  if (I.rNC.init > 0){
    I.rNC=as.data.frame(cbind("ID"=(total+1):(total+ I.rNC.init), "VL_rise" = round(runif(I.rNC.init, 1, runif(1,4.7,6.2))),
                              "removal.pd"=rep(NA, I.rNC.init),
                              #"VL"=abs(rnorm(I.rNC.init, 8, 1)*0.5), "VL_waning"=rep(8/runif(1, 15, 21),I.rNC.init),
                              #"VL"=runif(I.rNC.init,1, 8), "VL_waning"=rep(0.35,I.rNC.init),
                              "VL"=abs(rnorm(I.rNC.init, 8, 1)*0.5),"VL_waning"=rep(0.35,I.rNC.init),
                              "Inc.pd"=rep(0,I.rNC.init),
                              "Inf.pd"=round(runif(I.rNC.init, parms[["gamma1"]], parms[["gamma2"]])),
                              "Inf.days"=rep(0,I.rNC.init), # can change if we want to assume they are partway into inf period
                              "ID.days"=rep(0,I.rNC.init),
                              "Infectiousness"=rep(1,I.rNC.init)))
    
    I.rNC %>%
      mutate(
        V_doses = sample(c(3,2,0),
                         size=I.rNC.init,
                         prob=c(parms[["boost_vax_coverage_r"]],parms[["full_vax_coverage_r"]],parms[["refusal_r"]]),
                         replace=TRUE),
        V_t = case_when(V_doses==2 ~ round(runif(1,1,15)),
                        V_doses==3 ~ round(runif(1,1,9)),
                        V_doses==0 ~ 1000,
                        TRUE ~ V_doses),
        VE_i = case_when(V_doses==2 & V_t <=6 ~ VE_i1, 
                         V_doses==2 & V_t > 6 ~ VE_i2,
                         V_doses==3 & V_t <=6 ~ VE_i3,
                         V_doses==3 & V_t > 6 ~ VE_i4,
                         V_doses==0 ~ 0
        ),
        VE_s = case_when(V_doses==2 & V_t <=6 ~ VE_s1, # susceptibility (i.e. becoming infected)
                         V_doses==2 & V_t > 6 ~ VE_s2,
                         V_doses==3 & V_t <=6 ~ VE_s3,
                         V_doses==3 & V_t > 6 ~ VE_s4,
                         V_doses==0 ~ 0
        ),
        VE_p = case_when(V_doses==2 & V_t <=6 ~ VE_p1, # Progression to symptoms
                         V_doses==2 & V_t > 6 ~ VE_p2,
                         V_doses==3 & V_t <=6 ~ VE_p3,
                         V_doses==3 & V_t > 6 ~ VE_p4,
                         V_doses==0 ~ 0
        )
      ) -> I.rNC
    
    total <- total + I.rNC.init
  } else{
    I.rNC=as.data.frame(cbind("ID"=as.numeric(),"VL_rise"=as.numeric(), 
                              "removal.pd"=as.numeric(), "VL"=as.numeric(), "VL_waning"=as.numeric(), 
                              "V_doses" = as.numeric(),   
                              "V_t"=as.numeric(), #"V2_t"=as.numeric(),                        
                              "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric(),
                              "Inc.pd"=as.numeric(),
                              "Inf.pd"=as.numeric(),"Inf.days"=as.numeric(),
                              "ID.days"=as.numeric(),
                              "Infectiousness"=as.numeric()))
  }
  
  if (R.rNC.init > 0){
    R.rNC=as.data.frame(cbind("ID"=(total+1):(total + R.rNC.init),
                              "VL"=rep(0,R.rNC.init), 
                              "VL_waning"=rep(0,R.rNC.init),
                              "Inc.pd"=rep(0,R.rNC.init), 
                              "Inf.pd"=rep(0,R.rNC.init),
                              "Rec.days"=rep(0,R.rNC.init)))
    
    R.rNC %>%
      mutate(
        V_doses = sample(c(3,2,0),
                         size=R.rNC.init,
                         prob=c(parms[["boost_vax_coverage_r"]],parms[["full_vax_coverage_r"]],parms[["refusal_r"]]),
                         replace=TRUE),
        V_t = case_when(V_doses==2 ~ round(runif(1,1,15)),
                        V_doses==3 ~ round(runif(1,1,9)),
                        V_doses==0 ~ 1000,
                        TRUE ~ V_doses),
        VE_i = case_when(V_doses==2 & V_t <=6 ~ VE_i1, 
                         V_doses==2 & V_t > 6 ~ VE_i2,
                         V_doses==3 & V_t <=6 ~ VE_i3,
                         V_doses==3 & V_t > 6 ~ VE_i4,
                         V_doses==0 ~ 0
        ),
        VE_s = case_when(V_doses==2 & V_t <=6 ~ VE_s1, # susceptibility (i.e. becoming infected)
                         V_doses==2 & V_t > 6 ~ VE_s2,
                         V_doses==3 & V_t <=6 ~ VE_s3,
                         V_doses==3 & V_t > 6 ~ VE_s4,
                         V_doses==0 ~ 0
        ),
        VE_p = case_when(V_doses==2 & V_t <=6 ~ VE_p1, # Progression to symptoms
                         V_doses==2 & V_t > 6 ~ VE_p2,
                         V_doses==3 & V_t <=6 ~ VE_p3,
                         V_doses==3 & V_t > 6 ~ VE_p4,
                         V_doses==0 ~ 0
        )
      ) -> R.rNC
    
    total <- total + R.rNC.init
  } else{
    R.rNC=as.data.frame(cbind("ID"=as.numeric(),"VL"=as.numeric(), "VL_waning"=as.numeric(), 
                              #"arrival" = as.numeric(), "departure"= as.numeric(),  
                              "V_doses" = as.numeric(), #"V_refused" = as.numeric(),  
                              "V_t"=as.numeric(),# "V2_t"=as.numeric(),                        
                              "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric(),
                              "Inc.pd"=as.numeric(), "Inf.pd"=as.numeric(),
                              "Rec.days"=as.numeric()))
  }
  
  if (I.rC.init > 0){
    I.rC=as.data.frame(cbind("ID"=(total+1):(total+ I.rC.init),"VL_rise" = round(runif(I.rC.init, 1, runif(1,4.7,6.2))),
                             "removal.pd"=rep(NA, I.rC.init),
                             #"VL"=abs(rnorm(I.rC.init, 8, 1)*0.5), "VL_waning"=rep(8/runif(1, 15, 21),I.rC.init),
                             #"VL"=runif(I.rC.init,1, 8), "VL_waning"=rep(0.35,I.rC.init),
                             "VL"=abs(rnorm(I.rC.init, 8, 1)*0.5),"VL_waning"=rep(0.35,I.rC.init),
                             "Inc.pd"=rep(0,I.rC.init),
                             "Inf.pd"=round(runif(I.rC.init, parms[["gamma1"]], parms[["gamma2"]])),
                             "Inf.days"=rep(0,I.rC.init),
                             "ID.days"=rep(0,I.rC.init), # can change if we want to assume they are partway into inf period
                             "Infectiousness"=rep(1,I.rC.init), "asympt"=rbinom(I.rC.init, 1, parms[["alpha.r"]])))
    
    I.rC %>%
      mutate(
        V_doses = sample(c(3,2,0),
                         size=I.rC.init,
                         prob=c(parms[["boost_vax_coverage_r"]],parms[["full_vax_coverage_r"]],parms[["refusal_r"]]),
                         replace=TRUE),
        V_t = case_when(V_doses==2 ~ round(runif(1,1,15)),
                        V_doses==3 ~ round(runif(1,1,9)),
                        V_doses==0 ~ 1000,
                        TRUE ~ V_doses),
        VE_i = case_when(V_doses==2 & V_t <=6 ~ VE_i1, 
                         V_doses==2 & V_t > 6 ~ VE_i2,
                         V_doses==3 & V_t <=6 ~ VE_i3,
                         V_doses==3 & V_t > 6 ~ VE_i4,
                         V_doses==0 ~ 0
        ),
        VE_s = case_when(V_doses==2 & V_t <=6 ~ VE_s1, # susceptibility (i.e. becoming infected)
                         V_doses==2 & V_t > 6 ~ VE_s2,
                         V_doses==3 & V_t <=6 ~ VE_s3,
                         V_doses==3 & V_t > 6 ~ VE_s4,
                         V_doses==0 ~ 0
        ),
        VE_p = case_when(V_doses==2 & V_t <=6 ~ VE_p1, # Progression to symptoms
                         V_doses==2 & V_t > 6 ~ VE_p2,
                         V_doses==3 & V_t <=6 ~ VE_p3,
                         V_doses==3 & V_t > 6 ~ VE_p4,
                         V_doses==0 ~ 0
        )
      ) -> I.rC
    
    total <- total + I.rC.init
  } else{
    I.rC=as.data.frame(cbind("ID"=as.numeric(),"VL_rise" = as.numeric(),
                             "removal.pd"=as.numeric(),"VL"=as.numeric(), "VL_waning"=as.numeric(), 
                             "V_doses" = as.numeric(), #"V_refused" = as.numeric(),  
                             "V_t"=as.numeric(),# "V2_t"=as.numeric(),                        
                             "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric(),
                             "Inc.pd"=as.numeric(),
                             "Inf.pd"=as.numeric(),"Inf.days"=as.numeric(),
                             "ID.days"=as.numeric(),
                             "Infectiousness"=as.numeric(), "asympt"=as.numeric())) 
    
  }
  
  if (R.rC.init > 0){
    R.rC=as.data.frame(cbind("ID"=(total+1):(total + R.rC.init),"VL"=rep(0,R.rC.init), "VL_waning"=rep(0,R.rC.init), 
                             
                             "Inc.pd"=rep(0,R.rC.init), 
                             "Inf.pd"=rep(0,R.rC.init),
                             "Rec.days"=rep(0,R.rC.init)))
    
    R.rC %>%
      mutate(
        V_doses = sample(c(3,2,0),
                         size=R.rC.init,
                         prob=c(parms[["boost_vax_coverage_r"]],parms[["full_vax_coverage_r"]],parms[["refusal_r"]]),
                         replace=TRUE),
        V_t = case_when(V_doses==2 ~ round(runif(1,1,15)),
                        V_doses==3 ~ round(runif(1,1,9)),
                        V_doses==0 ~ 1000,
                        TRUE ~ V_doses),
        VE_i = case_when(V_doses==2 & V_t <=6 ~ VE_i1, 
                         V_doses==2 & V_t > 6 ~ VE_i2,
                         V_doses==3 & V_t <=6 ~ VE_i3,
                         V_doses==3 & V_t > 6 ~ VE_i4,
                         V_doses==0 ~ 0
        ),
        VE_s = case_when(V_doses==2 & V_t <=6 ~ VE_s1, # susceptibility (i.e. becoming infected)
                         V_doses==2 & V_t > 6 ~ VE_s2,
                         V_doses==3 & V_t <=6 ~ VE_s3,
                         V_doses==3 & V_t > 6 ~ VE_s4,
                         V_doses==0 ~ 0
        ),
        VE_p = case_when(V_doses==2 & V_t <=6 ~ VE_p1, # Progression to symptoms
                         V_doses==2 & V_t > 6 ~ VE_p2,
                         V_doses==3 & V_t <=6 ~ VE_p3,
                         V_doses==3 & V_t > 6 ~ VE_p4,
                         V_doses==0 ~ 0
        )
      ) -> R.rC
    
    total <- total + R.rC.init
  } else{
    R.rC=as.data.frame(cbind("ID"=as.numeric(),"VL"=as.numeric(),"VL_waning"=as.numeric(), 
                             "V_doses" = as.numeric(), #"V_refused" = as.numeric(),  
                             "V_t"=as.numeric(),# "V2_t"=as.numeric(),                        
                             "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric(),
                             "Inc.pd"=as.numeric(), "Inf.pd"=as.numeric(),
                             "Rec.days"=as.numeric())) 
  }
  
  rNC <- c(S.rNC$ID, E.rNC$ID, A.rNC$ID, I.rNC$ID,R.rNC$ID)
  
  
  total_hcw <- 5999 
  if (S.hcwNC.init > 0){
    S.hcwNC=as.data.frame(cbind("ID"=(total_hcw + 1):(total_hcw + S.hcwNC.init),"VL"=rep(NA,S.hcwNC.init)   ))
    
    S.hcwNC %>%
      mutate(
        V_doses = sample(c(3,2,0),
                         size=S.hcwNC.init,
                         prob=c(parms[["boost_vax_coverage_r"]],parms[["full_vax_coverage_r"]],parms[["refusal_r"]]),
                         replace=TRUE),
        V_t = case_when(V_doses==2 ~ round(runif(1,1,15)),
                        V_doses==3 ~ round(runif(1,1,9)),
                        V_doses==0 ~ 1000,
                        TRUE ~ V_doses),
        VE_i = case_when(V_doses==2 & V_t <=6 ~ VE_i1, 
                         V_doses==2 & V_t > 6 ~ VE_i2,
                         V_doses==3 & V_t <=6 ~ VE_i3,
                         V_doses==3 & V_t > 6 ~ VE_i4,
                         V_doses==0 ~ 0
        ),
        VE_s = case_when(V_doses==2 & V_t <=6 ~ VE_s1, # susceptibility (i.e. becoming infected)
                         V_doses==2 & V_t > 6 ~ VE_s2,
                         V_doses==3 & V_t <=6 ~ VE_s3,
                         V_doses==3 & V_t > 6 ~ VE_s4,
                         V_doses==0 ~ 0
        ),
        VE_p = case_when(V_doses==2 & V_t <=6 ~ VE_p1, # Progression to symptoms
                         V_doses==2 & V_t > 6 ~ VE_p2,
                         V_doses==3 & V_t <=6 ~ VE_p3,
                         V_doses==3 & V_t > 6 ~ VE_p4,
                         V_doses==0 ~ 0
        )
      ) -> S.hcwNC
    
    total_hcw <- total_hcw + S.hcwNC.init
  } else{
    S.hcwNC=as.data.frame(cbind("ID"=as.numeric(),"VL"=as.numeric(),  ## this used to say infectiousness
                                "V_doses" = as.numeric(), #"V_refused" = as.numeric(),  
                                "V_t"=as.numeric(), #"V2_t"=as.numeric()))                        
                                "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric()))
  }
  
  if (E.hcwNC.init > 0){
    E.hcwNC=as.data.frame(cbind("ID"=(total_hcw+1):(total_hcw + E.hcwNC.init),"VL"=rep(0,E.hcwNC.init),
                               
                                "Inc.pd"=round(runif(E.hcwNC.init, parms[["sigma1"]], parms[["sigma2"]])),"Days"=rep(0,E.hcwNC.init)))
    
    E.hcwNC %>%
      mutate(
        V_doses = sample(c(3,2,0),
                         size=E.hcwNC.init,
                         prob=c(parms[["boost_vax_coverage_r"]],parms[["full_vax_coverage_r"]],parms[["refusal_r"]]),
                         replace=TRUE),
        V_t = case_when(V_doses==2 ~ round(runif(1,1,15)),
                        V_doses==3 ~ round(runif(1,1,9)),
                        V_doses==0 ~ 1000,
                        TRUE ~ V_doses),
        VE_i = case_when(V_doses==2 & V_t <=6 ~ VE_i1, 
                         V_doses==2 & V_t > 6 ~ VE_i2,
                         V_doses==3 & V_t <=6 ~ VE_i3,
                         V_doses==3 & V_t > 6 ~ VE_i4,
                         V_doses==0 ~ 0
        ),
        VE_s = case_when(V_doses==2 & V_t <=6 ~ VE_s1, # susceptibility (i.e. becoming infected)
                         V_doses==2 & V_t > 6 ~ VE_s2,
                         V_doses==3 & V_t <=6 ~ VE_s3,
                         V_doses==3 & V_t > 6 ~ VE_s4,
                         V_doses==0 ~ 0
        ),
        VE_p = case_when(V_doses==2 & V_t <=6 ~ VE_p1, # Progression to symptoms
                         V_doses==2 & V_t > 6 ~ VE_p2,
                         V_doses==3 & V_t <=6 ~ VE_p3,
                         V_doses==3 & V_t > 6 ~ VE_p4,
                         V_doses==0 ~ 0
        )
      ) -> E.hcwNC
    
    total_hcw <- total_hcw + E.hcwNC.init
  } else{
    E.hcwNC=as.data.frame(cbind("ID"=as.numeric(),"VL"=as.numeric(),
                                "V_doses" = as.numeric(), #"V_refused" = as.numeric(),  
                                "V_t"=as.numeric(), #"V2_t"=as.numeric(),                        
                                "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric(), 
                                "Inc.pd"=as.numeric(),"Days"=as.numeric())) 
  }
  
  if (A.hcwNC.init > 0){
    A.hcwNC=as.data.frame(cbind("ID"=(total_hcw+1):(total_hcw+ A.hcwNC.init),"VL_rise" = round(runif(A.hcwNC.init, 1, runif(1,4.7,6.2))),
                                "removal.pd"=rep(NA, A.hcwNC.init),
                                "VL"=abs(rnorm(A.hcwNC.init, 8, 1)*0.5),"VL_waning"=rep(0.35,A.hcwNC.init),
                                "Inc.pd"=rep(0,A.hcwNC.init),
                                "Inf.pd"=round(runif(A.hcwNC.init, parms[["gamma1"]], parms[["gamma2"]])),
                                "Inf.days"=rep(0,A.hcwNC.init), # can change if we want to assume they are partway into inf period
                                "ID.days"=rep(0,A.hcwNC.init),
                                "Infectiousness"=rep(1,A.hcwNC.init)))
    
    A.hcwNC %>%
      mutate(
        V_doses = sample(c(3,2,0),
                         size=A.hcwNC.init,
                         prob=c(parms[["boost_vax_coverage_r"]],parms[["full_vax_coverage_r"]],parms[["refusal_r"]]),
                         replace=TRUE),
        V_t = case_when(V_doses==2 ~ round(runif(1,1,15)),
                        V_doses==3 ~ round(runif(1,1,9)),
                        V_doses==0 ~ 1000,
                        TRUE ~ V_doses),
        VE_i = case_when(V_doses==2 & V_t <=6 ~ VE_i1, 
                         V_doses==2 & V_t > 6 ~ VE_i2,
                         V_doses==3 & V_t <=6 ~ VE_i3,
                         V_doses==3 & V_t > 6 ~ VE_i4,
                         V_doses==0 ~ 0
        ),
        VE_s = case_when(V_doses==2 & V_t <=6 ~ VE_s1, # susceptibility (i.e. becoming infected)
                         V_doses==2 & V_t > 6 ~ VE_s2,
                         V_doses==3 & V_t <=6 ~ VE_s3,
                         V_doses==3 & V_t > 6 ~ VE_s4,
                         V_doses==0 ~ 0
        ),
        VE_p = case_when(V_doses==2 & V_t <=6 ~ VE_p1, # Progression to symptoms
                         V_doses==2 & V_t > 6 ~ VE_p2,
                         V_doses==3 & V_t <=6 ~ VE_p3,
                         V_doses==3 & V_t > 6 ~ VE_p4,
                         V_doses==0 ~ 0
        )
      ) -> A.hcwNC
    
    total_hcw <- total_hcw + A.hcwNC.init
  } else{
    A.hcwNC=as.data.frame(cbind("ID"=as.numeric(),"VL_rise" = as.numeric(),
                                "removal.pd"=as.numeric(),"VL"=as.numeric(), "VL_waning"=as.numeric(),
                                "V_doses" = as.numeric(), #"V_refused" = as.numeric(),  
                                "V_t"=as.numeric(), #"V2_t"=as.numeric(),                        
                                "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric(), 
                                "Inc.pd"=as.numeric(),
                                "Inf.pd"=as.numeric(),"Inf.days"=as.numeric(),
                                "ID.days"=as.numeric(),
                                "Infectiousness"=as.numeric())) 
  }
  
  if (I.hcwNC.init > 0){
    I.hcwNC=as.data.frame(cbind("ID"=(total_hcw+1):(total_hcw + I.hcwNC.init),"VL_rise" = round(runif(I.hcwNC.init, 1, runif(1,4.7,6.2))),
                                "removal.pd"=rep(NA, I.hcwNC.init),
                                "VL"=abs(rnorm(I.hcwNC.init, 8, 1)*0.5),"VL_waning"=rep(0.35,I.hcwNC.init),
                                "Inc.pd"=rep(0,I.hcwNC.init),
                                "Inf.pd"=round(runif(I.hcwNC.init, parms[["gamma1"]], parms[["gamma2"]])),
                                "Inf.days"=rep(0,I.hcwNC.init), # can change if we want to assume they are partway into inf period
                                "ID.days"=rep(0,I.hcwNC.init),
                                "Infectiousness"=rep(1,I.hcwNC.init)))
    
    I.hcwNC %>%
      mutate(
        V_doses = sample(c(3,2,0),
                         size=I.hcwNC.init,
                         prob=c(parms[["boost_vax_coverage_r"]],parms[["full_vax_coverage_r"]],parms[["refusal_r"]]),
                         replace=TRUE),
        V_t = case_when(V_doses==2 ~ round(runif(1,1,15)),
                        V_doses==3 ~ round(runif(1,1,9)),
                        V_doses==0 ~ 1000,
                        TRUE ~ V_doses),
        VE_i = case_when(V_doses==2 & V_t <=6 ~ VE_i1, 
                         V_doses==2 & V_t > 6 ~ VE_i2,
                         V_doses==3 & V_t <=6 ~ VE_i3,
                         V_doses==3 & V_t > 6 ~ VE_i4,
                         V_doses==0 ~ 0
        ),
        VE_s = case_when(V_doses==2 & V_t <=6 ~ VE_s1, # susceptibility (i.e. becoming infected)
                         V_doses==2 & V_t > 6 ~ VE_s2,
                         V_doses==3 & V_t <=6 ~ VE_s3,
                         V_doses==3 & V_t > 6 ~ VE_s4,
                         V_doses==0 ~ 0
        ),
        VE_p = case_when(V_doses==2 & V_t <=6 ~ VE_p1, # Progression to symptoms
                         V_doses==2 & V_t > 6 ~ VE_p2,
                         V_doses==3 & V_t <=6 ~ VE_p3,
                         V_doses==3 & V_t > 6 ~ VE_p4,
                         V_doses==0 ~ 0
        )
      ) -> I.hcwNC
    
    total_hcw <- total_hcw + I.hcwNC.init
  } else{
    I.hcwNC=as.data.frame(cbind("ID"=as.numeric(),"VL_rise" = as.numeric(),
                                "removal.pd"=as.numeric(),"VL"=as.numeric(), "VL_waning"=as.numeric(),
                                "V_doses" = as.numeric(), #"V_refused" = as.numeric(),  
                                "V_t"=as.numeric(), #"V2_t"=as.numeric(),                        
                                "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric(),
                                "Inc.pd"=as.numeric(),
                                "Inf.pd"=as.numeric(),"Inf.days"=as.numeric(),
                                "ID.days"=as.numeric(),
                                "Infectiousness"=as.numeric())) 
  }
  
  if (R.hcwNC.init > 0){
    R.hcwNC=as.data.frame(cbind("ID"=(total_hcw+1):(total_hcw + R.hcwNC.init),"VL"=rep(0,R.hcwNC.init), "VL_waning"=rep(0,R.hcwNC.init),
                               
                                "Inc.pd"=rep(0,R.hcwNC.init), "Inf.pd"=rep(0,R.hcwNC.init), #新增
                                "Rec.days"=rep(0,R.hcwNC.init)))
    
    R.hcwNC %>%
      mutate(
        V_doses = sample(c(3,2,0),
                         size=R.hcwNC.init,
                         prob=c(parms[["boost_vax_coverage_r"]],parms[["full_vax_coverage_r"]],parms[["refusal_r"]]),
                         replace=TRUE),
        V_t = case_when(V_doses==2 ~ round(runif(1,1,15)),
                        V_doses==3 ~ round(runif(1,1,9)),
                        V_doses==0 ~ 1000,
                        TRUE ~ V_doses),
        VE_i = case_when(V_doses==2 & V_t <=6 ~ VE_i1, 
                         V_doses==2 & V_t > 6 ~ VE_i2,
                         V_doses==3 & V_t <=6 ~ VE_i3,
                         V_doses==3 & V_t > 6 ~ VE_i4,
                         V_doses==0 ~ 0
        ),
        VE_s = case_when(V_doses==2 & V_t <=6 ~ VE_s1, # susceptibility (i.e. becoming infected)
                         V_doses==2 & V_t > 6 ~ VE_s2,
                         V_doses==3 & V_t <=6 ~ VE_s3,
                         V_doses==3 & V_t > 6 ~ VE_s4,
                         V_doses==0 ~ 0
        ),
        VE_p = case_when(V_doses==2 & V_t <=6 ~ VE_p1, # Progression to symptoms
                         V_doses==2 & V_t > 6 ~ VE_p2,
                         V_doses==3 & V_t <=6 ~ VE_p3,
                         V_doses==3 & V_t > 6 ~ VE_p4,
                         V_doses==0 ~ 0
        )
      ) -> R.hcwNC
    
    total_hcw <- total_hcw + R.hcwNC.init
  } else{
    R.hcwNC=as.data.frame(cbind("ID"=as.numeric(),"VL"=as.numeric(),"VL_waning"=as.numeric(),
                                "V_doses" = as.numeric(), #"V_refused" = as.numeric(),  
                                "V_t"=as.numeric(), #"V2_t"=as.numeric(),                        
                                "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric(),
                                "Inc.pd"=as.numeric(), "Inf.pd"=as.numeric(), 
                                "Rec.days"=as.numeric())) 
  }
  
  if (S.hcwC.init > 0){
    S.hcwC=as.data.frame(cbind("ID"=(total_hcw+1):(total_hcw + S.hcwC.init),"VL"=rep(NA,S.hcwC.init)  )) # can adjust if want to vary this
    
    S.hcwC %>%
      mutate(
        V_doses = sample(c(3,2,0),
                         size=S.hcwC.init,
                         prob=c(parms[["boost_vax_coverage_r"]],parms[["full_vax_coverage_r"]],parms[["refusal_r"]]),
                         replace=TRUE),
        V_t = case_when(V_doses==2 ~ round(runif(1,1,15)),
                        V_doses==3 ~ round(runif(1,1,9)),
                        V_doses==0 ~ 1000,
                        TRUE ~ V_doses),
        VE_i = case_when(V_doses==2 & V_t <=6 ~ VE_i1, 
                         V_doses==2 & V_t > 6 ~ VE_i2,
                         V_doses==3 & V_t <=6 ~ VE_i3,
                         V_doses==3 & V_t > 6 ~ VE_i4,
                         V_doses==0 ~ 0
        ),
        VE_s = case_when(V_doses==2 & V_t <=6 ~ VE_s1, # susceptibility (i.e. becoming infected)
                         V_doses==2 & V_t > 6 ~ VE_s2,
                         V_doses==3 & V_t <=6 ~ VE_s3,
                         V_doses==3 & V_t > 6 ~ VE_s4,
                         V_doses==0 ~ 0
        ),
        VE_p = case_when(V_doses==2 & V_t <=6 ~ VE_p1, # Progression to symptoms
                         V_doses==2 & V_t > 6 ~ VE_p2,
                         V_doses==3 & V_t <=6 ~ VE_p3,
                         V_doses==3 & V_t > 6 ~ VE_p4,
                         V_doses==0 ~ 0
        )
      ) -> S.hcwC
    
    total_hcw <- total_hcw + S.hcwC.init
  } else{
    S.hcwC=as.data.frame(cbind("ID"=as.numeric(),"VL"=as.numeric(), 
                               "V_doses" = as.numeric(), #"V_refused" = as.numeric(),  
                               "V_t"=as.numeric(), #"V2_t"=as.numeric()))                        
                               "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric())) 
  }
  
  if (E.hcwC.init > 0){
    E.hcwC=as.data.frame(cbind("ID"=(total_hcw+1):(total_hcw + E.hcwC.init),"VL"=rep(0,E.hcwC.init), 
                               
                               "Inc.pd"=round(runif(E.hcwC.init, parms[["sigma1"]], parms[["sigma2"]])),"Days"=rep(0,E.hcwC.init)))
    
    E.hcwC %>%
      mutate(
        V_doses = sample(c(3,2,0),
                         size=E.hcwC.init,
                         prob=c(parms[["boost_vax_coverage_r"]],parms[["full_vax_coverage_r"]],parms[["refusal_r"]]),
                         replace=TRUE),
        V_t = case_when(V_doses==2 ~ round(runif(1,1,15)),
                        V_doses==3 ~ round(runif(1,1,9)),
                        V_doses==0 ~ 1000,
                        TRUE ~ V_doses),
        VE_i = case_when(V_doses==2 & V_t <=6 ~ VE_i1, 
                         V_doses==2 & V_t > 6 ~ VE_i2,
                         V_doses==3 & V_t <=6 ~ VE_i3,
                         V_doses==3 & V_t > 6 ~ VE_i4,
                         V_doses==0 ~ 0
        ),
        VE_s = case_when(V_doses==2 & V_t <=6 ~ VE_s1, # susceptibility (i.e. becoming infected)
                         V_doses==2 & V_t > 6 ~ VE_s2,
                         V_doses==3 & V_t <=6 ~ VE_s3,
                         V_doses==3 & V_t > 6 ~ VE_s4,
                         V_doses==0 ~ 0
        ),
        VE_p = case_when(V_doses==2 & V_t <=6 ~ VE_p1, # Progression to symptoms
                         V_doses==2 & V_t > 6 ~ VE_p2,
                         V_doses==3 & V_t <=6 ~ VE_p3,
                         V_doses==3 & V_t > 6 ~ VE_p4,
                         V_doses==0 ~ 0
        )
      ) -> E.hcwC
    
    total_hcw <- total_hcw + E.hcwC.init
  } else{
    E.hcwC=as.data.frame(cbind("ID"=as.numeric(),"VL"=as.numeric(),
                               "V_doses" = as.numeric(), #"V_refused" = as.numeric(),  
                               "V_t"=as.numeric(), #"V2_t"=as.numeric(),                        
                               "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric(), 
                               "Inc.pd"=as.numeric(),"Days"=as.numeric())) 
  }
  
  if (A.hcwC.init > 0){
    A.hcwC=as.data.frame(cbind("ID"=(total_hcw+1):(total_hcw+ A.hcwC.init), "VL_rise" = round(runif(A.hcwC.init, 1, runif(1,4.7,6.2))),
                               "removal.pd"=rep(NA, A.hcwC.init),#"VL"=abs(rnorm(A.hcwC.init, 8, 1)*0.5), "VL_waning"=rep(8/runif(1, 15, 21),A.hcwC.init),
                               #"VL"=runif(A.hcwC.init,1, 8), "VL_waning"=rep(0.35,A.hcwC.init),
                               "VL"=abs(rnorm(A.hcwC.init, 8, 1)*0.5),"VL_waning"=rep(0.35,A.hcwC.init),
                               "Inc.pd"=rep(0,A.hcwC.init),
                               "Inf.pd"=round(runif(A.hcwC.init, parms[["gamma1"]], parms[["gamma2"]])),
                               "Inf.days"=rep(0,A.hcwC.init), # can change if we want to assume they are partway into inf period
                               "ID.days"=rep(0,A.hcwC.init),
                               "Infectiousness"=rep(1,A.hcwC.init)))
    A.hcwC %>%
      mutate(
        V_doses = sample(c(3,2,0),
                         size=A.hcwC.init,
                         prob=c(parms[["boost_vax_coverage_r"]],parms[["full_vax_coverage_r"]],parms[["refusal_r"]]),
                         replace=TRUE),
        V_t = case_when(V_doses==2 ~ round(runif(1,1,15)),
                        V_doses==3 ~ round(runif(1,1,9)),
                        V_doses==0 ~ 1000,
                        TRUE ~ V_doses),
        VE_i = case_when(V_doses==2 & V_t <=6 ~ VE_i1,
                         V_doses==2 & V_t > 6 ~ VE_i2,
                         V_doses==3 & V_t <=6 ~ VE_i3,
                         V_doses==3 & V_t > 6 ~ VE_i4,
                         V_doses==0 ~ 0
        ),
        VE_s = case_when(V_doses==2 & V_t <=6 ~ VE_s1, # susceptibility (i.e. becoming infected)
                         V_doses==2 & V_t > 6 ~ VE_s2,
                         V_doses==3 & V_t <=6 ~ VE_s3,
                         V_doses==3 & V_t > 6 ~ VE_s4,
                         V_doses==0 ~ 0
        ),
        VE_p = case_when(V_doses==2 & V_t <=6 ~ VE_p1, # Progression to symptoms
                         V_doses==2 & V_t > 6 ~ VE_p2,
                         V_doses==3 & V_t <=6 ~ VE_p3,
                         V_doses==3 & V_t > 6 ~ VE_p4,
                         V_doses==0 ~ 0
        )
      ) -> A.hcwC
    
    total_hcw <- total_hcw + A.hcwC.init
  } else{
    A.hcwC=as.data.frame(cbind("ID"=as.numeric(),"VL_rise" = as.numeric(),
                               "removal.pd"=as.numeric(),"VL"=as.numeric(), "VL_waning"=as.numeric(),
                               "V_doses" = as.numeric(), #"V_refused" = as.numeric(),  
                               "V_t"=as.numeric(), #"V2_t"=as.numeric(),                        
                               "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric(), 
                               "Inc.pd"=as.numeric(),
                               "Inf.pd"=as.numeric(),"Inf.days"=as.numeric(),
                               "ID.days"=as.numeric(),
                               "Infectiousness"=as.numeric())) 
  }
  
  if (I.hcwC.init > 0){
    I.hcwC=as.data.frame(cbind("ID"=(total_hcw+1):(total_hcw+ I.hcwC.init), "VL_rise" = round(runif(I.hcwC.init, 1, runif(1,4.7,6.2))),
                               "removal.pd"=rep(NA, I.hcwC.init),#"VL"=abs(rnorm(I.hcwC.init, 8, 1)*0.5), "VL_waning"=rep(8/runif(1, 15, 21),I.hcwC.init),
                               #"VL"=runif(I.hcwC.init,1, 8), "VL_waning"=rep(0.35,I.hcwC.init),
                               "VL"=abs(rnorm(I.hcwC.init, 8, 1)*0.5), "VL_waning"=rep(0.35,I.hcwC.init),
                               "Inc.pd"=rep(0,I.hcwC.init),
                               "Inf.pd"=round(runif(I.hcwC.init, parms[["gamma1"]], parms[["gamma2"]])),
                               "Inf.days"=rep(0,I.hcwC.init), # can change if we want to assume they are partway into inf period
                               "ID.days"=rep(0,I.hcwC.init),
                               "Infectiousness"=rep(1,I.hcwC.init)))
    I.hcwC %>%
      mutate(
        V_doses = sample(c(3,2,0),
                         size=I.hcwC.init,
                         prob=c(parms[["boost_vax_coverage_r"]],parms[["full_vax_coverage_r"]],parms[["refusal_r"]]),
                         replace=TRUE),
        V_t = case_when(V_doses==2 ~ round(runif(1,1,15)),
                        V_doses==3 ~ round(runif(1,1,9)),
                        V_doses==0 ~ 1000,
                        TRUE ~ V_doses),
        VE_i = case_when(V_doses==2 & V_t <=6 ~ VE_i1, 
                         V_doses==2 & V_t > 6 ~ VE_i2,
                         V_doses==3 & V_t <=6 ~ VE_i3,
                         V_doses==3 & V_t > 6 ~ VE_i4,
                         V_doses==0 ~ 0
        ),
        VE_s = case_when(V_doses==2 & V_t <=6 ~ VE_s1, # susceptibility (i.e. becoming infected)
                         V_doses==2 & V_t > 6 ~ VE_s2,
                         V_doses==3 & V_t <=6 ~ VE_s3,
                         V_doses==3 & V_t > 6 ~ VE_s4,
                         V_doses==0 ~ 0
        ),
        VE_p = case_when(V_doses==2 & V_t <=6 ~ VE_p1, # Progression to symptoms
                         V_doses==2 & V_t > 6 ~ VE_p2,
                         V_doses==3 & V_t <=6 ~ VE_p3,
                         V_doses==3 & V_t > 6 ~ VE_p4,
                         V_doses==0 ~ 0
        )
      ) -> I.hcwC
    
    total_hcw <- total_hcw + I.hcwC.init
  } else{
    I.hcwC=as.data.frame(cbind("ID"=as.numeric(),"VL_rise" = as.numeric(),"removal.pd"=as.numeric(),
                               "VL"=as.numeric(), "VL_waning"=as.numeric(),
                               "V_doses" = as.numeric(), #"V_refused" = as.numeric(),  
                               "V_t"=as.numeric(),# "V2_t"=as.numeric(),                        
                               "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric(), 
                               "Inc.pd"=as.numeric(),
                               "Inf.pd"=as.numeric(),"Inf.days"=as.numeric(),
                               "ID.days"=as.numeric(),
                               "Infectiousness"=as.numeric())) 
  }
  
  if (R.hcwC.init > 0){
    R.hcwC=as.data.frame(cbind("ID"=(total_hcw+1):(total_hcw + R.hcwC.init),"VL"=rep(0,R.hcwC.init),"VL_waning"=rep(0,R.hcwC.init),
                               "Inc.pd"=rep(0,R.hcwC.init),"Inf.pd"=rep(0,R.hcwC.init),
                               "Rec.days"=rep(0,R.hcwC.init)))
    
    R.hcwC %>%
      mutate(
        V_doses = sample(c(3,2,0),
                         size=R.hcwC.init,
                         prob=c(parms[["boost_vax_coverage_r"]],parms[["full_vax_coverage_r"]],parms[["refusal_r"]]),
                         replace=TRUE),
        V_t = case_when(V_doses==2 ~ round(runif(1,1,15)),
                        V_doses==3 ~ round(runif(1,1,9)),
                        V_doses==0 ~ 1000,
                        TRUE ~ V_doses),
        VE_i = case_when(V_doses==2 & V_t <=6 ~ VE_i1, 
                         V_doses==2 & V_t > 6 ~ VE_i2,
                         V_doses==3 & V_t <=6 ~ VE_i3,
                         V_doses==3 & V_t > 6 ~ VE_i4,
                         V_doses==0 ~ 0
        ),
        VE_s = case_when(V_doses==2 & V_t <=6 ~ VE_s1, # susceptibility (i.e. becoming infected)
                         V_doses==2 & V_t > 6 ~ VE_s2,
                         V_doses==3 & V_t <=6 ~ VE_s3,
                         V_doses==3 & V_t > 6 ~ VE_s4,
                         V_doses==0 ~ 0
        ),
        VE_p = case_when(V_doses==2 & V_t <=6 ~ VE_p1, # Progression to symptoms
                         V_doses==2 & V_t > 6 ~ VE_p2,
                         V_doses==3 & V_t <=6 ~ VE_p3,
                         V_doses==3 & V_t > 6 ~ VE_p4,
                         V_doses==0 ~ 0
        )
      ) -> R.hcwC
    
    total_hcw <- total_hcw + R.hcwC.init
  } else{
    R.hcwC=as.data.frame(cbind("ID"=as.numeric(),"VL"=as.numeric(),"VL_waning"=as.numeric(), 
                               "V_doses" = as.numeric(), #"V_refused" = as.numeric(),  
                               "V_t"=as.numeric(), #"V2_t"=as.numeric(),                        
                               "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric(), 
                               "Inc.pd"=as.numeric(),"Inf.pd"=as.numeric(),
                               "Rec.days"=as.numeric()))
  }
  
  if (I.hcwH.init > 0){
    I.hcwH=as.data.frame(cbind("ID"=(total_hcw+1):(total_hcw+ I.hcwH.init), "VL_rise" = round(runif(I.hcwH.init, 1, runif(1,4.7,6.2))),
                               "VL"=abs(rnorm(I.hcwH.init, 8, 1)*0.5), "VL_waning"=rep(0.35,I.hcwH.init),
                               "Inc.pd"=rep(0,I.hcwC.init),
                               "Rec.days"=rep(NA,I.hcwC.init),
                               "Inf.pd"=round(runif(I.hcwH.init, parms[["gamma1"]], parms[["gamma2"]])),
                               "Inf.days"=rep(0,I.hcwH.init),
                               "Home.pd"=round(runif(I.hcwH.init, parms[["gamma1"]], parms[["gamma2"]])),
                               "Home.days"=rep(0,I.hcwH.init), # can change if we want to assume they are partway into inf period
                               "Infectiousness"=rep(1,I.hcwC.init)))
    
    I.hcwH %>%
      mutate(
        V_doses = sample(c(3,2,0),
                         size=I.hcwH.init,
                         prob=c(parms[["boost_vax_coverage_r"]],parms[["full_vax_coverage_r"]],parms[["refusal_r"]]),
                         replace=TRUE),
        V_t = case_when(V_doses==2 ~ round(runif(1,1,15)),
                        V_doses==3 ~ round(runif(1,1,9)),
                        V_doses==0 ~ 1000,
                        TRUE ~ V_doses),
        VE_i = case_when(V_doses==2 & V_t <=6 ~ VE_i1, 
                         V_doses==2 & V_t > 6 ~ VE_i2,
                         V_doses==3 & V_t <=6 ~ VE_i3,
                         V_doses==3 & V_t > 6 ~ VE_i4,
                         V_doses==0 ~ 0
        ),
        VE_s = case_when(V_doses==2 & V_t <=6 ~ VE_s1, # susceptibility (i.e. becoming infected)
                         V_doses==2 & V_t > 6 ~ VE_s2,
                         V_doses==3 & V_t <=6 ~ VE_s3,
                         V_doses==3 & V_t > 6 ~ VE_s4,
                         V_doses==0 ~ 0
        ),
        VE_p = case_when(V_doses==2 & V_t <=6 ~ VE_p1, # Progression to symptoms
                         V_doses==2 & V_t > 6 ~ VE_p2,
                         V_doses==3 & V_t <=6 ~ VE_p3,
                         V_doses==3 & V_t > 6 ~ VE_p4,
                         V_doses==0 ~ 0
        )
      ) -> I.hcwH
    
    total_hcw <- total_hcw + I.hcwC.init
  } else{
    I.hcwH=as.data.frame(cbind("ID"=as.numeric(),"VL_rise" = as.numeric(),"VL"=as.numeric(), "VL_waning"=as.numeric(),
                               "V_doses" = as.numeric(), #"V_refused" = as.numeric(),  
                               "V_t"=as.numeric(), #"V2_t"=as.numeric(),                        
                               "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric(), 
                               "Inc.pd"=as.numeric(),
                               "Rec.days"=as.numeric(),
                               "Inf.pd"=as.numeric(),"Inf.days"=as.numeric(),
                               "Home.pd"=as.numeric(),"Home.days"=as.numeric(),
                               "Infectiousness"=as.numeric())) 
  }
  
  
  S.rNC3 = as.data.frame(cbind("ID"=as.numeric(),"VL"=as.numeric(), 
                               "V_doses" = as.numeric(),"V_t"=as.numeric(),                        
                               "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric(),
                               "Inc.pd"=as.numeric(), "Days"=as.numeric())) 
  S.hcwNC3 = as.data.frame(cbind("ID"=as.numeric(),"VL"=as.numeric(), 
                               "V_doses" = as.numeric(),"V_t"=as.numeric(),                        
                               "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric(),
                               "Inc.pd"=as.numeric(), "Days"=as.numeric())) 
  S.hcwC3 = as.data.frame(cbind("ID"=as.numeric(),"VL"=as.numeric(), 
                               "V_doses" = as.numeric(),"V_t"=as.numeric(),                        
                               "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric(),
                               "Inc.pd"=as.numeric(), "Days"=as.numeric())) 
  
  
  death.r=as.data.frame(cbind("ID"=as.numeric(),"VL"=as.numeric(), 
                              "V_doses" = as.numeric(),  
                              "V_t"=as.numeric(),                       
                              "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric())) 
  
  death.rC=as.data.frame(cbind("ID"=as.numeric(),"VL"=as.numeric(), 
                              "V_doses" = as.numeric(), 
                              "V_t"=as.numeric(),                    
                              "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric())) 
  
  
  rNC=as.data.frame(rNC)

  family1people=sample(rNC[,1],group_1p,replace=FALSE)
  family1Room <- seq(10000,10000+length(family1people)-1,1)
  family1=as.data.frame(cbind("ResID"=family1people, "Room" = family1Room))
  family = family1
  
  rNC %>%
    subset(!(rNC %in% family$ResID)) -> rNC 
  family2people1=sample(rNC[,1],group_2p,replace=FALSE)
  family2people2=sample(setdiff(rNC[,1],family2people1),group_2p,replace=FALSE)
  family2people=c(family2people1,family2people2)
  
  family2Room1 <- seq(20000,20000+length(family2people1)-1,1) 
  family2Room=c(family2Room1,family2Room1) #
  
  family2=as.data.frame(cbind("ResID"=family2people, "Room" = family2Room))

  family %>%
    bind_rows(family2) -> family 
  
  rNC %>%
    subset(!(rNC %in% family$ResID)) -> rNC 
  family3people1=sample(rNC[,1],group_3p,replace=FALSE)
  family3people=family3people1
  
  rNCleft=setdiff(rNC[,1],family3people) 
  family3people2=sample(rNCleft,group_3p,replace=FALSE)
  family3people=c(family3people1,family3people2)
  
  rNCleft=setdiff(rNC[,1],family3people) 
  family3people3=sample(rNCleft,group_3p,replace=FALSE)
  family3people=c(family3people1,family3people2,family3people3)
  
  family3Room1 <- seq(30000,30000+length(family3people1)-1,1) 
  family3Room=c(family3Room1,family3Room1,family3Room1) 
  
  family3=as.data.frame(cbind("ResID"=family3people, "Room" = family3Room))
  
  family %>%
    bind_rows(family3) -> family 
  
  rNC %>%
    subset(!(rNC %in% family$ResID)) -> rNC 庭
  family4people1=sample(rNC[,1],group_4p,replace=FALSE)
  family4people=family4people1
  
  rNCleft=setdiff(rNC[,1],family4people) 
  family4people2=sample(rNCleft,group_4p,replace=FALSE)
  family4people=c(family4people1,family4people2)
  
  rNCleft=setdiff(rNC[,1],family3people) 
  family4people3=sample(rNCleft,group_4p,replace=FALSE)
  family4people=c(family4people1,family4people2,family4people3)
  
  rNCleft=setdiff(rNC[,1],family3people) 
  family4people4=sample(rNCleft,group_4p,replace=FALSE)
  family4people=c(family4people1,family4people2,family4people3,family4people4)
  
  family4Room1 <- seq(40000,40000+length(family4people1)-1,1) 
  family4Room=c(family4Room1,family4Room1,family4Room1,family4Room1) 
  
  family4=as.data.frame(cbind("ResID"=family4people, "Room" = family4Room))
  
  family %>%
    bind_rows(family4) -> family 
  
  rNC %>%
    subset(!(rNC %in% family$ResID)) -> rNC 
  family5people1=sample(rNC[,1],group_5p,replace=FALSE)
  family5people=family5people1
  
  rNCleft=setdiff(rNC[,1],family5people) 
  family5people2=sample(rNCleft,group_5p,replace=FALSE)
  family5people=c(family5people1,family5people2)
  
  
  rNCleft=setdiff(rNC[,1],family5people) 
  family5people3=sample(rNCleft,group_5p,replace=FALSE)
  family5people=c(family5people1,family5people2,family5people3)
  
  
  rNCleft=setdiff(rNC[,1],family5people)
  family5people4=sample(rNCleft,group_5p,replace=FALSE)
  family5people=c(family5people1,family5people2,family5people3,family5people4)
  
  
  rNCleft=setdiff(rNC[,1],family5people) 
  family5people5=sample(rNCleft,group_5p,replace=FALSE)
  family5people=c(family5people1,family5people2,family5people3,family5people4,family5people5)
  
  family5Room1 <- seq(50000,50000+length(family5people1)-1,1) 
  family5Room=c(family5Room1,family5Room1,family5Room1,family5Room1,family5Room1) 
  
  family5=as.data.frame(cbind("ResID"=family5people, "Room" = family5Room))
  
  family %>%
    bind_rows(family5) -> family 
 
    ResidentFlag=rep(1,nrow(family)) 
    
    family %>%
    mutate(ResidentFlag=1) -> family  
    
      vax_count1_r <- nrow(S.rNC %>% subset(V_doses==2))+nrow(E.rNC %>% subset(V_doses==2))+nrow(A.rNC %>% subset(V_doses==2))+nrow(I.rNC %>% subset(V_doses==2))+nrow(R.rNC %>% subset(V_doses==2))
      +nrow(I.rC %>% subset(V_doses==2))+nrow(R.rC %>% subset(V_doses==2)) 
      vax_count2_r <- nrow(S.rNC %>% subset(V_doses==3))+nrow(E.rNC %>% subset(V_doses==3))+nrow(A.rNC %>% subset(V_doses==3))+nrow(I.rNC %>% subset(V_doses==3))+nrow(R.rNC %>% subset(V_doses==3))
      +nrow(I.rC %>% subset(V_doses==3))+nrow(R.rC %>% subset(V_doses==3))
      
      vax_count1_hcw <- nrow(S.hcwNC %>% subset(V_doses==2))+nrow(E.hcwNC %>% subset(V_doses==2))+nrow(A.hcwNC %>% subset(V_doses==2))+nrow(I.hcwNC %>% subset(V_doses==2))+nrow(R.hcwNC %>% subset(V_doses==2))
      +nrow(S.hcwC %>% subset(V_doses==2))+nrow(E.hcwC %>% subset(V_doses==2))+nrow(A.hcwC %>% subset(V_doses==2))+nrow(I.hcwC %>% subset(V_doses==2))+nrow(R.hcwC %>% subset(V_doses==2))
      +nrow(I.hcwH %>% subset(V_doses==2))
      vax_count2_hcw <- nrow(S.hcwNC %>% subset(V_doses==3))+nrow(E.hcwNC %>% subset(V_doses==3))+nrow(A.hcwNC %>% subset(V_doses==3))+nrow(I.hcwNC %>% subset(V_doses==3))+nrow(R.hcwNC %>% subset(V_doses==3))
      +nrow(S.hcwC %>% subset(V_doses==3))+nrow(E.hcwC %>% subset(V_doses==3))+nrow(A.hcwC %>% subset(V_doses==3))+nrow(I.hcwC %>% subset(V_doses==3))+nrow(R.hcwC %>% subset(V_doses==3))
      +nrow(I.hcwH %>% subset(V_doses==3))
    
    total_hcw <- total_hcw-5999
    
  Ns <- list("S.rNC"=S.rNC,
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
             "inc_r"=0,
             "inc_hcw"=0,
             "cum_inc_r"=0,
             "cum_inc_hcw"=0,
             "cum_inc_vax1_r"=0,
             "cum_inc_vax2_r"=0,
             "cum_inc_vax1_hcw"=0,
             "cum_inc_vax2_hcw"=0,
             "cum_inc_community"=0,
             "mortality"=0,
             "total"=total,
             "total_hcw"=total_hcw,
             "family"=family, 
             "vax_count1_r"=vax_count1_r,
             "vax_count2_r"=vax_count2_r,
             "vax_count1_hcw"=vax_count1_hcw,
             "vax_count2_hcw"=vax_count2_hcw,
             "asymptomatic_inc_r"=0,
             "symptomatic_inc_r"=0,
             "asymptomatic_inc_hcw"=0,
             "symptomatic_inc_hcw"=0)
  
  return(Ns)
}
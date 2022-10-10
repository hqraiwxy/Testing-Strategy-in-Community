recover_or_test_r <- function(df,parms,symptoms,t,TimeAndStrategy){
  
  df %>%
    subset(Inf.days==Inf.pd) %>% 
    mutate(Rec.days = Inf.pd) %>% 
    dplyr::select(ID,VL,VL_waning,Rec.days, 
                  V_doses,V_t,VE_i,VE_s,VE_p) -> recovered 
  df %>%
    subset(!(ID %in% recovered$ID)) -> df 
  
  df %>% 
    subset(removal.pd == ID.days) -> tested1 
  
  df %>% 
    subset(!(ID %in% tested1$ID) & !(ID %in% recovered$ID)) -> df 
  
  if(symptoms=="A"){    
    df %>%
      subset(((TimeAndStrategy[ID.days+1,1]==1 & VL>=parms[["VL.PCR_threshold"]]) | (TimeAndStrategy[ID.days+1,2]==1 & VL>=parms[["VL.Antigen_threshold"]] )) & is.na(removal.pd)) %>% #建立子集，检测当天病毒载量高于阈值并且没有被转移走？
      mutate(removal.pd = case_when(TimeAndStrategy[ID.days+1,1]==1 & VL>=parms[["VL.PCR_threshold"]] ~ ID.days + parms[["test_PCR_delay"]],
                                    TimeAndStrategy[ID.days+1,2]==1 & VL>=parms[["VL.Antigen_threshold"]]~ ID.days) 
             ) -> df1 #asymptomatic who are identified through testing
    
  }else{
    
    df %>%
      subset(ID.days == Inc.pd+parms[["id.I"]] | ((TimeAndStrategy[ID.days+1,1]==1 & VL>=parms[["VL.PCR_threshold"]]) | (TimeAndStrategy[ID.days+1,2]==1 & VL>=parms[["VL.Antigen_threshold"]] ))) %>% ## id.I=duration of pre-symptomatic transmission
      mutate(removal.pd = case_when(ID.days == Inc.pd+parms[["id.I"]] ~ ID.days, 
                                    TimeAndStrategy[ID.days+1,1]==1 & VL>=parms[["VL.PCR_threshold"]]~ ID.days + parms[["test_PCR_delay"]],
                                    TimeAndStrategy[ID.days+1,2]==1 & VL>=parms[["VL.Antigen_threshold"]]~ ID.days) 
             ) -> df1   
  }
  

  
  df %>%
    subset(((TimeAndStrategy[ID.days+1,1]==1 & VL<parms[["VL.PCR_threshold"]]) | (TimeAndStrategy[ID.days+1,2]==1 & VL<parms[["VL.Antigen_threshold"]] )) & is.na(removal.pd)  & !(ID %in% df1$ID)) %>%  # !(ID %in% df1$ID) excludes people identified through symptoms
    mutate(Inf.days=Inf.days + 1, 
           ID.days = ID.days + 1, 
           VL = case_when(Inf.days-Inc.pd <=VL_rise ~ VL * (Inf.days-Inc.pd+1)/(Inf.days-Inc.pd),
                          Inf.days-Inc.pd > VL_rise  ~ VL - VL_waning),
           VL = case_when(VL < 0 ~ 0, TRUE ~ VL),
           Infectiousness = case_when(VL<3 ~ 0, 
                                      VL>=3 & VL<7 ~ 0.5, 
                                      VL>=7 ~ 1)) -> df2   
  df %>%
    subset(((TimeAndStrategy[ID.days+1,1]==0 & TimeAndStrategy[ID.days+1,2]==0) | removal.pd > ID.days | is.na(ID.days)) & !(ID %in% df1$ID)) %>%
    mutate(Inf.days = Inf.days + 1,
           ID.days = ID.days + 1,
           VL = case_when(Inf.days-Inc.pd <=VL_rise ~ VL * (Inf.days-Inc.pd+1)/(Inf.days-Inc.pd),
                          Inf.days-Inc.pd > VL_rise ~ VL - VL_waning),
           VL = case_when(VL < 0 ~ 0, TRUE ~ VL),
           Infectiousness = case_when(VL<3 ~ 0, 
                                      VL>=3 & VL<7 ~ 0.5, 
                                      VL>=7 ~ 1)) -> df3   

  df1 %>% 
    subset(ID.days==removal.pd) %>%  #tested2
    mutate(ID.days = NA) -> tested2
  
  tested1 %>%
    bind_rows(tested2) -> tested
  
  df1 %>%
    subset(ID.days != removal.pd) %>%
    mutate(ID.days = ID.days + 1,
           Inf.days = Inf.days + 1,
           VL = case_when(Inf.days-Inc.pd <=VL_rise ~ VL * (Inf.days-Inc.pd+1)/(Inf.days-Inc.pd),
                          Inf.days-Inc.pd > VL_rise  ~ VL - VL_waning),
           VL = case_when(VL < 0 ~ 0, TRUE ~ VL),
           Infectiousness = case_when(VL<3 ~ 0, 
                                      VL>=3 & VL<7 ~ 0.5, 
                                      VL>=7 ~ 1)) %>%
    bind_rows(df2) %>%
    bind_rows(df3) -> df

  list("recovered"=recovered,  
       "tested"=tested, 
       "df"=df)  
}

recover_I.rC <- function(df,parms){
  
  df %>%
    subset(Inf.days==Inf.pd) %>%
    mutate(Rec.days = Inf.pd) %>%
    dplyr::select(ID,VL,VL_waning,Rec.days, 
                  V_doses,V_t,VE_i,VE_s,VE_p) -> recovered
  
  df %>%
    subset(!(ID %in% recovered$ID)) %>%
    mutate(Inf.days = Inf.days + 1,
           VL = case_when(Inf.days-Inc.pd <=VL_rise ~ VL * (Inf.days-Inc.pd+1)/(Inf.days-Inc.pd),
                          Inf.days-Inc.pd > VL_rise ~ VL - VL_waning),
           VL = case_when(VL < 0 ~ 0, TRUE ~ VL)) -> df  
  
  list("recovered"=recovered,
       "df"=df)
}
  

recover_or_test_hcw <- function(df,parms,total,symptoms,ratio,t,TimeAndStrategy){  # ratio of staff to residents 【ratio <- (N.hcwNC + N.hcwC)/(N.rNC + N.rC)】
  df %>%
    subset(Inf.days==Inf.pd) %>% 
    mutate(Rec.days = Inf.pd) %>%
    dplyr::select(ID,VL,VL_waning,Rec.days,Inc.pd,Inf.pd,
                  V_doses,V_t,VE_i,VE_s,VE_p) -> recovered 
  df %>%
    subset(!(ID %in% recovered$ID)) -> df 
  
  df %>%
    subset(removal.pd == ID.days) -> tested1 
  
  df %>% 
    subset(!(ID %in% tested1$ID) & !(ID %in% recovered$ID)) -> df 
  

  if(symptoms=="A"){    ## asymptomatics who are identified through testing
    df %>%
      subset(((TimeAndStrategy[ID.days+1,1]==1 & VL>=parms[["VL.PCR_threshold"]]) | (TimeAndStrategy[ID.days+1,2]==1 & VL>=parms[["VL.Antigen_threshold"]] )) & is.na(removal.pd)) %>%
      mutate(removal.pd = case_when(TimeAndStrategy[ID.days+1,1]==1 & VL>=parms[["VL.PCR_threshold"]] ~ ID.days + parms[["test_PCR_delay"]], #parms[["test_PCR_delay"]]
                                    TimeAndStrategy[ID.days+1,2]==1 & VL>=parms[["VL.Antigen_threshold"]]~ ID.days)) -> df1
  }else{
    
    df %>%
      subset(ID.days == Inc.pd+parms[["id.I"]] | ((TimeAndStrategy[ID.days+1,1]==1 & VL>=parms[["VL.PCR_threshold"]]) | (TimeAndStrategy[ID.days+1,2]==1 & VL>=parms[["VL.Antigen_threshold"]] )) & is.na(removal.pd)) %>%
      mutate(removal.pd = case_when(ID.days == Inc.pd+parms[["id.I"]] ~ ID.days, 
                                    TimeAndStrategy[ID.days,1]==1 & VL>=parms[["VL.PCR_threshold"]]~ ID.days + parms[["test_PCR_delay"]],#parms[["test_PCR_delay"]]
                                    TimeAndStrategy[ID.days,2]==1 & VL>=parms[["VL.Antigen_threshold"]]~ ID.days)) -> df1##  ## symptomatics who are identified through testing OR symptoms
  }


  df %>%
    subset(((TimeAndStrategy[ID.days+1,1]==1 & VL<parms[["VL.PCR_threshold"]]) | (TimeAndStrategy[ID.days+1,2]==1 & VL<parms[["VL.Antigen_threshold"]] )) & is.na(removal.pd) & !(ID %in% df1$ID)) %>%
    mutate(Inf.days=Inf.days + 1,
           ID.days=ID.days + 1,
           VL = case_when(Inf.days-Inc.pd <=VL_rise ~ VL * (Inf.days-Inc.pd+1)/(Inf.days-Inc.pd),
                          Inf.days-Inc.pd > VL_rise  ~ VL - VL_waning),
           VL = case_when(VL < 0 ~ 0, TRUE ~ VL),
           Infectiousness = case_when(VL<3 ~ 0, 
                                      VL>=3 & VL<7 ~ 0.5, 
                                      VL>=7 ~ 1)) -> df2   ## anyone who is tested but NOT identified 
  
  df %>%
    subset(((TimeAndStrategy[ID.days+1,1]==0 & TimeAndStrategy[ID.days+1,2]==0) | removal.pd > ID.days) & !(ID %in% df1$ID)) %>% ##没到检测周期定的时间or没到回复结果并被移除的时间【区别于res，工作人员的I.hcwC不需要VL】
    mutate(Inf.days = Inf.days + 1,
           ID.days = ID.days + 1,
           VL = case_when(Inf.days-Inc.pd <=VL_rise ~ VL * (Inf.days-Inc.pd+1)/(Inf.days-Inc.pd),
                          Inf.days-Inc.pd > VL_rise  ~ VL - VL_waning),
           VL = case_when(VL < 0 ~ 0, TRUE ~ VL),
           Infectiousness = case_when(VL<3 ~ 0, 
                                      VL>=3 & VL<7 ~ 0.5, 
                                      VL>=7 ~ 1)) -> df3  ## anyone who is not tested, does not have symptoms

    df1 %>% 
      subset(ID.days == removal.pd) %>%  
      mutate(ID.days = NA) -> tested2
    
    tested1 %>%
      bind_rows(tested2) %>%
      mutate(Home.pd=removal.pd+round(runif(1, parms[["gamma1"]], parms[["gamma2"]])),
             Home.days=removal.pd) %>%
      dplyr::select(-ID.days,-removal.pd) -> tested 
    
    df1 %>%
      subset(ID.days != removal.pd) %>%
      mutate(ID.days = ID.days + 1,
             Inf.days = Inf.days + 1,
             VL = case_when(Inf.days-Inc.pd <=VL_rise ~ VL * (Inf.days-Inc.pd+1)/(Inf.days-Inc.pd),
                            Inf.days-Inc.pd > VL_rise  ~ VL - VL_waning),
             VL = case_when(VL < 0 ~ 0, TRUE ~ VL),
             Infectiousness = case_when(VL<3 ~ 0, 
                                        VL>=3 & VL<7 ~ 0.5, 
                                        VL>=7 ~ 1))%>% 
      bind_rows(df2) %>%
      bind_rows(df3) -> df
    
    
  
  if (nrow(tested)>0){  
    
    if (isFALSE(parms[["shortage"]]) | ratio < parms[["shortage_threshold"]]){  
      new_E <- rbinom(1,nrow(tested),parms[["I.C"]]) 
      new_R <- rbinom(1,(nrow(tested)-new_E),parms[["prop_rhcwR"]])
      new_S <- nrow(tested) - new_E - new_R
      if (new_S >0){    
        new.hcwS <- as.data.frame(cbind("ID"=100000+c((total + 1):(total + new_S)),"VL"=NA))
        
        new.hcwS %>%
          mutate(
            V_doses = sample(c(3,2,0),
                             size=new_S,
                             prob=c(parms[["boost_vax_coverage_hcw"]],parms[["full_vax_coverage_hcw"]],parms[["refusal_hcw"]]),
                             replace=TRUE),
            V_t = case_when(V_doses==2 ~ round(runif(1,1,15)),
                            V_doses==3 ~ round(runif(1,1,9)),
                            V_doses==0 ~ 1000,
                            TRUE ~ V_doses),
            VE_i = case_when(V_doses==2 & V_t <=6 ~ parms[["VE_i1"]], 
                             V_doses==2 & V_t > 6 ~ parms[["VE_i2"]],
                             V_doses==3 & V_t <=6 ~ parms[["VE_i3"]],
                             V_doses==3 & V_t > 6 ~ parms[["VE_i4"]],
                             V_doses==0 ~ 0
                             ),
            VE_s = case_when(V_doses==2 & V_t <=6 ~ parms[["VE_s1"]], 
                             V_doses==2 & V_t > 6 ~ parms[["VE_s2"]],
                             V_doses==3 & V_t <=6 ~ parms[["VE_s3"]],
                             V_doses==3 & V_t > 6 ~ parms[["VE_s4"]],
                             V_doses==0 ~ 0
                             ),
            VE_p = case_when(V_doses==2 & V_t <=6 ~ parms[["VE_p1"]], 
                             V_doses==2 & V_t > 6 ~ parms[["VE_p2"]],
                             V_doses==3 & V_t <=6 ~ parms[["VE_p3"]],
                             V_doses==3 & V_t > 6 ~ parms[["VE_p4"]],
                             V_doses==0 ~ 0
                             )
          ) -> new.hcwS
      } else{
        new.hcwS <-  as.data.frame(cbind("ID"=as.character(),"VL"=as.numeric(),
                                         "V_doses" = as.numeric(), #"V_refused" = as.numeric(),  
                                         "V_t"=as.numeric(), #"V2_t"=as.numeric()))                        
                                         "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric()))
      }
      
      
     
      if (new_E >0){ 
        new.hcwE <- as.data.frame(cbind("ID"=100000+c((total + new_S + 1):(total+new_S + new_E)),"VL"=0,
                                   "Inc.pd"= t+round(runif(length(new_E), parms[["sigma1"]], parms[["sigma2"]])), "Days"=t))
         new.hcwE %>%
          mutate(
            V_doses = sample(c(3,2,0),
                             size=new_E,
                             prob=c(parms[["boost_vax_coverage_hcw"]],parms[["full_vax_coverage_hcw"]],parms[["refusal_hcw"]]),
                             replace=TRUE),
            V_t = case_when(V_doses==2 ~ round(runif(1,1,15)),
                            V_doses==3 ~ round(runif(1,1,9)),
                            V_doses==0 ~ 1000,
                            TRUE ~ V_doses
                            ),
            VE_i = case_when(V_doses==2 & V_t <=6 ~ parms[["VE_i1"]], 
                             V_doses==2 & V_t > 6 ~ parms[["VE_i2"]],
                             V_doses==3 & V_t <=6 ~ parms[["VE_i3"]],
                             V_doses==3 & V_t > 6 ~ parms[["VE_i4"]],
                             V_doses==0 ~ 0
                             ),
            VE_s = case_when(V_doses==2 & V_t <=6 ~ parms[["VE_s1"]],
                             V_doses==2 & V_t > 6 ~ parms[["VE_s2"]],
                             V_doses==3 & V_t <=6 ~ parms[["VE_s3"]],
                             V_doses==3 & V_t > 6 ~ parms[["VE_s4"]],
                             V_doses==0 ~ 0
                             ),
            VE_p = case_when(V_doses==2 & V_t <=6 ~ parms[["VE_p1"]], 
                             V_doses==2 & V_t > 6 ~ parms[["VE_p2"]],
                             V_doses==3 & V_t <=6 ~ parms[["VE_p3"]],
                             V_doses==3 & V_t > 6 ~ parms[["VE_p4"]],
                             V_doses==0 ~ 0
                             )
          ) -> new.hcwE
        
      } else{
        new.hcwE <- as.data.frame(cbind("ID"=as.character(),"VL"=as.numeric()),
                                  "V_doses" = as.numeric(), 
                                  "V_t"=as.numeric(),                     
                                  "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric(),
                                  "Inc.pd"=as.numeric(),"Days"=as.numeric())
      }
      
      
      if (new_R >0){   
        new.hcwR <- as.data.frame(cbind("ID"=100000+c((total + new_S + new_E + 1):(total+new_S + new_E + new_R)),"VL"=0,
                                        "VL_waning"=0, 
                                        "Inc.pd"=0, 
                                        "Inf.pd"=t, 
                                        "Rec.days"=t)) 
        new.hcwR %>%
          mutate(
            V_doses = sample(c(3,2,0),
                             size=new_R,
                             prob=c(parms[["boost_vax_coverage_hcw"]],parms[["full_vax_coverage_hcw"]],parms[["refusal_hcw"]]),
                             replace=TRUE),
            V_t = case_when(V_doses==2 ~ round(runif(1,1,15)),
                            V_doses==3 ~ round(runif(1,1,9)),
                            V_doses==0 ~ 1000,
                            TRUE ~ V_doses),
            VE_i = case_when(V_doses==2 & V_t <=6 ~ parms[["VE_i1"]], 
                             V_doses==2 & V_t > 6 ~ parms[["VE_i2"]],
                             V_doses==3 & V_t <=6 ~ parms[["VE_i3"]],
                             V_doses==3 & V_t > 6 ~ parms[["VE_i4"]],
                             V_doses==0 ~ 0
                             ),
            VE_s = case_when(V_doses==2 & V_t <=6 ~ parms[["VE_s1"]], 
                             V_doses==2 & V_t > 6 ~ parms[["VE_s2"]],
                             V_doses==3 & V_t <=6 ~ parms[["VE_s3"]],
                             V_doses==3 & V_t > 6 ~ parms[["VE_s4"]],
                             V_doses==0 ~ 0
                             ),
            VE_p = case_when(V_doses==2 & V_t <=6 ~ parms[["VE_p1"]], 
                             V_doses==2 & V_t > 6 ~ parms[["VE_p2"]],
                             V_doses==3 & V_t <=6 ~ parms[["VE_p3"]],
                             V_doses==3 & V_t > 6 ~ parms[["VE_p4"]],
                             V_doses==0 ~ 0
                             )
            
          ) -> new.hcwR
          
      } else{
        new.hcwR <-  as.data.frame(cbind("ID"=as.character(),"VL"=as.numeric(),
                                         "V_doses" = as.numeric(),  
                                         "V_t"=as.numeric(),                       
                                         "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric(),
                                         "VL_waning"=as.numeric(),
                                         "Inc.pd"=as.numeric(), "Inf.pd"=as.numeric(), 
                                         "Rec.days"=as.numeric()))
      }
    } else{   
      new.hcwS <-  as.data.frame(cbind("ID"=as.character(),"VL"=as.numeric(),
                                       "V_doses" = as.numeric(),  
                                       "V_t"=as.numeric(),                      
                                       "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric()))
      new.hcwE <- as.data.frame(cbind("ID"=as.character(),"VL"=as.numeric(),
                                      "V_doses" = as.numeric(),  
                                      "V_t"=as.numeric(),                       
                                      "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric(),
                                      "Inc.pd"=as.numeric(),"Days"=as.numeric()))
      new.hcwR <-  as.data.frame(cbind("ID"=as.character(),"VL"=as.numeric(),
                                       "V_doses" = as.numeric(), 
                                       "V_t"=as.numeric(),                        
                                       "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric(),
                                       "VL_waning"=as.numeric(),
                                       "Inc.pd"=as.numeric(),"Inf.pd"=as.numeric(),
                                       "Rec.days"=as.numeric()))
    }
  } else{   
    new.hcwS <-  as.data.frame(cbind("ID"=as.character(),"VL"=as.numeric(),
                                     "V_doses" = as.numeric(),
                                     "V_t"=as.numeric(),                   
                                     "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric()))
    new.hcwE <- as.data.frame(cbind("ID"=as.character(),"VL"=as.numeric(),
                                    "V_doses" = as.numeric(), 
                                    "V_t"=as.numeric(),                        
                                    "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric(),
                                    "Inc.pd"=as.numeric(),"Days"=as.numeric()))
    new.hcwR <-  as.data.frame(cbind("ID"=as.character(),"VL"=as.numeric(),
                                     "V_doses" = as.numeric(), 
                                     "V_t"=as.numeric(),                     
                                     "VE_i"=as.numeric(), "VE_s"=as.numeric(), "VE_p"=as.numeric(),
                                     "VL_waning"=as.numeric(),
                                     "Inc.pd"=as.numeric(),"Inf.pd"=as.numeric(),
                                     "Rec.days"=as.numeric()))
  }
  
  tested %>%
   subset(!(ID>100000)) -> tested   
  
  list("recovered"=recovered,
       "tested"=tested,
       "df"=df,
       "new.hcwE"=new.hcwE,
       "new.hcwS"=new.hcwS,
       "new.hcwR"=new.hcwR)
  
}

recover_home_hcw <- function(I.hcwH,parms){      
  I.hcwH %>%
    subset(Home.pd==Home.days) %>%  
    mutate(Rec.days=Home.pd) %>%   
    dplyr::select(ID,VL,VL_waning,Rec.days,
                  V_doses,V_t,VE_i,VE_s,VE_p,Inc.pd,Inf.pd) -> recovered   
                                                             
  I.hcwH %>%
    subset(!(ID %in% recovered$ID)) %>%
    mutate(Inf.days = Inf.days + 1,
           Home.days = Home.days + 1,
           VL = case_when(Inf.days-Inc.pd <=VL_rise ~ VL * (Inf.days-Inc.pd+1)/(Inf.days-Inc.pd),
                          Inf.days-Inc.pd > VL_rise  ~ VL - VL_waning),
           VL = case_when(VL < 0 ~ 0, TRUE ~ VL),
           Infectiousness = case_when(VL<3 ~ 0, 
                                      VL>=3 & VL<7 ~ 0.5, 
                                      VL>=7 ~ 1)) -> I.hcwH      
  
  list("recovered"=recovered,
       "I.hcwH"=I.hcwH)  
  
}


E_to_I_r <- function(df,parms){
  
  df %>% 
    subset(Days==Inc.pd) %>%
    dplyr::select(ID,VL,V_doses,V_t,VE_i,VE_s,VE_p,Inc.pd,Days) %>% 
    mutate(asympt = rbinom(length(ID),1,(1-(1-parms[["alpha.r"]])*(1-VE_p))),  
           Inf.pd=Inc.pd+round(runif(1, parms[["gamma1"]], parms[["gamma2"]])),
           Inf.days=Inc.pd,
           ID.days=Inc.pd,
           removal.pd = NA) -> infected

  infected %>%
    subset(asympt==1) %>%
    mutate(
           VL_rise = round(runif(sum(asympt==1), 1, runif(1,4.7,6.2))),
           VL=abs(rnorm(sum(asympt==1), 8, 1)*0.5),
           VL_waning=rep(0.35,sum(asympt==1)),
           Infectiousness = case_when(VL<3 ~ 0, 
                                      VL>=3 & VL<7 ~ 0.5, 
                                      VL>=7 ~ 1))-> infected_asympt
  #保证输出的A结构不变
  infected_asympt %>%
    dplyr::select(ID,VL,V_doses,V_t,VE_i,VE_s,VE_p,Inc.pd,Inf.pd,Inf.days,ID.days,removal.pd,VL_rise,VL_waning,Infectiousness) -> infected_asympt
  
  infected %>%
    subset(asympt==0) %>%
    mutate(
           VL_rise = round(runif(sum(asympt==0), 1, runif(1,4.7,6.2))),
           VL=abs(rnorm(sum(asympt==0), 8, 1)*0.5),
           VL_waning=rep(0.35,sum(asympt==0)),
           Infectiousness = case_when(VL<3 ~ 0, 
                                      VL>=3 & VL<7 ~ 0.5, 
                                      VL>=7 ~ 1)) -> infected_sympt
  
 
  infected_sympt %>%
    dplyr::select(ID,VL,V_doses,V_t,VE_i,VE_s,VE_p,Inc.pd,Inf.pd,Inf.days,ID.days,removal.pd,VL_rise,VL_waning,Infectiousness) -> infected_sympt

  df %>%
    subset(Inc.pd != Days) %>%
    mutate(Days=Days + 1) -> df
  
  list("infected_asympt"=infected_asympt,
       "infected_sympt"=infected_sympt,
       "df"=df)
}

E_to_I_hcw <- function(df,parms){
  
  df %>%
    subset(Days==Inc.pd) %>%
    dplyr::select(ID,VL,V_doses,V_t,VE_i,VE_s,VE_p,Inc.pd,Days) %>%
    mutate(asympt = rbinom(length(ID),1,(1-(1-parms[["alpha.hcw"]])*(1-VE_p))),
           Inf.pd=Inc.pd+round(runif(1, parms[["gamma1"]], parms[["gamma2"]])),
           Inf.days=Inc.pd,
           ID.days=Inc.pd, 
           removal.pd = NA) -> infected

  infected %>%
    subset(asympt==1) %>%
    mutate(
           VL_rise = round(runif(sum(asympt==1), 1, runif(1,4.7,6.2))),
           VL=abs(rnorm(sum(asympt==1), 8, 1)*0.5),
           VL_waning=rep(0.35,sum(asympt==1)),
           Infectiousness = case_when(VL<3 ~ 0, 
                                      VL>=3 & VL<7 ~ 0.5, 
                                      VL>=7 ~ 1)) %>%
    dplyr::select(-asympt) -> infected_asympt
  
  infected %>%
    subset(asympt==0) %>%
    mutate(
           VL_rise = round(runif(sum(asympt==0), 1, runif(1,4.7,6.2))),
           VL=abs(rnorm(sum(asympt==0), 8, 1)*0.5),
           VL_waning=rep(0.35,sum(asympt==0)),
           Infectiousness = case_when(VL<3 ~ 0, 
                                      VL>=3 & VL<7 ~ 0.5, 
                                      VL>=7 ~ 1)) %>%
    dplyr::select(-asympt) -> infected_sympt
  
  df %>%
    subset(Inc.pd != Days) %>%
    mutate(Days=Days + 1) -> df
  
  list("infected_asympt"=infected_asympt,
       "infected_sympt"=infected_sympt,
       "df"=df)
}

move_hcw <- function(S.hcwNC, 
                     E.hcwNC, 
                     A.hcwNC,
                     I.hcwNC, 
                     R.hcwNC,
                     S.hcwC, 
                     E.hcwC, 
                     A.hcwC,
                     I.hcwC, 
                     R.hcwC,
                     prop.rNC, prop.hcwNC){
  
  SEAI.hcwNC <- c(S.hcwNC$ID,E.hcwNC$ID,A.hcwNC$ID,I.hcwNC$ID)
  SEAI.hcwC <- c(S.hcwC$ID,E.hcwC$ID,A.hcwC$ID,I.hcwC$ID)
  
  if(prop.rNC > prop.hcwNC){ # need to move hcw from C to NC 
    
    move <- round((prop.rNC - prop.hcwNC)*(length(SEAI.hcwC) + length(SEAI.hcwNC) + nrow(R.hcwNC) + nrow(R.hcwC)))  ##需要移动的人数并取整
    left <- length(SEAI.hcwC) + nrow(R.hcwC) - move  
    if (left==0){  
      move <- move-1  
    }
    
    if (move > length(SEAI.hcwC)){ 
      R.move <- move - length(SEAI.hcwC)  
      SEAI.move <- length(SEAI.hcwC) 
    } else{
      R.move <- 0  
      SEAI.move <- move  
    }   
    
    SEAI.hcwC %>%
      as.data.frame() %>%
      setNames("ID") %>%
      sample_n(SEAI.move,replace=FALSE) -> ID_move
    
    S.hcwC %>%
      subset(ID %in% ID_move$ID) %>%
      bind_rows(S.hcwNC) -> S.hcwNC  
    
    S.hcwC %>%
      subset(!(ID %in% ID_move$ID)) -> S.hcwC  
    
    E.hcwC %>%
      subset(ID %in% ID_move$ID) %>%
      bind_rows(E.hcwNC) -> E.hcwNC
    
    E.hcwC %>%
      subset(!(ID %in% ID_move$ID)) -> E.hcwC
    
    A.hcwC %>%
      subset(ID %in% ID_move$ID) %>%
      bind_rows(A.hcwNC) -> A.hcwNC
    
    A.hcwC %>%
      subset(!(ID %in% ID_move$ID)) -> A.hcwC
    
    I.hcwC %>%
      subset(ID %in% ID_move$ID) %>%
      bind_rows(I.hcwNC) -> I.hcwNC
    
    I.hcwC %>%
      subset(!(ID %in% ID_move$ID)) -> I.hcwC
    
    if (R.move > 0){   
      R.hcwC %>%
        sample_n(R.move,replace=FALSE) -> ID_move.R
      
      R.hcwC %>%
        subset(ID %in% ID_move.R$ID) %>%
        bind_rows(R.hcwNC) -> R.hcwNC
      
      R.hcwC %>%
        subset(!(ID %in% ID_move.R$ID)) -> R.hcwC
    }
    
    
  } else if (prop.rNC < prop.hcwNC){ # need to move hcw from NC to C
    
    move <- round((prop.hcwNC-prop.rNC)*(length(SEAI.hcwC) + length(SEAI.hcwNC) + nrow(R.hcwNC) + nrow(R.hcwC)))
    left <- length(SEAI.hcwNC) + nrow(R.hcwNC) - move
    if (left==0){
      move <- move-1
    }
    
    if (move > length(SEAI.hcwNC)){
      R.move <- move - length(SEAI.hcwNC)
      SEAI.move <- length(SEAI.hcwNC)
    } else{
      R.move <- 0
      SEAI.move <- move
    }
    
    SEAI.hcwNC %>%
      as.data.frame() %>%
      setNames("ID") %>%
      sample_n(SEAI.move,replace=FALSE) -> ID_move
    
    S.hcwNC %>%
      subset(ID %in% ID_move$ID) %>%
      bind_rows(S.hcwC) -> S.hcwC
    
    S.hcwNC %>%
      subset(!(ID %in% ID_move$ID)) -> S.hcwNC
    
    E.hcwNC %>%
      subset(ID %in% ID_move$ID) %>%
      bind_rows(E.hcwC) -> E.hcwC
    
    E.hcwNC %>%
      subset(!(ID %in% ID_move$ID)) -> E.hcwNC
    
    A.hcwNC %>%
      subset(ID %in% ID_move$ID) %>%
      bind_rows(A.hcwC) -> A.hcwC
    
    A.hcwNC %>%
      subset(!(ID %in% ID_move$ID)) -> A.hcwNC
    
    I.hcwNC %>%
      subset(ID %in% ID_move$ID) %>%
      bind_rows(I.hcwC) -> I.hcwC
    
    I.hcwNC %>%
      subset(!(ID %in% ID_move$ID)) -> I.hcwNC
    
    if (R.move > 0){
      R.hcwNC %>%
        sample_n(R.move,replace=FALSE) -> ID_move.R
      
      R.hcwNC %>%
        subset(ID %in% ID_move.R$ID) %>%
        bind_rows(R.hcwC) -> R.hcwC
      
      R.hcwNC %>%
        subset(!(ID %in% ID_move.R$ID)) -> R.hcwNC
    }
  }
  
  list(S.hcwNC,E.hcwNC,A.hcwNC,I.hcwNC,R.hcwNC,S.hcwC,E.hcwC,A.hcwC,I.hcwC,R.hcwC)
  
}

death <- function(df,mu,parms,total,t,dt){ # departure and death function
  new_dead <- data.frame()
  if (nrow(df)>0){
    for (i in 1:nrow(df)){ # loop through each person to see if they die
      death_stat <- rbinom(1,1,mu)
      if (death_stat == 1){
        new_dead <- bind_rows(new_dead,df %>% subset(ID==i))
        
      }
    }
  }

  
  df %>%
    subset(!(ID %in% new_dead$ID)) -> df 
  
  if(nrow(new_dead)>0){
  new_dead %>%
    dplyr::select(ID,VL,V_doses,V_t,VE_i,VE_s,VE_p) -> new_dead
  }
  
  list(df,nrow(new_dead),new_dead)  
  
}


covid_death <- function(df,mu,parms,total,t,dt){ 
  new_dead <- data.frame()
  if (nrow(df)>0){
    for (i in 1:nrow(df)){ # loop through each person to see if they die
      death_stat <- rbinom(1,1,mu) 
      if (death_stat == 1){
        new_dead <- bind_rows(new_dead,df %>% subset(ID==i))
      }
    }
  }
  df %>%
    subset(!(ID %in% new_dead$ID)) -> df
  
  if(nrow(new_dead)>0){
    new_dead %>%
      dplyr::select(ID,VL,V_doses,V_t,VE_i,VE_s,VE_p) -> new_dead
  }
  
  list(df,nrow(new_dead),new_dead) 
}



get_VL <- function(list, day){
  
  VL <- NULL
  
  for(j in 1:length(list)){
    if("ID" %in% names(list[[j]])){
      
      as.data.frame(list[[j]]) %>%
        dplyr::select(ID,VL) %>%
        mutate(Sim.day=day, State=names(list)[j]) -> new_VL
      
      VL <- rbind(VL, new_VL)}
  }
  
  return(VL)
  
}


assign_rooms <- function(family,R.room,int.room,t,parms){
  
  if (length(R.room)>0){
    if (int.room == 1 & t> parms[["int_time"]]){
      
      family$ResidentFlag <- ifelse(family$ResID %in% R.room,1,family$ResidentFlag)
    }
  }
    return(family)
}


expose_roommate <- function(family,S_ID,A.rNC,I.rNC){
  if (nrow(A.rNC) + nrow(I.rNC) > 0){  
    family %>%   
      subset(family$ResID %in% c(A.rNC$ID,I.rNC$ID)) -> familyAI
    family %>%  
      subset(family$Room %in% familyAI$Room) -> familyS
    
    exp <- ifelse(S_ID %in% familyS$ResID,1,0) 
  } else{
    exp <- 0
  }
  
  return(exp)
  
}

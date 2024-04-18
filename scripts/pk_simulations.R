library(tidyverse)
library(mrgsolve)
library(xgxr)
library(here)


belva.model <- mread(here('data',"pk_models", "belva_model.cpp"))
cobi.model <- mread(here('data',"pk_models", "cobi_model.cpp"))

### Belva Simulations 

n=500

cohort_4 = as_data_set(ev(ID=1:n, CMT = 1, ii = 24, until = 30*24, amt = 50)) %>% 
  realize_addl() %>% 
  mutate(TRT01P = "Belvarafenib 50mg QD")

cohort_5 = as_data_set(ev(ID=1:n, CMT = 1, ii = 12, until = 30*24, amt = 100)) %>% 
  realize_addl() %>% 
  mutate(TRT01P = "Belvarafenib 100mg BID")

cohort_6 = as_data_set(ev(ID=1:n, CMT = 1, ii = 24, until = 30*24, amt = 200)) %>% 
  realize_addl() %>% 
  mutate(TRT01P = "Belvarafenib 200mg QD")

cohort_7 = as_data_set(ev(ID=1:n, CMT = 1, ii = 12, until = 30*24, amt = 200)) %>% 
  realize_addl() %>% 
  mutate(TRT01P = "Belvarafenib 200mg BID")

cohort_8 = as_data_set(ev(ID=1:n, CMT = 1, ii = 12, until = 30*24, amt = 300)) %>% 
  realize_addl() %>% 
  mutate(TRT01P = "Belvarafenib 300mg BID")


cohort_9 = as_data_set(ev(ID=1:n, CMT = 1, ii = 12, until = 30*24, amt = 400)) %>% 
  realize_addl() %>% 
  mutate(TRT01P = "Belvarafenib 400mg BID")


df.dose = bind_rows(cohort_4,cohort_5,cohort_6,cohort_7,cohort_8,cohort_9) %>%
  group_by(ID, TRT01P) %>%
  mutate(ID = cur_group_id()) %>% ungroup() %>% 
  arrange(ID,time) 

set.seed(12345)
df.sim.out <- belva.model %>% 
  data_set(df.dose) %>% 
  carry_out(ID,amt,cmt,evid) %>% 
  mrgsim(end = 30*24, delta = 0.5, tad = TRUE, recover = c("TRT01P")) %>% 
  as_tibble()


df.sim.out %>% 
  filter(evid == 0) %>% 
  filter(time > 21*24, time < 30*24) %>% 
  write_csv(here("data/PK_variability/belva_sim_pk_data.csv"))


### Cobi Simulations
n=500

cohort_1 = as_data_set(ev(ID=1:n, cmt = 1, ii = 1, until = 30, amt = 20)) %>% 
  realize_addl() %>% 
  mutate(TRT01P = "Cobi 20mg QD")

cohort_2 = as_data_set(ev(ID=1:n, cmt = 1, ii = 2, until = 30, amt = 20)) %>% 
  realize_addl() %>% 
  mutate(TRT01P = "Cobi 20mg QOD")

cohort_3 = as_data_set(ev(ID=1:n, cmt = 1, ii = 1, until = 30, amt = 40)) %>% 
  realize_addl() %>% 
  mutate(TRT01P = "Cobi 40mg QD")

cohort_4 = as_data_set(ev(ID=1:n, cmt = 1, ii = 1, until = 30, amt = 60)) %>% 
  realize_addl() %>% 
  mutate(TRT01P = "Cobi 60mg QD")


cohort_5 <- tibble(time = c(0,48,96,168,216,264,336,384,432,504,552,600)/24) %>% 
  mutate(evid=1, cmt = 1, amt = 40) %>% 
  expand_grid(ID = 1:n) %>% 
  mutate(TRT01P = "Cobi 40mg TIW") %>% 
  arrange(ID,time)




df.dose = bind_rows(cohort_1,cohort_2,cohort_3,cohort_4, cohort_5) %>%
  group_by(ID, TRT01P) %>%
  mutate(ID = cur_group_id()) %>% ungroup() %>% 
  arrange(ID,time) 

set.seed(12345)
df.sim.out <- cobi.model %>% 
  data_set(df.dose) %>% 
  carry_out(ID,amt,cmt,evid) %>% 
  mrgsim(end = 30, delta = 1/24, tad = TRUE, recover = c("TRT01P")) %>% 
  as_tibble()

df.sim.out %>% 
  filter(evid == 0) %>% 
  filter(time > 21, time < 30) %>% 
  write_csv(here("data/PK_variability/cobi_sim_pk_data.csv"))

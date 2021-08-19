
library(deSolve)
library(tidyverse)
library(RColorBrewer)



model <- function(t, state, parms)
{ with(as.list(c(state,parms)),
       { 
         dPa <-  r*Pa*(1-(Pa+Pb)/K)  - epsA*phi*Va*Pa - g*Pa
         dIa <-  epsA*phi*Va*Pa - etaA*Ia - g*Ia
         dVa <-  etaA*Ia*betaA - m*Va  - phi*Va*Pa 
         dPb <-  r*Pb*(1-(Pb+Ia)/K)  - epsA*phi*Vb*Pb -g*Pb
         dIb <-  epsA*phi*Vb*Pb - etaA*Ib - g*Ib
         dVb <-  etaA*Ib*betaA - m*Vb  - phi*Vb*Pb 
         list(c(dPa, dVa, dIa,dPb, dVb, dIb))
       }) }



dff <- data.frame(Pa=as.numeric(),
                  Pb=as.numeric(),
                  Ia=as.numeric(),
                  Ib=as.numeric(),
                  Va=as.numeric(),
                  Vb=as.numeric(),
                  time=as.numeric(),
                  epsA=as.numeric(),
                  epsB=as.numeric(),
                  r=as.numeric(),
                  K=as.numeric(),
                  m=as.numeric(),
                  g=as.numeric(),
                  conditions = as.character())

r <- c(0.06,0.06,0.06,0.06,0.09)
k <- c(1e6,1e6,1e6,5e5,1e6)
m <- c(0.1,0.1,0.5,0.1,0.1)
g <- c(0,0.02,0,0,0)

df <- data.frame(r,k,m,g) 
df <- df %>% mutate(conditions = paste0("r = ", r, "; k = ", k, ";m = ", m, "; g = ",g))

for (i in 1:nrow(df)) {

wdf <- df[i,]  

p <- c(r = wdf$r ,
       betaB=35.8,
       m=wdf$m,
       phi=1e-7,
       K=wdf$k, 
       epsB = 0.3,
       etaB = 1/7.7,
       betaA=104,
       epsA = 0.64,
       etaA = 1/3.63,
       g = wdf$g
)



Pa = as.numeric(p['m']/(p['betaA']*p['epsA']-1)/p['phi'])
Pb = as.numeric(p['m']/(p['betaB']*p['epsB']-1)/p['phi'])
Ia = as.numeric(p['r']*Pa/p['etaA']*(1-(Pa+Pb)/p['K']))
Ib = as.numeric(p['r']*Pb/p['etaB']*(1-(Pa+Pb)/p['K']))
Va = as.numeric(p['r']/(p['epsA']*p['phi'])*(1-(Pa+Pb)/p['K']))
Vb = as.numeric(p['r']/(p['epsB']*p['phi'])*(1-(Pa+Pb)/p['K']))

s <- c(Pa=Pa*1.05,
       Va=Va*1.05,
       Ia = Ia*1.05,
       Pb=Pb*1.05,
       Vb=Vb*1.05,
       Ib = Ib*1.05)

times <- seq(0, 200000, by = 0.01)

out <- ode(y = s, times = times, func = model, parms = p)
out <- as.data.frame(out)
out <- out[seq(19000000,20000000), ]%>% 
  mutate(r = wdf$r,
         m=wdf$m,
         K=wdf$k, 
         g = wdf$g,
         conditions = wdf$conditions)
dff <- dff %>% add_row(out)
}


f_data<-dff %>% filter(conditions != "r = 0.06; k = 5e+05;m = 0.5; g = 0.02")  


f_data$conditions <- factor(f_data$conditions, ordered = TRUE, 
                          levels = c('r = 0.06; k = 1e+06;m = 0.1; g = 0',
                                     'r = 0.06; k = 1e+06;m = 0.1; g = 0.02',
                                     'r = 0.06; k = 1e+06;m = 0.5; g = 0',
                                     'r = 0.06; k = 5e+05;m = 0.1; g = 0',
                                     'r = 0.09; k = 1e+06;m = 0.1; g = 0',
                                     'r = 0.06; k = 5e5 ;m = 0.5; g = 0.02'))




dat <- dat %>% mutate(Vb = NA_real_,
                      Ib = NA_real_,
                      Pb = NA_real_) %>% filter(conditions == "r = 0.06; k = 5e5 ;m = 0.5; g = 0.02") %>% 
  select(-time,-clade)

f_data <- dff %>% mutate(time_days = time/24) %>% 
  select(-time)
#%>% 
f_data <- f_data %>%   add_row(dat)

save(f_data, file = 'c:/Users/lindellab/Google Drive/PhD/prezentatons and posters/second paper/phase paln.RData')

f_data_A <- f_data %>% 
  select(Va, Pa, conditions, time_days) %>% filter(time_days>8250)

#x <- seq(8250,8333,0.001)
f_data_A <- f_data_A %>% filter(time_days >8250 )



f_data_A <- f_data_A %>%   add_row(wdat)

f_data_A$conditions <- factor(f_data_A$conditions, ordered = TRUE, 
                              levels = c('r = 0.06; k = 1e+06;m = 0.1; g = 0',
                                         'r = 0.06; k = 1e+06;m = 0.1; g = 0.02',
                                         'r = 0.06; k = 1e+06;m = 0.5; g = 0',
                                         'r = 0.06; k = 5e+05;m = 0.1; g = 0',
                                         'r = 0.09; k = 1e+06;m = 0.1; g = 0',
                                         'r = 0.06; k = 5e5 ;m = 0.5; g = 0.02'))




cladeA_plot <- ggplot(f_data_A) +
  geom_point(aes(x = Pa, y = Va, color = conditions))+
  theme_bw()+
  scale_colour_brewer(palette = "Set1")+
  scale_x_log10()+
  scale_y_log10()

cladeB_plot <-ggplot(f_data_A)+ 
  geom_point(aes(x = Pb, y = Vb, color = conditions))+
  theme_bw()+
  scale_colour_brewer(palette = "Set1")+
  scale_x_log10()+
  scale_y_log10()

full_plot <- cladeA_plot / cladeB_plot


ggsave(cladeB_plot,height = 12, width = 24, file = 'c:/Users/lindellab/Google Drive/PhD/prezentatons and posters/second paper/supp fig6 clade B2.jpeg')
ggsave(cladeA_plot,height = 12, width = 24, file = 'c:/Users/lindellab/Google Drive/PhD/prezentatons and posters/second paper/supp fig6 clade A.svg')


library(tidyverse)
library(RColorBrewer)

#### change in stedy state clade B populations as function of virulence
r = 0.03
m=0.1
fie=1e-7
K=1e6 
betaB = 35.8
betaA = 103.78
eps = 0.64








eps <- seq(0.001,1,0.001)

df <-  as.data.frame(eps)
df <- df %>%  mutate(B = (r/(eps*fie))*(1-(m/((eps*betaB-1)*fie*K))),
                     A = (r/(eps*fie))*(1-(m/((eps*betaA-1)*fie*K))),
                     Sa = m/((betaA*eps - 1)*fie),
                     Sb = m/((betaB*eps - 1)*fie)) %>% 
              pivot_longer(c(A,B), names_to = "clade", values_to ="phage_density")
df2 <- df %>%  mutate(
                     Ha = m/((betaA*eps - 1)*fie),
                     Hb = m/((betaB*eps - 1)*fie)) %>% 
  pivot_longer(c(Ha,Hb), names_to = "Host", values_to ="cell ml")

df <- df %>% mutate(limits = if_else(clade == 'A', betaA*eps - 1, betaB*eps - 1  )) %>% 
  mutate(phage_density = if_else(limits <0 | phage_density < 0,NA_real_,phage_density))

df2 <- df2 %>% mutate(`cell ml` = if_else(`cell ml` <0 |`cell ml` > K  , K, `cell ml`)) 

full_fig <- df %>% 
  ggplot()+
  geom_line(aes(y = phage_density, x = eps, color = clade)#, size = 1
            )+
  #scale_color_manual(values =c('red', 'blue'))+
 # geom_area(mapping = aes(x = eps ,y = ifelse(eps>0.057 & eps <= 0.11 , phage_density, NA_real_)), fill ='#E41A1C', na.rm = T, alpha = 0.3)+
  #geom_area(mapping = aes(x = eps ,y = ifelse(eps>0.168 & eps<= 0.321 & clade == 'B' , phage_density, NA_real_)), fill ='#377EB8', na.rm = T, alpha = 0.3)+
 # geom_area(mapping = aes(x = eps ,y = ifelse(eps>0.019 & eps <= 0.033 , phage_density, NA_real_)), fill ='#E41A1C', na.rm = T, alpha = 0.3)+
  #geom_area(mapping = aes(x = eps ,y = ifelse(eps>0.056 & eps<= 0.095 & clade == 'B' , phage_density, NA_real_)), fill ='#377EB8', na.rm = T, alpha = 0.3)+
  
   geom_point(aes(x= 0.3, y = r/(0.3*fie)*(1-(m/((0.3*betaB-1)*fie*K )))), color = '#377EB8',
             #size = 4, 
             shape = 1)+
  geom_point(aes(x= 0.64,y = r/(0.64*fie)*(1-(m/((0.64*betaA-1)*fie*K )))), color = '#E41A1C',
             #size = 4, 
             shape = 1)+
  theme_bw()+
  #ggthemes::theme_few()+
  scale_colour_brewer(palette = "Set1")+
  scale_x_log10()+
  scale_y_log10()#+
  coord_cartesian(ylim=c(1e1,3e6))#+
  #labs(captionh="X  mesured mean virulence")+
  #theme(legend.position="bottom",
   #     plot.caption = element_text(hjust = 0.5))
  
full_fig
ggsave(full_fig,width =3 ,height =2 ,units = 'in',file ='c:/Users/imaid/Google Drive/PhD/prezentatons and posters/second paper/dens vs vir 003 05.svg',
)

cell_fig <- df2 %>% 
  ggplot()+
  geom_line(aes(y = `cell ml`, x = eps, color = Host)#, size = 1
  )+
  #scale_color_manual(values =c('red', 'blue'))+
  # geom_area(mapping = aes(x = eps ,y = ifelse(eps>0.057 & eps <= 0.11 , phage_density, NA_real_)), fill ='#E41A1C', na.rm = T, alpha = 0.3)+
  #geom_area(mapping = aes(x = eps ,y = ifelse(eps>0.168 & eps<= 0.321 & clade == 'B' , phage_density, NA_real_)), fill ='#377EB8', na.rm = T, alpha = 0.3)+
  # geom_area(mapping = aes(x = eps ,y = ifelse(eps>0.019 & eps <= 0.033 , phage_density, NA_real_)), fill ='#E41A1C', na.rm = T, alpha = 0.3)+
  #geom_area(mapping = aes(x = eps ,y = ifelse(eps>0.056 & eps<= 0.095 & clade == 'B' , phage_density, NA_real_)), fill ='#377EB8', na.rm = T, alpha = 0.3)+
  
  geom_point(aes(x= 0.3, y = m/((betaB*0.3 - 1)*fie)), color = '#377EB8',
             #size = 4, 
             shape = 1)+
  geom_point(aes(x= 0.64,y = m/((betaA*0.64 - 1)*fie)), color = '#E41A1C',
             #size = 4, 
             shape = 1)+
  theme_bw()+
  #ggthemes::theme_few()+
  scale_colour_brewer(palette = "Set1")+
  scale_x_log10()+
  scale_y_log10()#+
coord_cartesian(ylim=c(1e1,3e6))#+
#labs(captionh="X  mesured mean virulence")+
#theme(legend.position="bottom",
#     plot.caption = element_text(hjust = 0.5))

cell_fig
ggsave(cell_fig,width =3 ,height =2 ,units = 'in',file ='c:/Users/imaid/Google Drive/PhD/prezentatons and posters/second paper/cells vs vir3.svg',
)


#


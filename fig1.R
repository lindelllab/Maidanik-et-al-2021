library(tidyverse)
library(RColorBrewer)
library(scatterplot3d)
#### change in stedy state clade B populations as function of virulence
setwd("C:/Users/lindellab/Google Drive/PhD/prezentatons and posters/second paper/")
r = 0.06
m=0.1
fie=1e-7
K=1e6 
betaB = 35.8
betaA = 103.78
eps = 0.64


epsA <- seq(0.001,1,0.001)
epsB <- seq(0.001,1,0.001)



df <-  as.data.frame(crossing(epsA, epsB)) %>% filter(epsA>epsB)

df <- df %>%  mutate(B = (r/(epsB*fie))*(1-(m/((epsA*betaA-1)*fie*K))-(m/((epsB*betaB-1)*fie*K))),
                     A = (r/(epsA*fie))*(1-(m/((epsA*betaA-1)*fie*K))-(m/((epsB*betaB-1)*fie*K))),
                     Sa = m/((betaA*epsA - 1)*fie),
                     Sb = m/((betaB*epsB - 1)*fie)) %>%
              mutate(A = if_else(betaA*epsA - 1 <=0 | A<=0 ,NA_real_,A  ),
                    B = if_else(betaB*epsB - 1 <=0 | B<=0 ,NA_real_,B  ),
                    Z = B/A,
                    S_sum = Sa+Sb) %>% 
              mutate(biger_then_K = if_else(S_sum>K,1,0))

df


df2 <- df %>%  mutate(B = (r/(epsB*fie))*(1-(m/((epsA*betaA-1)*fie*K))-(m/((epsB*betaB-1)*fie*K))),
                     A = (r/(epsA*fie))*(1-(m/((epsA*betaA-1)*fie*K))-(m/((epsB*betaB-1)*fie*K))),
                     Sa = m/((betaA*epsA - 1)*fie),
                     Sb = m/((betaB*epsB - 1)*fie)) %>%
  mutate(A = if_else(betaA*epsA - 1 <=0  ,NA_real_,A  ),
         B = if_else(betaB*epsB - 1 <=0  ,NA_real_,B  ))


df3 <- df2 %>% filter(A<0)
df4 <- df2 %>% filter(B>0)
           



# pivot_longer(c(A,B), names_to = "clade", values_to ="phage_density")

# df <- df %>% mutate(A = if_else(betaA*epsA - 1 <=0 | A<=0 ,NA_real_,A  ),
#                     B = if_else(betaB*epsB - 1 <=0 | B<=0 ,NA_real_,B  ),
#                     Z = A/B) %>% 
#   mutate(phage_density = if_else(limits <0 | phage_density < 0,NA_real_,phage_density)#,
#          #Z = if_else(limits <0 | Z < 0,NA_real_,Z)
#          )




dfB <- df %>% filter(clade =='B')
#####3d----
dfB_plot <- dfB[seq(1, nrow(dfB),100),]
scatterplot3d(x=log10(dfB$epsB),y=log10(dfB$epsA),z = log10(dfB$phage_density), 
              xlab = "B virulence log10",
              ylab = "A virulence log10",
              zlab = "log10(phage/ml)",
              highlight.3d = TRUE#,
              #angle = 120
              )

library(plotly)
library(htmlwidgets)

p <- plot_ly() %>%
  add_trace( x=~log10(dfB_plot$epsB),y=~log10(dfB_plot$epsA),z = ~log10(dfB_plot$phage_density), type="scatter3d",
             marker = list(color = ~log10(dfB_plot$phage_density), colorscale = c('#FFE1A1', '#683531'), showscale = TRUE)) %>%
  layout(scene = list(
    aspectmode = "manual", aspectratio = list(x=1, y=1, z=1),
    xaxis = list(title = "B virulence log10"),
    yaxis = list(title = "A virulence log10)",
    zaxis = list(title = "log10(phage/ml)")
    )))

saveWidget(p, "3d clade B vir vs dens.html", selfcontained = F, libdir = "lib")

###2d----


df <- na.omit(df)

a <- ifelse(data$category == 0.061, "'#377EB8'", "black")

p1 <-   ggplot(df)+
  geom_point(aes(x = epsA, y=epsB, color = log10(Z)), size=0.1)+
  geom_point(aes(x = 0.64, y=0.3),color="red")+
  geom_hline(yintercept=0.061, linetype="dashed", color = "blue")+
  geom_vline(xintercept=0.021, linetype="dashed", color = "red")+
  #geom_point(aes(x = df4$epsA, y=df4$epsB), color ="grey",size =1)+
  scale_y_continuous(breaks = c(0.061,0.25,0.5,0.75,1))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.001))+
  theme(axis.text.y = element_text(color = c('blue', 'black','black', 'black', 'black')))+
  scale_color_continuous(type = "viridis",name = "Vb/Va", direction =-1, breaks = c(0.3,0.6,0.9,1.2))+theme_bw()#+scale_x_log10()+scale_y_log10()

ggsave(p1, filename = "steady state ratios bw log minimal virulence2.jpeg",width = 4.2,height = 3.5, units = 'in')  

p2 <- df %>% ggplot(aes(x = epsB, y=epsA, z = Z,colour =Z))+
  stat_contour() +
  scale_color_continuous()





#fig 4c----



epsB <- seq(0.001,1,0.001)
epsA <-seq(0.001,1,0.001)
df<-  data.frame(epsB,epsA)

df3 <- df3 %>%  mutate(#virA_0.4 = (r/(epsB*fie))*(1-(m/((0.4*betaA-1)*fie*K))-(m/((epsB*betaB-1)*fie*K))),
                       virA_0.6 = (r/(epsB*fie))*(1-(m/((0.6*betaA-1)*fie*K))-(m/((epsB*betaB-1)*fie*K))),
                       #virA_0.8 = (r/(epsB*fie))*(1-(m/((0.8*betaA-1)*fie*K))-(m/((epsB*betaB-1)*fie*K)))
                       virB_0.3 = (r/(epsA*fie))*(1-(m/((epsA*betaA-1)*fie*K))-(m/((0.3*betaB-1)*fie*K)))) %>%
              mutate(
                
                virA_0.6 = if_else(betaB*epsB - 1 <=0 | virA_0.6<=0 ,NA_real_,virA_0.6  ),
                virB_0.3 = if_else(betaA*epsA - 1 <=0 | virB_0.3<=0 ,NA_real_,virB_0.3  ))
                


df3 <- df3 %>% pivot_longer(cols = c(virA_0.4,virA_0.6,virA_0.8), names_to = 'vir_A', values_to ='phage_ml' )

p4 <- ggplot(df3)+
  geom_line(aes(x = epsB,
                y =  virA_0.6), color = '#377EB8')+
  geom_line(aes(x = epsA,
                y = virB_0.3),color = '#E41A1C')+ 
geom_point(aes(x= 0.3, y = r/(0.3*fie)*(1-(m/((0.3*betaB-1)*fie*K ))-(m/((0.64*betaA-1)*fie*K )))), color = '#377EB8',
           size = 2.2, 
           shape = 1)+
  geom_point(aes(x= 0.64,y = r/(0.64*fie)*(1-(m/((0.64*betaA-1)*fie*K ))-(m/((0.3*betaB-1)*fie*K )))), color = '#E41A1C',
             size = 2.2, 
             shape = 1)+scale_x_log10()+scale_y_log10()+ theme_bw()+
  xlab(element_blank())+
  ylab(element_blank())+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave(p4,width =3 ,height =3 ,units = 'in',file ='c:/Users/lindellab/Google Drive/PhD/prezentatons and posters/second paper/phage vs vir3.svg')

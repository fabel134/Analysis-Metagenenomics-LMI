#########################################DATOS DE San Juan Raya A. marmorata################################
#import library
library(plyr)
library(Rmisc)
library(readxl)
library(ggplot2)
library(viridis)
library(dplyr)
library(forcats)

#Choose all setting to plots
tema=theme(axis.text.x = element_text(color="black",size=12, angle=0,hjust=0.5,vjust=1.5, family = "sans" ),
           #axis.text.x = element_text(color="black",size=12, angle=90,hjust=0.5,vjust=1.5),
           axis.text.y = element_text(color="black",size=12, vjust = 1.5, family = "sans"),
           axis.title = element_text(color="black",size=12, face = "bold", family = "sans"),
           axis.title.x.bottom = element_blank(),
           panel.border =element_rect(color = "black", fill = NA),#element_blank(),
           strip.text.x = element_text(size=12, color="black",face="bold", family = "sans"),
           strip.text.y = element_text(size=12, color="black",face="bold", family = "sans"),
           strip.placement = "outside", strip.background = element_rect(fill = "white"), 
           panel.background = element_rect(fill = "white",colour = "white",size = 0.8, linetype = "solid"),
           panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),
           legend.position = "right", legend.text = element_text(color = "black",size=12, family = "sans"), 
           legend.direction = "vertical", legend.title = element_text(color = "black",size=12, face = "bold", family = "sans"),
           legend.key.size = unit(0.4,"cm"))


#Import data
setwd("~/LIM/CF_AM")

df <- read_excel("All_time_A_marmorata_T-0-6-12-18.xlsx")
#View(df)
levels(df$Tratamiento)
df$Tratamiento <- factor(df$Tratamiento, levels = c("MOCK","EV","CAMP","EV+CAMP","PFC","PFCS")) # Reordering group factor levels

k <- df %>%
  mutate(Muestreo = fct_relevel(Muestreo, 
                                "T0","T6","T12","T18")) ## reorder the time

##Estadistica descriptiva por tratamientos (Tratamiento)

##Altura
tapply(df$`Altura (cm)`,df$Tratamiento,summary)
##sd
tapply(df$`Altura (cm)`,df$Tratamiento,sd)

##Grosor
tapply(df$`Grosor (cm)`,df$Tratamiento,summary)
##sd
tapply(df$`Grosor (cm)`,df$Tratamiento,sd)

##n° de hojas 
tapply(df$`No. de Hojas`,df$Tratamiento,summary)
##sd
tapply(df$`No. de Hojas`,df$Tratamiento,sd)

###obtenemos datos por medición

#####################Altura######################################
altura=summarySE(df,measurevar=c("Altura (cm)"),groupvars=c("Tratamiento"))
altura#n=6, ya sabemos que tenemos replicas

#Evaluamos la distribución y homocedasticidad en todos los tratamientos
#### Realizamos el test de normalidad
#Ho: no hay diferencia entre la distribucion normal y esta distribucion
#Ho se rechaza con p<0.05

tapply(df$`Altura (cm)`,df$Tratamiento,shapiro.test)
##p es menor a 0.05, se rechaza la Ho (no es distribucion normal)
###Usamos prueba no paramétrica

kruskal.test(`Altura (cm)` ~ Tratamiento, data = df)##muy Significativo
#kruskal.test(`Altura (cm)` ~ Muestreo, data = df)##muy Significativo

pairwise.wilcox.test(df$`Altura (cm)`, df$Tratamiento,
                     p.adjust.method = "BH")

#pairwise.wilcox.test(df$`Altura (cm)`, df$Muestreo,
#                     p.adjust.method = "BH")

############################# Plot largo hoja ################################
##usamos reorder para ordenar ascendentemente
ggplot(k, aes(x=reorder(Tratamiento,`Altura (cm)`), y=`Altura (cm)`,fill=reorder(Tratamiento,`Altura (cm)`))) + 
  geom_boxplot()+
  #geom_dotplot(binaxis='y', stackdir='center'
   #                           ,stackratio=0.7, dotsize=0.4,position=position_dodge())+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  labs(title = substitute(bold(paste('Foliar area of six months ',italic('A. marmorata'),
                                     ' plants at greenhouse conditions.'))),
       subtitle = "Kruskal-Wallis chi-squared = 26.713, df = 5, p-value = 6.488e-05", x = "Treatment",
       y="Length (cm²)",fill="Treatment")+
  theme(text = element_text(size = 11, face = "plain",family = "Gill_Sans_MT"),##Todo el texto
                                              strip.text.x = element_text(size = 10, face = "bold",hjust = 0.5), ##Texto de titulos de graficos
                                              plot.title = element_text(hjust = 0, size = 12),
                                              plot.title.position = "plot",
                                              plot.subtitle = element_text(hjust = 0, size = 11))+
  scale_fill_viridis(alpha = 1, begin = 0, end = 0.92, direction = -1,
                     discrete = T, option = "D")+
  facet_wrap(~ Muestreo)

############################# Plot largo hoja_ solo tratamientos ################################
##usamos reorder para ordenar ascendentemente
ggplot(k, aes(x=reorder(Tratamiento,`Altura (cm)`), y=`Altura (cm)`,fill=reorder(Tratamiento,`Altura (cm)`))) + 
  geom_boxplot(outlier.size = 0, aes(fill=Tratamiento))+
  scale_fill_manual(values = c("#000000","#EB73B3","#F64971","#F8B620","#21B087","#2CB5C0"))+
  #geom_boxplot(outlier.size = 0)+ ##delete outliers
  #geom_dotplot(binaxis='y', stackdir='center'
  #                           ,stackratio=0.7, dotsize=0.4,position=position_dodge())+
  labs(title = substitute(bold(paste('Foliar area of six months ',italic('A. marmorata'),
                                     ' plants all time and conditions.'))),
       subtitle = "Kruskal-Wallis chi-squared = 39.082, df = 5, p-value = 2.287e-07", x = "Treatment",
       y="Length (cm²)",fill="Treatment")+
  theme(text = element_text(size = 11, face = "plain",family = "Gill_Sans_MT"),##Todo el texto
        strip.text.x = element_text(size = 10, face = "bold",hjust = 0.5), ##Texto de titulos de graficos
        plot.title = element_text(hjust = 0, size = 12),
        plot.title.position = "plot",
        plot.subtitle = element_text(hjust = 0, size = 11))+
    #scale_fill_viridis(alpha = 0.7, begin = 0, end = 0.9, direction = 1,
    #                 discrete = T, option = "H")+
  geom_jitter(aes(colour = Muestreo),shape=16, position=position_jitter(0.2))+
   scale_color_manual(values = c("T0"="gray80","T6"="gray60","T12"="gray40","T18"="gray30"))+ tema
 


########################### Plot Muestreo ###############################

ggplot(k, aes(x = Muestreo, y = `Altura (cm)`, fill = Muestreo, `Altura (cm)`)) + ##usamos reorder para ordenar ascendentemente
  geom_boxplot()+
  #geom_dotplot(binaxis='y', stackdir='center'
  #                           ,stackratio=0.7, dotsize=0.4,position=position_dodge())+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  labs(title = substitute(bold(paste('Foliar area of six months ',italic('A. marmorata'),
                                     ' plants at greenhouse conditions.'))),
       #subtitle = "Kruskal-Wallis chi-squared = 26.713, df = 5, p-value = 6.488e-05", x = "Treatment",
       y="Lenght (cm²)",fill="Treatment")+
  theme(text = element_text(size = 11, face = "plain",family = "Gill_Sans_MT"),##Todo el texto
        strip.text.x = element_text(size = 10, face = "bold",hjust = 0.5), ##Texto de titulos de graficos
        plot.title = element_text(hjust = 0, size = 12),
        plot.title.position = "plot",
        plot.subtitle = element_text(hjust = 0, size = 11))+
  scale_fill_viridis(alpha = 1, begin = 0, end = 0.92, direction = -1,
                     discrete = T, option = "D")+
  facet_wrap( ~ df$Tratamiento)



#####################Grosor######################################

Grosor=summarySE(df,measurevar=c("Grosor (cm)"),groupvars=c("Tratamiento"))
Grosor#n=6, ya sabemos que tenemos replicas

#Evaluamos la distribución y homocedasticidad en todos los tratamientos
#### Realizamos el test de normalidad
#Ho: no hay diferencia entre la distribucion normal y esta distribucion
#Ho se rechaza con p<0.05
tapply(df$`Grosor (cm)`,df$Tratamiento,shapiro.test)
#NO TODOS SON NORMALES  


kruskal.test(`Grosor (cm)` ~ Tratamiento, data = df)

##muy Significativo

pairwise.wilcox.test(df$`Grosor (cm)`, df$Tratamiento,
                     p.adjust.method = "BH")



##usamos reorder para ordenar ascendentemente
ggplot(k, aes(x=reorder(Tratamiento,`Grosor (cm)`), y=`Grosor (cm)`,fill=reorder(Tratamiento,`Grosor (cm)`))) + 
  geom_boxplot()+
  #geom_dotplot(binaxis='y', stackdir='center'
  #                           ,stackratio=0.7, dotsize=0.4,position=position_dodge())+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  labs(title = substitute(bold(paste('Foliar area of six months ',italic('A. marmorata'),
                                     ' plants at greenhouse conditions.'))),
       subtitle = "Kruskal-Wallis chi-squared = 43.831, df = 5, p-value = 2.507e-08", x = "Treatment",
       y="Leaf width (cm)",fill="Treatment")+
  theme(text = element_text(size = 11, face = "plain",family = "Gill_Sans_MT"),##Todo el texto
        strip.text.x = element_text(size = 10, face = "bold",hjust = 0.5), ##Texto de titulos de graficos
        plot.title = element_text(hjust = 0, size = 12),
        plot.title.position = "plot",
        plot.subtitle = element_text(hjust = 0, size = 11))+
  scale_fill_viridis(alpha = 1, begin = 0, end = 0.92, direction = -1,
                     discrete = T, option = "D")+
  facet_wrap(~ Muestreo)

############################# Plot grosor hoja_ solo tratamientos ################################
##usamos reorder para ordenar ascendentemente
ggplot(k, aes(x=reorder(Tratamiento,`Grosor (cm)`), y=`Grosor (cm)`,fill=reorder(Tratamiento,`Grosor (cm)`))) + 
  geom_boxplot(outlier.size = 0, aes(fill=Tratamiento))+
  scale_fill_manual(values = c("#000000","#EB73B3","#F64971","#F8B620","#21B087","#2CB5C0"))+
  #geom_boxplot(outlier.size = 0)+ ##delete outliers
  #geom_dotplot(binaxis='y', stackdir='center'
  #                           ,stackratio=0.7, dotsize=0.4,position=position_dodge())+
  labs(title = substitute(bold(paste('Foliar area of six months ',italic('A. marmorata'),
                                     ' plants all time and conditions.'))),
       subtitle = "Kruskal-Wallis chi-squared = 21.107, df = 5, p-value = 0.0007732", x = "Treatment",
       y="Leaf width (cm)",fill="Treatment")+
  theme(text = element_text(size = 11, face = "plain",family = "Gill_Sans_MT"),##Todo el texto
        strip.text.x = element_text(size = 10, face = "bold",hjust = 0.5), ##Texto de titulos de graficos
        plot.title = element_text(hjust = 0, size = 12),
        plot.title.position = "plot",
        plot.subtitle = element_text(hjust = 0, size = 11))+
  #scale_fill_viridis(alpha = 0.7, begin = 0, end = 0.9, direction = 1,
  #                 discrete = T, option = "H")+
  geom_jitter(aes(colour = Muestreo),shape=16, position=position_jitter(0.2))+
  scale_color_manual(values = c("T0"="gray80","T6"="gray60","T12"="gray40","T18"="gray30"))+ tema


##### Grosor por muestreo ####

ggplot(k, aes(x = Muestreo, y = `Grosor (cm)`, fill = Muestreo, `Grosor (cm)`)) + ##usamos reorder para ordenar ascendentemente
  geom_boxplot()+
  #geom_dotplot(binaxis='y', stackdir='center'
  #                           ,stackratio=0.7, dotsize=0.4,position=position_dodge())+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  labs(title = substitute(bold(paste('Foliar area of six months ',italic('A. marmorata'),
                                     ' plants at greenhouse conditions.'))),
       #subtitle = "Kruskal-Wallis chi-squared = 26.713, df = 5, p-value = 6.488e-05", x = "Treatment",
       y="width cm",fill="Treatment")+
  theme(text = element_text(size = 11, face = "plain",family = "Gill_Sans_MT"),##Todo el texto
        strip.text.x = element_text(size = 10, face = "bold",hjust = 0.5), ##Texto de titulos de graficos
        plot.title = element_text(hjust = 0, size = 12),
        plot.title.position = "plot",
        plot.subtitle = element_text(hjust = 0, size = 11))+
  scale_fill_viridis(alpha = 1, begin = 0, end = 0.92, direction = -1,
                     discrete = T, option = "D")+
  facet_wrap( ~ df$Tratamiento)

#####################No. de Hojas######################################
Nhojas=summarySE(df,measurevar=c("No. de Hojas"),groupvars=c("Tratamiento"))
Nhojas#n=6, ya sabemos que tenemos replicas

#Evaluamos la distribución y homocedasticidad en todos los tratamientos
#### Realizamos el test de normalidad
#Ho: no hay diferencia entre la distribucion normal y esta distribucion
#Ho se rechaza con p<0.05
tapply(df$`No. de Hojas`,df$Tratamiento,shapiro.test)
#NO TODOS SON NORMALES  

kruskal.test(`No. de Hojas` ~ Tratamiento, data = df)
##muy Significativo

pairwise.wilcox.test(df$`No. de Hojas`, df$Tratamiento,
                     p.adjust.method = "BH")
#distrib
##usamos reorder para ordenar ascendentemente
ggplot(k, aes(x=reorder(Tratamiento,`No. de Hojas`), y=`No. de Hojas`,fill=reorder(Tratamiento,`No. de Hojas`))) + 
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center'
                             ,stackratio=0.7, dotsize=0.4,position=position_dodge())+
  #geom_jitter(shape=16, position=position_jitter(0.2))+
  labs(title = substitute(bold(paste('Foliar area of six months ',italic('A. marmorata'),
                                     ' plants at greenhouse conditions.'))),
       subtitle = "Kruskal-Wallis chi-squared = 15.787, df = 5, p-value = 0.007479", x = "Treatment",
       y="Number of leaves",fill="Treatment")+
  theme(text = element_text(size = 11, face = "plain",family = "Gill_Sans_MT"),##Todo el texto
        strip.text.x = element_text(size = 10, face = "bold",hjust = 0.5), ##Texto de titulos de graficos
        plot.title = element_text(hjust = 0, size = 12),
        plot.title.position = "plot",
        plot.subtitle = element_text(hjust = 0, size = 11))+
  scale_fill_viridis(alpha = 1, begin = 0, end = 0.92, direction = -1,
                     discrete = T, option = "D")+
  facet_wrap (~ df$Muestreo)

############################# Plot No. de hojas solo tratamientos ################################
##usamos reorder para ordenar ascendentemente

ggplot(k, aes(x=reorder(Tratamiento,`No. de Hojas`), y= `No. de Hojas`,fill=reorder(Tratamiento,`No. de Hojas`))) + 
  geom_boxplot(outlier.size = 0, aes(fill=Tratamiento))+
  scale_fill_manual(values = c("#000000","#EB73B3","#F64971","#F8B620","#21B087","#2CB5C0"))+
  #geom_boxplot(outlier.size = 0)+ ##delete outliers
  #geom_dotplot(binaxis='y', stackdir='center'
  #                           ,stackratio=0.7, dotsize=0.4,position=position_dodge())+
  labs(title = substitute(bold(paste('Foliar area of six months ',italic('A. marmorata'),
                                     ' plants all time and conditions.'))),
       subtitle = "Kruskal-Wallis chi-squared = 11.833, df = 5, p-value = 0.03715", x = "Treatment",
       y="Number of leaves",fill="Treatment")+
  theme(text = element_text(size = 11, face = "plain",family = "Gill_Sans_MT"),##Todo el texto
        strip.text.x = element_text(size = 10, face = "bold",hjust = 0.5), ##Texto de titulos de graficos
        plot.title = element_text(hjust = 0, size = 12),
        plot.title.position = "plot",
        plot.subtitle = element_text(hjust = 0, size = 11))+
  #scale_fill_viridis(alpha = 0.7, begin = 0, end = 0.9, direction = 1,
  #                 discrete = T, option = "H")+
  geom_jitter(aes(colour = Muestreo),shape=16, position=position_jitter(0.2))+
  scale_color_manual(values = c("T0"="gray80","T6"="gray60","T12"="gray40","T18"="gray30"))+tema

### numero de hojas por muestreo ##

ggplot(k, aes(x = Muestreo, y = `No. de Hojas`, fill = Muestreo, `No. de Hojas`)) + ##usamos reorder para ordenar ascendentemente
  geom_boxplot()+
  #geom_dotplot(binaxis='y', stackdir='center'
  #                           ,stackratio=0.7, dotsize=0.4,position=position_dodge())+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  labs(title = substitute(bold(paste('Foliar area of six months ',italic('A. marmorata'),
                                     ' plants at greenhouse conditions.'))),
       #subtitle = "Kruskal-Wallis chi-squared = 26.713, df = 5, p-value = 6.488e-05", x = "Treatment",
       y="Number of leaves",fill="Treatment")+
  theme(text = element_text(size = 11, face = "plain",family = "Gill_Sans_MT"),##Todo el texto
        strip.text.x = element_text(size = 10, face = "bold",hjust = 0.5), ##Texto de titulos de graficos
        plot.title = element_text(hjust = 0, size = 12),
        plot.title.position = "plot",
        plot.subtitle = element_text(hjust = 0, size = 11))+
  scale_fill_viridis(alpha = 1, begin = 0, end = 0.92, direction = -1,
                     discrete = T, option = "D")+
  facet_wrap( ~ df$Tratamiento)

#import data

#########################################DATOS DE San Juan Raya################################
library(readxl)
setwd("~/LIM/R-analisis/statics/scripts_Irving")

datos_San_Juan_Raya_mayo_2023 <- read_excel("datos San Juan Raya mayo 2023.xlsx")
#View(datos_San_Juan_Raya_mayo_2023)
levels(datos_San_Juan_Raya_mayo_2023$Tratamiento)
datos_San_Juan_Raya_mayo_2023$Tratamiento <- factor(datos_San_Juan_Raya_mayo_2023$Tratamiento,      # Reordering group factor levels
                                                    levels = c("MOCK","EV","CAMP","EV + CAMP","PFC","PFCS"))
##Estadistica descriptiva por tratamientos (Tratamiento)

##Altura
tapply(datos_San_Juan_Raya_mayo_2023$`Altura (cm)`,datos_San_Juan_Raya_mayo_2023$Tratamiento,summary)
##sd
tapply(datos_San_Juan_Raya_mayo_2023$`Altura (cm)`,datos_San_Juan_Raya_mayo_2023$Tratamiento,sd)

##Grosor
tapply(datos_San_Juan_Raya_mayo_2023$`Grosor (cm)`,datos_San_Juan_Raya_mayo_2023$Tratamiento,summary)
##sd
tapply(datos_San_Juan_Raya_mayo_2023$`Grosor (cm)`,datos_San_Juan_Raya_mayo_2023$Tratamiento,sd)

##n° de hojas 
tapply(datos_San_Juan_Raya_mayo_2023$`No. de Hojas`,datos_San_Juan_Raya_mayo_2023$Tratamiento,summary)
##sd
tapply(datos_San_Juan_Raya_mayo_2023$`No. de Hojas`,datos_San_Juan_Raya_mayo_2023$Tratamiento,sd)

###obtenemos datos por medición
library(Rmisc)

#####################Altura######################################
altura=summarySE(datos_San_Juan_Raya_mayo_2023,measurevar=c("Altura (cm)"),groupvars=c("Tratamiento"))
altura#n=6, ya sabemos que tenemos replicas

#Evaluamos la distribución y homocedasticidad en todos los tratamientos
#### Realizamos el test de normalidad
#Ho: no hay diferencia entre la distribucion normal y esta distribucion
#Ho se rechaza con p<0.05

tapply(datos_San_Juan_Raya_mayo_2023$`Altura (cm)`,datos_San_Juan_Raya_mayo_2023$Tratamiento,shapiro.test)
##p es menor a 0.05, se rechaza la Ho (no es distribucion normal)
###Usamos prueba no paramétrica

kruskal.test(`Altura (cm)` ~ Tratamiento, data = datos_San_Juan_Raya_mayo_2023)##muy Significativo

pairwise.wilcox.test(datos_San_Juan_Raya_mayo_2023$`Altura (cm)`, datos_San_Juan_Raya_mayo_2023$Tratamiento,
                     p.adjust.method = "BH")

library(ggplot2)
library(viridis)
#library(ggsignif)
##usamos reorder para ordenar ascendentemente
ggplot(datos_San_Juan_Raya_mayo_2023, aes(x=reorder(Tratamiento,`Altura (cm)`), y=`Altura (cm)`,fill=reorder(Tratamiento,`Altura (cm)`))) + 
  geom_boxplot()+
  #geom_dotplot(binaxis='y', stackdir='center'
   #                           ,stackratio=0.7, dotsize=0.4,position=position_dodge())+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  labs(title = substitute(bold(paste('Foliar area of six months ',italic('A. marmorata'),
                                     ' plants at greenhouse conditions.'))),
       subtitle = "Kruskal-Wallis chi-squared = 26.713, df = 5, p-value = 6.488e-05", x = "Treatment",
       y="area (cm²)",fill="Treatment")+
  theme(text = element_text(size = 11, face = "plain",family = "Gill_Sans_MT"),##Todo el texto
                                              strip.text.x = element_text(size = 10, face = "bold",hjust = 0.5), ##Texto de titulos de graficos
                                              plot.title = element_text(hjust = 0, size = 12),
                                              plot.title.position = "plot",
                                              plot.subtitle = element_text(hjust = 0, size = 11))+
  scale_fill_viridis(alpha = 1, begin = 0, end = 0.92, direction = -1,
                     discrete = T, option = "D")



#####################Grosor######################################
library(Rmisc)
Grosor=summarySE(datos_San_Juan_Raya_mayo_2023,measurevar=c("Grosor (cm)"),groupvars=c("Tratamiento"))
Grosor#n=6, ya sabemos que tenemos replicas


#Evaluamos la distribución y homocedasticidad en todos los tratamientos
#### Realizamos el test de normalidad
#Ho: no hay diferencia entre la distribucion normal y esta distribucion
#Ho se rechaza con p<0.05
tapply(datos_San_Juan_Raya_mayo_2023$`Grosor (cm)`,datos_San_Juan_Raya_mayo_2023$Tratamiento,shapiro.test)
#NO TODOS SON NORMALES  


kruskal.test(`Grosor (cm)` ~ Tratamiento, data = datos_San_Juan_Raya_mayo_2023)

##muy Significativo

pairwise.wilcox.test(datos_San_Juan_Raya_mayo_2023$`Grosor (cm)`, datos_San_Juan_Raya_mayo_2023$Tratamiento,
                     p.adjust.method = "BH")



##usamos reorder para ordenar ascendentemente
ggplot(datos_San_Juan_Raya_mayo_2023, aes(x=reorder(Tratamiento,`Grosor (cm)`), y=`Grosor (cm)`,fill=reorder(Tratamiento,`Grosor (cm)`))) + 
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
                     discrete = T, option = "D")

#####################No. de Hojas######################################
Nhojas=summarySE(datos_San_Juan_Raya_mayo_2023,measurevar=c("No. de Hojas"),groupvars=c("Tratamiento"))
Nhojas#n=6, ya sabemos que tenemos replicas



#Evaluamos la distribución y homocedasticidad en todos los tratamientos
#### Realizamos el test de normalidad
#Ho: no hay diferencia entre la distribucion normal y esta distribucion
#Ho se rechaza con p<0.05
tapply(datos_San_Juan_Raya_mayo_2023$`No. de Hojas`,datos_San_Juan_Raya_mayo_2023$Tratamiento,shapiro.test)
#NO TODOS SON NORMALES  


kruskal.test(`No. de Hojas` ~ Tratamiento, data = datos_San_Juan_Raya_mayo_2023)
##muy Significativo

pairwise.wilcox.test(datos_San_Juan_Raya_mayo_2023$`No. de Hojas`, datos_San_Juan_Raya_mayo_2023$Tratamiento,
                     p.adjust.method = "BH")



#distrib
library(ggplot2)
library(viridis)
#library(ggsignif)
##usamos reorder para ordenar ascendentemente
ggplot(datos_San_Juan_Raya_mayo_2023, aes(x=reorder(Tratamiento,`No. de Hojas`), y=`No. de Hojas`,fill=reorder(Tratamiento,`No. de Hojas`))) + 
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
                     discrete = T, option = "D")

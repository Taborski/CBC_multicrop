setwd("~/Desktop/stage3A/R/traitement_donnees")
source("packages.R") # charge les packages indiqués dans le fichier packages.R
rm(list=ls()) # nettoyer l'espace de travail au cas où
source("funs.R") # charge les fonctions indiquées dans le fichier funs.R (ici uniquement la fonction de chargement des données)
##### ouverture des fichiers #####
setwd("~/Desktop/Data_Stage")
ouverture_calcul_patch('variation_agregation') -> variation_agregation
ouverture_calcul_patch('variation_foraging') -> variation_foraging
ouverture_calcul_patch('variation_infection') %>% 
  filter(infection.rate !=2)-> variation_infection
ouverture_calcul_patch('variation_overwintering') -> variation_overwintering
ouverture_calcul_patch('variation_mortality') -> variation_mortality
ouverture_calcul_patch('variation_dispersal')-> variation_dispersal

setwd("~/Desktop/stage3A/R/traitement_donnees")

config_palette <-  c("#A6CEE3" ,"#1F78B4", "#B2DF8A","#33A02C" , "#FDBF6F", "#FF7F00","#CAB2D6","#6A3D9A","#2D004B") # palette bleu vert orange
limite_delai <- c(0,43)

##### AS facteurs #####
sensi_infection <- function(data.spa){
  data.spa %>% 
    filter(infection.rate == 4) %>% 
    mutate(infection=if_else(is.na(inf.date) == T,0,1)) %>% #permet de calculer cb au total de patch on été infectés
    mutate(visite=if_else(is.na(time.for.crop.loss) == T,0,1)) %>% #permet de calculer cb au total de patch on étés régulés
    
    group_by(filename.spatial, infection.rate, land.cover) %>% 
    summarise(tx_visite= sum(visite)/sum(infection)) %>% 
    ungroup() %>% 
    
    group_by(infection.rate,land.cover) %>% 
    summarise(moy=mean(tx_visite),low=moy - sd(tx_visite)/sqrt(n()), upper= moy + sd(tx_visite)/sqrt(n())) -> points1
  
  data.spa %>% # tibble temporaire pour récupérer notre valeur de simulation de base
    filter(infection.rate == 4) %>% 
    mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% #!! mettre ça partout ! 
    group_by(filename.spatial, infection.rate,land.cover) %>% 
    filter(time.for.crop.loss<151) %>% 
    
    summarise(var1= mean(time.for.crop.loss)) %>% 
    ungroup() %>% 
    
    group_by(infection.rate, land.cover) %>% 
    summarise(moy=mean(var1),low=moy - sd(var1)/sqrt(n()), upper= moy + sd(var1)/sqrt(n()) ) -> points2
  
  data.spa %>% 
    mutate(infection=if_else(is.na(inf.date) == T,0,1)) %>% #permet de calculer cb au total de patch on été infectés
    mutate(visite=if_else(is.na(time.for.crop.loss) == T,0,1)) %>% #permet de calculer cb au total de patch on étés régulés
    
    group_by(filename.spatial,infection.rate,land.cover) %>% 
    summarise(tx_visite= sum(visite)/sum(infection)) %>% 
    ungroup() %>% 
    
    group_by(infection.rate,land.cover) %>% 
    summarise(moy=mean(tx_visite),low=moy - sd(tx_visite)/sqrt(n()), upper= moy + sd(tx_visite)/sqrt(n())) %>% 
    
    ungroup() %>%
    ggplot() + 
    geom_point(mapping=aes(x=infection.rate, y=moy))+
    geom_line(mapping = aes(x =infection.rate, y = moy)) +
    geom_errorbar(mapping=aes(x=infection.rate, ymin=low, ymax=upper),width=0.5)+
    facet_wrap(.~land.cover)+
    geom_point(points1,mapping=aes(x=infection.rate,y=moy),color="red")+
    scale_colour_manual(values=config_palette)+
    scale_x_continuous(breaks = 1:9)+
    scale_y_continuous(name="Visitation rate (%)",breaks=c(0,0.25,0.5,0.75,1),labels = c(0,25,50,75,100),limits = c(0,1) )+
    theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))+
    xlab("I (%)")-> visit
  
  data.spa %>% # tibble temporaire calcule la moyenne totale pour chaque expérience puis on calcule l'écart à la moyenne
    mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% #!! mettre ça partout ! 
    group_by(filename.spatial,infection.rate,land.cover) %>% 
    filter(time.for.crop.loss<151) %>% 
    
    summarise(var1= mean(time.for.crop.loss)) %>% 
    ungroup() %>% 
    
    group_by(infection.rate,land.cover) %>% 
    summarise(moy=mean(var1),low=moy - sd(var1)/sqrt(n()), upper= moy + sd(var1)/sqrt(n()) ) %>% 
    
    ggplot() + 
    geom_point(mapping=aes(x=infection.rate, y=moy))+
    geom_line(mapping = aes(x =infection.rate, y = moy)) +
    geom_errorbar(mapping=aes(x=infection.rate, ymin=low, ymax=upper),width=0.5)+
    facet_wrap(.~land.cover)+
    geom_point(points2,mapping=aes(x=infection.rate,y=moy),color="red")+
    
    scale_colour_manual(values=config_palette)+
    scale_x_continuous(breaks = 1:9)+
    scale_y_continuous(name="Delay (days)", limits=limite_delai )+
    theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))+
    xlab("I (%)")-> delay
  
  return(arrangeGrob(visit,delay,
                      ncol = 2,
                      top=textGrob("Infection rate", gp=gpar(fontsize=15,font=8))))
}
#grid.arrange(sensi_infection(variation_infection))
sensi_overwintering <- function(data.spa){
  data.spa <- mutate(data.spa, overwintering=overwintering * 10)
  data.spa %>% 
    filter(overwintering == 1) %>% 
    mutate(infection=if_else(is.na(inf.date) == T,0,1)) %>% #permet de calculer cb au total de patch on été infectés
    mutate(visite=if_else(is.na(time.for.crop.loss) == T,0,1)) %>% #permet de calculer cb au total de patch on étés régulés
    
    group_by(filename.spatial, overwintering, land.cover) %>% 
    summarise(tx_visite= sum(visite)/sum(infection)) %>% 
    ungroup() %>% 
    
    group_by( overwintering,land.cover) %>% 
    summarise(moy=mean(tx_visite),low=moy - sd(tx_visite)/sqrt(n()), upper= moy + sd(tx_visite)/sqrt(n())) -> points1
  
  data.spa %>% # tibble temporaire pour récupérer notre valeur de simulation de base
    filter(overwintering == 1) %>% 
    mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% #!! mettre ça partout ! 
    group_by(filename.spatial, overwintering,land.cover) %>% 
    filter(time.for.crop.loss<151) %>% 
    
    summarise(var1= mean(time.for.crop.loss)) %>% 
    ungroup() %>% 
    
    group_by( overwintering, land.cover) %>% 
    summarise(moy=mean(var1),low=moy - sd(var1)/sqrt(n()), upper= moy + sd(var1)/sqrt(n()) ) -> points2
  
  data.spa %>% 
    mutate(infection=if_else(is.na(inf.date) == T,0,1)) %>% #permet de calculer cb au total de patch on été infectés
    mutate(visite=if_else(is.na(time.for.crop.loss) == T,0,1)) %>% #permet de calculer cb au total de patch on étés régulés
    
    group_by(filename.spatial,overwintering,land.cover) %>% 
    summarise(tx_visite= sum(visite)/sum(infection)) %>% 
    ungroup() %>% 
    
    group_by(overwintering, land.cover) %>% 
    summarise(moy=mean(tx_visite),low=moy - sd(tx_visite)/sqrt(n()), upper= moy + sd(tx_visite)/sqrt(n())) %>% 
    
    ungroup() %>%
    ggplot() + 
    geom_point(mapping=aes(x=overwintering, y=moy))+
    geom_line(mapping = aes(x =overwintering, y = moy)) +
    geom_errorbar(mapping=aes(x=overwintering, ymin=low, ymax=upper),width=0.5)+
    geom_point(points1, mapping=aes(x=overwintering, y=moy),color="red")+
    facet_wrap(.~land.cover)+
    scale_colour_manual(values=config_palette)+
    scale_y_continuous(name="Visitation rate (%)",breaks=c(0,0.25,0.5,0.75,1),labels = c(0,25,50,75,100),limits = c(0,1) )+
    scale_x_continuous(breaks = 1:10)+
    theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))+
    xlab("C")-> visit
  
  data.spa %>% # tibble temporaire calcule la moyenne totale pour chaque expérience puis on calcule l'écart à la moyenne
    mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% #!! mettre ça partout ! 
    group_by(filename.spatial,overwintering,land.cover) %>% 
    filter(time.for.crop.loss<151) %>% 
    
    summarise(var1= mean(time.for.crop.loss)) %>% 
    ungroup() %>% 
    
    group_by(overwintering,land.cover) %>% 
    summarise(moy=mean(var1),low=moy - sd(var1)/sqrt(n()), upper= moy + sd(var1)/sqrt(n()) ) %>% 
    
    ggplot() + 
    geom_point(mapping=aes(x=overwintering, y=moy))+
    geom_line(mapping = aes(x =overwintering, y = moy)) +
    geom_errorbar(mapping=aes(x=overwintering, ymin=low, ymax=upper),width=0.5)+
    geom_point(points2, mapping=aes(x=overwintering, y=moy),color="red")+
    facet_wrap(.~land.cover)+
    scale_x_continuous(breaks = 1:10)+
    scale_colour_manual(values=config_palette)+
    
    scale_y_continuous(name="Delay (days)", limits=limite_delai)+
    theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))+
    xlab("C")-> delay
  
  return(arrangeGrob(visit,delay,
                      ncol = 2,
                      top=textGrob("Carrying capacity", gp=gpar(fontsize=15,font=8))))
}
#grid.arrange(sensi_overwintering(variation_overwintering))
sensi_agregation <- function(data.spa){
  data.spa %>% 
    filter(agregation ==7) %>% 
    mutate(infection=if_else(is.na(inf.date) == T,0,1)) %>% #permet de calculer cb au total de patch on été infectés
    mutate(visite=if_else(is.na(time.for.crop.loss) == T,0,1)) %>% #permet de calculer cb au total de patch on étés régulés
    
    group_by(filename.spatial,agregation, land.cover) %>% 
    summarise(tx_visite= sum(visite)/sum(infection)) %>% 
    ungroup() %>% 
    
    group_by(agregation,land.cover) %>% 
    summarise(moy=mean(tx_visite),low=moy - sd(tx_visite)/sqrt(n()), upper= moy + sd(tx_visite)/sqrt(n())) -> points1
  
  data.spa %>% # tibble temporaire pour récupérer notre valeur de simulation de base
    filter(agregation == 7) %>% 
    mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% #!! mettre ça partout ! 
    group_by(filename.spatial,agregation, land.cover) %>% 
    filter(time.for.crop.loss<151) %>% 
    
    summarise(var1= mean(time.for.crop.loss)) %>% 
    ungroup() %>% 
    
    group_by(agregation, land.cover) %>% 
    summarise(moy=mean(var1),low=moy - sd(var1)/sqrt(n()), upper= moy + sd(var1)/sqrt(n()) ) -> points2
  
  data.spa %>% 
    mutate(infection=if_else(is.na(inf.date) == T,0,1)) %>% #permet de calculer cb au total de patch on été infectés
    mutate(visite=if_else(is.na(time.for.crop.loss) == T,0,1)) %>% #permet de calculer cb au total de patch on étés régulés
    
    group_by(filename.spatial,agregation,land.cover) %>% 
    summarise(tx_visite= sum(visite)/sum(infection)) %>% 
    ungroup() %>% 
    
    group_by(agregation,land.cover) %>% 
    summarise(moy=mean(tx_visite),low=moy - sd(tx_visite)/sqrt(n()), upper= moy + sd(tx_visite)/sqrt(n())) %>% 
    
    ungroup() %>%
    ggplot() + 
    geom_point(mapping=aes(x=agregation, y=moy))+
    geom_line(mapping = aes(x =agregation, y = moy)) +
    geom_errorbar(mapping=aes(x=agregation, ymin=low, ymax=upper),width=0.5)+
    geom_point(points1,mapping=aes(x=agregation,y=moy), color= "red")+
    scale_colour_manual(values=config_palette)+
    scale_x_continuous(breaks = 1:9)+
    facet_wrap(.~land.cover)+
    scale_y_continuous(name="Visitation rate (%)",breaks=c(0,0.25,0.5,0.75,1),labels = c(0,25,50,75,100),limits = c(0,1) )+
    theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))+
    xlab("g")-> visit
  
  data.spa %>% # tibble temporaire calcule la moyenne totale pour chaque expérience puis on calcule l'écart à la moyenne
    mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% #!! mettre ça partout ! 
    group_by(filename.spatial,agregation,land.cover) %>% 
    filter(time.for.crop.loss<151) %>% 
    
    summarise(var1= mean(time.for.crop.loss)) %>% 
    ungroup() %>% 
    
    group_by(agregation,land.cover) %>% 
    summarise(moy=mean(var1),low=moy - sd(var1)/sqrt(n()), upper= moy + sd(var1)/sqrt(n()) ) %>% 
    
    ggplot() + 
    geom_point(mapping=aes(x=agregation, y=moy))+
    geom_line(mapping = aes(x =agregation, y = moy)) +
    geom_errorbar(mapping=aes(x=agregation, ymin=low, ymax=upper),width=0.5)+
    geom_point(points2,mapping=aes(x=agregation,y=moy), color= "red")+
    scale_colour_manual(values=config_palette)+
    scale_x_continuous(breaks = 1:9)+
    scale_y_continuous(name="Delay (days)", limits=limite_delai)+
    facet_wrap(.~land.cover)+
    theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))+
    xlab("g")-> delay
  
  return(arrangeGrob(visit,delay,
                      ncol = 2,
                      top=textGrob("Aggregation", gp=gpar(fontsize=15,font=8))))
}
#grid.arrange(sensi_agregation(variation_agregation))
sensi_mortality <- function(data.spa){
  data.spa %>% 
    filter(mortality == 0.06) %>% 
    mutate(infection=if_else(is.na(inf.date) == T,0,1)) %>% #permet de calculer cb au total de patch on été infectés
    mutate(visite=if_else(is.na(time.for.crop.loss) == T,0,1)) %>% #permet de calculer cb au total de patch on étés régulés
    
    group_by(filename.spatial, mortality, land.cover) %>% 
    summarise(tx_visite= sum(visite)/sum(infection)) %>% 
    ungroup() %>% 
    
    group_by(mortality,land.cover) %>% 
    summarise(moy=mean(tx_visite),low=moy - sd(tx_visite)/sqrt(n()), upper= moy + sd(tx_visite)/sqrt(n())) -> points1
  
  data.spa %>% # tibble temporaire pour récupérer notre valeur de simulation de base
    filter(mortality == 0.06) %>% 
    mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% #!! mettre ça partout ! 
    group_by(filename.spatial, mortality, land.cover) %>% 
    filter(time.for.crop.loss<151) %>% 
    
    summarise(var1= mean(time.for.crop.loss)) %>% 
    ungroup() %>% 
    
    group_by( mortality, land.cover) %>% 
    summarise(moy=mean(var1),low=moy - sd(var1)/sqrt(n()), upper= moy + sd(var1)/sqrt(n()) ) -> points2
  
  data.spa %>% 
    mutate(infection=if_else(is.na(inf.date) == T,0,1)) %>% #permet de calculer cb au total de patch on été infectés
    mutate(visite=if_else(is.na(time.for.crop.loss) == T,0,1)) %>% #permet de calculer cb au total de patch on étés régulés
    
    group_by(filename.spatial,mortality,land.cover) %>% 
    summarise(tx_visite= sum(visite)/sum(infection)) %>% 
    ungroup() %>% 
    
    group_by(mortality, land.cover) %>% 
    summarise(moy=mean(tx_visite),low=moy - sd(tx_visite)/sqrt(n()), upper= moy + sd(tx_visite)/sqrt(n())) %>% 
    
    ungroup() %>%
    ggplot() + 
    geom_point(mapping=aes(x=mortality, y=moy))+
    geom_line(mapping = aes(x =mortality, y = moy)) +
    geom_errorbar(mapping=aes(x=mortality, ymin=low, ymax=upper),width=0.01)+
    geom_point(points1, mapping=aes(x= mortality,y=moy),color="red")+
    facet_wrap(.~land.cover)+
    scale_colour_manual(values=config_palette)+
    scale_y_continuous(name="Visitation rate (%)",breaks=c(0,0.25,0.5,0.75,1),labels = c(0,25,50,75,100),limits = c(0,1) )+
    theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))+
    xlab("m")-> visit
  
  data.spa %>% # tibble temporaire calcule la moyenne totale pour chaque expérience puis on calcule l'écart à la moyenne
    mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% #!! mettre ça partout ! 
    group_by(filename.spatial,mortality, land.cover) %>% 
    filter(time.for.crop.loss<151) %>% 
    
    summarise(var1= mean(time.for.crop.loss)) %>% 
    ungroup() %>% 
    
    group_by(mortality,land.cover) %>% 
    summarise(moy=mean(var1),low=moy - sd(var1)/sqrt(n()), upper= moy + sd(var1)/sqrt(n()) ) %>% 
    
    ggplot() + 
    geom_point(mapping=aes(x=mortality, y=moy))+
    geom_line(mapping = aes(x =mortality, y = moy)) +
    geom_errorbar(mapping=aes(x=mortality, ymin=low, ymax=upper),width=0.01)+
    geom_point(points2, mapping=aes(x= mortality,y=moy),color="red")+
    facet_wrap(.~land.cover)+
    scale_colour_manual(values=config_palette)+
    scale_y_continuous(name="Delay (days)", limits=limite_delai)+
    theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))+
    xlab("m")-> delay
  
  return(arrangeGrob(visit,delay,
                      ncol = 2,
                      top=textGrob(" NE Mortality", gp=gpar(fontsize=15,font=8))))
}
#grid.arrange(sensi_mortality(variation_mortality))
sensi_foraging.freq <- function(data.spa){
  data.spa %>% 
    filter(foraging.freq ==1) %>% 
    mutate(infection=if_else(is.na(inf.date) == T,0,1)) %>% #permet de calculer cb au total de patch on été infectés
    mutate(visite=if_else(is.na(time.for.crop.loss) == T,0,1)) %>% #permet de calculer cb au total de patch on étés régulés
    
    group_by(filename.spatial,foraging.freq, land.cover) %>% 
    summarise(tx_visite= sum(visite)/sum(infection)) %>% 
    ungroup() %>% 
    
    group_by(foraging.freq,land.cover) %>% 
    summarise(moy=mean(tx_visite),low=moy - sd(tx_visite)/sqrt(n()), upper= moy + sd(tx_visite)/sqrt(n())) -> points1
  
  data.spa %>% 
    mutate(infection=if_else(is.na(inf.date) == T,0,1)) %>% #permet de calculer cb au total de patch on été infectés
    mutate(visite=if_else(is.na(time.for.crop.loss) == T,0,1)) %>% #permet de calculer cb au total de patch on étés régulés
    
    group_by(filename.spatial,foraging.freq, land.cover) %>% 
    summarise(tx_visite= sum(visite)/sum(infection)) %>% 
    ungroup() %>% 
    
    group_by(foraging.freq,land.cover) %>% 
    summarise(moy=mean(tx_visite),low=moy - sd(tx_visite)/sqrt(n()), upper= moy + sd(tx_visite)/sqrt(n())) %>% 
    ungroup() %>%
    ggplot() + 
    geom_point(mapping=aes(x=foraging.freq, y=moy))+
    geom_line(mapping = aes(x =foraging.freq, y = moy)) +
    geom_errorbar(mapping=aes(x=foraging.freq, ymin=low, ymax=upper),width=0.5)+
    geom_point(points1,mapping=aes(x=foraging.freq, y=moy),color="red")+
    scale_x_continuous(breaks=1:12)+
    facet_wrap(.~land.cover)+
    scale_y_continuous(name="Visitation rate (%)",breaks=c(0,0.25,0.5,0.75,1),labels = c(0,25,50,75,100),limits = c(0,1) )+
    theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))+
    xlab("f")-> visit
  
  data.spa %>% # tibble temporaire pour récupérer notre valeur de simulation de base
    filter(foraging.freq == 1) %>% 
    mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% #!! mettre ça partout ! 
    group_by(filename.spatial,foraging.freq, land.cover) %>% 
    filter(time.for.crop.loss<151) %>% 
    
    summarise(var1= mean(time.for.crop.loss)) %>% 
    ungroup() %>% 
    
    group_by(foraging.freq, land.cover) %>% 
    summarise(moy=mean(var1),low=moy - sd(var1)/sqrt(n()), upper= moy + sd(var1)/sqrt(n()) ) -> points2
  
  data.spa %>% #delai
    mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% #!! mettre ça partout ! 
    group_by(filename.spatial,foraging.freq, land.cover) %>% 
    filter(time.for.crop.loss<151) %>% 
    
    summarise(var1= mean(time.for.crop.loss)) %>% 
    ungroup() %>% 
    
    group_by(foraging.freq, land.cover) %>% 
    summarise(moy=mean(var1),low=moy - sd(var1)/sqrt(n()), upper= moy + sd(var1)/sqrt(n()) ) %>% 
    
    ggplot() + 
    geom_point(mapping=aes(x=foraging.freq, y=moy))+
    geom_line(mapping = aes(x =foraging.freq, y = moy)) +
    geom_errorbar(mapping=aes(x=foraging.freq, ymin=low, ymax=upper),width=0.5)+
    geom_point(points2,mapping=aes(x=foraging.freq, y=moy),color="red")+
    facet_wrap(.~land.cover)+
    scale_y_continuous(name="Delay (days)", limits=limite_delai)+
    scale_x_continuous(breaks=1:12)+
    theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))+
    xlab("f")-> delay
  
  return(arrangeGrob(visit,delay,
                      ncol = 2,
                      top=textGrob("NE movement frequency ", gp=gpar(fontsize=15,font=8))))
}
#grid.arrange(sensi_foraging.freq(variation_foraging))
sensi_dispersal <- function(data.spa){
  data.spa %>% 
    mutate(pattern= ifelse(pattern == "Oneof",paste(pattern, ability),"Random")) -> data.spa

  data.spa %>% #taux de visite
    filter(pattern =="Oneof 9") %>% 
    mutate(infection=if_else(is.na(inf.date) == T,0,1)) %>% #permet de calculer cb au total de patch on été infectés
    mutate(visite=if_else(is.na(time.for.crop.loss) == T,0,1)) %>% #permet de calculer cb au total de patch on étés régulés
    
    group_by(filename.spatial,pattern, land.cover) %>% 
    summarise(tx_visite= sum(visite)/sum(infection)) %>% 
    ungroup() %>% 
    
    group_by(pattern,land.cover) %>% 
    summarise(moy=mean(tx_visite),low=moy - sd(tx_visite)/sqrt(n()), upper= moy + sd(tx_visite)/sqrt(n())) -> points1
  
  data.spa %>% #taux de visite
    mutate(infection=if_else(is.na(inf.date) == T,0,1)) %>% #permet de calculer cb au total de patch on été infectés
    mutate(visite=if_else(is.na(time.for.crop.loss) == T,0,1)) %>% #permet de calculer cb au total de patch on étés régulés
    
    group_by(filename.spatial,pattern, land.cover) %>% 
    summarise(tx_visite= sum(visite)/sum(infection)) %>% 
    ungroup() %>% 
    
    group_by(pattern,land.cover) %>% 
    summarise(moy=mean(tx_visite),low=moy - sd(tx_visite)/sqrt(n()), upper= moy + sd(tx_visite)/sqrt(n())) %>% 
    ungroup() %>%
    ggplot() + 
    geom_point(mapping=aes(x=pattern, y=moy))+
    geom_errorbar(mapping=aes(x=pattern, ymin=low, ymax=upper),width=0.5)+
    geom_point(points1,mapping=aes(x=pattern, y=moy),color="red")+
    facet_wrap(.~land.cover)+
    scale_y_continuous(name="Visitation rate (%)",breaks=c(0,0.25,0.5,0.75,1),labels = c(0,25,50,75,100),limits = c(0,1) )+
    scale_x_discrete(limits=c("Random", "Oneof 9", "Oneof 25", "Oneof 49"))+
    theme(strip.text.x = element_text(size = 14),
          axis.text.y = element_text(size=14),
          axis.text.x = element_text(size = 14, angle = 45, hjust=1, vjust=1),
          axis.title = element_text(size=14))+
    xlab("b")-> visit
  
  data.spa %>% # delai: tibble temporaire pour récupérer notre valeur de simulation de base
    filter(pattern == "Oneof 9") %>% 
    mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% #!! mettre ça partout ! 
    group_by(filename.spatial,pattern, land.cover) %>% 
    filter(time.for.crop.loss<151) %>% 
    
    summarise(var1= mean(time.for.crop.loss)) %>% 
    ungroup() %>% 
    
    group_by(pattern, land.cover) %>% 
    summarise(moy=mean(var1),low=moy - sd(var1)/sqrt(n()), upper= moy + sd(var1)/sqrt(n()) ) -> points2
  
  data.spa %>% #delai
    mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% #!! mettre ça partout ! 
    group_by(filename.spatial,pattern, land.cover) %>% 
    filter(time.for.crop.loss<151) %>% 
    
    summarise(var1= mean(time.for.crop.loss)) %>% 
    ungroup() %>% 
    
    group_by(pattern, land.cover) %>% 
    summarise(moy=mean(var1),low=moy - sd(var1)/sqrt(n()), upper= moy + sd(var1)/sqrt(n()) ) %>% 
    
    ggplot() + 
    geom_point(mapping=aes(x=pattern, y=moy))+
    geom_errorbar(mapping=aes(x=pattern, ymin=low, ymax=upper),width=0.5)+
    geom_point(points2,mapping=aes(x=pattern, y=moy),color="red")+
    facet_wrap(.~land.cover)+
    scale_y_continuous(name="Delay (days)", limits=limite_delai)+
    scale_x_discrete(limits=c("Random", "Oneof 9", "Oneof 25", "Oneof 49"))+
    theme(strip.text.x = element_text(size = 14),
          axis.text.y = element_text(size=14),
          axis.text.x = element_text(size = 14, angle = 45, hjust=1, vjust=1),
          axis.title = element_text(size=14))+
    xlab("b")-> delay
  
  return(arrangeGrob(visit,delay,
                     ncol = 2,
                     top=textGrob("NE behaviour ", gp=gpar(fontsize=15,font=8))))
}
#grid.arrange(sensi_dispersal(variation_dispersal))


x11()
analyse_sensi_par_crop <- function(config='all'){
  if (config != "all") {
    sensi_foraging.freq(filter(variation_foraging, crop.order ==config)) -> foraging
    sensi_agregation(filter(variation_agregation, crop.order ==config)) -> agregation
    sensi_mortality(filter(variation_mortality, crop.order ==config)) -> mortality
    sensi_overwintering(filter(variation_overwintering, crop.order ==config)) -> overwintering
    sensi_infection(filter(variation_infection, crop.order ==config))-> infection
    sensi_dispersal(filter(variation_dispersal, crop.order == config)) -> dispersion
    texte <- config
  }
  else {
    sensi_foraging.freq(variation_foraging) -> foraging
    sensi_agregation(variation_agregation) -> agregation
    sensi_mortality(variation_mortality) -> mortality
    sensi_overwintering(variation_overwintering) -> overwintering
    sensi_infection(variation_infection)-> infection
    sensi_dispersal(variation_dispersal) -> dispersion
    texte <- "toutes les config"
  }
  grid.arrange(overwintering,infection,agregation,mortality,foraging, dispersion,ncol=2,nrow=3,top=texte)
}
analyse_sensi_par_crop()

x11()
sensi_foraging.freq_crop.order <- function(data.spa){
  data.spa %>% #visit
    mutate(infection=if_else(is.na(inf.date) == T,0,1)) %>% #permet de calculer cb au total de patch on été infectés
    mutate(visite=if_else(is.na(time.for.crop.loss) == T,0,1)) %>% #permet de calculer cb au total de patch on étés régulés
    
    group_by(filename.spatial,foraging.freq, land.cover,crop.order) %>% 
    summarise(tx_visite= sum(visite)/sum(infection)) %>% 
    ungroup() %>% 
    
    group_by(foraging.freq,land.cover, crop.order) %>% 
    summarise(moy=mean(tx_visite),low=moy - sd(tx_visite)/sqrt(n()), upper= moy + sd(tx_visite)/sqrt(n())) %>% 
    
    ungroup() %>%
    ggplot() + 
    geom_vline(aes(xintercept=1),color="red",alpha=0.3)+
    geom_point(mapping=aes(x=foraging.freq, y=moy, color=crop.order))+
    geom_line(mapping = aes(x =foraging.freq, y = moy, color=crop.order)) +
    geom_errorbar(mapping=aes(x=foraging.freq, ymin=low, ymax=upper, color=crop.order),width=0.5)+
    scale_x_continuous(breaks=1:12)+
    facet_wrap(.~land.cover)+
    scale_y_continuous(name="Visitation rate", limits=c(0,1))+
    scale_colour_manual(values=config_palette, limits = c("1-2-3","1-3-2","2-1-3","2-3-1","3-1-2","3-2-1","random","random9","random25"))+
    theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))+
    xlab("foraging frequency")-> visit
  
  data.spa %>% #delai
    mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% 
    group_by(filename.spatial,foraging.freq, land.cover,crop.order) %>% 
    filter(time.for.crop.loss<151) %>% 
    
    summarise(var1= mean(time.for.crop.loss)) %>% 
    ungroup() %>% 
    
    group_by(foraging.freq, land.cover,crop.order) %>% 
    summarise(moy=mean(var1),low=moy - sd(var1)/sqrt(n()), upper= moy + sd(var1)/sqrt(n()) ) %>% 
    
    ggplot() + 
    geom_vline(aes(xintercept=1),color="red",alpha=0.3)+
    geom_point(mapping=aes(x=foraging.freq, y=moy, color=crop.order))+
    geom_line(mapping = aes(x =foraging.freq, y = moy, color=crop.order)) +
    geom_errorbar(mapping=aes(x=foraging.freq, ymin=low, ymax=upper, color=crop.order),width=0.5)+
    facet_wrap(.~land.cover)+
    scale_y_continuous(name="Delay (days)", limits=limite_delai)+
    scale_colour_manual(values=config_palette, limits = c("1-2-3","1-3-2","2-1-3","2-3-1","3-1-2","3-2-1","random","random9","random25"))+
    scale_x_continuous(breaks=1:12)+
    theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))+
    xlab("foraging frequency")-> delay
  
  return(arrangeGrob(visit,delay,
                     ncol = 2,
                     top=textGrob("FORAGING FREQUENCY ", gp=gpar(fontsize=15,font=8))))
}
grid.arrange(sensi_foraging.freq_crop.order(variation_foraging),top="par config")

x11()
dev.print(device = jpeg, file = "images/AS_tot_all.jpg",width=1920, height=1080)
print('finito')

##### calcul de la croploss sans regul pour la méthode 2 #####
variation_infection %>% 
  select(infection.rate,land.cover,crop.order,inf.date,sensibility.d1,croploss) %>% 
  filter(infection.rate == 4) %>% 
  mutate(croploss_16_2.2=croploss_function(inf.date,sensibility.d1,2,16,2.2,2.2*720/3.2)) %>% 
  mutate(croploss_25_2.2=croploss_function(inf.date,sensibility.d1,2,25,2.2,2.2*720/3.2)) %>%
  mutate(croploss_33_2.2=croploss_function(inf.date,sensibility.d1,2,33,2.2,2.2*720/3.2)) %>%
  
  mutate(croploss_16_3.2=croploss_function(inf.date,sensibility.d1,2,16,3.2,3.2*720/3.2)) %>% 
  mutate(croploss_25_3.2=croploss_function(inf.date,sensibility.d1,2,25,3.2,3.2*720/3.2)) %>%
  mutate(croploss_33_3.2=croploss_function(inf.date,sensibility.d1,2,33,3.2,3.2*720/3.2)) %>%
  
  mutate(croploss_16_4.2=croploss_function(inf.date,sensibility.d1,2,16,3.2,4.2*720/3.2)) %>% 
  mutate(croploss_25_4.2=croploss_function(inf.date,sensibility.d1,2,25,3.2,4.2*720/3.2)) %>% 
  mutate(croploss_33_4.2=croploss_function(inf.date,sensibility.d1,2,33,3.2,4.2*720/3.2)) %>% 
  gather("croploss_16_2.2","croploss_16_3.2","croploss_16_4.2",
         "croploss_25_2.2","croploss_25_3.2","croploss_25_4.2",
         "croploss_33_2.2","croploss_33_3.2","croploss_33_4.2",
         "croploss",key=method,value=croploss_final) %>% 
  mutate(method=factor(method,levels=c("croploss_16_2.2","croploss_16_3.2","croploss_16_4.2",
                                       "croploss_25_2.2","croploss_25_3.2","croploss_25_4.2",
                                       "croploss_33_2.2","croploss_33_3.2","croploss_33_4.2",
                                       "croploss"))) %>% 
  group_by(infection.rate,method) %>% 
  summarise(moy=mean(croploss_final)) -> selection

variation_infection %>% 
  select(infection.rate,land.cover,crop.order,inf.date,sensibility.d1,croploss) %>% 
  mutate(croploss_16_2.2=croploss_function(inf.date,sensibility.d1,2,16,2.2,2.2*720/3.2)) %>% 
  mutate(croploss_25_2.2=croploss_function(inf.date,sensibility.d1,2,25,2.2,2.2*720/3.2)) %>%
  mutate(croploss_33_2.2=croploss_function(inf.date,sensibility.d1,2,33,2.2,2.2*720/3.2)) %>%
  
  mutate(croploss_16_3.2=croploss_function(inf.date,sensibility.d1,2,16,3.2,3.2*720/3.2)) %>% 
  mutate(croploss_25_3.2=croploss_function(inf.date,sensibility.d1,2,25,3.2,3.2*720/3.2)) %>%
  mutate(croploss_33_3.2=croploss_function(inf.date,sensibility.d1,2,33,3.2,3.2*720/3.2)) %>%
  
  mutate(croploss_16_4.2=croploss_function(inf.date,sensibility.d1,2,25,3.2,4.2*720/3.2)) %>% 
  mutate(croploss_25_4.2=croploss_function(inf.date,sensibility.d1,2,25,3.2,4.2*720/3.2)) %>% 
  mutate(croploss_33_4.2=croploss_function(inf.date,sensibility.d1,2,33,3.2,4.2*720/3.2)) %>% 
  gather("croploss_16_2.2","croploss_16_3.2","croploss_16_4.2",
         "croploss_25_2.2","croploss_25_3.2","croploss_25_4.2",
         "croploss_33_2.2","croploss_33_3.2","croploss_33_4.2",
         "croploss",key=method,value=croploss_final) %>% 
  mutate(method=factor(method,levels=c("croploss_16_2.2","croploss_16_3.2","croploss_16_4.2",
                                       "croploss_25_2.2","croploss_25_3.2","croploss_25_4.2",
                                       "croploss_33_2.2","croploss_33_3.2","croploss_33_4.2",
                                       "croploss"))) %>% 
  group_by(infection.rate,method) %>% 
  summarise(moy=mean(croploss_final)) %>% 
  
  ggplot()+
  geom_col(mapping=aes(x=infection.rate, y=moy))+
  geom_col(selection,mapping=aes(x=infection.rate, y=moy),fill="darkred")+
  facet_wrap(.~method,ncol=3)+
  ylab("moyenne de patch croploss")+
  scale_x_continuous(breaks=3:9)


##### AS methode croploss #####
AS_croploss_method <- function(data.spa){
  tibble(inf.date=c(1:50)) %>% 
    mutate(method1_90 =croploss_function(inf.date,1,1,90)) %>%
    mutate(method1_30 =croploss_function(inf.date,1,1,30)) %>% 
    
    mutate(method2_16_3.2=croploss_function(inf.date,1,2,16,3.2,3.2*720/3.2)) %>% 
    mutate(method2_25_3.2=croploss_function(inf.date,1,2,25,3.2,3.2*720/3.2)) %>%
    
    
    gather('method1_90','method1_30',
           'method2_16_3.2',
           'method2_25_3.2',key="method", value='croploss') %>% 
    ggplot()+
    geom_vline(xintercept = 26,colour="darkred", alpha = 0.5,size = 2, linetype= "dashed")+
    geom_vline(xintercept = 13,colour="darkred", alpha = 0.5,size = 2, linetype= "dashed")+
    geom_text(x=16,y=102,label="I= 8 %",colour="darkred",size=6)+
    geom_text(x=29,y=102,label="I= 4 %",colour="darkred",size=6)+
    geom_line(mapping=aes(x=inf.date, y=croploss,colour = method),size = 2)+
    theme()-> courbes 
  
  data.spa %>% 
    select(filename.spatial,land.cover,proportionSNH,crop.order,inf.date,sensibility.d1, time.for.crop.loss) %>% 
    mutate(method1_90= croploss_function(inf.date,sensibility.d1,1,90)) %>% 
    #mutate(method1_60= croploss_function(inf.date,sensibility.d1,1,60)) %>% 
    mutate(method1_30= croploss_function(inf.date,sensibility.d1,1,30)) %>% 
    
    mutate(method2_16_3.2=croploss_function(inf.date,sensibility.d1,2,16,3.2,3.2*720/3.2)) %>% 
    mutate(method2_25_3.2=croploss_function(inf.date,sensibility.d1,2,25,3.2,3.2*720/3.2)) %>%
    #mutate(method2_33_3.2=croploss_function(inf.date,sensibility.d1,2,33,3.2,3.2*720/3.2)) %>%
    
    gather('method1_90','method1_30',
           'method2_16_3.2',
           'method2_25_3.2',
           #'method2_33_3.2','method1_60',
           key="method", value='croploss') %>% 
    mutate(regulation= if_else(time.for.crop.loss < 51, 90 * (exp((- 1 / 30) * ( time.for.crop.loss))),12.48 * exp((50 - time.for.crop.loss)/40)+ 4.5)) %>% 
    mutate(regulation = if_else(is.na(regulation)==T,0, regulation)) %>% 
    mutate(yield = 100 * (1 - croploss /100 * (1 - regulation / 100))) -> stock
  
  stock %>% 
    group_by(filename.spatial,proportionSNH, land.cover, method) %>% 
    summarise( var1=mean(croploss * (1 - regulation/100))) %>%
    ungroup() %>% 
    
    group_by(proportionSNH,land.cover,method) %>% 
    summarise(moy=mean(var1), low= moy - sd(var1)/sqrt(n()), upper=moy+ sd(var1)/sqrt(n()), deviation= sd(var1)) %>% 
    ungroup() -> landCover_loss
  
  stock %>% 
    group_by(filename.spatial,proportionSNH,method) %>% 
    summarise( var1=mean(croploss * (1 - regulation/100))) %>%
    mutate(land.cover = ' Landscape') %>% 
    ungroup() %>% 
    
    group_by(proportionSNH,land.cover,method) %>% 
    summarise(moy=mean(var1), low= moy - sd(var1)/sqrt(n()), upper=moy+ sd(var1)/sqrt(n()), deviation= sd(var1)) %>% 
    ungroup() -> tot_loss
  
  
  ggplot(bind_rows(landCover_loss, tot_loss))+
    theme_bw()+
    geom_point(mapping = aes(x = proportionSNH, y = moy , colour = method)) +
    geom_line(mapping = aes(x = proportionSNH, y = moy , colour = method)) +
    geom_errorbar(mapping = aes( x = proportionSNH, ymin = low ,  ymax = upper,  colour = method),width=largeur) +
    scale_y_continuous(name = "Patch crop loss (%)") +
    xlab('') +
    #scale_colour_manual(values=config_palette, limits = c("1-2-3","1-3-2","2-1-3","2-3-1","3-1-2","3-2-1","random","random9","random25"))+
    facet_wrap(.~land.cover, ncol=4)+
    labs(fill = "Crop configuration") +
    theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none")-> graphe
  
  grid.arrange(courbes,graphe,ncol=2)
}

AS_croploss_method(data.spatial2)

print('finito')

##### pour les 3 configs #####
petite_palette <- c("#3fa9e1" , "#87ca4b","#6A3D9A")
### infection ###
variation_infection %>% 
  filter(crop.order=="1-2-3"| crop.order == "2-1-3"| crop.order == "random9") -> data.spa

data.spa %>% 
  mutate(infection=if_else(is.na(inf.date) == T,0,1)) %>% #permet de calculer cb au total de patch on été infectés
  mutate(visite=if_else(is.na(time.for.crop.loss) == T,0,1)) %>% #permet de calculer cb au total de patch on étés régulés
  
  group_by(filename.spatial,infection.rate,crop.order,land.cover) %>% 
  summarise(tx_visite= sum(visite)/sum(infection)) %>% 
  ungroup() %>% 
  
  group_by(infection.rate,crop.order,land.cover) %>% 
  summarise(moy=mean(tx_visite),low=moy - sd(tx_visite)/sqrt(n()), upper= moy + sd(tx_visite)/sqrt(n())) %>% 
  
  ungroup() %>%
  ggplot() + 
  geom_vline(xintercept = 4, colour = "darkred")+
  
  geom_point(mapping=aes(x=infection.rate, y=moy,colour=crop.order))+
  geom_line(mapping = aes(x =infection.rate, y = moy, colour=crop.order)) +
  geom_errorbar(mapping=aes(x=infection.rate, ymin=low, ymax=upper,colour=crop.order),width=0.5)+
  facet_wrap(.~land.cover)+
  #geom_point(points1,mapping=aes(x=infection.rate,y=moy),color="red")+
  
  scale_colour_manual(values=petite_palette)+
  scale_x_continuous(breaks = 1:9)+
  scale_y_continuous(name="Visitation rate", limits=c(0,1))+
  theme(strip.text.x = element_text(size = 14),
        axis.text = element_text(size=14),
        axis.title = element_text(size=14))+
  xlab("Infection rate (%)")-> visit

data.spa %>% # tibble temporaire calcule la moyenne totale pour chaque expérience puis on calcule l'écart à la moyenne
  mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% #!! mettre ça partout ! 
  group_by(filename.spatial,infection.rate,crop.order,land.cover) %>% 
  filter(time.for.crop.loss<151) %>% 
  
  summarise(var1= mean(time.for.crop.loss)) %>% 
  ungroup() %>% 
  
  group_by(infection.rate,land.cover,crop.order) %>% 
  summarise(moy=mean(var1),low=moy - sd(var1)/sqrt(n()), upper= moy + sd(var1)/sqrt(n()) ) %>% 
  
  ggplot() + 
  geom_vline(xintercept = 4, colour = "darkred")+
  
  geom_point(mapping=aes(x=infection.rate, y=moy,colour=crop.order))+
  geom_line(mapping = aes(x =infection.rate, y = moy, colour=crop.order)) +
  geom_errorbar(mapping=aes(x=infection.rate, ymin=low, ymax=upper,colour=crop.order),width=0.5)+
  facet_wrap(.~land.cover)+
  #geom_point(points2,mapping=aes(x=infection.rate,y=moy),color="red")+
  
  scale_x_continuous(breaks = 1:9)+
  scale_y_continuous(name="Delay (days)", limits=limite_delai)+
  scale_colour_manual(values=petite_palette, limits = c("1-2-3","2-1-3","random9"))+
  theme(strip.text.x = element_text(size = 14),
        axis.text = element_text(size=14),
        axis.title = element_text(size=14))+
  theme(legend.position="none")+
  xlab("Infection rate (%)")-> delay

legend <- get_legend(visit)

arrangeGrob(visit+ theme(legend.position = "none"),delay,
            ncol = 2,
            top=textGrob("INFECTION RATE", gp=gpar(fontsize=15,font=8))) ->infection

### agregation ###
variation_agregation %>% 
  filter(crop.order=="1-2-3"| crop.order == "2-1-3"| crop.order == "random9") -> data.spa

data.spa %>% 
  mutate(infection=if_else(is.na(inf.date) == T,0,1)) %>% #permet de calculer cb au total de patch on été infectés
  mutate(visite=if_else(is.na(time.for.crop.loss) == T,0,1)) %>% #permet de calculer cb au total de patch on étés régulés
  
  group_by(filename.spatial,agregation,crop.order,land.cover) %>% 
  summarise(tx_visite= sum(visite)/sum(infection)) %>% 
  ungroup() %>% 
  
  group_by(agregation,crop.order,land.cover) %>% 
  summarise(moy=mean(tx_visite),low=moy - sd(tx_visite)/sqrt(n()), upper= moy + sd(tx_visite)/sqrt(n())) %>% 
  
  ungroup() %>%
  ggplot() + 
  geom_vline(xintercept = 7, colour = "darkred")+
  
  geom_point(mapping=aes(x=agregation, y=moy,colour=crop.order))+
  geom_line(mapping = aes(x =agregation, y = moy, colour=crop.order)) +
  geom_errorbar(mapping=aes(x=agregation, ymin=low, ymax=upper,colour=crop.order),width=0.5)+
  facet_wrap(.~land.cover)+
  #geom_point(points1,mapping=aes(x=agregation,y=moy),color="red")+
  
  scale_colour_manual(values=petite_palette)+
  scale_x_continuous(breaks = 1:9)+
  scale_y_continuous(name="Visitation rate", limits=c(0,1))+
  theme(strip.text.x = element_text(size = 14),
        axis.text = element_text(size=14),
        axis.title = element_text(size=14))+
  xlab("Aggregation")-> visit

data.spa %>% # tibble temporaire calcule la moyenne totale pour chaque expérience puis on calcule l'écart à la moyenne
  mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% #!! mettre ça partout ! 
  group_by(filename.spatial,agregation,crop.order,land.cover) %>% 
  filter(time.for.crop.loss<151) %>% 
  
  summarise(var1= mean(time.for.crop.loss)) %>% 
  ungroup() %>% 
  
  group_by(agregation,land.cover,crop.order) %>% 
  summarise(moy=mean(var1),low=moy - sd(var1)/sqrt(n()), upper= moy + sd(var1)/sqrt(n()) ) %>% 
  
  ggplot() + 
  geom_vline(xintercept = 7, colour = "darkred")+
  
  geom_point(mapping=aes(x=agregation, y=moy,colour=crop.order))+
  geom_line(mapping = aes(x =agregation, y = moy, colour=crop.order)) +
  geom_errorbar(mapping=aes(x=agregation, ymin=low, ymax=upper,colour=crop.order),width=0.5)+
  facet_wrap(.~land.cover)+
  #geom_point(points2,mapping=aes(x=agregation,y=moy),color="red")+
  
  scale_x_continuous(breaks = 1:9)+
  scale_y_continuous(name="Delay (days)", limits=limite_delai)+
  scale_colour_manual(values=petite_palette, limits = c("1-2-3","2-1-3","random9"))+
  theme(strip.text.x = element_text(size = 14),
        axis.text = element_text(size=14),
        axis.title = element_text(size=14))+
  theme(legend.position="none")+
  xlab("Aggregation")-> delay

arrangeGrob(visit+ theme(legend.position = "none"),delay,
            ncol = 2,
            top=textGrob("AGGREGATION", gp=gpar(fontsize=15,font=8))) ->agregation

### mortality ###
variation_mortality %>% 
  filter(crop.order=="1-2-3"| crop.order == "2-1-3"| crop.order == "random9") -> data.spa

data.spa %>% 
  mutate(infection=if_else(is.na(inf.date) == T,0,1)) %>% #permet de calculer cb au total de patch on été infectés
  mutate(visite=if_else(is.na(time.for.crop.loss) == T,0,1)) %>% #permet de calculer cb au total de patch on étés régulés
  
  group_by(filename.spatial,mortality,crop.order,land.cover) %>% 
  summarise(tx_visite= sum(visite)/sum(infection)) %>% 
  ungroup() %>% 
  
  group_by(mortality,crop.order,land.cover) %>% 
  summarise(moy=mean(tx_visite),low=moy - sd(tx_visite)/sqrt(n()), upper= moy + sd(tx_visite)/sqrt(n())) %>% 
  
  ungroup() %>%
  ggplot() + 
  geom_vline(xintercept =0.06, colour = "darkred")+
  
  geom_point(mapping=aes(x=mortality, y=moy,colour=crop.order))+
  geom_line(mapping = aes(x =mortality, y = moy, colour=crop.order)) +
  geom_errorbar(mapping=aes(x=mortality, ymin=low, ymax=upper,colour=crop.order),width=0.01)+
  facet_wrap(.~land.cover)+
  #geom_point(points1,mapping=aes(x=mortality,y=moy),color="red")+
  
  scale_colour_manual(values=petite_palette)+
  scale_y_continuous(name="Visitation rate", limits=c(0,1))+
  theme(strip.text.x = element_text(size = 14),
        axis.text = element_text(size=14),
        axis.title = element_text(size=14))+
  xlab("Mortality")-> visit

data.spa %>% # tibble temporaire calcule la moyenne totale pour chaque expérience puis on calcule l'écart à la moyenne
  mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% #!! mettre ça partout ! 
  group_by(filename.spatial,mortality,crop.order,land.cover) %>% 
  filter(time.for.crop.loss<151) %>% 
  
  summarise(var1= mean(time.for.crop.loss)) %>% 
  ungroup() %>% 
  
  group_by(mortality,land.cover,crop.order) %>% 
  summarise(moy=mean(var1),low=moy - sd(var1)/sqrt(n()), upper= moy + sd(var1)/sqrt(n()) ) %>% 
  
  ggplot() + 
  geom_vline(xintercept = 0.06, colour = "darkred")+
  
  geom_point(mapping=aes(x=mortality, y=moy,colour=crop.order))+
  geom_line(mapping = aes(x =mortality, y = moy, colour=crop.order)) +
  geom_errorbar(mapping=aes(x=mortality, ymin=low, ymax=upper,colour=crop.order),width=0.01)+
  facet_wrap(.~land.cover)+
  #geom_point(points2,mapping=aes(x=mortality,y=moy),color="red")+
  
  scale_y_continuous(name="Delay (days)", limits=limite_delai)+
  scale_colour_manual(values=petite_palette, limits = c("1-2-3","2-1-3","random9"))+
  theme(strip.text.x = element_text(size = 14),
        axis.text = element_text(size=14),
        axis.title = element_text(size=14))+
  theme(legend.position="none")+
  xlab("Mortality")-> delay

arrangeGrob(visit+ theme(legend.position = "none"),delay,
            ncol = 2,
            top=textGrob("MORTALITY", gp=gpar(fontsize=15,font=8))) ->mortality

### overwintering ###
variation_overwintering %>% 
  mutate(overwintering=overwintering*10) %>% 
  filter(crop.order=="1-2-3"| crop.order == "2-1-3"| crop.order == "random9") -> data.spa

data.spa %>% 
  mutate(infection=if_else(is.na(inf.date) == T,0,1)) %>% #permet de calculer cb au total de patch on été infectés
  mutate(visite=if_else(is.na(time.for.crop.loss) == T,0,1)) %>% #permet de calculer cb au total de patch on étés régulés
  
  group_by(filename.spatial,overwintering,crop.order,land.cover) %>% 
  summarise(tx_visite= sum(visite)/sum(infection)) %>% 
  ungroup() %>% 
  
  group_by(overwintering,crop.order,land.cover) %>% 
  summarise(moy=mean(tx_visite),low=moy - sd(tx_visite)/sqrt(n()), upper= moy + sd(tx_visite)/sqrt(n())) %>% 
  
  ungroup() %>%
  ggplot() + 
  geom_vline(xintercept =1, colour = "darkred")+
  
  geom_point(mapping=aes(x=overwintering, y=moy,colour=crop.order))+
  geom_line(mapping = aes(x =overwintering, y = moy, colour=crop.order)) +
  geom_errorbar(mapping=aes(x=overwintering, ymin=low, ymax=upper,colour=crop.order),width=0.5)+
  facet_wrap(.~land.cover)+
  #geom_point(points1,mapping=aes(x=overwintering,y=moy),color="red")+
  
  scale_colour_manual(values=petite_palette)+
  scale_y_continuous(name="Visitation rate", limits=c(0,1))+
  theme(strip.text.x = element_text(size = 14),
        axis.text = element_text(size=14),
        axis.title = element_text(size=14))+
  xlab("Overwintering")-> visit

data.spa %>% # tibble temporaire calcule la moyenne totale pour chaque expérience puis on calcule l'écart à la moyenne
  mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% #!! mettre ça partout ! 
  group_by(filename.spatial,overwintering,crop.order,land.cover) %>% 
  filter(time.for.crop.loss<151) %>% 
  
  summarise(var1= mean(time.for.crop.loss)) %>% 
  ungroup() %>% 
  
  group_by(overwintering,land.cover,crop.order) %>% 
  summarise(moy=mean(var1),low=moy - sd(var1)/sqrt(n()), upper= moy + sd(var1)/sqrt(n()) ) %>% 
  
  ggplot() + 
  geom_vline(xintercept = 1, colour = "darkred")+
  
  geom_point(mapping=aes(x=overwintering, y=moy,colour=crop.order))+
  geom_line(mapping = aes(x =overwintering, y = moy, colour=crop.order)) +
  geom_errorbar(mapping=aes(x=overwintering, ymin=low, ymax=upper,colour=crop.order),width=0.5)+
  facet_wrap(.~land.cover)+
  #geom_point(points2,mapping=aes(x=mortality,y=moy),color="red")+
  
  scale_y_continuous(name="Delay (days)", limits=limite_delai)+
  scale_colour_manual(values=petite_palette, limits = c("1-2-3","2-1-3","random9"))+
  theme(strip.text.x = element_text(size = 14),
        axis.text = element_text(size=14),
        axis.title = element_text(size=14))+
  theme(legend.position="none")+
  xlab("Overwintering")-> delay

arrangeGrob(visit+ theme(legend.position = "none"),delay,
            ncol = 2,
            top=textGrob("OVERWINTERING", gp=gpar(fontsize=15,font=8))) ->overwintering
### foraging ###
variation_foraging %>% 
  filter(crop.order=="1-2-3"| crop.order == "2-1-3"| crop.order == "random9") -> data.spa

data.spa %>% 
  mutate(infection=if_else(is.na(inf.date) == T,0,1)) %>% #permet de calculer cb au total de patch on été infectés
  mutate(visite=if_else(is.na(time.for.crop.loss) == T,0,1)) %>% #permet de calculer cb au total de patch on étés régulés
  
  group_by(filename.spatial,foraging.freq,crop.order,land.cover) %>% 
  summarise(tx_visite= sum(visite)/sum(infection)) %>% 
  ungroup() %>% 
  
  group_by(foraging.freq,crop.order,land.cover) %>% 
  summarise(moy=mean(tx_visite),low=moy - sd(tx_visite)/sqrt(n()), upper= moy + sd(tx_visite)/sqrt(n())) %>% 
  
  ungroup() %>%
  ggplot() + 
  geom_vline(xintercept =1, colour = "darkred")+
  
  geom_point(mapping=aes(x=foraging.freq, y=moy,colour=crop.order))+
  geom_line(mapping = aes(x =foraging.freq, y = moy, colour=crop.order)) +
  geom_errorbar(mapping=aes(x=foraging.freq, ymin=low, ymax=upper,colour=crop.order),width=0.5)+
  facet_wrap(.~land.cover)+
  #geom_point(points1,mapping=aes(x=foraging.freq,y=moy),color="red")+
  
  scale_colour_manual(values=petite_palette)+
  scale_y_continuous(name="Visitation rate", limits=c(0,1))+
  theme(strip.text.x = element_text(size = 14),
        axis.text = element_text(size=14),
        axis.title = element_text(size=14))+
  xlab("Foraging frequency")-> visit

data.spa %>% # tibble temporaire calcule la moyenne totale pour chaque expérience puis on calcule l'écart à la moyenne
  mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% #!! mettre ça partout ! 
  group_by(filename.spatial,foraging.freq,crop.order,land.cover) %>% 
  filter(time.for.crop.loss<151) %>% 
  
  summarise(var1= mean(time.for.crop.loss)) %>% 
  ungroup() %>% 
  
  group_by(foraging.freq,land.cover,crop.order) %>% 
  summarise(moy=mean(var1),low=moy - sd(var1)/sqrt(n()), upper= moy + sd(var1)/sqrt(n()) ) %>% 
  
  ggplot() + 
  geom_vline(xintercept = 1, colour = "darkred")+
  
  geom_point(mapping=aes(x=foraging.freq, y=moy,colour=crop.order))+
  geom_line(mapping = aes(x =foraging.freq, y = moy, colour=crop.order)) +
  geom_errorbar(mapping=aes(x=foraging.freq, ymin=low, ymax=upper,colour=crop.order),width=0.5)+
  facet_wrap(.~land.cover)+
  #geom_point(points2,mapping=aes(x=mortality,y=moy),color="red")+
  
  scale_y_continuous(name="Delay (days)", limits=limite_delai)+
  scale_colour_manual(values=petite_palette, limits = c("1-2-3","2-1-3","random9"))+
  theme(strip.text.x = element_text(size = 14),
        axis.text = element_text(size=14),
        axis.title = element_text(size=14))+
  theme(legend.position="none")+
  xlab("Foraging frequency")-> delay

arrangeGrob(visit+ theme(legend.position = "none"),delay,
            ncol = 2,
            top=textGrob("NE FORAGING FREQUENCY", gp=gpar(fontsize=15,font=8))) ->foraging

### dispersion ###
variation_dispersal %>% 
  mutate(pattern= ifelse(pattern == "Oneof",paste(pattern, ability),"Random")) %>% 
  filter(crop.order=="1-2-3"| crop.order == "2-1-3"| crop.order == "random9") -> data.spa

data.spa %>% 
  mutate(infection=if_else(is.na(inf.date) == T,0,1)) %>% #permet de calculer cb au total de patch on été infectés
  mutate(visite=if_else(is.na(time.for.crop.loss) == T,0,1)) %>% #permet de calculer cb au total de patch on étés régulés
  
  group_by(filename.spatial,pattern,crop.order,land.cover) %>% 
  summarise(tx_visite= sum(visite)/sum(infection)) %>% 
  ungroup() %>% 
  
  group_by(pattern,crop.order,land.cover) %>% 
  summarise(moy=mean(tx_visite),low=moy - sd(tx_visite)/sqrt(n()), upper= moy + sd(tx_visite)/sqrt(n())) %>% 
  
  ungroup() %>% 
  ggplot() + 
  
  geom_point(mapping=aes(x=pattern, y=moy,color = crop.order))+
  geom_errorbar(mapping=aes(x=pattern, ymin=low, ymax=upper, color = crop.order),width=0.5)+
  
  facet_wrap(.~land.cover)+
  scale_y_continuous(name="Visitation rate", limits=c(0,1))+
  scale_x_discrete(limits=c("Random", "Oneof 9", "Oneof 25", "Oneof 49"))+
  scale_colour_manual(values=petite_palette, limits = c("1-2-3","2-1-3","random9"))+
  theme(strip.text.x = element_text(size = 14),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size = 14, angle = 45, hjust=1, vjust=1),
        axis.title = element_text(size=14))+
  xlab("Dispersal pattern")-> visit

data.spa %>% # tibble temporaire calcule la moyenne totale pour chaque expérience puis on calcule l'écart à la moyenne
  mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% #!! mettre ça partout ! 
  group_by(filename.spatial,pattern,crop.order,land.cover) %>% 
  filter(time.for.crop.loss<151) %>% 
  
  summarise(var1= mean(time.for.crop.loss)) %>% 
  ungroup() %>% 
  
  group_by(pattern,land.cover,crop.order) %>% 
  summarise(moy=mean(var1),low=moy - sd(var1)/sqrt(n()), upper= moy + sd(var1)/sqrt(n()) ) %>% 
  
  ggplot() + 
  
  geom_point(mapping=aes(x=pattern, y=moy,color = crop.order))+
  geom_errorbar(mapping=aes(x=pattern, ymin=low, ymax=upper, color = crop.order),width=0.5)+
  
  facet_wrap(.~land.cover)+
  scale_y_continuous(name="Delay (days)", limits=limite_delai)+
  scale_x_discrete(limits=c("Random", "Oneof 9", "Oneof 25", "Oneof 49"))+
  scale_colour_manual(values=petite_palette, limits = c("1-2-3","2-1-3","random9"))+
  theme(strip.text.x = element_text(size = 14),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size = 14, angle = 45, hjust=1, vjust=1),
        axis.title = element_text(size=14))+
  theme(legend.position = "none")+
  xlab("Dispersal pattern")-> delay

arrangeGrob(visit+ theme(legend.position = "none"),delay,
            ncol = 2,
            top=textGrob("FORAGING PATTERN", gp=gpar(fontsize=15,font=8))) ->dispersal

### assemblage ###
grid.arrange(overwintering,infection,agregation,mortality,foraging, dispersal,ncol=2,nrow=3)

##### vrac #####
variation_foraging%>% #violin plot
  filter(land.cover=="crop 3") %>% 
  mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>%
  filter(time.for.crop.loss< 152) %>% 
  group_by(foraging.freq,crop.order) %>% 
  ggplot(mapping=aes(x=as.factor(foraging.freq), y=time.for.crop.loss))+
  geom_violin(scale="count")+
  stat_summary(fun.y=mean, geom="point",size=2)+
  facet_wrap(.~crop.order)+
  scale_y_continuous(name="delay")

variation_foraging_1crop %>% #histogramme des délais
  #filter(crop.order == "1-2-3") %>% 
  mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>%
  filter(time.for.crop.loss< 152) %>% 
  group_by(foraging.freq,time.for.crop.loss) %>% 
  summarise(nombre=n()/10) %>% 
  ggplot(mapping=aes(x=time.for.crop.loss, y=nombre))+
  scale_colour_manual(values=config_palette)+
  geom_col()+
  facet_wrap(~foraging.freq)+
  scale_y_continuous(name="nombre")+
  scale_x_continuous(breaks=seq(0,44,by=2))


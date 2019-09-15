source("packages.R") # charge les packages indiqués dans le fichier packages.R
rm(list=ls()) # nettoyer l'espace de travail au cas où
source("funs.R") # charge les fonctions indiquées dans le fichier funs.R 

##### enregistrement données au format rds ##### 
#setwd("~/Desktop/Data_Stage")
#ouverture_calcul_patch('tests')->data.spatial  #il faut le data/
#  saveRDS( file = "sim1.rds")
#setwd("~/Desktop/stage3A/R/traitement_donnees")

##### ouverture fichier rds #####
#crop1 <- ouverture_calcul_patch('data/1crop')
data.spatial_1 <- readRDS(file='simul1.rds')
#data.spatial_2 <- readRDS(file='simul2.rds')

data.spatial_1 %>% 
  filter(proportionSNH <=30) -> data.spatial2

config_palette <- c("#A6CEE3" ,"#1F78B4", "#B2DF8A","#33A02C" , "#FDBF6F", "#FF7F00","#CAB2D6","#2D004B","#6A3D9A") # palette des configurations
config_order <- c("1-2-3","1-3-2","2-1-3","2-3-1","3-1-2","3-2-1","random","random9","random25")
largeur=0.5 # utilisé pour définir la largeur générique des barres des écart types des graphes

##### parametres de l'expe ######
data.spatial2 %>% # verif du nombre de répétition
  group_by(proportionSNH,crop.order) %>% 
  summarise(iteration= round(n()/round(1089*(1-unique(proportionSNH/100))))) %>% 
  ungroup() %>% 
  mutate(rows=row_number()) %>% 
  filter(rows==1) %>% 
  select(iteration)->repet

parametres <- paste("I=",unique(data.spatial2$infection.rate),"  ",
                    "A=",unique(data.spatial2$agregation),"  ",
                    "pattern=",unique(data.spatial2$pattern),"  ",
                    "mortality=",unique(data.spatial2$mortality),"  ",
                    "repet=",repet) # String avec les différents paramètres de l'expé

##### variation SNH #####
graphe_SNH <- function(data.spa,lignes=F){
  # la fonction crée le graphe en calculant chacun des indicateurs. 
  #data.spa correspond au tableau en entrée et lignes à un booléen pour savoir si on trace la droite verticale SNH=7
  
  # l'affichage de chaue indicateur se ait en trois étapes, le calcul de celui-ci pour le paysage (1er pipe) puis de même pour chaque culture (2eme pipe). 
  # Enfin on concaténe les lignes des deux tibbles et on crée le graphe 
  
  ## visit rate ##
  data.spa %>% 
    mutate(infection=if_else(is.na(inf.date) == T,0,1)) %>% #permet de calculer cb au total de patch on été infectés
    mutate(visite=if_else(is.na(time.for.crop.loss) == T,0,1)) %>% #permet de calculer cb au total de patch on étés régulés
    group_by(filename.spatial,proportionSNH,crop.order) %>% 
    summarise(tx_visite= sum(visite)/sum(infection)) %>% 
    ungroup() %>% 
    group_by(proportionSNH, crop.order) %>% 
    summarise(moy=mean(tx_visite),low=moy - sd(tx_visite)/sqrt(n())/3, upper= moy + sd(tx_visite)/sqrt(n())/3) %>% 
    mutate(land.cover=" Landscape ") %>% 
    ungroup()-> total_visit
  
  data.spa %>% 
    mutate(infection=if_else(is.na(inf.date) == T,0,1)) %>% #permet de calculer cb au total de patch on été infectés
    mutate(visite=if_else(is.na(time.for.crop.loss) == T,0,1)) %>% #permet de calculer cb au total de patch on étés régulés
    group_by(filename.spatial,proportionSNH,crop.order,land.cover) %>% 
    summarise(tx_visite= sum(visite)/sum(infection)) %>% 
    ungroup() %>% 
    group_by(crop.order,proportionSNH,land.cover) %>% 
    summarise(moy=mean(tx_visite),
              low=moy - sd(tx_visite)/sqrt(n()),
              upper= if_else(moy + sd(tx_visite)/sqrt(n())>1,1,moy + sd(tx_visite)/sqrt(n()))) %>% 
    ungroup()-> landCover_visit
  
  bind_rows(total_visit, landCover_visit) %>% 
    arrange(land.cover,moy) %>% 
    mutate(order=row_number()) %>% 
    ggplot() +
    geom_point(mapping=aes(x=proportionSNH, y=moy,colour=crop.order))+
    geom_line(mapping = aes(x =proportionSNH, y = moy, colour = crop.order)) +
    geom_errorbar(mapping=aes(x=proportionSNH, ymin=low, ymax=upper, colour=crop.order),width=largeur)+
    scale_y_continuous(name="Visitation rate (%)",breaks=seq(0,1,by=0.25),labels = seq(0,100,by=25) )+
    scale_colour_manual(values=config_palette, limits = c("1-2-3","1-3-2","2-1-3","2-3-1","3-1-2","3-2-1","random","random9","random25"))+
    facet_wrap(.~land.cover, 
               scales="free_x", 
               ncol=4)+
    theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))+
    xlab("")-> visit_para
  
  ## delai moyen ##
  data.spa %>% # tibble temporaire calcule la moyenne totale pour chaque expérience puis on calcule l'écart à la moyenne
    mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>%  
    group_by(filename.spatial,proportionSNH, crop.order) %>% 
    filter(time.for.crop.loss<151) %>% 
    summarise(var1= mean(time.for.crop.loss)) %>% 
    ungroup() %>% 
    group_by(proportionSNH, crop.order) %>% 
    summarise(moy=mean(var1),low=moy - sd(var1)/sqrt(n()), upper= moy + sd(var1)/sqrt(n()) ) %>% 
    mutate(land.cover=" Landscape ")-> total_delai_mean
  
  data.spa %>% # calcule la moyenne par landcover pour chaque expérience puis on calcule l'écart à la moyenne
    mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% 
    group_by(filename.spatial,proportionSNH,crop.order, land.cover) %>% 
    filter(time.for.crop.loss<151) %>% 
    summarise(var1= mean(time.for.crop.loss)) %>% 
    ungroup() %>% 
    group_by(proportionSNH, crop.order, land.cover) %>% 
    summarise(moy=mean(var1),low=moy - sd(var1)/sqrt(n()), upper= moy + sd(var1)/sqrt(n()) )-> landCover_delai_mean
  
  bind_rows(total_delai_mean,landCover_delai_mean) %>% 
    ggplot() +
    geom_point(mapping=aes(y=moy, x=proportionSNH, colour=crop.order))+
    geom_line(mapping = aes(x = proportionSNH, y = moy, colour = crop.order)) +
    geom_errorbar(mapping=aes(x=proportionSNH, ymin=low, ymax=upper, colour=crop.order),show.legend=F,width=largeur)+
    scale_colour_manual(values=config_palette, limits = c("1-2-3","1-3-2","2-1-3","2-3-1","3-1-2","3-2-1","random","random9","random25"))+
    facet_wrap(.~land.cover, scales="free_x", ncol=4)+
    
    # titre et axes
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("Delay (days)") + 
    theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))+
    xlab("") -> delai_mean_para
  
  ## avg crop loss ##
  data.spa %>% 
    mutate(inf.date=ifelse(is.na(inf.date)==T,-10,inf.date)) %>% 
    filter(inf.date>0) %>% 
    group_by(filename.spatial,proportionSNH,crop.order) %>% 
    summarise( var1=mean(croploss * (1 - regulation/100))) %>% 
    ungroup() %>% 
    group_by(proportionSNH, crop.order) %>% 
    summarise(moy=mean(var1), low= moy - sd(var1)/sqrt(n())/3, upper=moy+ sd(var1)/sqrt(n())/3, deviation= sd(var1)) %>% 
    mutate(land.cover = " Landscape ") %>% 
    ungroup() -> total_loss
  
  data.spa %>% 
    mutate(inf.date=ifelse(is.na(inf.date)==T,-10,inf.date)) %>% 
    filter(inf.date>0) %>% 
    group_by(filename.spatial,proportionSNH,crop.order,land.cover) %>% 
    summarise( var1=mean(croploss * (1 - regulation/100))) %>%
    ungroup() %>% 
    group_by(proportionSNH, crop.order,land.cover) %>% 
    summarise(moy=mean(var1), low= moy - sd(var1)/sqrt(n()), upper=moy+ sd(var1)/sqrt(n()), deviation= sd(var1)) %>% 
    ungroup() -> landCover_loss
  
  ggplot(bind_rows(total_loss, landCover_loss)) +
    geom_point(mapping = aes(x = proportionSNH, y = moy , colour = crop.order)) +
    geom_line(mapping = aes(x = proportionSNH, y = moy , colour = crop.order)) +
    geom_errorbar(mapping = aes( x = proportionSNH, ymin = low ,  ymax = upper,  colour = crop.order),width=largeur) +
    scale_y_continuous(name = "Patch crop loss (%)") +
    xlab('') +
    scale_colour_manual(values=config_palette, limits = c("1-2-3","1-3-2","2-1-3","2-3-1","3-1-2","3-2-1","random","random9","random25"))+
    facet_wrap(.~land.cover, ncol=4)+
    labs(fill = "Crop configuration") +
    theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))+
    theme(plot.title = element_text(hjust = 0.5)) -> croploss_para
  
  ## rendement ##
  data.spa %>% 
    group_by(filename.spatial,proportionSNH,crop.order) %>% 
    summarise( var1=sum(yield)) %>% 
    ungroup() %>% 
    group_by(proportionSNH, crop.order) %>% 
    summarise(moy=mean(var1)/3, low= moy - sd(var1)/sqrt(n())/3, upper=moy+ sd(var1)/sqrt(n())/3, deviation= sd(var1)) %>% 
    mutate(land.cover = " Landscape (/3)") %>% 
    ungroup() -> total_yield
  
  data.spa %>% 
    group_by(filename.spatial,proportionSNH,crop.order,land.cover) %>% 
    summarise(var1=sum(yield)) %>% 
    ungroup() %>% 
    group_by(proportionSNH, crop.order,land.cover) %>% 
    summarise(moy=mean(var1), low= moy - sd(var1)/sqrt(n()), upper=moy+ sd(var1)/sqrt(n()), deviation= sd(var1)) %>% 
    ungroup() -> landCover_yield
  
  ggplot(bind_rows(total_yield, landCover_yield)) +
    geom_point(mapping = aes(x = proportionSNH, y = moy / 10000, colour = crop.order)) +
    geom_line(mapping = aes(x = proportionSNH, y = moy / 10000, colour = crop.order)) +
    geom_errorbar(mapping = aes( x = proportionSNH, ymin = low / 10000,  ymax = upper / 10000,  colour = crop.order),width=largeur) +
    scale_y_continuous(name = "Yield (10⁴)",breaks = seq(0.5,3.5,by=0.5)) +
    xlab('') +
    scale_colour_manual(values=config_palette, limits = c("1-2-3","1-3-2","2-1-3","2-3-1","3-1-2","3-2-1","random","random9","random25"))+
    facet_wrap(.~land.cover, ncol=4)+
    labs(colour = "Crop configuration") +
    theme(legend.title=element_text(size=14),
          legend.text=element_text(size=14))+
    theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))+
    theme(plot.title = element_text(hjust = 0.5)) -> yield_para
  
  ## graphe final ##
  x11()
  legend <-get_legend(yield_para)
  
  if(lignes == F){# pas de ligne
    gamme <- c(0,5,seq(10,max(unique(data.spa$proportionSNH)), by = 5))
    grid.arrange(arrangeGrob(visit_para+theme(legend.position="none")+scale_x_continuous(breaks = gamme),
                             delai_mean_para+theme(legend.position="none")+scale_x_continuous(breaks = gamme),
                             croploss_para+theme(legend.position="none")+scale_x_continuous(breaks = gamme),
                             yield_para+theme(legend.position="none")+scale_x_continuous(breaks = gamme),
                             ncol=1,nrow=4,bottom='Semi Natural Habitats (%)'),
                 legend,ncol=2,
                 widths=c(20, 2))
  }else{
    gamme <- c(0,5,7,seq(10,max(unique(data.spa$proportionSNH)), by = 5))
    grid.arrange(arrangeGrob(visit_para+theme(legend.position="none")+geom_vline(xintercept=7, alpha=0.3)+scale_x_continuous(breaks = gamme),
                             delai_mean_para+theme(legend.position="none")+geom_vline(xintercept=7, alpha=0.3)+scale_x_continuous(breaks = gamme),
                             croploss_para+theme(legend.position="none")+geom_vline(xintercept=7, alpha=0.3)+scale_x_continuous(breaks = gamme),
                             yield_para+theme(legend.position="none")+geom_vline(xintercept=7, alpha=0.3)+scale_x_continuous(breaks = gamme),
                             ncol=1,nrow=4),
                 legend,ncol=2,
                 widths=c(20, 3))
  }
}

graphe_SNH(data.spatial2,T)

#dev.print(device = jpeg, file = "images/var_SNH_0-30_oneof_croploss1.jpg",width=1920,height=1080)

print("finito")

##### pour un certain %  d'ESN ####
visit_land.cover_crop.order=function(objet){
  objet %>% 
    mutate(infection=if_else(is.na(inf.date) == T,0,1)) %>% #permet de calculer cb au total de patch on été infectés (les 1)
    mutate(visite=if_else(is.na(time.for.crop.loss) == T,0,1)) %>% #permet de calculer cb au total de patch on étés régulés (les 1)
    group_by(filename.spatial,crop.order) %>% 
    summarise(tx_visite= sum(visite)/sum(infection)) %>% # par fichier on a le taux de visite pour tous les patches infectés
    ungroup() %>% 
    group_by(crop.order) %>% 
    summarise(moy=mean(tx_visite), # on a le taux de visite par configuration 
              low=moy - sd(tx_visite)/sqrt(n()), # la variation est due au nombre de fichier et non au nombre de patches donc on divise pas par 3 !
              upper= ifelse(moy + sd(tx_visite)/sqrt(n())>1,1,moy + sd(tx_visite)/sqrt(n()))) %>% 
    mutate(land.cover=" landscape ") %>% 
    ungroup()-> total
  
  objet %>% 
    mutate(infection=if_else(is.na(inf.date) == T,0,1)) %>% #permet de calculer cb au total de patch on été infectés (les 1)
    mutate(visite=if_else(is.na(time.for.crop.loss) == T,0,1)) %>% #permet de calculer cb au total de patch on étés régulés (les 1)
    group_by(filename.spatial,crop.order,land.cover) %>%  
    summarise(tx_visite= sum(visite)/sum(infection)) %>%  # taux de visite pour tous les patches infectés par crop et par fichiers
    ungroup() %>% 
    group_by(crop.order,land.cover) %>% 
    summarise(moy=mean(tx_visite),# taux de visite par configuration et landCover
              low=moy - sd(tx_visite)/sqrt(n()),
              upper= if_else(moy + sd(tx_visite)/sqrt(n())>1,1,moy + sd(tx_visite)/sqrt(n()))) %>% 
    ungroup() -> landCover
  
  bind_rows(total, landCover) %>% 
    arrange(land.cover,moy) %>% 
    mutate(order=row_number()) %>% 
    ggplot() + 
    geom_col( mapping=aes(x=reorder(order,-moy), y=moy*100,fill=crop.order))+
    geom_errorbar(mapping=aes(x=reorder(order,-moy), ymin=low*100, ymax=upper*100), width=0.3)+
    scale_y_continuous(name="Visit rate (%)")+
    scale_fill_manual(values=config_palette)+
    facet_wrap(.~land.cover, scales="free_x", ncol=4)+
    xlab("")+
    theme(axis.text.x = element_blank())+
    theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))+
    theme(plot.title = element_text(hjust = 0.5))->affich
  
  return(affich)
}
delai_land.cover_crop.order_median=function(objet){
  objet %>% # tibble temporaire calcule la médiane totale pour chaque expérience puis on calcule l'écart à la moyenne
    mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% #!! 180 pour pas régulé
    filter(time.for.crop.loss<150) %>% 
    group_by(filename.spatial,crop.order, time.for.crop.loss) %>% 
    summarise(nb_delta= n()) %>% 
    mutate(nb_delta=nb_delta/sum(nb_delta)*100) %>% 
    mutate(cum_nb_delta=cumsum(nb_delta)) %>% 
    filter(cum_nb_delta>=50) %>% 
    slice(1:1) %>% 
    ungroup() %>% 
    group_by(crop.order) %>% 
    summarise(moy=mean(time.for.crop.loss),low=moy - sd(time.for.crop.loss)/sqrt(n()), upper= moy + sd(time.for.crop.loss)/sqrt(n()) ) %>% 
    mutate(land.cover=" landscape ")-> total
  
  objet %>% # calcule la médiane par landcover pour chaque expérience puis on calcule l'écart à la moyenne
    mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% #!! 180 pour pas régulé
    filter(time.for.crop.loss<150) %>% # ceux qui sont régulés après 150 ne le sont pas. 
    group_by(filename.spatial,crop.order,land.cover, time.for.crop.loss) %>% 
    summarise(nb_delta= n()) %>% 
    mutate(nb_delta=nb_delta/sum(nb_delta)*100) %>% 
    mutate(cum_nb_delta=cumsum(nb_delta)) %>% 
    filter(cum_nb_delta>=50) %>% 
    slice(1:1) %>% 
    ungroup() %>% 
    group_by(crop.order, land.cover) %>% 
    summarise(moy=mean(time.for.crop.loss),low=moy - sd(time.for.crop.loss)/sqrt(n()), upper= moy + sd(time.for.crop.loss)/sqrt(n()) ) -> landCover 
  
  
  bind_rows(total,landCover) %>% 
    arrange(land.cover,moy) %>% 
    mutate(order=row_number()) %>% 
    ungroup() %>% 
    ggplot() + 
    geom_col( mapping=aes(x=reorder(order,moy), y=moy,fill=crop.order))+
    geom_errorbar(mapping=aes(x=reorder(order,moy), ymin=low, ymax=upper), width=0.3)+
    scale_fill_manual(values=config_palette)+
    facet_wrap(.~land.cover, scales="free_x", ncol=4)+
    # titre et s axes
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("Delay (days)") + 
    xlab("") +
    theme(axis.text.x = element_blank())+
    theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))+
    theme(plot.title = element_text(hjust = 0.5))->affich
  
  return(affich)
}
delai_land.cover_crop.order_mean=function(objet){
  objet %>% # tibble temporaire calcule la médiane totale pour chaque expérience puis on calcule l'écart à la moyenne
    mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% #!! 180 pour pas régulé
    filter(time.for.crop.loss<150) %>% 
    group_by(filename.spatial,crop.order) %>% 
    summarise(var1= mean(time.for.crop.loss)) %>% 
    ungroup() %>% 
    group_by(crop.order) %>% 
    summarise(moy=mean(var1),low=moy - sd(var1)/sqrt(n()), upper= moy + sd(var1)/sqrt(n()) ) %>% 
    mutate(land.cover=" landscape ")-> total
  
  objet %>% # calcule la médiane par landcover pour chaque expérience puis on calcule l'écart à la moyenne
    mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% #!! 180 pour pas régulé
    filter(time.for.crop.loss<150) %>% # ceux qui sont régulés après 150 ne le sont pas. 
    group_by(filename.spatial,crop.order,land.cover) %>% 
    summarise(var1= mean(time.for.crop.loss)) %>% 
    ungroup() %>% 
    group_by(crop.order, land.cover) %>% 
    summarise(moy=mean(var1),low=moy - sd(var1)/sqrt(n()), upper= moy + sd(var1)/sqrt(n()) )  -> landCover 
  
  bind_rows(total,landCover) %>% 
    arrange(land.cover,moy) %>% 
    mutate(order=row_number()) %>% 
    ungroup() %>% 
    ggplot() + 
    geom_col( mapping=aes(x=reorder(order,moy), y=moy,fill=crop.order))+
    geom_errorbar(mapping=aes(x=reorder(order,moy), ymin=low, ymax=upper), width=0.3)+
    scale_fill_manual(values=config_palette)+
    facet_wrap(.~land.cover, scales="free_x", ncol=4)+
    
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("Delay (days)") + 
    xlab("") +
    theme(axis.text.x = element_blank())+
    theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))+
    theme(plot.title = element_text(hjust = 0.5))-> affich
  
  return(affich)
}
croploss_land.cover_crop.order=function(objet){
  objet %>% 
    mutate(inf.date=ifelse(is.na(inf.date)==T,-10,inf.date)) %>% 
    filter(inf.date>0) %>% 
    group_by(filename.spatial,crop.order) %>% 
    summarise( var1=mean(croploss * (1 - regulation/100))) %>% 
    ungroup() %>% 
    group_by(crop.order) %>% 
    summarise(moy=mean(var1), low= moy - sd(var1)/sqrt(n()), upper=moy+ sd(var1)/sqrt(n()), deviation= sd(var1)) %>% 
    mutate(land.cover = " landscape") %>% 
    ungroup() -> total
  
  objet %>% 
    mutate(inf.date=ifelse(is.na(inf.date)==T,-10,inf.date)) %>% 
    filter(inf.date>0) %>% 
    group_by(filename.spatial,crop.order,land.cover) %>% 
    summarise( var1=mean(croploss * (1 - regulation/100))) %>% 
    ungroup() %>% 
    group_by(crop.order,land.cover) %>% 
    summarise(moy=mean(var1), low= moy - sd(var1)/sqrt(n()), upper=moy+ sd(var1)/sqrt(n()), deviation= sd(var1)) %>% 
    ungroup() -> landCover
  
  bind_rows(total, landCover) %>% 
    arrange(land.cover,moy) %>% 
    mutate(order=row_number()) %>% 
    ggplot() + 
    geom_col( mapping=aes(x=reorder(order,moy), y=moy,fill=crop.order))+
    geom_errorbar(mapping=aes(x=reorder(order,moy), ymin=low, ymax=upper), width=0.3)+
    scale_y_continuous(name="Patch crop loss (%)")+
    scale_fill_manual(values=config_palette)+
    facet_wrap(.~land.cover, scales="free_x", ncol=4)+
    
    xlab("")+
    theme(axis.text.x = element_blank())+
    theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))+
    theme(plot.title = element_text(hjust = 0.5))->affich
  return(affich)
}
yield_land.cover_crop.order=function(objet){
  objet %>% 
    group_by(filename.spatial,crop.order) %>% 
    summarise( var1=sum(yield)) %>% 
    ungroup() %>% 
    group_by(crop.order) %>% 
    summarise(moy=mean(var1)/3, 
              low= moy - sd(var1)/sqrt(n())/3, # division par trois pour mettre à l'echelle l'intervalle
              upper=moy+ sd(var1)/sqrt(n())/3, 
              deviation= sd(var1)) %>% 
    mutate(land.cover = " landscape ") %>% 
    ungroup() -> total
  
  objet %>% 
    group_by(filename.spatial,crop.order,land.cover) %>% 
    summarise(var1= sum(yield)) %>% 
    ungroup() %>% 
    group_by(crop.order,land.cover) %>% 
    summarise(moy=mean(var1),
              low=moy - sd(var1)/sqrt(n()),
              upper=moy + sd(var1)/sqrt(n())) %>% 
    ungroup() -> landCover
  
  bind_rows(total, landCover) %>% 
    arrange(land.cover,moy) %>% 
    mutate(order=row_number()) %>% 
    ggplot() + 
    geom_col( mapping=aes(x=reorder(order,-moy), y=moy/10000,fill=crop.order))+
    geom_errorbar(mapping=aes(x=reorder(order,-moy), ymin=low/10000, ymax=upper/10000), width=0.3)+
    scale_y_continuous(name="Yield (10⁴)")+
    xlab('')+
    labs(fill='crop configuration')+
    scale_fill_manual(values=config_palette, limits= config_order)+
    facet_wrap(.~land.cover, scales="free_x", ncol=4)+
    
    theme(axis.text.x = element_blank()) +
    theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=18))+
    theme(legend.title=element_text(size=14),
          legend.text=element_text(size=14))+
    theme(plot.title = element_text(hjust = 0.5))->affich
  
  return(affich)
}
graphe_1SNH <- function(data.spa,propSNH){
  tempo <- filter(data.spa, proportionSNH == propSNH)
  visit <- visit_land.cover_crop.order(tempo)
  delai <- delai_land.cover_crop.order_mean(tempo)
  yield <- yield_land.cover_crop.order(tempo)
  croploss <- croploss_land.cover_crop.order(tempo)
  legend <- get_legend(yield)
  parametres <- paste(
    "SNH=",unique(tempo$proportionSNH),"  ",
    "I=",unique(tempo$infection.rate),"  ",
    "A=",unique(tempo$agregation),"  ",
    "NE behaviour=","OneOf  ",
    "mortality=",unique(tempo$mortality),
    "repetitions=",repet)
  x11()
  grid.arrange(arrangeGrob(visit+theme(legend.position="none"),
                           delai+theme(legend.position="none"),
                           croploss+theme(legend.position="none"),
                           yield+theme(legend.position="none"),
                           ncol=1,nrow=4),
               legend,ncol=2,
               widths=c(20, 2),
               top=parametres)
}

graphe_1SNH(data.spatial2,7)

dev.print(device = jpeg, file = "images/figzrz.jpg",width=1920,height=1080)

print('finito')

###### optimisation (commenté) ##### 
visit_opti<- function(data.spa,value){
  # fonction qui calcule la proportion d'ESN minimale pour avoir un taux de taux de visite >= value. Sa sortie étant un plot
  data.spa %>% 
    mutate(infection=if_else(is.na(inf.date) == T,0,1)) %>% #permet de calculer cb au total de patch on été infectés
    mutate(visite=if_else(is.na(time.for.crop.loss) == T,0,1)) %>% #permet de calculer cb au total de patch on étés régulés
    group_by(filename.spatial,proportionSNH,crop.order) %>% 
    summarise(tx_visite= sum(visite)/sum(infection)) %>% 
    ungroup() %>% 
    group_by(proportionSNH, crop.order) %>% 
    summarise(moy=mean(tx_visite),low=moy - sd(tx_visite)/sqrt(n())/3, upper= moy + sd(tx_visite)/sqrt(n())/3) %>% 
    mutate(land.cover=" landscape ") %>% 
    ungroup()-> total_visit_tempo
  
  data.spa %>% 
    mutate(infection=if_else(is.na(inf.date) == T,0,1)) %>% #permet de calculer cb au total de patch on été infectés
    mutate(visite=if_else(is.na(time.for.crop.loss) == T,0,1)) %>% #permet de calculer cb au total de patch on étés régulés
    group_by(filename.spatial,proportionSNH,crop.order,land.cover) %>% 
    summarise(tx_visite= sum(visite)/sum(infection)) %>% 
    ungroup() %>% 
    group_by(proportionSNH,crop.order,land.cover) %>% 
    summarise(moy=mean(tx_visite),
              low=moy - sd(tx_visite)/sqrt(n()),
              upper= if_else(moy + sd(tx_visite)/sqrt(n())>1,1,moy + sd(tx_visite)/sqrt(n()))) %>% 
    ungroup() -> landCover_visit_tempo
  
  bind_rows(total_visit_tempo, landCover_visit_tempo) %>% 
    group_by(land.cover,crop.order) %>%  
    filter(moy >= value) %>% 
    arrange(crop.order,land.cover, proportionSNH) %>% 
    slice(1:1) %>% 
    ungroup() %>% 
    arrange(land.cover, proportionSNH) %>% 
    mutate(order=row_number()) %>% 
    ggplot() +
      geom_col(mapping=aes(x=order, y=proportionSNH, fill= crop.order), width=0.5)+
      scale_fill_manual(values=config_palette)+
      ylim(0,35)+
      xlab('') +
      theme(axis.text.x = element_blank())+
      theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))+
      ylab('')+
      ggtitle(paste('Visit rate =',value * 100,"%"))+
      facet_wrap(.~land.cover, scales="free_x", ncol=4)  -> optimize_visit
  
  return(optimize_visit)
} 
delay_opti <- function(data.spa,value){
  # fonction qui calcule la proportion d'ESN minimale pour avoir un délai <= value. Sa sortie étant un plot
  data.spa %>% 
    mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% #!! mettre ça partout ! 
    group_by(filename.spatial,proportionSNH, crop.order) %>% 
    filter(time.for.crop.loss<151) %>% 
    summarise(var1= mean(time.for.crop.loss)) %>% 
    ungroup() %>% 
    group_by(proportionSNH, crop.order) %>% 
    summarise(moy=mean(var1),low=moy - sd(var1)/sqrt(n()), upper= moy + sd(var1)/sqrt(n()) ) %>% 
    mutate(land.cover=" landscape ")-> total_delay_mean_tempo
  
  data.spa %>% 
    mutate(time.for.crop.loss=if_else(is.na(time.for.crop.loss) == T,180,as.numeric(time.for.crop.loss))) %>% #!! mettre ça partout ! 
    group_by(filename.spatial,proportionSNH, crop.order, land.cover) %>% 
    filter(time.for.crop.loss<151) %>% 
    summarise(var1= mean(time.for.crop.loss)) %>% 
    ungroup() %>% 
    group_by(proportionSNH, crop.order, land.cover) %>% 
    summarise(moy=mean(var1),low=moy - sd(var1)/sqrt(n()), upper= moy + sd(var1)/sqrt(n()) )  -> landCover_delay_mean_tempo
  
  bind_rows(total_delay_mean_tempo, landCover_delay_mean_tempo) %>% 
    filter(moy <= value) %>% 
    group_by(crop.order,land.cover) %>% 
    arrange(crop.order,land.cover, proportionSNH) %>% 
    slice(1:1) %>% 
    ungroup() %>% 
    arrange(land.cover, proportionSNH) %>% 
    mutate(order=row_number()) %>%
    ggplot(mapping=aes(x=order, y=proportionSNH, fill=crop.order))+
    geom_col(width=0.5)+
    xlab('')+
    theme(axis.text.x = element_blank())+
    theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))+
    ylab('')+
    ylim(0,35)+
    scale_fill_manual(values=config_palette)+
    facet_wrap(.~land.cover, ncol=4, scales='free_x')+
    ggtitle(paste('Delay =',value)) -> optimize_delay
  
  return(optimize_delay)
}
croploss_opti <- function(data.spa, value){
  # fonction qui calcule la proportion d'ESN minimale pour avoir une croploss <= value. Sa sortie étant un plot
  data.spa %>% 
    mutate(inf.date=ifelse(is.na(inf.date)==T,-10,inf.date)) %>% 
    filter(inf.date>0) %>% 
    group_by(filename.spatial,proportionSNH,crop.order) %>% 
    summarise( var1=mean(croploss * (1 - regulation/100))) %>% 
    ungroup() %>% 
    group_by(proportionSNH,crop.order) %>% 
    summarise(moy=mean(var1), low= moy - sd(var1)/sqrt(n())/3, upper=moy+ sd(var1)/sqrt(n())/3, deviation= sd(var1)) %>% 
    mutate(land.cover = " landscape ") %>% 
    ungroup() -> total_loss_tempo
  
  data.spa %>% 
    mutate(inf.date=ifelse(is.na(inf.date)==T,-10,inf.date)) %>% 
    filter(inf.date>0) %>% 
    group_by(filename.spatial,proportionSNH,crop.order,land.cover) %>% 
    summarise( var1=mean(croploss * (1 - regulation/100))) %>% 
    ungroup() %>% 
    group_by(proportionSNH,crop.order,land.cover) %>% 
    summarise(moy=mean(var1), low= moy - sd(var1)/sqrt(n()), upper=moy+ sd(var1)/sqrt(n()), deviation= sd(var1)) %>% 
    ungroup()-> landCover_loss_tempo
  
  bind_rows(total_loss_tempo, landCover_loss_tempo) %>% 
    group_by(crop.order,land.cover) %>% 
    arrange(crop.order,land.cover, -moy) %>% 
    filter(moy <= value) %>% 
    slice(1:1) %>% 
    ungroup() %>% 
    arrange(land.cover, proportionSNH) %>% 
    mutate(order=row_number()) %>%
    ggplot(mapping=aes(x=order, y=proportionSNH, fill=crop.order))+
    geom_col(width=0.5)+
    xlab('')+
    theme(axis.text.x = element_blank())+
    theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))+
    ylab('')+
    ylim(0,35)+
    scale_fill_manual(values=config_palette)+
    labs(fill='Crop order')+
    ggtitle(paste('Patch crop loss =',value))+
    facet_wrap(.~land.cover, ncol=4, scales='free_x') -> optimize_croploss
  
  return(optimize_croploss)
}
yield_opti <- function(data.spa){
  # fonction qui détermine la proportion d'ESN pour lequel le rendement est maximal. La fonction retourne un ggplot.
  data.spa %>% 
    group_by(filename.spatial,proportionSNH,crop.order) %>% 
    summarise( var1=sum(yield)) %>% 
    ungroup() %>% 
    group_by(proportionSNH,crop.order) %>% 
    summarise(moy=mean(var1)/3, low= moy - sd(var1)/sqrt(n())/3, upper=moy+ sd(var1)/sqrt(n())/3, deviation= sd(var1)) %>% 
    mutate(land.cover = " landscape (/3)") %>% 
    ungroup() -> total_yield_tempo
  
  data.spa %>% 
    group_by(filename.spatial,proportionSNH,crop.order,land.cover) %>% 
    summarise( var1=sum(yield)) %>% 
    ungroup() %>% 
    group_by(proportionSNH,crop.order,land.cover) %>% 
    summarise(moy=mean(var1), low= moy - sd(var1)/sqrt(n()), upper=moy+ sd(var1)/sqrt(n()), deviation= sd(var1)) %>% 
    ungroup() -> landCover_yield_tempo
  
  bind_rows(total_yield_tempo, landCover_yield_tempo) %>% 
    group_by(land.cover,crop.order) %>% 
    arrange(crop.order,land.cover, -moy) %>% 
    slice(1:1) %>% 
    ungroup() %>% 
    arrange(land.cover, proportionSNH) %>% 
    mutate(order=row_number()) %>%
    ggplot() +
    geom_col(mapping=aes(x=order, y=proportionSNH, fill= crop.order), width=0.5)+
    scale_colour_manual(values=config_palette) +
    scale_fill_manual(values=config_palette)+
    xlab('') +
    theme(axis.text.x = element_blank())+
    theme(legend.title=element_text(size=14),
          legend.text=element_text(size=14))+
    theme(strip.text.x = element_text(size = 14),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))+
    ylab('')+
    ylim(0,35)+
    ggtitle('maximum yield ')+
    labs(fill = "Crop configuration") +
    facet_wrap(.~land.cover, ncol=4, scales='free_x')  -> optimize_yield
  return(optimize_yield)
}
graphe_opti <- function(visit=0.85,delay=10,croploss=25){
  # Cette fonction a pour but de créer le graphe en fonction des valeurs d'optimisations souhaitées. 
  optimize_visit <- visit_opti(data.spatial_1, visit)
  optimize_delay <- delay_opti(data.spatial_1, delay)
  optimize_croploss <- croploss_opti(data.spatial_1, croploss)
  optimize_yield <- yield_opti(data.spatial_1)
  x11()
  legend <-get_legend(optimize_yield)

  grid.arrange(arrangeGrob(optimize_visit+theme(legend.position="none"),
                           optimize_delay+theme(legend.position="none"),
                          optimize_croploss+theme(legend.position="none"),
                           optimize_yield+theme(legend.position="none"),
                          ncol=1,nrow=4),
               legend, 
               ncol=2,
               top=parametres,
               left='Semi Natural Habitat proportion (%)',
               widths=c(20, 2))
}

graphe_opti()
#dev.print(device = jpeg, file = "images/opti_croploss1.jpg",width=1920,height=1080)    

plot('finito')
## récupérer la légende
get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

## fonction croploss
croploss_function <- function(x,d1,type=1,a=25,b=3.2,c=720,mini=10){
  # fonction général du calcule de la croploss
  
  if (type == 1){#premier type de croploss,  exponentielle décroissante
    results <- a * exp(-(1/90) * (x - d1)) 
    
  }else{#deuxième type de croploss, pic à une certaine date (forme d'une densité de probabilité de loi normale)
    results <- mini + c/(b*sqrt(2*pi))*exp(-((x - a - d1)^2)/(2*b^2))
  }
  
  return(results) # un nombre
}

##ouverture du fichier spatial
fun.load.data.spatial <- function(data_path){
  
  # Random # .[0-9]{2}
  filename.spatial <- list.files(data_path, 
                                 pattern="^(spatial)")# on récupère tous les noms de fichiers qui nous intéressent
  print(filename.spatial)

  column.names.spatial <- c(# liste des noms des futures colonnes
    "year", 
    "date",
    # spatial
    "pxcor",
    "pycor",
    "land.cover",
    "sensibility.d1",
    "sensibility.d2",
    "infection.rate",
    # outputs
    "closest.snh",
    "nb.visits",
    "nb.infection",
    "nb.curation",
    "inf.date",
    "time.for.crop.loss",
    "birthYear")
  
  col.types <- cols(year = "i", # types des futures colonnes
                    date = "i",
                    pxcor = "i",
                    pycor = "i",
                    land.cover = "f",
                    sensibility.d1="i",
                    sensibility.d2="i",
                    infection.rate="n",
                    closest.snh = "i",
                    nb.visits = "i",
                    nb.infection = "i",
                    nb.curation = "i",
                    inf.date = "i",
                    time.for.crop.loss = "i",
                    birthYear = "i")
  
  
  ## load 
  tibble(filename.spatial) %>% 
    mutate(file_contents = map(filename.spatial, 
                               ~ read_delim(file.path(data_path, .), 
                                            delim = " ", 
                                            col_names = column.names.spatial, 
                                            col_types = col.types))) %>% 
    unnest() %>% 
    
    mutate(crop.proportion=str_extract_all(filename.spatial,"(C[:digit:]{1,3}-){1,5}")) %>%
    mutate(crop.proportion=str_replace_all(crop.proportion,"-",""))%>%
    mutate(crop.proportion=str_replace(crop.proportion,"C",""))%>%
    separate(crop.proportion, c("C1","C2"), sep = "C", remove = FALSE)%>%
    
    
    mutate(crop.order=str_extract(filename.spatial,"(R[:digit:]{1,2}-){1,}")) %>% 
    mutate(crop.order= str_replace_all(crop.order,"R","")) %>% 
    mutate(crop.order= str_sub(crop.order,1,5)) %>%
    
    mutate(nb.crops=str_count(filename.spatial,"-C")) %>% # compte le nombre de crops
    mutate(name2=str_replace_all(filename.spatial,"-[CISR][:digit:]{1,3}","")) %>%
    mutate(name2=str_replace_all(name2,"_[:digit:]{1,3}","")) %>%
    
    ### transformation
    # extraction des différentes composantes de filename.spatial pour en faire des var. à part entières
    separate(name2, c("spatial", "rand1", "proportionSNH","agregation","crop.distribution", "overwintering", "nb.init", "mortality", "date.to.flee", "foraging.freq", "ability", "pattern","productionNE", "rand2" ), sep = "-", remove = FALSE) %>% 
    mutate(productionNE=str_replace(productionNE,"J","")) %>% 
    
    # on retire les compossantes de filename.spatial qu'on ne veut pas conserver
    select(- spatial, - rand1, - rand2, - crop.proportion, - name2) %>% 
    
    # pour inf.rate et mortality il faut gérer le fait qu'on a enlevé les . décimales lors de la création de filename
    mutate(agregation = str_replace(agregation, "^0", "0."), 
           mortality = str_replace(mortality, "0", "0.")) %>% 
            
    # les var. sont chr donc il faut les passer en numeric
    mutate(proportionSNH = as.numeric(proportionSNH),
           agregation = as.numeric(agregation),
           overwintering = as.numeric(overwintering)/10,
           nb.init = as.numeric(nb.init),
           mortality = as.numeric(mortality),
           date.to.flee = as.numeric(date.to.flee),
           foraging.freq = as.numeric(foraging.freq),
           ability = as.numeric(ability)) -> data.spatial.plus.var
  
  return(data.spatial.plus.var)
}

## mise en forme et calcul du tibble avec regulation, crop loss et yield
ouverture_calcul_patch=function(chemin,fonction){
  data.spatial <- fun.load.data.spatial(chemin) 
  print(data.spatial)
  data.spatial %>% 
    mutate(crop.order= if_else(crop.distribution=="randomly","random", crop.order)) %>% 
    mutate(crop.order= if_else(crop.distribution=="randomly9","random9", crop.order)) %>% 
    mutate(crop.order= if_else(crop.distribution=="randomly22","random25", crop.order)) %>% 
    mutate(crop.order= if_else(crop.distribution=="randomly25","random25", crop.order)) %>% 
    
    mutate(type_croploss = 1) %>% 
    mutate(croploss = if_else(inf.date==180,0,croploss_function(inf.date,sensibility.d1,1))) %>% 
    
    mutate(regulation= if_else(time.for.crop.loss < 51, 90 * (exp((- 1 / 30) * ( time.for.crop.loss))),12.48 * exp((50 - time.for.crop.loss)/40)+ 4.5)) %>% 
    mutate(regulation = if_else(is.na(regulation)==T,0, regulation)) %>% 
    mutate(yield = 100 * (1 - croploss /100 * (1 - regulation / 100))) %>% 
    select(-date, -year, - pxcor, -pycor, - crop.distribution,-birthYear, - date.to.flee, - nb.init) %>% 
    mutate( land.cover = ifelse(land.cover==1,'crop 1',ifelse(land.cover==2,'crop 2','crop 3')))->  result
  return(result)}



### à virer du Rproj
##fonction qui sors une palette de couleur pour que chaque groupe soit d'une couleur différente. 
hsd.color=function(hsd.object){
  bar_col=brewer.pal(length(unique(hsd.object$groups$groups)),'Set3')
  bar_colf=bar_col[1]
  tempo=1
  for (i in 2:length(hsd.object$groups$groups)) {
    if (hsd.object$groups$groups[i-1] == hsd.object$groups$groups[i])
    {
      bar_colf = c(bar_colf, list(bar_col[tempo]))
    }else {
      tempo= tempo +1
      bar_colf = c(bar_colf, list(bar_col[tempo]))
    }}
  return(paste(bar_colf))
  
  bar.group(distribution.hsd$groups,ylim=c(0, 500000 ), col=hsd.color(distribution.hsd), ylab= "yield")
}

## fonction pour les groupes sur les graphes
grouper <- function(objet_trie){
  groupes=rep('', times=length(objet_trie$low))
  tempo <- 0
  
  for(indice in 1: length(groupes)){
    if (groupes[indice]==''){
      tempo <-  tempo + 1 
      letter <- letters[tempo]
      
      for (i in 1:length(objet_trie$crop.order)){
        if (objet_trie$low[indice] <= objet_trie$upper[i] ){
          if(objet_trie$upper[indice] >= objet_trie$low[i]){
            groupes[i] = paste(groupes[i], letter)
          }
        }
      }
    }
  }
  return(bind_cols(objet_trie,tibble(groups=groupes)))
}

##partie déplacement ennemis
##ouverture pas finie 
fun.load.data.enemies <- function(data_path){
  
  # Random # .[0-9]{2}
  filename.enemies <- list.files(data_path, 
                                 pattern="^(adult_predators_movements)")# on récupère tous les noms de fichiers qui nous intéressent
  print(filename.enemies)
  
  column.names.enemies <- c(# liste des noms des futures colonnes
    "ID",
    "birthYear",
    "birthDate",
    "year", 
    "date",
    "landCover")
  
  col.types <- cols(year = "i", # types des futures colonnes
                    date = "i",
                    landCover = "f",
                    ID =  "c",
                    birthYear = "i",
                    birthDate = "i")
  
  ## load 
  tibble(filename.enemies) %>% 
    mutate(file_contents = map(filename.enemies, 
                               ~ read_delim(file.path(data_path, .), 
                                            delim = " ", 
                                            col_names = column.names.enemies, 
                                            col_types = col.types))) %>% 
    unnest() %>% 
    
    mutate(crop.proportion=str_extract_all(filename.enemies,"(C[:digit:]{1,3}-){1,}")) %>% 
    mutate(crop.proportion=str_replace_all(crop.proportion,"-",""))%>%
    mutate(crop.proportion=str_replace(crop.proportion,"C",""))%>%
    separate(crop.proportion, c("C1","C2"), sep = "C", remove = FALSE, convert=TRUE)%>%
    
    mutate(crop.order=str_extract(filename.enemies,"(R[:digit:]{1,2}-){1,}")) %>% 
    mutate(crop.order= str_replace_all(crop.order,"R","")) %>% 
    mutate(crop.order= str_sub(crop.order,1,5)) %>% 
    
    mutate(name2=str_replace_all(filename.enemies,"-[CISR][:digit:]{1,3}","")) %>%
    mutate(name2=str_replace_all(name2,"_[:digit:]{1,3}","")) %>%
    ### transformation
    # extraction des différentes composantes de filename.enemies pour en faire des var. à part entières
    separate(name2, c("adult_predators_movements", "rand1", "proportionSNH","agregation","crop.distribution", "overwintering", "nb.init", "mortality", "date.to.flee", "foraging.freq", "ability", "pattern","productionNE", "rand2" ), sep = "-", remove = FALSE) %>% 
    mutate(productionNE=str_replace(productionNE,"J","")) %>% 
    
    # on retire les compossantes de filename.enemies qu'on ne veut pas conserver
    select(- adult_predators_movements, - rand1, - rand2, - crop.proportion, - name2) %>% 
    
    # pour inf.rate et mortality il faut gérer le fait qu'on a enlevé les . décimales lors de la création de filename
    mutate(agregation = str_replace(agregation, "^0", "0."), 
           mortality = str_replace(mortality, "0", "0.")) %>% 
    
    # les var. sont chr donc il faut les passer en numeric
    mutate(proportionSNH = as.numeric(proportionSNH),
           agregation = as.numeric(agregation),
           overwintering = as.numeric(overwintering)/10,
           nb.init = as.numeric(nb.init),
           mortality = as.numeric(mortality),
           date.to.flee = as.numeric(date.to.flee),
           foraging.freq = as.numeric(foraging.freq),
           ability = as.numeric(ability)) -> data.enemies.plus.var
  
  return(data.enemies.plus.var)
}

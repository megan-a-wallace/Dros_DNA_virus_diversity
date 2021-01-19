#################################################################
## INLA Models for patterns of Drosophila DNA virus prevalence ##
#################################################################

##Written by Megan A. Wallace & Gregory F. Albery
##08/2020-01/2021

#Script to use R-INLA, and functions from package ggregplot to model variation in DNA virus prevalence across space and time
##Creates maps of spatial fields, effect size plots, and calculates range of spatial autocorrelation

##Installing 
install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
install.packages("remotes")
remotes::install_github("gfalbery/ggregplot")

##Loading packages
library(tidyverse); library(magrittr); library(ggregplot); library(INLA); library(colorspace);library(ggplot2); library(RColorBrewer); library(grid); library(inlabru);library("rnaturalearth"); library("rnaturalearthdata"); library(cowplot); library(patchwork); library(rockchalk);library(colorRamps);library(lemon)

# Establishing Themes and Palettes ####

Darren_palette<-colorRampPalette(c("firebrick3","darkgoldenrod1","steelblue","white"))(10)
Darren_scale<-colorRampPalette(c("white","steelblue","darkgoldenrod1","firebrick3"))

##pallettes and themes from ggregplot
ParasitePalettes<-c("PuRd","PuBu","BuGn","Purples","Oranges")
ParasiteColours<-c("#DD1c77","#2B8CBE","#2CA25F",brewer.pal(5,"Purples")[4],brewer.pal(5,"Oranges")[4])

AlberPalettes <- c("YlGnBu","Reds","BuPu", "PiYG")
AlberColours <- sapply(AlberPalettes, function(a) RColorBrewer::brewer.pal(5, a)[4])
AlberColours[length(AlberColours)+1:2] <- RColorBrewer::brewer.pal(11, AlberPalettes[[4]])[c(2,10)]

AlberTheme <- theme_bw() +
  theme(axis.title.x = element_text(vjust = -0.35,
                                    size = 12,
                                    colour = "black"),
        axis.title.y = element_text(vjust = 1.2,
                                    size = 12,
                                    colour = "black"),
        plot.title = element_text(size=16),
        strip.background = element_rect(fill = "white", colour = "dark grey"))

theme_set(AlberTheme)

##input data
dat <- read.delim("Kallithea_diversity_mapping/Supplementary_File_S4_justGLMM.txt")

VirusVars <- c("Kallithea", "Esparto", "Linvill_Road", "Viltain", "Vesanto") %>% paste0("_virus") %>% c("Relative_Galbut")

dat %>%
  select(all_of(VirusVars)) %>%
  mutate_at(VirusVars[1:5], ~as.numeric(.x>0.01)) %>%
  mutate_at(VirusVars[6], ~as.numeric(.x>0.001)) %>%
  rename_all(~paste0("Binary", .x) %>% str_remove("_virus$") %>% str_remove("Relative_")) %>%
  bind_cols(dat, .) -> dat

Covar <- c("Year", "Perc_Dsim", "Wolbachia", "Season")

dat$Season <- factor(dat$Season,levels = c("S","F"))
dat$Year <- factor(dat$Year,levels = c("Y2014","Y2015","2016"))

dat %>% select(contains("Binary")) %>%
  names %>% setdiff(c("BinaryEsparto", "BinaryViltain")) ->
  Resps
#removing esparto and viltain because there aren't enough samples for a reasonable model

##Evaluating the addition of a spatial mesh, and a separate, uncorrelated spatial mesh per season to Binomial INLA models, with all fixed effects and the location random effect, and its impact on DIC

IMList <- list()
for(r in seq_along(Resps)){
  
  print(Resps[r])
  
  IM1 <- INLAModelAdd(Resps[r], dat,
                      Explanatory = Covar,
                      Add = NULL,
                      Random = "Loc", RandomModel = "iid",
                      Family = "binomial",
                      AddSpatial = T, Coordinates = c("Lon", "Lat"),
                      Groups = T,GroupVar = "Season", GroupModel = "Rep")
  
  IMList[[Resps[r]]] <- IM1
  
}

##Support for the addition of a spatial random effect to the models in the case of Kallithea, Linvill Road and for Galbut, a spatiotemporal effect with a separate mesh for each season
Resps %>% map(~IMList[[.x]]$FinalModel %>% list(IMList[[.x]]$Spatial$Model,IMList[[.x]]$Spatial$SpatiotemporalModel) %>% INLADICFig() + ggtitle(.x)) %>%
  arrange_ggplot2(ncol = 2)

##And looking at the model summaries for the final models, including a spatial component for Kallithea and Linvill Road, and a spatiotemporal component for Galbut
IMList[[Resps[1]]]$Spatial$Model %>% summary()

IMList[[Resps[2]]]$Spatial$Model %>% summary()

IMList[[Resps[3]]]$FinalModel %>% summary()

IMList[[Resps[4]]]$Spatial$SpatiotemporalModel %>% summary()

####
#Looking at the proportion of variance explained by the random effects, both location and the spatial effect
IMList$BinaryKallithea$Spatial$Model %>% INLARep(Family = "binomial",SPDEModel = IMList$BinaryKallithea$Spatial$SPDE)
IMList$BinaryLinvill_Road$Spatial$Model %>% INLARep(Family = "binomial",SPDEModel = IMList$BinaryLinvill_Road$Spatial$SPDE)
IMList$BinaryLinvill_Road$FinalModel %>% INLARep(Family = "binomial")
IMList$BinaryGalbut$Spatial$SpatiotemporalModel %>% INLARep(Family = "binomial",SPDEModel = IMList$BinaryGalbut$Spatial$SPDE)

#And looking at the range over which spatial autocorrelation fades
cols<-c("#E69F00", "#56B4E9", "#0072B2")

MaxRange = 50
Priors <- MaxRange/2
PriorProbabilities <- 0.5
ModelList=list(IMList$BinaryKallithea$Spatial$Model,IMList$BinaryLinvill_Road$Spatial$Model,IMList$BinaryGalbut$Spatial$SpatiotemporalModel)
MeshList=list(IMList$BinaryKallithea$Spatial$Mesh,IMList$BinaryLinvill_Road$Spatial$Mesh,IMList$BinaryGalbut$Spatial$Mesh)
Priors <- rep(Priors, length(ModelList))
PriorProbabilities <- rep(PriorProbabilities, length(ModelList))
WNames = c("w","w","wTemporal") 
Resolution = 100
ModelNames = c("Kallithea virus","Linvill Road virus", "Galbut EVE")

SpFi.w = 1:length(ModelList) %>% lapply(function(j){
  inla.spde2.result(inla = ModelList[[j]],
                    name = WNames[[j]],
                    spde = inla.spde2.pcmatern(mesh = MeshList[[j]],
                                               prior.range = c(Priors[[j]], PriorProbabilities[[j]]),
                                               prior.sigma = c(.5, .5)),
                    do.transfer = TRUE)
})

Kappa <- lapply(SpFi.w,  function(j)
  inla.emarginal(function(x) x,
                 j$marginals.kappa[[1]] ))

d.vec <- seq(0, MaxRange, length = Resolution)

Cor<-lapply(Kappa,function(f){
  Cor.M <- as.numeric((f * d.vec) * besselK(f * d.vec, 1))
  Cor.M[1] <- 1
  return(data.frame(d.vec,Cor.M))
})

Cor <- dplyr::bind_rows(Cor)

Cor$Model <- as.factor(rep(1:length(Kappa), each = Resolution))
levels(Cor$Model) <- ModelNames

ReturnList <- list()

ReturnList$Figure <- ggplot(Cor, aes(d.vec, Cor.M, colour = Model, lty = Model)) +
  geom_line(size = 1) + coord_fixed(ratio = MaxRange) +
  labs(colour ="Model",x = "Distance (degrees)", y = "Correlation")
ReturnList$Data <- Cor
ReturnList$Kappa <- unlist(Kappa)

range_data_Kallithea<-ReturnList$Data[ReturnList$Data$Model=="Kallithea virus",]
range_data_Linvill_Road<-ReturnList$Data[ReturnList$Data$Model=="Linvill Road virus",]
range_data_Galbut<-ReturnList$Data[ReturnList$Data$Model=="Galbut EVE",]

#Halving range
range_data_Kallithea$d.vec[which(abs(range_data_Kallithea$Cor.M-0.5)==min(abs(range_data_Kallithea$Cor.M-0.5)))]->halving_distance_Kallithea
range_data_Linvill_Road$d.vec[which(abs(range_data_Linvill_Road$Cor.M-0.5)==min(abs(range_data_Linvill_Road$Cor.M-0.5)))]->halving_distance_Linvill_Road
range_data_Galbut$d.vec[which(abs(range_data_Galbut$Cor.M-0.5)==min(abs(range_data_Galbut$Cor.M-0.5)))]->halving_distance_Galbut

#Kappa (inverse range)
kappa_Kallithea<-ReturnList$Kappa[1]
kappa_Linvill_Road<-ReturnList$Kappa[2]
kappa_Galbut<-ReturnList$Kappa[3]

###Figure showing range over which the spatial autocorrelation fades in each of the models
ReturnList$Figure + scale_color_manual(name = "Model", values = cols, labels = c("Kallithea virus","Linvill Road virus","Galbut EVE"),breaks=c("Kallithea virus", "Linvill Road virus","Galbut EVE")) + theme(legend.position = c(.9, .9),legend.justification = c("right", "top"))

##Looking at how adding the static spatial field affects the fixed effects effect sizes
##Plotting the effect sizes of the models with and without the spatial components

###Runnning the internal bits of Efxplot to make what I want
graphlist_Kallithea<-list()
ModelList_Kallithea <-list(IMList$BinaryKallithea$FinalModel,IMList$BinaryKallithea$Spatial$Model,IMList$BinaryKallithea$Spatial$SpatiotemporalModel)

for(i in 1:length(ModelList_Kallithea)){
  
  model_Kallithea<-ModelList_Kallithea[[i]]
  
  graph_Kallithea<-as.data.frame(summary(model_Kallithea)$fixed)
  colnames(graph_Kallithea)[which(colnames(graph_Kallithea)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
  colnames(graph_Kallithea)[which(colnames(graph_Kallithea)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
  colnames(graph_Kallithea)[which(colnames(graph_Kallithea)%in%c("mean"))]<-c("Estimate")
  
  graph_Kallithea$Model<-i
  graph_Kallithea$Factor<-rownames(graph_Kallithea)
  
  graphlist_Kallithea[[i]]<-graph_Kallithea
  
}

graph_Kallithea <- bind_rows(graphlist_Kallithea)
graph_Kallithea$Sig <- with(graph_Kallithea, ifelse(Lower*Upper>0, "*", ""))
graph_Kallithea$Model <- as.factor(graph_Kallithea$Model)

levels(graph_Kallithea$Model)<-c("Wolb + Year + Season + %Dsim + f(Loc)","Wolb + Year + Season + %Dsim + f(Loc) + w","Wolb + Year + Season + %Dsim + f(Loc) + wTemporal")
position <- ifelse(length(unique(graph_Kallithea$Model))  ==  1, "none", "right")

graph_Kallithea$Factor <- factor(graph_Kallithea$Factor, levels = c("(Intercept)","Intercept","Perc_Dsim","SeasonF","Year2016","YearY2015","Wolbachia"))
combineLevels(graph_Kallithea$Factor,levs = c("Intercept", "(Intercept)"), newLabel = c("Intercept"))->graph_Kallithea$Factor
levels(graph_Kallithea$Factor) <- c("% Dsim","Season:Late","Year:2016","Year:2015","Wolbachia","Intercept")

#to remove intercept
#VarNames <- VarNames[!str_detect(VarNames, "ntercept")]
#graph <- graph %>% filter(Factor %in% VarNames)

min_Kallithea<-min(graph_Kallithea$Lower,na.rm = T)
max_Kallithea<-max(graph_Kallithea$Upper,na.rm = T)
graph_Kallithea$starloc <- 17.7

cols<-c("#E69F00", "#56B4E9", "#0072B2")
Kallithea_efxplot <- ggplot(as.data.frame(graph_Kallithea),
                            aes(x = as.factor(Factor),
                                y = Estimate,
                                group = Model,
                                colour = Model))+
  scale_y_continuous(breaks = seq(-8, 18, by = 4),limits = c(-8,18))+
  geom_point(position = position_dodge(w = 0.5), size = 2.1) +
  geom_errorbar(position = position_dodge(w = 0.5),
                aes(ymin = Lower, ymax = Upper), size = 0.7,
                width = 0.3) +
  scale_color_manual(values=cols)+
  geom_hline(aes(yintercept = 0),lty = 2) + labs(x = NULL) + coord_flip() +
  theme(legend.position = "none") + #=position if you want the legend 
  geom_text(aes(label = Sig, y = starloc),
            position = position_dodge(w = 0.5),
            show.legend = F,size=6) +
  ggtitle(label = "Kallithea virus")

##Linvill Road
graphlist_Linvill_Road<-list()
ModelList_Linvill_Road <-list(IMList$BinaryLinvill_Road$FinalModel,IMList$BinaryLinvill_Road$Spatial$Model,IMList$BinaryLinvill_Road$Spatial$SpatiotemporalModel)

for(i in 1:length(ModelList_Linvill_Road)){
  
  model_Linvill_Road<-ModelList_Linvill_Road[[i]]
  
  graph_Linvill_Road<-as.data.frame(summary(model_Linvill_Road)$fixed)
  colnames(graph_Linvill_Road)[which(colnames(graph_Linvill_Road)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
  colnames(graph_Linvill_Road)[which(colnames(graph_Linvill_Road)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
  colnames(graph_Linvill_Road)[which(colnames(graph_Linvill_Road)%in%c("mean"))]<-c("Estimate")
  
  graph_Linvill_Road$Model<-i
  graph_Linvill_Road$Factor<-rownames(graph_Linvill_Road)
  
  graphlist_Linvill_Road[[i]]<-graph_Linvill_Road
  
}

graph_Linvill_Road <- bind_rows(graphlist_Linvill_Road)
graph_Linvill_Road$Sig <- with(graph_Linvill_Road, ifelse(Lower*Upper>0, "*", ""))
graph_Linvill_Road$Model <- as.factor(graph_Linvill_Road$Model)

levels(graph_Linvill_Road$Model)<-c("Wolb + Year + Season + %Dsim + f(Loc)","Wolb + Year + Season + %Dsim + f(Loc) + w","Wolb + Year + Season + %Dsim + f(Loc) + wTemporal")
position <- ifelse(length(unique(graph_Linvill_Road$Model))  ==  1, "none", "right")

graph_Linvill_Road$Factor <- factor(graph_Linvill_Road$Factor, levels = c("(Intercept)","Intercept","Perc_Dsim","SeasonF","Year2016","YearY2015","Wolbachia"))
combineLevels(graph_Linvill_Road$Factor,levs = c("Intercept", "(Intercept)"), newLabel = c("Intercept"))->graph_Linvill_Road$Factor
levels(graph_Linvill_Road$Factor) <- c("% Dsim","Season:Late","Year:2016","Year:2015","Wolbachia","Intercept")

#to remove intercept
#VarNames <- VarNames[!str_detect(VarNames, "ntercept")]
#graph <- graph %>% filter(Factor %in% VarNames)

min_Linvill_Road<-min(graph_Linvill_Road$Lower,na.rm = T)
max_Linvill_Road<-max(graph_Linvill_Road$Upper,na.rm = T)
graph_Linvill_Road$starloc <- 17.7

Linvill_Road_efxplot <- ggplot(as.data.frame(graph_Linvill_Road),
                               aes(x = as.factor(Factor),
                                   y = Estimate,
                                   group = Model,
                                   colour = Model)) +
  scale_y_continuous(breaks = seq(-8, 18, by = 4),limits = c(-8,18))+
  scale_x_discrete(labels = rep("           ", 6), breaks = 1:6)+
  geom_point(position = position_dodge(w = 0.5), size = 2.1) +
  geom_errorbar(position = position_dodge(w = 0.5),
                aes(ymin = Lower, ymax = Upper), size = 0.7,
                width = 0.3) +
  scale_color_manual(values=cols)+
  geom_hline(aes(yintercept = 0),lty = 2) + labs(x = NULL) + coord_flip() +
  theme(legend.position = "none") + #=position if you want the legend 
  geom_text(aes(label = Sig, y = starloc),
            position = position_dodge(w = 0.5),
            show.legend = F,size=6) +
  ggtitle(label = "Linvill Road virus")

##Vesanto
graphlist_Vesanto<-list()
ModelList_Vesanto <-list(IMList$BinaryVesanto$FinalModel,IMList$BinaryVesanto$Spatial$Model,IMList$BinaryVesanto$Spatial$SpatiotemporalModel)

for(i in 1:length(ModelList_Vesanto)){
  
  model_Vesanto<-ModelList_Vesanto[[i]]
  
  graph_Vesanto<-as.data.frame(summary(model_Vesanto)$fixed)
  colnames(graph_Vesanto)[which(colnames(graph_Vesanto)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
  colnames(graph_Vesanto)[which(colnames(graph_Vesanto)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
  colnames(graph_Vesanto)[which(colnames(graph_Vesanto)%in%c("mean"))]<-c("Estimate")
  
  graph_Vesanto$Model<-i
  graph_Vesanto$Factor<-rownames(graph_Vesanto)
  
  graphlist_Vesanto[[i]]<-graph_Vesanto
  
}

graph_Vesanto <- bind_rows(graphlist_Vesanto)
graph_Vesanto$Sig <- with(graph_Vesanto, ifelse(Lower*Upper>0, "*", ""))
graph_Vesanto$Model <- as.factor(graph_Vesanto$Model)

levels(graph_Vesanto$Model)<-c("Wolb + Year + Season + %Dsim + f(Loc)","Wolb + Year + Season + %Dsim + f(Loc) + w","Wolb + Year + Season + %Dsim + f(Loc) + wTemporal")
position <- ifelse(length(unique(graph_Vesanto$Model))  ==  1, "none", "right")

graph_Vesanto$Factor <- factor(graph_Vesanto$Factor, levels = c("(Intercept)","Intercept","Perc_Dsim","SeasonF","Year2016","YearY2015","Wolbachia"))
combineLevels(graph_Vesanto$Factor,levs = c("Intercept", "(Intercept)"), newLabel = c("Intercept"))->graph_Vesanto$Factor
levels(graph_Vesanto$Factor) <- c("% Dsim","Season:Late","Year:2016","Year:2015","Wolbachia","Intercept")

#to remove intercept
#VarNames <- VarNames[!str_detect(VarNames, "ntercept")]
#graph <- graph %>% filter(Factor %in% VarNames)

min_Vesanto<-min(graph_Vesanto$Lower,na.rm = T)
max_Vesanto<-max(graph_Vesanto$Upper,na.rm = T)
graph_Vesanto$starloc <- 17.7

Vesanto_efxplot <- ggplot(as.data.frame(graph_Vesanto),
                          aes(x = as.factor(Factor),
                              y = Estimate,
                              group = Model,
                              colour = Model)) +
  scale_y_continuous(breaks = seq(-8, 18, by = 4),limits = c(-8,18))+
  scale_x_discrete(labels = rep("           ", 6), breaks = 1:6)+
  geom_point(position = position_dodge(w = 0.5), size = 2.1) +
  geom_errorbar(position = position_dodge(w = 0.5),
                aes(ymin = Lower, ymax = Upper), size = 0.7,
                width = 0.3) +
  scale_color_manual(values=cols)+
  geom_hline(aes(yintercept = 0),lty = 2) + labs(x = NULL) + coord_flip() +
  theme(legend.position = "none") + #=position if you want the legend 
  geom_text(aes(label = Sig, y = starloc),
            position = position_dodge(w = 0.5),
            show.legend = F,size = 6) +
  ggtitle(label = "Vesanto virus")

##Galbut EVE
graphlist_Galbut<-list()
ModelList_Galbut <-list(IMList$BinaryGalbut$FinalModel,IMList$BinaryGalbut$Spatial$Model,IMList$BinaryGalbut$Spatial$SpatiotemporalModel)

for(i in 1:length(ModelList_Galbut)){
  
  model_Galbut<-ModelList_Galbut[[i]]
  
  graph_Galbut<-as.data.frame(summary(model_Galbut)$fixed)
  colnames(graph_Galbut)[which(colnames(graph_Galbut)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
  colnames(graph_Galbut)[which(colnames(graph_Galbut)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
  colnames(graph_Galbut)[which(colnames(graph_Galbut)%in%c("mean"))]<-c("Estimate")
  
  graph_Galbut$Model<-i
  graph_Galbut$Factor<-rownames(graph_Galbut)
  
  graphlist_Galbut[[i]]<-graph_Galbut
  
}

graph_Galbut <- bind_rows(graphlist_Galbut)
graph_Galbut$Sig <- with(graph_Galbut, ifelse(Lower*Upper>0, "*", ""))
graph_Galbut$Model <- as.factor(graph_Galbut$Model)

levels(graph_Galbut$Model)<-c("Non-spatial","Spatial","Spatiotemporal")
position <- ifelse(length(unique(graph_Galbut$Model))  ==  1, "none", "right")

graph_Galbut$Factor <- factor(graph_Galbut$Factor, levels = c("(Intercept)","Intercept","Perc_Dsim","SeasonF","Year2016","YearY2015","Wolbachia"))
combineLevels(graph_Galbut$Factor,levs = c("Intercept", "(Intercept)"), newLabel = c("Intercept"))->graph_Galbut$Factor
levels(graph_Galbut$Factor) <- c("% Dsim","Season:Late","Year:2016","Year:2015","Wolbachia","Intercept")

#to remove intercept
#VarNames <- VarNames[!str_detect(VarNames, "ntercept")]
#graph <- graph %>% filter(Factor %in% VarNames)

min_Galbut<-min(graph_Galbut$Lower,na.rm = T)
max_Galbut<-max(graph_Galbut$Upper,na.rm = T)
graph_Galbut$starloc <- 17.7

Galbut_efxplot <- ggplot(as.data.frame(graph_Galbut),
                         aes(x = as.factor(Factor),
                             y = Estimate,
                             group = Model,
                             colour = Model)) +
  scale_y_continuous(breaks = seq(-8, 18, by = 4),limits = c(-8,18))+
  scale_x_discrete(labels = rep("           ", 6), breaks = 1:6)+
  geom_point(position = position_dodge(w = 0.5), size = 2.1) +
  geom_errorbar(position = position_dodge(w = 0.5),
                aes(ymin = Lower, ymax = Upper), size = 0.7,
                width = 0.3) +
  scale_color_manual(values=cols)+#using my usual cbp pallette
  geom_hline(aes(yintercept = 0),lty = 2) + labs(x = NULL) + coord_flip() +
  theme(legend.position = position) + #=position if you want the legend 
  geom_text(aes(label = Sig, y = starloc),
            position = position_dodge(w = 0.5),
            show.legend = F,size = 6) +
  ggtitle(label = "Galbut EVE")

ggdraw() +
  draw_plot(Kallithea_efxplot, x = 0, y = 0, width = 0.28, height = 1) +
  draw_plot(Linvill_Road_efxplot, x = 0.27, y = 0, width = 0.2, height = 1) +
  draw_plot(Vesanto_efxplot, x = 0.46, y = 0, width = 0.2, height = 1) +
  draw_plot(Galbut_efxplot, x = 0.65, y = 0, width = 0.35, height = 1)

#########
##Plotting the supported spatial and spatiotemporal fields
#########

##Altering ggfield a bit to get fly level prevalence on the legend, all scales the same and Latitude and Longitude as axis labels, just running the code from the function separately

# ggplot for INLA fields ####
###Kallithea virus
Projection_Kallithea <- inla.mesh.projector(IMList$BinaryKallithea$Spatial$Mesh,dims = c(300, 300))
Full.Projection_Kallithea <- expand.grid(x = Projection_Kallithea$x, y = Projection_Kallithea$y)
Dim1 <- nrow(Full.Projection_Kallithea)

WName_Kallithea <- IMList$BinaryKallithea$Spatial$Model$summary.hyperpar %>% rownames
WName_Kallithea[str_detect(WName_Kallithea, "Range for")] %>%
  str_split(" ") %>% map_chr(last) ->
  WName_Kallithea

Full.Projection_Kallithea$value <- c(inla.mesh.project(Projection_Kallithea,IMList$BinaryKallithea$Spatial$Model$summary.random$w$mean))

Full.Projection_Kallithea$Fill <- Full.Projection_Kallithea$value

###altering so that the underlying plot is the fly level prevalence
Full.Projection_Kallithea %<>% mutate_at("Fill", ~logistic(as.numeric(as.character(.x))))
Full.Projection_Kallithea$Fill->q_Kallithea #this is the likelihood of an infected pool (site-level prevalence)
##Now calculating the fly level prevalence from the rearranged max likelihood function, assuming 40 flies in all pools 
Full.Projection_Kallithea$Fill<-(1-(exp((log(1-q_Kallithea))/40)))*100

Full.Projection_Kallithea <- na.omit(Full.Projection_Kallithea)

FieldPlot_Kallithea <- ggplot(Full.Projection_Kallithea,aes(x, y))

FieldPlot_Kallithea <- FieldPlot_Kallithea +
  # geom_tile(colour = "black") +
  geom_tile(colour = NA, aes(fill = Full.Projection_Kallithea$Fill)) +
  coord_fixed() + labs(fill = "Prev. (%)") +
  labs(x = "Longitude", y = "Latitude")

###Linvill Road virus 
Projection_Linvill <- inla.mesh.projector(IMList$BinaryLinvill_Road$Spatial$Mesh,dims = c(300, 300))
Full.Projection_Linvill <- expand.grid(x = Projection_Linvill$x, y = Projection_Linvill$y)
Dim1 <- nrow(Full.Projection_Linvill)

WName_Linvill <- IMList$BinaryLinvill_Road$Spatial$Model$summary.hyperpar %>% rownames
WName_Linvill[str_detect(WName_Linvill, "Range for")] %>%
  str_split(" ") %>% map_chr(last) ->
  WName_Linvill

Full.Projection_Linvill$value <- c(inla.mesh.project(Projection_Linvill,IMList$BinaryLinvill_Road$Spatial$Model$summary.random$w$mean))

Full.Projection_Linvill$Fill <- Full.Projection_Linvill$value

Full.Projection_Linvill %<>% mutate_at("Fill", ~logistic(as.numeric(as.character(.x))))
Full.Projection_Linvill$Fill->q_Linvill #this is the likelihood of an infected pool (site-level prevalence)
##Now calculating the fly level prevalence from the rearranged max likelihood function, assuming 40 flies in all pools this time, as you've got a full projection rather than specific locations with pool sizes attached...and we've used a simple binomial model 
Full.Projection_Linvill$Fill<-(1-(exp((log(1-q_Linvill))/40)))*100 

Full.Projection_Linvill <- na.omit(Full.Projection_Linvill)

FieldPlot_Linvill <- ggplot(Full.Projection_Linvill,aes(x, y))

FieldPlot_Linvill <- FieldPlot_Linvill +
  geom_tile(colour = NA, aes(fill = Full.Projection_Linvill$Fill)) +
  coord_fixed() + labs(fill = "Prev. (%)") +
  labs(x = "Longitude", y = "Latitude")

###And now the Galbut EVE with Early and Late season 
Projection_Galbut <- inla.mesh.projector(IMList$BinaryGalbut$Spatial$Mesh,dims = c(300, 300))
Full.Projection_Galbut <- expand.grid(x = Projection_Galbut$x, y = Projection_Galbut$y)
Dim1 <- nrow(Full.Projection_Galbut)

WName_Galbut<- IMList$BinaryGalbut$Spatial$SpatiotemporalModel$summary.hyperpar %>% rownames
WName_Galbut[str_detect(WName_Galbut, "Range for")] %>%
  str_split(" ") %>% map_chr(last) ->
  WName_Galbut

Full.Projection_Galbut[,paste0("Group",1:2)] <-apply(matrix(IMList$BinaryGalbut$Spatial$SpatiotemporalModel$summary.random[[WName_Galbut]]$mean, ncol = 2), 2,function(x) c(inla.mesh.project(Projection_Galbut, x)))
Full.Projection_Galbut <-
  reshape2::melt(Full.Projection_Galbut,
                 id.vars = c(names(Full.Projection_Galbut)[-which(names(Full.Projection_Galbut)%in%paste0("Group",
                                                                                                          1:2))]))
Full.Projection_Galbut$Fill <- Full.Projection_Galbut$value

Full.Projection_Galbut %<>% mutate_at("Fill", ~logistic(as.numeric(as.character(.x))))
Full.Projection_Galbut$Fill->q_Galbut #this is the likelihood of an infected pool (site-level prevalence)
##Now calculating the fly level prevalence from the rearranged max likelihood function, assuming 40 flies in all pools this time
Full.Projection_Galbut$Fill<-(1-(exp((log(1-q_Galbut))/40)))*100
Full.Projection_Galbut$Group <- factor(rep(1:2, each = Dim1),levels = c(1,2),labels = c("Galbut EVE - Early Season","Galbut EVE - Late Season"))
Full.Projection_Galbut <- na.omit(Full.Projection_Galbut)

FieldPlot_Galbut <- ggplot(Full.Projection_Galbut,aes(x, y))

FieldPlot_Galbut <- FieldPlot_Galbut +
  geom_tile(colour = NA, aes(fill = Full.Projection_Galbut$Fill)) +
  coord_fixed() + 
  labs(fill = "Prev. (%)") +
  labs(x = "Longitude", y = "Latitude")

FieldPlot_Galbut <- FieldPlot_Galbut + facet_rep_wrap( ~ Group, repeat.tick.labels = TRUE) +
  theme(strip.background = element_rect(fill = "white")) 

###Making the plots
world <- ne_countries(scale = "medium", returnclass = "sf")

theme_set(theme_cowplot())

Kallithea_Linvill_plot<-FieldPlot_Kallithea +
  scale_fill_continuous_sequential(AlberPalettes[[1]],breaks = seq(0,3.5,by=0.5)) + 
  geom_sf(data = world,inherit.aes = F, fill = NA, colour = "grey", size = 1) +
  coord_sf(xlim = c(min(dat$Lon), max(dat$Lon)), ylim = c(min(dat$Lat), max(dat$Lat))) +
  geom_point(data = dat, inherit.aes = F, aes(Lon, Lat), size = 3) +
  geom_point(data = dat, inherit.aes = F, aes(Lon, Lat), size = 2, colour = "white") +
  labs(title = "Kallithea virus") +
  FieldPlot_Linvill +
  scale_fill_continuous_sequential(AlberPalettes[[1]],breaks = seq(1,7,by=1.5)) + 
  geom_sf(data = world,inherit.aes = F, fill = NA, colour = "grey", size = 1) +
  coord_sf(xlim = c(min(dat$Lon), max(dat$Lon)), ylim = c(min(dat$Lat), max(dat$Lat))) +
  geom_point(data = dat, inherit.aes = F, aes(Lon, Lat), size = 3) +
  geom_point(data = dat, inherit.aes = F, aes(Lon, Lat), size = 2, colour = "white") +
  labs(title = "Linvill Road virus") 

Galbut_plot<-FieldPlot_Galbut +
  scale_fill_continuous_sequential(AlberPalettes[[1]],breaks = seq(1,9,by=2)) + 
  geom_sf(data = world,inherit.aes = F, fill = NA, colour = "grey", size = 1) +
  coord_sf(xlim = c(min(dat$Lon), max(dat$Lon)), ylim = c(min(dat$Lat), max(dat$Lat))) +
  geom_point(data = dat, inherit.aes = F, aes(Lon, Lat), size = 3) +
  geom_point(data = dat, inherit.aes = F, aes(Lon, Lat), size = 2, colour = "white") +
  labs(title = "Galbut EVE") + theme(panel.spacing = unit(3,"lines"))

inlabru::multiplot(plotlist = list(Kallithea_Linvill_plot,Galbut_plot),cols = 1)

##Now with the Darren scale

Kallithea_Linvill_plot<-FieldPlot_Kallithea +
  scale_fill_gradientn(colours = Darren_scale(11)[-(11)],breaks = seq(0.5,3.5,by=0.5)) + 
  geom_sf(data = world,inherit.aes = F, fill = NA, colour = "grey", size = 1) +
  coord_sf(xlim = c(min(dat$Lon), max(dat$Lon)), ylim = c(min(dat$Lat), max(dat$Lat))) +
  geom_point(data = dat, inherit.aes = F, aes(Lon, Lat), size = 3) +
  geom_point(data = dat, inherit.aes = F, aes(Lon, Lat), size = 2, colour = "white") +
  labs(title = "Kallithea virus") + theme(plot.title = element_text(hjust = 0.5,size = 16, face = "plain")) +
  FieldPlot_Linvill +
  scale_fill_gradientn(colours = Darren_scale(11)[-(11)],breaks = seq(1,7,by=1.5)) +  
  geom_sf(data = world,inherit.aes = F, fill = NA, colour = "grey", size = 1) +
  coord_sf(xlim = c(min(dat$Lon), max(dat$Lon)), ylim = c(min(dat$Lat), max(dat$Lat))) +
  geom_point(data = dat, inherit.aes = F, aes(Lon, Lat), size = 3) +
  geom_point(data = dat, inherit.aes = F, aes(Lon, Lat), size = 2, colour = "white") +
  labs(title = "Linvill Road virus") + theme(plot.title = element_text(hjust = 0.5,size = 16, face = "plain"))

Galbut_plot<-FieldPlot_Galbut +
  scale_fill_gradientn(colours = Darren_scale(11)[-(11)],breaks = seq(1,9,by=2)) + 
  geom_sf(data = world,inherit.aes = F, fill = NA, colour = "grey", size = 1) +
  coord_sf(xlim = c(min(dat$Lon), max(dat$Lon)), ylim = c(min(dat$Lat), max(dat$Lat))) +
  geom_point(data = dat, inherit.aes = F, aes(Lon, Lat), size = 3) +
  geom_point(data = dat, inherit.aes = F, aes(Lon, Lat), size = 2, colour = "white")+
  #labs(title = "Galbut EVE") +
  theme(panel.spacing = unit(4,"lines"),strip.text.x = element_text(size = 16, face = "plain",vjust = 1.5),
        plot.title = element_text(hjust = 0.5,size = 16, face = "plain"))

inlabru::multiplot(plotlist = list(Kallithea_Linvill_plot,Galbut_plot),cols = 1)

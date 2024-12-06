
# full code for BIOL806 final project

# load packages 
library(sf) #for GIS data
library(tidyverse) 
library(ggplot2)
library(ggspatial) #plot maps
library(rnaturalearth) 
library(rnaturalearthdata)
library(ggmap)
library(patchwork) # for combining plots (creating an inset map)
library(basemaps) 
library(mapedit) # allows for interactive map selection in basemaps
library(maps) # basic maps package
library(dunn.test) #Dunn test for MP percent type
library(FSA) #for dunn test
library(car)
library(multcomp)
library(lme4) # for fitting glmms 
library(sjPlot) # for visualizing mixed models 
library(stringr)



# import primary dataset 

#load HSE MP data - this is main dataset with info on sites and MP counts 
MP <- read_csv('MP_no_Blackwater.csv')

#Correct column types 
MP$Site<-as.factor(MP$Site)
MP$`Sample Name`<-as.factor(MP$`Sample Name`)
MP$Estuary_Position<-as.factor(MP$Estuary_Position)
MP$Marsh_Position<-as.factor(MP$Marsh_Position)
MP$Depth<-as.factor(MP$Depth)
MP$Total.MP<-as.numeric(MP$Total.MP)


# SITE MAP CHUNK

#load shp file with site locations 
MP_geo <- st_read("/Users/gabriellejarrett/Desktop/HSEgeodata/MP_Sediment_Cores-point.shp") %>% 
  mutate(Name = gsub("Browns River", "Brown's", Name), # change names to match between MP and MP_geo
         Name = gsub("Taylor River", "Taylor's", Name))
st_crs(MP_geo) # retrieve reference coordinate system - lat lon here 
# MP_geo crs = 4326

# join geo data to MP - get coordinates for each site
MP_geo_combined <- MP %>% 
  left_join(MP_geo, by = c("Sample Name" = "Name"))

# separate lat and lon into separate columns 
MP_geo_sep <- MP_geo_combined %>%
  separate(geometry, into = c("lon", "lat"), sep = ", ") %>%
  mutate(lon=as.numeric(gsub("c\\(","",lon)),
         lat=as.numeric(gsub("\\)","",lat)))

# pull in Google maps satellite imagery of HSE 
register_google(key = "AIzaSyC969CwmD5ID53tIQsq1SGMeHqZkuBh3Uc")  # create API key for google maps
location <- c(left = -70.9, bottom = 42.85, right = -70.8, top = 42.95) #defining coordinate boundaries for HSE 

#pull in satellite map for HSE 
estuary_map <- get_googlemap(center = c(lon = -70.82, lat = 42.90),
                             zoom = 13, maptype = "satellite")
ggmap(estuary_map)

# add sample site coord points onto HSE map: 
site_map <- ggmap(estuary_map) +
  geom_point(data = MP_geo_sep %>% 
               filter(Marsh_Position == 'High'), # remove low marsh - only have one point per site on map
             aes(x = lon, y = lat, shape = Estuary_Position, color = Estuary_Position), size = 4) +
  #  aes(x = lon, y = lat, shape = Estuary_Position, color = Estuary_Position, size = Total.MP)) + #look at how many MP found in diff parts of estuary visually - don't love
  theme_void() +
  scale_shape_manual(values = c("Upper" = 16, "Lower" = 17)) +  # Customize shapes - diff numbers
  scale_color_manual(values = c("Lower" = 'lightsteelblue2', "Upper" = 'steelblue')) +
  labs(title = "Hampton-Seabrook Estuary",
       x = "Longitude", y = "Latitude", shape = "Estuary Position",
       color = "Estuary Position") + 
  theme(plot.title = element_text(hjust = 0.5, size = 15))  #center title 


site_map

# Adding inset in corner showing larger region with box around HSE:
# get map of east coast from google maps 
nh_sat <- get_googlemap(center = c(lon = -70.82, lat = 42.90),
                        zoom = 8, maptype = "satellite") 

# define coordinates for box around HSE (roughly)
xmin <- -70.9  #western boundary
xmax <- -70.8  #eastern boundary
ymin <- 42.84  #southern boundary
ymax <- 43.0   #northern boundary

# graph map with box on top
nh_sat2 <- ggmap(nh_sat) +
  annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
           fill = NA, color = "red", linetype = "solid", size = 1) +
  theme_void() + 
  theme(panel.border = element_rect(color = "black", size = 2, fill = NA)) #add black border

# insert nh map into corner 
site_map2 <- site_map + inset_element(nh_sat2,  
                                      left = 0.41,  # decrease 'left' to shift the inset further left
                                      bottom = 0.01,  # keep 'bottom' constant or lower for alignment
                                      right = 1.2,  # increase 'right' to expand width
                                      top = 0.4,  # increase 'top' to expand height
                                      align_to = "panel") 
site_map2
# site map - DONE! 




# MP present in all samples figure: 
site_MP <- ggplot(MP %>% 
                    group_by(Site, Estuary_Position) %>% 
                    summarize(mean_MP = mean(Total.MP)) %>% 
                    mutate(site_number = Site)
                  , aes(x=Site, y=mean_MP)) +
  geom_bar(stat = 'identity', fill = "lightsteelblue2") +
  theme(strip.text.x = element_text(size=11),
        strip.background = element_rect(fill="beige"), legend.title = element_blank()) +
  theme(axis.text.x = element_text(size = 11)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text = element_text(size=12),
        axis.text.x = element_text(angle = 40, vjust = 1.1, hjust = 1), 
        axis.title = element_text(size=14)) +  
  labs(x ="Site", y = "Mean Microplastic Abundance")




# MP count by estuary and marsh figure
# MP abundance with all 3 variables of interest (sediment depth, marsh habitat, estuary position) ! 

# get SE 
MP_SE <- MP %>%
  group_by(Depth, Marsh_Position, Estuary_Position) %>%
  summarize(n = n(), Mean.Total.MP = mean(Total.MP, na.rm = TRUE), SD = sd(Total.MP)) %>%
  mutate(SE = SD/sqrt(n))

MP.Graph.Estuary <- ggplot(MP_SE, aes(x = Mean.Total.MP, y = Depth, fill = Marsh_Position)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(xmin = Mean.Total.MP-SE, xmax = Mean.Total.MP+SE), width = .2,
                position =position_dodge(.9)) +
  facet_wrap(~Estuary_Position, 
             labeller = as_labeller(c("Lower" = "Lower Estuary", 
                                      "Upper" = "Upper Estuary"))) + #rename lower and upper to include "estuary" 
  scale_fill_manual(values = c("lightsteelblue2", "steelblue"), labels = c("High Marsh", "Low Marsh")) +
  scale_y_discrete(limits = rev) +
  theme(legend.text = element_text(size = 11), legend.title = element_text(size = 12)) +
  labs(fill = "Marsh Position") + #rename legend title 
  theme(strip.text.x = element_text(size = 11),
        strip.background = element_rect(fill = "beige")) + #USED TO BE ELEMENT_BLANK()
  theme(axis.text.x = element_text(size = 11)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14)) +  
  labs(x ="Mean Microplastic Abundance", y = "Sediment Depth (cm)") 





# MP count by sediment depth figure:
#get SE 
MP_Means_Depth <- MP %>%
  group_by(Depth) %>%
  summarize(n = n(), Mean.Total.MP = mean(Total.MP, na.rm = TRUE), SD_Depth = sd(Total.MP)) %>%
  mutate(SE_Depth = SD_Depth/sqrt(n)) 


MP.Depth <- ggplot(MP_Means_Depth, aes(x = Mean.Total.MP, y = Depth)) +
  geom_col(position = "dodge", fill = "lightsteelblue2") + 
  scale_y_discrete(limits = rev) + #switch order of y axis - want 0-2 at top to mimic sediment depth profile 
  theme(strip.text.x = element_text(size=11),
        strip.background = element_rect(fill="beige"), legend.title = element_blank()) +
  theme(axis.text.x = element_text(size = 11)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14)) +  
  labs(x ="Mean Microplastic Abundance", y = "Sediment Depth (cm)") +
  geom_errorbar(aes(xmin = Mean.Total.MP-SE_Depth, xmax = Mean.Total.MP+SE_Depth), width = .2,
                position = position_dodge(.9)) 





# MP type percent 

# take average of each MP type by depth 
MP_type_avg <- MP %>% 
  group_by(Depth) %>% 
  summarize(mean_total_MP = mean(Total.MP), mean_fiber = mean(`TOTAL Fibers`), mean_fiber_bundle = mean(`TOTAL Fiber Bundles`), mean_frag = mean(`Total Fragments`), mean_pellet = mean(`TOTAL Pellets`), mean_film = mean(`TOTAL Films`), mean_foam = mean(`TOTAL Foams`), mean_styrofoam = mean(`TOTAL Styrofoam`))

# calculate percent of total for each MP type 
percent_avg <- MP_type_avg %>% 
  mutate(percent_fiber = (mean_fiber/mean_total_MP), percent_fiber_bundle = (mean_fiber_bundle/mean_total_MP), percent_frag = (mean_frag/mean_total_MP), percent_pellet = (mean_pellet/mean_total_MP), percent_film = (mean_film/mean_total_MP), percent_foam = (mean_foam/mean_total_MP), percent_styrofoam = (mean_styrofoam/mean_total_MP))

# pivot percent_avg longer 
percent_avg_longer <- percent_avg %>% 
  pivot_longer(cols = c(percent_fiber, percent_fiber_bundle, percent_frag, percent_pellet, percent_film, percent_foam, percent_styrofoam), 
               names_to = "Type",
               names_prefix = "percent.",
               values_to = "Percent",
               values_drop_na = TRUE) %>% 
  mutate(Percent2 = (Percent*100)) #percents out of 100 not 1

# graph percent_avg_longer 
MP_type_percent_depth <- ggplot(percent_avg_longer %>% 
                                  mutate(Type = as.factor(Type)) %>% 
                                  mutate(Type = fct_reorder(Type, Percent2))
                                , aes(x = Percent2, y = Depth, fill = Type)) +
  geom_bar(stat = 'identity') +
  scale_y_discrete(limits=rev) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab('Percent Microplastic Type') + ylab('Sediment Depth (cm)') +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_fill_manual(name = 'Microplastic Type', 
                    labels = c('Percent Styrofoams', 'Percent Foams', 'Percent Pellets', 'Percent Films', 'Percent Fiber Bundles', 'Percent Fragments', 'Percent Fibers'),
                    values = c("#fdbf6f", "#e31b1b", "#fb9a99", "#34a02d", "#b2df8a", "#1f78b4", "#a6cee3"))




#stats 

#using Kruskal-Wallis for nonparametric one-way
MP_type_kruskal <- kruskal.test(Percent2 ~ Type, data = percent_avg_longer)
MP_type_kruskal

dunn_test <- dunnTest(percent_avg_longer$Percent2 ~ percent_avg_longer$Type, method = "bonferroni")


# Extract results and format as a data frame
dunn_results <- as.data.frame(dunn_test$res)

# Beautify the table with knitr::kable
dunn_table <- dunn_results %>%
  mutate(
    Z = round(Z, 3),          # Round Z-scores
    P.adj = formatC(P.adj, format = "e", digits = 3) # Format adjusted p-values
  ) %>%
  knitr::kable(
    caption = "Results of Dunn's Test with Adjusted P-values",
    col.names = c("Comparison", "Z-Score", "P-Value", "Adjusted P-Value"),
    align = c("l", "c", "c", "c"),
    format = "html"
  ) %>%
  kableExtra::kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE
  )



# MAPS 

gdb_path <- "/Users/gabriellejarrett/Desktop/MPAnthro.gdb"  
st_layers(gdb_path)  # List layers to identify exact names



#pollution discharge sites layer:
pollution_discharge_sites <- st_read(gdb_path, layer = "PollutionDischargeSites") 

pollution_discharge_sites <- pollution_discharge_sites %>% 
  st_transform(crs = 4326) %>%  #change UTM to lat long  
  mutate(lon = st_coordinates(.)[, 1],  #extract lon from POINT
         lat = st_coordinates(.)[, 2])  #extract lon from POINT

# map pollution discharge sites onto HSE map - pull in satellite imagery for HSE again
estuary_map2 <- get_googlemap(center = c(lon = -70.82, lat = 42.90),
                              zoom = 13, maptype = "satellite")
poll_source_map <- ggmap(estuary_map2) +
  geom_point(data = pollution_discharge_sites %>% 
               filter(CATEGORY == 'WW'), #only interested in WW (wastewater treatment effluenct sites) in HSE - waste water effluent known source of MP 
             aes(x = lon, y = lat, color = 'hotpink'), size = 3) +
  labs(title = "Waste Water Effluent Discharge Sites") +
  theme(plot.title = element_text(hjust = 0.5, size = 15))  #center title 
poll_source_map

# okay slay made the map ! What now lol ??  How do I analyze this? What does this mean in relation to my data/site?? 


#map WWE pollution sites and HSE sample sites on same map! 
poll_and_sites <- ggmap(estuary_map2) + 
  # Sample site points (High marsh only)
  geom_point(data = MP_geo_sep %>%  
               filter(Marsh_Position == 'High'), # Filter for High marsh
             aes(x = lon, y = lat, color = Marsh_Position, shape = Marsh_Position), 
             size = 3) +
  # Pollution discharge sites (Wastewater effluent only)
  geom_point(data = pollution_discharge_sites %>%  
               filter(CATEGORY == 'WW'), # Filter for wastewater effluent
             aes(x = lon, y = lat, color = CATEGORY, shape = CATEGORY), 
             size = 3) +
  scale_shape_manual(values = c("High" = 19, "WW" = 17),  
                     labels = c("High" = "Sample Sites", "WW" = "Pollution Discharge Sites")) + #customize legend labels 
  scale_color_manual(values = c("High" = "steelblue", "WW" = "red3"), 
                     labels = c("High" = "Sample Sites", "WW" = "Pollution Discharge Sites")) +
  labs(shape = "Site Type", color = "Site Type") +
  theme_void() +
  labs(title = "Pollution Discharge Sites in HSE") +
  theme(plot.title = element_text(hjust = 0.5, size = 15))





#FLD hazard map
flood_haz <- st_read(gdb_path, layer = "FLD_HAZ_NH_HSE")

flood_haz <- flood_haz %>% 
  st_transform(crs = 4326) #change UTM to lat long  

# map flood haz onto HSE estuary map
estuary_map3 <- get_googlemap(center = c(lon = -70.82, lat = 42.90),
                              zoom = 13, maptype = "satellite")

# add flood_haz layer on top of HSE 
fld_haz_map <- ggmap(estuary_map3) +
  geom_sf(data = flood_haz, aes(fill = FLD_ZONE), inherit.aes = FALSE, alpha = 0.35, color = NA) +
  scale_fill_manual(name = "Flood Zone Codes",  # change legend title 
                    values = c("AE" = "red3", 
                               "X" = "lightgreen", 
                               "A" = "yellow", 
                               "AO" = "orange", 
                               "VE" = "hotpink"
                    )) +
  theme_minimal() +
  labs(title = "HSE Flood Hazard Zones") +
  theme(plot.title = element_text(hjust = 0.5, size = 15))  #center title 


fld_haz_map
# Flood zone codes: starts w A or V = high risk flood area ; X = low risk flood area 
# estuary is high risk flood area - coastal wetland 
# flooding/extreme weather events known to be significant events to transport MP - potential factor affecting MP deposition in salt marsh sediment in HSE - HSE is high flood risk, flood risk only going to increase w SLR (Eastern US coast hot spot for SLR)

flood_haz$FloodRisk <- ifelse(flood_haz$FLD_ZONE %in% c("AE", "A", "AO", "VE"), "High Flood Risk", "Low Flood Risk")


fld_and_sites <- ggmap(estuary_map3) +
  geom_sf(data = flood_haz, aes(fill = FloodRisk), inherit.aes = FALSE, alpha = 0.4, color = NA) +
  scale_fill_manual(name = "Flood Zones",  # change legend title 
                    values = c("High Flood Risk" = "red3", 
                               "Low Flood Risk" = "lightgreen")) +
  theme_void() +
  labs(title = "HSE Flood Hazard Zones") +
  theme(plot.title = element_text(hjust = 0.5, size = 15)) +  #center title 
  geom_point(data = MP_geo_sep %>%   #sample site points: 
               filter(Marsh_Position == 'High'), # remove low marsh - only have one point per site on map
             aes(x = lon, y = lat, shape = Estuary_Position, color = Estuary_Position), size = 3)  +
  scale_shape_manual(values = c("Upper" = 16, "Lower" = 17)) +  # Customize shapes - diff numbers
  scale_color_manual(values = c("Lower" = 'lightsteelblue2', "Upper" = 'steelblue')) +
  labs(shape = "Estuary Position", color = "Estuary Position")


fld_and_sites 




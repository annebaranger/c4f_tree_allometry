## spatial cross-validation
library(targets)
library(sf)
library(dplyr)

tar_load(mod.data)
civ=read_sf(dsn="data/civ/",layer="COTE D'IVE") |> 
  st_transform(crs="epsg:4326")

plot.data<-mod.data |> 
  dplyr::select(database,id_plot,system,lat,long) |> 
  unique() 

plot.data |> 
  ggplot()+
  geom_point(aes(x=long,y=lat,color=system))+
  geom_sf(data=civ,fill=NA)
  

box<-st_bbox(civ)
breaks_y=seq(box[["ymin"]],box[["ymax"]],length.out=4)
breaks_x=seq(box[["xmin"]],box[["xmax"]],length.out=4)

plot.data |> 
  mutate(long.cat=cut(long,breaks = breaks_x),
         lat.cat=cut(lat,breaks = breaks_y),
         subdataset=as.numeric(as.factor(paste0(long.cat,lat.cat)))) |> 
  ggplot()+
  geom_point(aes(x=long,y=lat,color=as.factor(subdataset)))+
  geom_sf(data=civ,fill=NA)

km=kmeans(as.matrix(plot.data[c("long","lat")]),centers=7)
km.center=as.data.frame(km$centers) |> tibble::rownames_to_column(var="cls") |> 
  mutate(cluster_near=NA,
         cls=as.numeric(cls))
for (i in 1:dim(km.center)[1]){
  print(i)
  if(km.center$lat[i]>quantile(mod.data$lat,probs = 0.99)[[1]]){
    ind=which.min(st_distance(st_as_sf(km.center[i,],
                         coords=c("long","lat")),
                st_as_sf(km.center[-i,],
                         coords=c("long","lat"))
    ))
  km.center$cluster_near[i]=as.numeric(km.center[-i,][ind,"cls"])
  }
}

cbind(plot.data,
      cluster=km$cluster) |> 
  left_join(km.center[,c("cls","cluster_near")],by=c("cluster"="cls")) |> 
  mutate(cluster=case_when(!is.na(cluster_near)~cluster_near,
                           TRUE~cluster)) |> 
  group_by(cluster) |> 
  mutate(n_clust=n()) |> 
  ungroup() |> 
  ggplot()+
  geom_point(aes(x=long,y=lat,color=as.factor(n_clust)))+
  # geom_point(data=km.center,aes(x=long,y=lat),color="red")+
  geom_sf(data=civ,fill=NA)+
  geom_hline(yintercept =   quantile(mod.data$lat,probs = 0.99)[[1]])


cbind(plot.data,
      cluster=km$cluster) |> 
  left_join(km.center[,c("cls","cluster_near")],by=c("cluster"="cls")) |> 
  mutate(cluster=case_when(!is.na(cluster_near)~cluster_near,
                           TRUE~cluster)) |> 
  group_by(cluster) |> 
  mutate(n_clust=n()) |> 
  ungroup() |> 
  ggplot(aes(as.factor(cluster),fill=system))+
  geom_bar()


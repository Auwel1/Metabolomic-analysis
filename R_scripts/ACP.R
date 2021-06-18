# Chargement des librairies 
library(BBmisc)
library(ade4)
library(FactoMineR)
library(factoextra)

#Chargement de données 
setwd("/run/media/aurelien/ACOFFE/Stage/pipeline_python/pythonProject/data_frame_results/")
area = read.csv("geom_mean_area_means.csv")

setwd("/run/media/aurelien/ACOFFE/Stage/pipeline_python/pythonProject/data_frames_results/")
enr = read.csv("result_enrichement.csv")



# test jeu de donnée avec acides aminés séparés
AA = c("Glutamate","Alanine","Aspartate","Glutamine",
       "Glycine","Serine")

areanoaa = area
enrnoaa = enr

for (aa in AA){
  areanoaa  <- subset(areanoaa, areanoaa$X != aa)
  enrnoaa <- subset(enrnoaa, enrnoaa$X != aa)
}

noaa_a_means = areanoaa[grep('.*_mean', colnames(areanoaa))]
rownames(noaa_a_means) <- areanoaa$X

# T48AB_area = areanoaa$T48_AB_mean
# 
# T48Cont_area = areanoaa$T48_Cont_mean
# 
# T48AB_enr = enrnoaa$T48_AB_mean
# T48Cont_enr = enrnoaa$T48_Cont_mean 
# 
# plot(T48AB_area,T48AB_enr, xlim = c(0, 55))
# text(T48AB_area, T48AB_enr,  labels = areanoaa$X, cex = 0.7, pos = 3)


# extraction des sous tableaux de moyenne
area_means = area[grep('.*_mean', colnames(area))]
rownames(area_means) <- area$X


#transposition 
transpose = data.frame(t(area_means))
colnames(transpose) <- rownames(area_means)

# ACP
global_PCA = dudi.pca(transpose, scannf = FALSE, nf = 2)
pca2 = PCA(transpose)
pca2$var$contrib

fviz_eig(global_PCA, addlabels = TRUE)

fviz_pca_var(global_PCA)
fviz_cos2(global_PCA, choice = "var", axes = 1:2)

fviz_contrib(res.pca, choice = "ind", axes = 1:2)
var <- get_pca_var(global_PCA)

# visu ACP par lignées 
group = rownames(transpose)
group[grepl('.*_A_.*', rownames(transpose))] <- 'A'
group[grepl('.*_B_.*', rownames(transpose))] <- 'B'
group[grepl('.*_AB_.*', rownames(transpose))] <- 'AB'
group[grepl('.*_Cont_.*', rownames(transpose))] <- 'Cont'
group = c(group)


group_PCA = fviz_pca_ind(X = global_PCA,
             col.ind = group, palette = c("khaki2", "green3", "red", "deepskyblue4"),
             repel = TRUE, legend.title = "Cell line", mean.point = FALSE,
             geom = c("point", "text"),pointshape = 19 , )
ggpubr::ggpar(group_PCA, title = "Average metabolite area - PCA")


# visu ACP par points de temps 
time_point = rownames(transpose)
time_point[grepl('T0_.*', rownames(transpose))] <- 'T0'
time_point[grepl('T24_.*', rownames(transpose))] <- 'T24'
time_point[grepl('T48_.*', rownames(transpose))] <- 'T48'

TP_PCA = fviz_pca_ind(X = global_PCA,
             col.ind = time_point, palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE, legend.title = "Time Point", mean.point = FALSE )
ggpubr::ggpar(TP_PCA, title = "Average metabolite area - PCA")

# Graphique des contributions pour chaque dimension

fviz_contrib(global_PCA, choice = "var", axes = 1)

fviz_contrib(global_PCA, choice = "var", axes = 2)

# Plot ACP par point de temps et lignées 
uwu = data.frame("PC1" = global_PCA$li$Axis1, "PC2" = global_PCA$li$Axis2)

ggplot(uwu, aes( x = PC1, y = PC2, color = group, shape = time_point ))+
  geom_point(size = 5)+ 
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2)+ 
  xlab("PC 1 (42.4%)") +
  ylab("PC 2 (30.6%)") +
  scale_color_manual(breaks = group, values = c("red","khaki2", "green3", "deepskyblue4"),
                     name = "Cell line")+
  ggtitle("PCA- Mean area values") +
  labs(shape = "Time point")

# ACP sur les métabolites 
# Sans les acides aminés
PCA2 = dudi.pca(noaa_a_means, scannf = FALSE, nf =2)
fviz_eig(PCA2, addlabels = TRUE)

fviz_pca_var(PCA2)
fviz_contrib(PCA2, choice = "var", axes = 1)

fviz_contrib(PCA2, choice = "var", axes = 2)

fviz_pca_ind(X = PCA2)
fviz_pca_biplot(PCA2)


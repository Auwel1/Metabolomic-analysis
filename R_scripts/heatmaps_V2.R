setwd('/run/media/aurelien/ACOFFE/Stage/pipeline_python/pythonProject/data_frame_results/')
# Charger les librairies
library(pheatmap)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(psych)
library(dplyr)

#Charger les données 
area = read.csv("geom_mean_area_means.csv", row.names = 1)

setwd('/run/media/aurelien/ACOFFE/Stage/pipeline_python/pythonProject/data_frames_results/')
enr = read.csv("result_enrichement.csv", row.names = 1)

setwd('/run/media/aurelien/ACOFFE/Stage/jeu_données_brut/')

all_data = readxl::read_xlsx('dec2020_ALL_DATA_res_isocor.xlsx')

names = list('T0_Cont', 'T24_Cont', 'T48_Cont',
             'T0_A', 'T24_A', 'T48_A', 
             'T0_B', 'T24_B','T48_B',
             'T0_AB','T24_AB', 'T48_AB')

cellnumber = data.frame(2E6, 3.1E6,4.4E6,
                        2E6,2.8E6,3.6E6,
                        2E6,2.9E6,3.3E6,
                        2E6,2E6,6E5)
colnames(cellnumber) <- names
metabolites = rownames(area)
# applique la correction vis à vis du nombre de  cellules



# Séparer les colonnes de moyennes et d'écart types 
pattern_mean = '.*_mean'
pattern_sd = '.*_sd'
samples = c('A', 'B', 'AB')
states = c('marked_only', 'non_marked', 'total')

area_means = area[grep(pattern_mean, names(area))]

times = c('T0', 'T24', 'T48')
control = c('T0_Cont', 'T24_Cont', 'T48_Cont')


area_sd = area[grep(pattern_sd, names(area))]

enr_means = enr[grep(pattern_mean, names(enr))]
non_enr_means = apply(enr_means,2, function(x) 100 -x) 

area_marked_only = area_means*(enr_means/100)
area_non_marked_only = area_means*(non_enr_means/100)

for(n in names){
  print(n)
  if(length(grep('.*_A',n)) !=0){
    all_data[grep(paste0(n,'[^B]'), all_data$sample),]$corrected_area <-
      unlist(sapply(all_data[grep(paste0(n,'[^B]'), all_data$sample),]$corrected_area, 
                    function(x) x/cellnumber[n], simplify = "array"))
    
  }
  else{
    all_data[grep(n, all_data$sample),]$corrected_area <-
      unlist(sapply(all_data[grep(n, all_data$sample),]$corrected_area, 
                    function(x) x/cellnumber[n], simplify = "array"))
    
  }}


enr_means =enr_means[order(rownames(enr_means)),]

point = 0
table_list = list()
for (subname in names){
  
  #divise les données en tableaux par échantillons
  tmp = assign(subname, subset(all_data, grepl(paste0(subname,'_[[:digit:]]'),all_data$sample)))
  
  #fait la somme des isotopologues pour chaque réplicat d'échantillon 
  t_en = assign(paste0(subname,'_group'),
                aggregate(tmp$corrected_area, by = list(tmp$sample, tmp$metabolite) , FUN = sum))
  

  if(point ==0 ){
    # Créer un tableau par condition (total, enrichie, non enrichie) 
    mean_no_mark_replicates = data.frame("metabolites" = t_en$Group.2)
    mean_marked_replicates = data.frame("metabolites" = t_en$Group.2)
    mean_area_replicates = data.frame("metabolites" = t_en$Group.2)
    
  }
  # Remplis chaque tableau
  mean_area_replicates[subname] = t_en$x
  mean_no_mark_replicates[subname] = t_en$x
  mean_marked_replicates[subname] = t_en$x
  
  for(m in metabolites){
    t_mnm = subset(mean_no_mark_replicates, metabolites == m)
    t_mnm[,subname] <- t_mnm[,subname]*(non_enr_means[m,paste0(subname,'_mean')]/100) 
    mean_no_mark_replicates[mean_no_mark_replicates["metabolites"] == m,] <- t_mnm
    t_mo = subset(mean_marked_replicates, metabolites == m)
    t_mo[,subname] <- t_mo[,subname]*(enr_means[m,paste0(subname,'_mean')]/100)
    mean_marked_replicates[mean_marked_replicates["metabolites"] == m,] <- t_mo
  }
  # Enregistre les tableaux complets puis réinitialise
  if( point == 2){
    table_list[paste0('m_a_r_', subname)] <- list(assign(paste0('m_a_r_', subname), 
                                                         data.frame(mean_area_replicates)))
    
    table_list[paste0('m_a_r_', subname,'_m_only')] <- list(assign(paste0('m_a_r_', subname,'_m_only'), 
                                                                   data.frame(mean_marked_replicates)))
    
    table_list[paste0('m_a_r_', subname,'_nm_only')] <- list(assign(paste0('m_a_r_', subname,'_nm_only'), 
                                                                    data.frame(mean_no_mark_replicates)))
    
    
    point <- 0
  }
  else{
    point = point + 1
  }
}



pt = 0
for(n in ls()){
  if(length(grep("_group", n))!= 0 & length(grep("Cont_group", n))== 0 & length(grep("std", n)) == 0 &
     length(grep("FC", n)) == 0){
    
    pt = pt+1
    
    nam = n
    print(nam)
    forrich = as.list(strsplit(nam, '_')[[1]])
    tmp = get(n)
    
    a = tmp[grep('_1', tmp$Group.1 ), ]
    rownames(a) = a$Group.2
    a = arrange(a)
    am = a; anm = a
    am['x'] <- a$x*(enr_means[grep(paste0(forrich[1],'_',forrich[2],'_'), colnames(enr_means))]/100)
    anm['x'] <- a$x*((100 - enr_means[grep(paste0(forrich[1],'_',forrich[2]), colnames(enr_means))])/100)
    
    b = tmp[grep('_2', tmp$Group.1 ), ]; bnm = b; bm = b
    bm['x'] <- b$x*(enr_means[grep(paste0(forrich[1],'_',forrich[2],'_'), colnames(enr_means))]/100)
    bnm['x'] <- b$x*((100 - enr_means[grep(paste0(forrich[1],'_',forrich[2]), colnames(enr_means))])/100)
    
    c = tmp[grep('_3', tmp$Group.1 ), ]; cm = c; cnm =c 
    cm['x'] <- c$x*(enr_means[grep(paste0(forrich[1],'_',forrich[2],'_'), colnames(enr_means))]/100)
    cnm['x'] <- c$x*((100 - enr_means[grep(paste0(forrich[1],'_',forrich[2]), colnames(enr_means))])/100)
    
    
    for(t in times){
      if(length(grep(t, n)) !=0){
        u = t
      }
    }
    
    cont = control[grep(u, control)]
    
    a2 = subset(get(paste0(cont,'_group')), get(paste0(cont,'_group'))$Group.1 == paste0(u,'_Cont_1'))
    a2m = a2; a2nm = a2
    a2m['x'] <- a2$x*(enr_means[grep(cont, colnames(enr_means))]/100)
    a2nm['x'] <- a2$x*((100 - enr_means[grep(cont, colnames(enr_means))])/100)
    
    b2 = subset(get(paste0(cont,'_group')), get(paste0(cont,'_group'))$Group.1 == paste0(u,'_Cont_2'))
    b2m = b2; b2nm = b2
    b2m['x'] <- b2$x*(enr_means[grep(cont, colnames(enr_means))]/100)
    b2nm['x'] <- b2$x*((100 - enr_means[grep(cont, colnames(enr_means))])/100)
    
    c2 = subset(get(paste0(cont,'_group')), get(paste0(cont,'_group'))$Group.1 == paste0(u,'_Cont_3'))
    c2m = c2; c2nm = c2
    c2m['x'] <- c2$x*(enr_means[grep(cont, colnames(enr_means))]/100)
    c2nm['x'] <- c2$x*((100 - enr_means[grep(cont, colnames(enr_means))])/100)
    
    compare = cbind(a$x,b$x,c$x,a2$x,b2$x,c2$x)
    comparem = cbind(am$x,bm$x,cm$x,a2m$x,b2m$x,c2m$x)
    comparenm = cbind(anm$x,bnm$x,cnm$x,a2nm$x,b2nm$x,c2nm$x)
    row.names(compare) = a$Group.2; row.names(comparem) = a$Group.2; row.names(comparenm) = a$Group.2
    compare = data.frame(compare);comparem = data.frame(comparem); comparenm = data.frame(comparenm)
    compare[compare == 0] <- NaN; comparem[comparem == 0] <- NaN; comparenm[comparenm == 0] <- NaN
    compare['std'] = apply(compare, 1, function(x)  sd(x, na.rm = TRUE))
    comparem['std'] = apply(comparem, 1, function(x)  sd(x, na.rm = TRUE))
    comparenm['std'] = apply(comparenm, 1, function(x)  sd(x, na.rm = TRUE))
    
    
    
    
    assign(paste0('std_',nam,'_vs_',cont), compare); comparebis = compare/compare$std
    assign(paste0('std_',nam,'_vs_',cont,'_mark'), comparem) ; comparembis = comparem/comparem$std
    assign(paste0('std_',nam,'_vs_',cont,'_nomark'), comparenm); comparenmbis = comparenm/comparenm$std
    
    comparebis[c(1:6)][is.na(comparebis[c(1:6)])] <- 1E-5
    comparembis[c(1:6)][is.na(comparembis[c(1:6)])] <- 1E-5
    comparenmbis[c(1:6)][is.na(comparenmbis[c(1:6)])] <- 1E-5
    
    
    
    if(pt ==1){
      final_datframe_tot = data.frame(row.names = rownames(compare))
      finaldatframe_m = data.frame(row.names = rownames(compare)) 
      finaldatframe_nm = data.frame(row.names = rownames(compare))
    }
    
    
    
    final_datframe_tot[paste0('FC_',nam,'/',cont,'_tot')] = apply(comparebis[c(1,2,3)],
                                                                  1, function(x) exp(mean(log(x))) )/apply(comparebis[c(4,5,6)],1,
                                                                                        function(x) exp(mean(log(x))) )
   
    finaldatframe_m[paste0('FC_',nam,'/',cont,'_m')] = apply(comparembis[c(1,2,3)],
                                                             1, function(x) exp(mean(log(x))) )/apply(comparembis[c(4,5,6)],1, 
                                                                                   function(x) exp(mean(log(x))) )

    finaldatframe_nm[paste0('FC_',nam,'/',cont,'_nm')] = apply(comparenmbis[c(1,2,3)],
                                                               1, function(x) exp(mean(log(x))) )/apply(comparenmbis[c(4,5,6)],1, 
                                                                                     function(x) exp(mean(log(x))) )
    
    
    
  }
}

colnames(final_datframe_tot)
final_datframe_tot = subset(final_datframe_tot, select = c(1,4,7,3,6,9,2,5,8))
finaldatframe_m = subset(finaldatframe_m, select = c(1,4,7,3,6,9,2,5,8))
finaldatframe_nm = subset(finaldatframe_nm, select = c(1,4,7,3,6,9,2,5,8))

arrange(final_datframe_tot)

log_tot = log2(final_datframe_tot)
log_m = log2(finaldatframe_m)
log_nm = log2(finaldatframe_nm)


#Paramètres graphiques 
paletteFunc <- colorRampPalette(c('blue','purple4', 'plum4', 'white','lemonchiffon1','orange','orangered', 'red'))
palette = paletteFunc(48)
scale = seq(-5,5,0.25)

setwd('/run/media/aurelien/ACOFFE/Stage/essais/')
plot_list = list()
FC_list = list()
lFC_list = list()
# Boucle de construction de heatmaps
for(s in samples){
  
  pat = paste0('_',s , '_')
  

  tmp_tot = log_tot[grep(pat, colnames(log_tot))]
  tmp_m = log_m[grep(pat, colnames(log_tot))]
  tmp_nm = log_nm[grep(pat, colnames(log_tot))]
  
  #Création des heatmaps 
  h_t = pheatmap(na.omit(tmp_tot), legend_labels = row.names(tmp_tot), cluster_cols=F,
                 main = paste0(s,'_total'), breaks = scale, color = palette,
                 cluster_rows = F)
  
  h_nm = pheatmap(na.omit(tmp_nm), legend_labels = row.names(tmp_nm), cluster_cols=F,
                  main = paste0(s, '_non_marked'), breaks = scale, color = palette,
                  cluster_rows = F)
  
  h_m = pheatmap(na.omit(tmp_m), legend_labels = row.names(tmp_m), cluster_cols=F, 
                 main = paste0(s, '_marked'), breaks = scale, color = palette,
                 cluster_rows = F)
  #intégration en une liste d'objets 
  plot_list[['a']] = h_t[[4]]
  plot_list[['b']] = h_nm[[4]]
  plot_list[['c']] = h_m[[4]]
  
  #Sortie graphique
  g<-do.call(grid.arrange,plot_list )
  ggsave(paste0(s, '_heatmaps_new.pdf'),g, height = 35, width = 15, device = "pdf")
  plot_list = list()
}


a = T0_Cont_group[grep('_1', T0_Cont_group$Group.1),]
rownames(a) = a$Group.2

am = a; anm = a
am['x'] <- a$x*(enr_means[grep(paste0('T0_Cont','_'), colnames(enr_means))]/100)
anm['x'] <- a$x*((100 - enr_means[grep(paste0('T0_Cont','_'), colnames(enr_means))])/100)

b =  T0_Cont_group[grep('_2', T0_Cont_group$Group.1),]; bnm = b; bm = b
bm['x'] <- b$x*(enr_means[grep(paste0('T0_Cont','_'), colnames(enr_means))]/100)
bnm['x'] <- b$x*((100 - enr_means[grep(paste0('T0_Cont','_'), colnames(enr_means))])/100)

c =  T0_Cont_group[grep('_3', T0_Cont_group$Group.1),]; cm = c; cnm =c 
cm['x'] <- c$x*(enr_means[grep(paste0('T0_Cont','_'), colnames(enr_means))]/100)
cnm['x'] <- c$x*((100 - enr_means[grep(paste0('T0_Cont','_'), colnames(enr_means))])/100)


cont = "T48_Cont"

a2 = subset(T48_Cont_group, T48_Cont_group$Group.1 == paste0("T48",'_Cont_1'))
a2m = a2; a2nm = a2
a2m['x'] <- a2$x*(enr_means[grep(cont, colnames(enr_means))]/100)
a2nm['x'] <- a2$x*((100 - enr_means[grep(cont, colnames(enr_means))])/100)

b2 = subset(T48_Cont_group, T48_Cont_group$Group.1 == paste0("T48",'_Cont_2'))
b2m = b2; b2nm = b2
b2m['x'] <- b2$x*(enr_means[grep(cont, colnames(enr_means))]/100)
b2nm['x'] <- b2$x*((100 - enr_means[grep(cont, colnames(enr_means))])/100)

c2 = subset(T48_Cont_group, T48_Cont_group$Group.1 == paste0("T48",'_Cont_3'))
c2m = c2; c2nm = c2
c2m['x'] <- c2$x*(enr_means[grep(cont, colnames(enr_means))]/100)
c2nm['x'] <- c2$x*((100 - enr_means[grep(cont, colnames(enr_means))])/100)

compare = cbind(a$x,b$x,c$x,a2$x,b2$x,c2$x)
comparem = cbind(am$x,bm$x,cm$x,a2m$x,b2m$x,c2m$x)
comparenm = cbind(anm$x,bnm$x,cnm$x,a2nm$x,b2nm$x,c2nm$x)
row.names(compare) = a$Group.2; row.names(comparem) = a$Group.2; row.names(comparenm) = a$Group.2
compare = data.frame(compare);comparem = data.frame(comparem); comparenm = data.frame(comparenm)
compare['std'] = apply(compare, 1, function(x)  sd(x))
comparem['std'] = apply(comparem, 1, function(x)  sd(x))
comparenm['std'] = apply(comparenm, 1, function(x)  sd(x))


compare[compare == 0] <- 1E-5
comparem[comparem == 0] <- 1E-5
comparenm[comparenm == 0] <- 1E-5


assign(paste0('std_',"T0Cont",'_vs_',cont), compare); comparebis = compare/compare$std
assign(paste0('std_',"T0Cont",'_vs_',cont,'_mark'), comparem) ; comparembis = comparem/comparem$std
assign(paste0('std_',"T0Cont",'_vs_',cont,'_nomark'), comparenm); comparenmbis = comparenm/comparenm$std

comparebis[c(1:6)][comparebis[c(1:6)] == 0] <- 1E-5
comparembis[c(1:6)][comparembis[c(1:6)] == 0] <- 1E-5
comparenmbis[c(1:6)][comparenmbis[c(1:6)] == 0] <- 1E-5


Cont_datframe_tot = data.frame(row.names = a$Group.2)
Cont_datframe_m = data.frame(row.names = a$Group.2) 
Cont_datframe_nm = data.frame(row.names = a$Group.2)



Cont_datframe_tot[paste0('FC_',"T0Cont",'/',cont,'_tot')] = apply(comparebis[c(1,2,3)],
                                                              1, function(x) exp(mean(log(x))) )/apply(comparebis[c(4,5,6)],1,
                                                                                                       function(x) exp(mean(log(x))) )

Cont_datframe_m[paste0('FC_',"T0Cont",'/',cont,'_m')] = apply(comparembis[c(1,2,3)],
                                                         1, function(x) exp(mean(log(x))) )/apply(comparembis[c(4,5,6)],1, 
                                                                                                  function(x) exp(mean(log(x))) )

Cont_datframe_nm[paste0('FC_',"T0Cont",'/',cont,'_nm')] = apply(comparenmbis[c(1,2,3)],
                                                           1, function(x) exp(mean(log(x))) )/apply(comparenmbis[c(4,5,6)],1, 
                                                                                                    function(x) exp(mean(log(x))) )
Comparison_Controls = cbind(Cont_datframe_tot,Cont_datframe_nm,Cont_datframe_m)



# write.csv(Comparison_Controls, file = "/run/media/aurelien/ACOFFE/Stage/essais/heatmapv2/FC_Cont0vsCont48")
# write.csv(final_datframe_tot, file = "/run/media/aurelien/ACOFFE/Stage/essais/heatmapv2/FC_Cont_tot.csv")
# write.csv(finaldatframe_m, file = "/run/media/aurelien/ACOFFE/Stage/essais/heatmapv2/FC_Cont_mark.csv")
# write.csv(finaldatframe_nm, file = "/run/media/aurelien/ACOFFE/Stage/essais/heatmapv2/FC_Cont_nomark.csv")


setwd(dir = "/run/media/aurelien/ACOFFE/Stage/essais/heatmapv2/")
tmp = read.csv(file = "FC_Cont_tot.csv", row.names = 1)
write.csv(log2(tmp), file = "/run/media/aurelien/ACOFFE/Stage/essais/heatmapv2/logFC_Cont_tot.csv")

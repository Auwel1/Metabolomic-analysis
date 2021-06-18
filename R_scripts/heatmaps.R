setwd('/run/media/aurelien/ACOFFE/Stage/pipeline_python/pythonProject/data_frame_results/')
# Charger les librairies
library(pheatmap)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(psych)

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

# applique la correction vis à vis du nombre de  cellules



# Séparer les colonnes de moyennes et d'écart types 
pattern_mean = '.*_mean'
pattern_sd = '.*_sd'


area_means = area[grep(pattern_mean, names(area))]




area_sd = area[grep(pattern_sd, names(area))]

enr_means = enr[grep(pattern_mean, names(enr))]

# Calculs coefficients de variation
var_coeff_tab = area_sd/area_means 
length(var_coeff_tab[var_coeff_tab < 0.5])
# Pourcentage de coefficients de variations > 0.5
(table(var_coeff_tab > 0.5)["TRUE"]/table(var_coeff_tab < 0.5)["TRUE"])*100



#Créer les tableaux d'aires enrichies en C13 et non enrichies 
non_enr_means = apply(enr_means,2, function(x) 100 -x) 

area_marked_only = area_means*(enr_means/100)
area_non_marked_only = area_means*(non_enr_means/100)

samples = c('A', 'B', 'AB')
states = c('marked_only', 'non_marked', 'total')

for(n in names){
  print(n)
  if(length(grep('.*_A',n)) !=0){
    all_data[grep(paste0(n,'[^B]'), all_data$sample),]$corrected_area <-
      unlist(sapply(all_data[grep(paste0(n,'[^B]'), all_data$sample),]$corrected_area, 
                    function(x) x/cellnumber[n], simplify = "array"))
    all_data$corrected_area[all_data[grep(n, all_data$sample),]$corrected_area == 0] <- 1E-5
    area_means[grep(n, colnames(area_means))][area_means[grep(n, colnames(area_means))]== 0] <- 1E-5
    area_marked_only[grep(n, colnames(area_marked_only))][area_marked_only[grep(n, colnames(area_marked_only))]== 0] <- 1E-5
    area_non_marked_only[grep(n, colnames(area_non_marked_only))][area_non_marked_only[grep(n, colnames(area_non_marked_only))]== 0] <- 1E-5
    
  }
  else{
    all_data[grep(n, all_data$sample),]$corrected_area <-
      unlist(sapply(all_data[grep(n, all_data$sample),]$corrected_area, 
                    function(x) x/cellnumber[n], simplify = "array"))
    all_data$corrected_area[all_data[grep(n, all_data$sample),]$corrected_area == 0] <- 1E-5
    area_means[grep(n, colnames(area_means))][area_means[grep(n, colnames(area_means))]== 0] <- 1E-5
    area_marked_only[grep(n, colnames(area_marked_only))][area_marked_only[grep(n, colnames(area_marked_only))]== 0] <- 1E-5
    area_non_marked_only[grep(n, colnames(area_non_marked_only))][area_non_marked_only[grep(n, colnames(area_non_marked_only))]== 0] <- 1E-5
    
  }
  
} 


#Créer les tableau de conditions contrôle
control_tot = area_means[grep('.*Cont.*', names(area_means))]
control_marked_only = area_marked_only[grep('.*Cont.*', names(area_means))]
control_non_marked_only = area_non_marked_only[grep('.*Cont.*', names(area_means))]

#Fonction permettant d'enlever les valeurs infinies 
noInf <- function(z){
  z <- do.call(data.frame,
               lapply(z,
                      function(x) replace(x, is.infinite(x), NA)))
  return(z)
}



#jeu test
# AB = area_marked_only[grep('.*AB.*', names(area_means))]
# FC_AB = log2(AB/control_marked_only)
# FC_AB <- noInf(FC_AB)
# rownames(FC_AB) = row.names(AB)
# pheatmap(na.omit(FC_AB), legend_labels = row.names(FC_AB), cluster_cols=F)



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
  
  pat = paste0('.*_',s , '_.*')
  # Calcul des FC totaux

  tmp_tot = assign(paste0('FC_tot_', s), (area_means[grep(pat, names(area_means))]/control_tot))
  l_tot = assign(paste0('logFC_tot_',s), log2(tmp_tot))
  rownames(l_tot) = rownames(area)
  rownames(tmp_tot) = rownames(area)
  rownames(l_tot) = rownames(area)
  FC_list[paste0('FC_tot_', s)] <- list(tmp_tot)
  lFC_list[paste0('logFC_tot_',s)] <- list(l_tot)
  
  # FC aires non marquées
  tmp_nm = assign(paste0('FC_notmarked_', s), (area_non_marked_only[grep(
  pat, names(area_means))]/control_non_marked_only))
  l_nm = assign(paste0('logFC_nm_',s), log2(tmp_nm))
  rownames(l_nm) = rownames(area)
  rownames(tmp_nm) = rownames(area)
  FC_list[paste0('FC_notmarked_', s)] <- list(tmp_nm)
  lFC_list[paste0('logFC_nm_',s)] <- list(l_nm)
  
  #FC aires marquées 

  tmp_m = assign(paste0('FC_marked_', s), 
                 (area_marked_only[grep(pat, names(area_means))]/(control_marked_only)))
  l_m = assign(paste0('logFC_m_',s), log2(tmp_m))
  rownames(l_m) = rownames(area)
  rownames(tmp_m) = rownames(area)
  FC_list[paste0('FC_marked_', s)] <- list(tmp_m)
  lFC_list[paste0('logFC_m_',s)] <- list(l_m)
  
  #Création des heatmaps 
  h_t = pheatmap(na.omit(l_tot), legend_labels = row.names(tmp_tot), cluster_cols=F,
                 main = paste0(s,'_total'), breaks = scale, color = palette,
                 cluster_rows = F)
  
  h_nm = pheatmap(na.omit(l_nm), legend_labels = row.names(tmp_nm), cluster_cols=F,
                  main = paste0(s, '_non_marked'), breaks = scale, color = palette,
                  cluster_rows = F)
  
  h_m = pheatmap(na.omit(l_m), legend_labels = row.names(tmp_m), cluster_cols=F, 
           main = paste0(s, '_marked'), breaks = scale, color = palette,
           cluster_rows = F)
  #intégration en une liste d'objets 
  plot_list[['a']] = h_t[[4]]
  plot_list[['b']] = h_nm[[4]]
  plot_list[['c']] = h_m[[4]]
  
  #Sortie graphique
  g<-do.call(grid.arrange,plot_list )
  ggsave(paste0(s, '_heatmaps.pdf'),g, height = 35, width = 15, device = "pdf")
  plot_list = list()
}


  metabolites = rownames(area)

point = 0
table_list = list()

# Boucle de construction des tableaux pour les calculs de p-values (t-test)
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


times = c('T0', 'T24', 'T48')
control = c('T0_Cont', 'T24_Cont', 'T48_Cont')


# Fonction de calcul de p-values 
p_val_calculator <- function(dataframe, control_dat){
  point <- 0
  flag = 0
  times = c('T0', 'T24', 'T48')
  pval_list = c()
  pval_matrix = data.frame(row.names = times)
  for(m in metabolites){
    flag = flag +1
    for( t in times){
      point = point +1
      tut = dataframe[m == dataframe$metabolites,
      grep(paste0('.*', t, '.*'), colnames(dataframe))]
      cont = control_dat[m == control_dat$metabolites,
                            grep(paste0('.*', t, '.*'), colnames(control_dat))]
      agglo = append(cont, unlist(tut))
      std = sd(agglo)
      if(length(log2(tut) == Inf) == 3 | length(log2(cont)) == 3){
        ttst =  wilcox.test(tut/std, cont/std, alternative = "two.sided")
      }
      else{
      ttst =  wilcox.test(noInf(log2(tut/std)), noInf(log2(cont/std)), alternative = "two.sided")
      }
      pval = p.adjust(ttst$p.value, "hochberg")
      pval_list = append(pval_list, pval)
    }
    pval_matrix[m] = pval_list
    pval_list = c()
    
  }
  return(pval_matrix)
}

# Boucle de calcul des p-values

iter = 1
for(a in names(table_list)){
  # Evite les calculs des contrôles contre eux même
  if(length(grep(pattern = '.*Cont.*', names(table_list[[a]][2]))) != 0){
    print('control')
  }
  else{
    if(iter == 1){
    assign(paste0('pvalue_',a),p_val_calculator(table_list[[a]] , 
                                                table_list[[1]]))}
    if(iter == 2){
      assign(paste0('pvalue_',a), p_val_calculator(table_list[[a]] , 
                                                 table_list[[2]]))}
    if(iter == 3){
      assign(paste0('pvalue_',a),p_val_calculator(table_list[[a]] , 
                                                  table_list[[3]]))
      iter <- 0}
    iter = iter +1
  }
  
} 
library(ggrepel)
#volcano 
volc = data.frame(cbind((FC_marked_AB$T48_AB_mean),unlist(pvalue_m_a_r_T48_AB_m_only['T48',])))

ggplot(data = volc , aes(x=log2(X1), y=-log10(X2), label = rownames(volc))) +
  geom_point() +
    geom_point(x= log2(volc$X1), y = -log10(volc$X2)) +
      xlim(c(-4,8)) +
        geom_text(label = rownames(volc))

# Calculs de Z-scores + p-values :


pvals = data.frame(row.names = rownames(data.frame(FC_list[1])))
for(l in names(FC_list)){
  tmp = data.frame(FC_list[l])
  for(i in colnames(tmp)){
    # Plot des densités
    pvalues_list = c()
    actlist = unlist(tmp[i])
    actlist[is.na(actlist)] = 0
    names(actlist) <- rownames(tmp)
    densl2 = density(na.omit(log2(actlist)))
    dens = density(na.omit(actlist))
    plot(dens, main = paste0(l, '_',i))
    zscores = (actlist - mean(actlist))/sd(actlist)
    names(zscores) <- rownames(tmp)
    for(v in 0:length(actlist)){
      r_pval = p.adjust(pnorm(q = zscores[v], mean = mean(zscores), sd = sd(zscores),
                     lower.tail = FALSE ), "bonferroni")
      l_pval = p.adjust(pnorm(q = zscores[v], mean = mean(zscores), sd = sd(zscores),
                              lower.tail = TRUE ), "bonferroni")
      if(min(c(r_pval, l_pval))!= Inf){
      pvalues_list <- c(pvalues_list, min(c(r_pval,l_pval))) 
      }
    }

    pvals[paste0(i)] =  pvalues_list
  }
}




all_FC = data.frame(row.names = rownames(data.frame(FC_list[l])))
for(l in names(FC_list)){
  for(i in names(FC_list[l])){
    all_FC[paste0(l,'_',i)] = FC_list[l][i]
  }
}

for(i in ls()){
  if(length(grep("control_", i) !=0 )){
    co = get(i)
    for(c in control){
      pat = paste0('.*_',c , '_.*')
      dv_co = cbind(area_means[grep(pat, names(area_means))], control_tot)
      dv_tco['std'] = apply(dv_tot, 1, function(x)  sd(x))
      assign(paste0('FC_',i,'0vs48h'), co$T48_Cont_mean)
    }
  }
}










# a = subset(T48_AB_group, T48_AB_group$Group.1 == 'T48_AB_1')
# b = subset(T48_AB_group, T48_AB_group$Group.1 == 'T48_AB_2')
# c = subset(T48_AB_group, T48_AB_group$Group.1 == 'T48_AB_3')
# 
# a2 = subset(T48_Cont_group, T48_Cont_group$Group.1 == 'T48_Cont_1')
# b2 = subset(T48_Cont_group, T48_Cont_group$Group.1 == 'T48_Cont_2')
# c2 = subset(T48_Cont_group, T48_Cont_group$Group.1 == 'T48_Cont_3')
# 
# spec_df = data.frame(row.names = a$Group.2, 'T48_AB_1' = a$x, 'T48_AB_2' = b$x, 'T48_AB_2' = c$x,
#                      'T48_Cont_group_1' = a2$x, 'T48_Cont_group_2' = b2$x,'T48_Cont_group_3' = c2$x )
# 
# 
# owo = compare/compare$std 
# apply(owo[c(1,2,3)], 1, function(x) mean(x))/apply(owo[c(4,5,6)], 1, function(x) mean(x))


# 
# m_aFC = log2(all_FC[,grep('_marked_', colnames(all_FC))])
# nm_aFC = log2(all_FC[,grep('_notmarked_', colnames(all_FC))])
# tot_aFC = log2(all_FC[,grep('_tot_', colnames(all_FC))])
# write.csv(x = all_FC, file = "/run/media/aurelien/ACOFFE/Stage/essais/all_FC_V2.csv")
# 
# write.csv(x = m_aFC, file = "/run/media/aurelien/ACOFFE/Stage/essais/all_FC_m.csv")
# write.csv(x = nm_aFC, file = "/run/media/aurelien/ACOFFE/Stage/essais/all_FC_nm.csv")
# write.csv(x = tot_aFC, file = "/run/media/aurelien/ACOFFE/Stage/essais/all_FC_tot.csv")
# 
# log2_all_FC = apply(all_FC, MARGIN = 2, FUN = function(x) log2(x))
# write.csv(x = log2_all_FC, file = "/run/media/aurelien/ACOFFE/Stage/essais/log_all_FC.csv")
# all_FC_log2 = log2(all_FC)
# write.csv(all_FC_log2, file = "/run/media/aurelien/ACOFFE/Stage/essais/all_FC_log2.csv")

# #
# a = subset(T48_AB_group, T48_AB_group$Group.1 == 'T48_AB_1')
# b = subset(T48_AB_group, T48_AB_group$Group.1 == 'T48_AB_2')
# c = subset(T48_AB_group, T48_AB_group$Group.1 == 'T48_AB_3')
# 
# a2 = subset(T48_Cont_group, T48_Cont_group$Group.1 == 'T48_Cont_1')
# b2 = subset(T48_Cont_group, T48_Cont_group$Group.1 == 'T48_Cont_2')
# c2 = subset(T48_Cont_group, T48_Cont_group$Group.1 == 'T48_Cont_3')
# 
# spec_df = data.frame(row.names = a$Group.2, 'T48_AB_1' = a$x, 'T48_AB_2' = b$x, 'T48_AB_2' = c$x,
#                      'T48_Cont_group_1' = a2$x, 'T48_Cont_group_2' = b2$x,'T48_Cont_group_3' = c2$x )
# 
# 
# 
# write.csv(x = spec_df, file = "/run/media/aurelien/ACOFFE/Stage/essais/samples2cond.csv" )


# write.csv(pvals, "/run/media/aurelien/ACOFFE/Stage/essais/normix_pval.csv")
# write.csv(x = T48_AB, "/run/media/aurelien/ACOFFE/Stage/essais/T48_AB.csv")
# write.csv(x = T48_Cont , "/run/media/aurelien/ACOFFE/Stage/essais/T48_Cont.csv" )





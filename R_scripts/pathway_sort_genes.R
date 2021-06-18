library(readr)
setwd(dir = "/run/media/aurelien/ACOFFE/Stage/integration_job/")



amino_acids = read_csv("ensembl_gene_list/alanine_aspartate_glutamate_metabolism.csv")
# aspartate = read_csv("ensembl_gene_list/aspartate.csv")
# alanine = read_csv("ensembl_gene_list/alanine.csv")
# glutamate = read_csv("ensembl_gene_list/glutamate.csv")
# glycolyse = read_csv("ensembl_gene_list/glycolysis.csv")
# pentose_phosphate <- read_csv("ensembl_gene_list/pentose_phosphate.csv")
# purine = read_csv("ensembl_gene_list/Purine.csv")
# pyrimidine = read_csv("ensembl_gene_list/Pyrimidine.csv")
# pyruvate = read_csv("ensembl_gene_list/Pyruvate.csv")
# tca = read_csv("ensembl_gene_list/TCA.csv")

sublist = ls()


dataset = (read.csv("new_result_part3/new_result_part3/control_normoxie_ldha_ldhb_normoxie/table/DEG_0_05.tsv", 
                                         header = TRUE, sep = '\t'))

for(t in sublist){
  tmp = get(t)
  sub = assign(paste0("sublist_",t), subset(dataset, 
                                      substr(dataset$external_gene_name, 
                                             1,nchar(dataset$external_gene_name)) 
                                      %in% tmp$name))
  write.csv(sub, paste0("ensembl_gene_list/","sub_ContvsAB_normo_",t,".csv"))
}

# substr(dataset$ensembl_gene_id_version, 1,nchar(dataset$ensembl_gene_id_version)-3)

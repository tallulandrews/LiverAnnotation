require("biomaRt")
ensembl = useMart("ensembl", mirror='useast')
ensg = useDataset("hsapiens_gene_ensembl", ensembl)
musg = useDataset("mmusculus_gene_ensembl", ensembl)
rat = useDataset("rnorvegicus_gene_ensembl", ensembl)

attributes = listAttributes(musg)


ensg2musg <- getBM(attributes=c('ensembl_gene_id', "mmusculus_homolog_ensembl_gene"), mart = ensg)
ensg2musg <- ensg2musg[ensg2musg[,2] != "",]

ensg2rat <- getBM(attributes=c('ensembl_gene_id', "rnorvegicus_homolog_ensembl_gene"), mart = ensg)
ensg2rat <- ensg2rat[ensg2rat[,2] != "",]

musg2rat <- getBM(attributes=c('ensembl_gene_id', "rnorvegicus_homolog_ensembl_gene"), mart = musg)
musg2rat <- musg2rat[musg2rat[,2] != "",]

ensg2symbol <- getBM(attributes=c('ensembl_gene_id', "hgnc_symbol"), mart = ensg)
ensg2symbol <- ensg2symbol[ensg2symbol[,2] != "",]


musg2symbol <- getBM(attributes=c('ensembl_gene_id', "mgi_symbol"), mart = musg)
musg2symbol <- musg2symbol[musg2symbol[,2] != "",]

rat2symbol <- getBM(attributes=c('ensembl_gene_id', "rgd_symbol"), mart = rat)
rat2symbol <- rat2symbol[rat2symbol[,2] != "",]



saveRDS(list(ensg2musg=ensg2musg, ensg2symbol=ensg2symbol, musg2symbol=musg2symbol, ensg2rat=ensg2rat, rat2symbol=rat2symbol, musg2rat=musg2rat), "Ensembl_name_maps.rds")

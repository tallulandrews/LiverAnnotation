### Currently only supports Human Data
# Automatically generates various plots to help with manual annotation of clusters.
# Includes: 
#  Dotplot / heatmaps of key marker genes
#  TCR and BCR constant chain expression
#  Correlation of hepatocytes with Halpern spatial stratification reference
#  Correlation of T-cells with Zheng sorted reference

option_list <- list(
    make_option(c("-i", "--input_rds"),
        help="RDS of a normalized, scaled, and clustered Seurat Object"),
    make_option(c("-o", "--out_prefix"), default="Manual_anno",
        help="prefix for output files [default %default]"),
    make_option(c("--organism"), default="Human",
        help="Organism, currently supported: [Human, Mouse, Rat]"),
    make_option(c("-c", "--clusters"), default="seurat_clusters",
        help="Name of the metadata column containing the clusters to use."),
    make_option(c("-a", "--autoannotation"), default="marker_labs",
        help="Name of the autoannotation labels to use."),
    )

OPTS <- parse_args(OptionParser(option_list=option_list))

OPTS$org <- "Hsap"
if (OPTS$organism == "Mouse") {
	OPTS$org <- "Mmus"
}
if (OPTS$organism == "Rat") {
	OPTS$org <- "Rat"
}

print(OPTS)


# Read in prereqs & reference data
require("Seurat")
require("scmap")
require("SingleCellExperiment")
require(proxy)
require(gplots)
require(ggplot2)
require(RColorBrewer)
source("~/scripts/LiverMap2.0/Colour_Scheme.R")
source("~/Annotation_Package/Generic_functions.R")
source("~/R-Scripts/Ensembl_Stuff.R")
Immune_Profiles <- readRDS("~/Annotation_Package/Zheng_immune_profiles.rds");
Hepatocyte_Human_spatial_profiles <- readRDS("~/Annotation_Package/Halpern_Human_Ortho_sig_profiles.rds");
Hepatocyte_Mouse_spatial_profiles <- readRDS("~/Annotation_Package/Halpern_Mouse_sig_profiles.rds");
soupX <- read.table("~/Annotation_Package/SoupXGenesets_markers.csv", sep=",", header=T)
markers <- readRDS("~/Annotation_Package/markers.rds")

# Read input 
obj <- readRDS(OPTS$input_rds)

clusters <- obj@meta.data[,OPTS$clusters]

# For subsetting & specific analyses.
autoanno <- obj@meta.data[,OPTS$autoannotation]

tab <- table(clusters, autoanno); tab <- tab/rowSums(tab);
auto_anno <- colnames(tab)[apply(tab, 1, function(x){which(x==max(x))})]

coarse_anno <- as.character(auto_anno)
coarse_anno[grep("Hep", coarse_anno)] <- "Hep"
coarse_anno[grep("NK*cell", coarse_anno)] <- "Lympho"
coarse_anno[grep("T*cell", coarse_anno)] <- "Lympho"
coarse_anno[grep("B*cell", coarse_anno)] <- "Lympho"

Lympho_clusters <- which(coarse_anno == "Lympho")
Hepato_clusters <- which(coarse_anno == "Hep")
ambiguous <- which(coarse_anno == "ambiguous")
Lympho_clusters <- c(ambiguous, Lympho_clusters)
Hepato_clusters <- c(ambiguous, Hepato_clusters)

gene_cluster_means <- group_rowmeans(obj@assays[[1]]@data, clusters)
scale_cluster_means <- group_rowmeans(obj@assays[[1]]@scale.data, clusters)

# Make mouse-gene named expression matrix;
mouse_genes1 <- General_Map(rownames(scale_cluster_means), in.org="Rat", out.org="Mmus", in.name="symbol", out.name="symbol")
mouse_genes2 <- General_Map(rownames(scale_cluster_means), in.org="Hsap", out.org="Mmus", in.name="symbol", out.name="symbol")
mouse_genes <- rownames(scale_cluster_means);
mouse_genes[!(mouse_genes1 == "")] <- mouse_genes1[!(mouse_genes1 == "")]
mouse_genes[!(mouse_genes2 == "")] <- mouse_genes2[!(mouse_genes2 == "")]

tmp <- Rename_Rows(scale_cluster_means, mouse_genes)
scale_cluster_means_mouse <- tmp$data_matrix;

# Hepatocyte Annotation
sync_profiles <- function(to_anno_means, ref_profiles, anno.org=c("Hsap", "Mmus", "Rat"), ref.org=c("Hsap", "Mmus", "Rat")) {
	if (FALSE) {
		remap <- rownames(to_anno_means);
		remapped <- General_Map(remap, in.org=anno.org, in.name="symbol", out.org=ref.org, out.name="symbol")
		keep <- remapped != "" & !duplicated(remapped);	
		to_anno_means <- to_anno_means[keep,]
		remapped <- remapped[keep]
		rownames(to_anno_means) <- remapped
	}
	common_genes <- sort(rownames(to_anno_means)[rownames(to_anno_means) %in% rownames(ref_profiles)])
	to_anno_means <- to_anno_means[match(common_genes, rownames(to_anno_means)),]
	ref_profiles <- ref_profiles[match(common_genes, rownames(ref_profiles)),]

	tmp <- colnames(to_anno_means)
	to_anno_means <- t(apply(to_anno_means, 1, scale));
	colnames(to_anno_means) <- tmp

	tmp <- colnames(ref_profiles)
	ref_profiles <- t(apply(ref_profiles, 1, scale));
	colnames(ref_profiles) <- tmp
	return(list(anno=to_anno_means, ref=ref_profiles));
}


if (OPTS$org == "Hsap") {
	synced <- sync_profiles(scale_cluster_means[,Hepato_clusters], Hepatocyte_Human_spatial_profiles, anno.org=OPTS$org, ref.org="Hsap");
} else {
	synced <- sync_profiles(scale_cluster_means_mouse[,Hepato_clusters], Hepatocyte_Mouse_spatial_profiles, anno.org=OPTS$org, ref.org="Mmus");
}
	
cluster_profiles <- synced$anno
hep_ref <- synced$ref;

cors <- simil(t(cluster_profiles), t(hep_ref))
cos <- simil(t(cluster_profiles), t(hep_ref), method="cosine")

heat_col <- rev(c(rev(brewer.pal(4, "Reds")), "white", brewer.pal(4,"Blues")))

png(paste(OPTS$out_prefix, "Hep_vs_Halpern_cosine.png", sep="_"), width=8, height=8, units="in", res=300)
heatmap.2(cos, trace="none", Rowv=TRUE, Colv=FALSE, dendrogram=c("none"), symm=TRUE, scale="none", symbreaks=TRUE, col=heat_col)
dev.off()
png(paste(OPTS$out_prefix, "Hep_vs_Halpern_cors.png", sep="_"), width=8, height=8, units="in", res=300)
heatmap.2(cors, trace="none", Rowv=TRUE, Colv=FALSE, dendrogram=c("none"), symm=TRUE, scale="none", symbreaks=TRUE, col=heat_col)
dev.off()


# Immune Annotation

immune_ref <- Immune_Profiles$profiles[Immune_Profiles$specificity > 5,]

if (length(Lympho_clusters) > 1) {
	synced <- sync_profiles(scale_cluster_means[,Lympho_clusters], immune_ref, anno.org=OPTS$org, ref.org="Hsap")
	require(proxy)
	cors <- simil(t(synced$anno), t(synced$ref))
	cos <- simil(t(synced$anno), t(synced$ref), method="cosine")

	require(gplots)
	require(RColorBrewer)
	heat_col <- rev(c(rev(brewer.pal(4, "Reds")), "white", brewer.pal(4,"Blues")))

	png(paste(OPTS$out_prefix, "Immune_vs_Zheng_cosine.png", sep="_"), width=8, height=8, units="in", res=300)
	heatmap.2(cos, trace="none", Rowv=TRUE, Colv=FALSE, dendrogram=c("none"), symm=TRUE, scale="none", symbreaks=TRUE, col=heat_col)
	dev.off()
	png(paste(OPTS$out_prefix, "Immune_vs_Zheng_cors.png", sep="_"), width=8, height=8, units="in", res=300)
	heatmap.2(cors, trace="none", Rowv=TRUE, Colv=FALSE, dendrogram=c("none"), symm=TRUE, scale="none", symbreaks=TRUE, col=heat_col)
	dev.off()
}

#Immune + scamp

tmp <- colnames(immune_ref);
immune_scmap <- t(apply(immune_ref, 1, scale))
colnames(immune_scmap) <- tmp;
	

immune_cells <- obj[,clusters %in% Lympho_clusters]
immune_cells_sce <- as.SingleCellExperiment(immune_cells)
immune_cells_sce <- immune_cells_sce[rownames(immune_cells_sce) %in% rownames(immune_cells@assays[[1]]@scale.data),]

immune_cells_sce <- immune_cells_sce[match(rownames(immune_cells@assays[[1]]@scale.data), rownames(immune_cells_sce)),]

assays(immune_cells_sce)[["logcounts"]] <- immune_cells@assays[[1]]@scale.data;
rowData(immune_cells_sce)$feature_symbol <- rownames(immune_cells@assays[[1]]@scale.data);
if (OPTS$org != "Hsap") {
#	rowData(immune_cells_sce)$feature_symbol <- General_Map(rownames(immune_cells@assays[[1]]@scale.data), anno.org=OPTS$org, out.org="Hsap", in.name="symbol", out.name="symbol");
#	immune_cells_sce <- immune_cells_sce[!duplicated(rowData(immune_cells_sce)$feature_symbol),]
}

cell_level <- scmapCluster(projection=immune_cells_sce, index_list = list(zheng=immune_scmap), threshold=0)

cell_immune_lab <- cell_level$combined_labs;
names(cell_immune_lab) <- colnames(immune_cells_sce)

obj@meta.data$Immune_scmap <- cell_immune_lab[match(rownames(obj@meta.data), names(cell_immune_lab))]
obj@meta.data$Immune_scmap[is.na(obj@meta.data$Immune_scmap)] <- "Non-Lympho"

png(paste(OPTS$out_prefix, "scmap_immune.png", sep="_"), width=9, height=7, units="in", res=300)
DimPlot(obj, reduction="umap", group.by="Immune_scmap", label=TRUE, pt.size=0.5)
dev.off()



#Immune marker gene plots.
name_col <- "Hgene"
#if (OPTS$org == "Mmus") {
#	name_col <- "Mgene"
#}
#if (OPTS$org == "Rat") {
#	name_col <- "Rgene"
#}

zheng_immune_markers=markers$immune_other[markers$immune_other[,name_col] %in% rownames(obj),]
zheng_empiric_markers=markers$immune_zheng[markers$immune_zheng[,name_col] %in% rownames(obj),]

pdf(paste(OPTS$out_prefix, "zheng_markers.pdf", sep="_"), width=20, height=8)
Seurat::DotPlot(obj, features=zheng_immune_markers[,name_col], group.by=OPTS$clusters)+ theme(axis.text.x = element_text(angle = 90))
dev.off()
pdf(paste(OPTS$out_prefix, "zheng_emp_markers.pdf", sep="_"), width=20, height=8)
Seurat::DotPlot(obj,features=zheng_empiric_markers[,name_col], group.by=OPTS$clusters)+ theme(axis.text.x = element_text(angle = 90))
dev.off()


# Other marker gene plots



cytokines <- c(rownames(obj)[grep("^CXC",rownames(obj), ignore.case=TRUE)],  rownames(obj)[grep("^IL",rownames(obj), ignore.case=TRUE)])
TCR_BCR <- c("TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2", "IGKC", "IGHE", "IGHM", "IGLC1", "IGLC3", "IGLC2", "JCHAIN", "IGKV1-12", "IGKV4-1", "IGLV3-1", "IGLV6-57", "IGLL5", "IGLC7")
Other_Surface_Signal <- c("TNF", rownames(obj)[grep("^HLA",rownames(obj), ignore.case=TRUE)], "PTPRC", "CD8A", "CD3E", "MARCO", "CD79B", "NK67", "GNLY")

#immune_phenotype <- c(as.character(soupX[,2]),"TNF", "JCHAIN", "IGKV1-12", "IGKV4-1", "IGLV3-1", "IGLV6-57", "IGLL5", "IGLC7", rownames(obj)[grep("^CXC",rownames(obj), ignore.case=TRUE)], rownames(obj)[grep("^IL",rownames(obj), ignore.case=TRUE)], rownames(obj)[grep("^HLA",rownames(obj), ignore.case=TRUE)])

if (OPTS$org != "Hsap") {
#	immune_phenotype2 <- General_Map(immune_phenotype, in.org="Hsap", out.org=OPTS$org, in.name="symbol", out.name="symbol");
#	immune_phenotype2[immune_pheontype2 == ""] <- immune_pheontype[immune_pheontype2 == ""]
#	immune_pheontype <- immune_phenotype2
}

pdf(paste(OPTS$out_prefix, "Immune_pheno_cytokines.pdf", sep="_"), width=20, height=8)
Seurat::DotPlot(immune_cells, features=cytokines, group.by=OPTS$clusters)+ theme(axis.text.x = element_text(angle = 90))
dev.off()

pdf(paste(OPTS$out_prefix, "Immune_pheno_TCR_BCR.pdf", sep="_"), width=15, height=8)
Seurat::DotPlot(immune_cells, features=TCR_BCR, group.by=OPTS$clusters)+ theme(axis.text.x = element_text(angle = 90))
dev.off()

pdf(paste(OPTS$out_prefix, "Immune_pheno_Other.pdf", sep="_"), width=15, height=8)
Seurat::DotPlot(immune_cells, features=Other_Surface_Signal, group.by=OPTS$clusters)+ theme(axis.text.x = element_text(angle = 90))
dev.off()

#Heatmap
markers$best <- markers$best[markers$best[,name_col] %in% rownames(scale_cluster_means),]
row_colours <- Cell_type_colours[match(map_cell_types(markers$best[,1]), Cell_type_colours[,1]),2]

heatdata <- scale_cluster_means[match(markers$best[,name_col], rownames(scale_cluster_means)),]

png(paste(OPTS$out_prefix, "best_markers.png", sep="_"), width=8, height=8, units="in", res=300)
heatmap.2(heatdata, RowSideColors=row_colours, scale="row", trace="none")
dev.off()

out <- as.matrix(table(clusters, obj$Phase))
out <- cbind(out, auto_anno)

write.table(out, row.names=T, col.names=T, file=paste(OPTS$out_prefix, "CCPhase.txt", sep="_"), sep=",", quote=FALSE)


################################################################
##                         color                              ##
################################################################
{
	library(RColorBrewer)
	color1=paste0(colorRampPalette(brewer.pal(n=8,name="Greens"))(25),"3F")
	color2=paste0(colorRampPalette(brewer.pal(n=8,name="Blues"))(8),"3F")
	color3=paste0(colorRampPalette(brewer.pal(n=8,name="Reds"))(5),"3F")
	color4=paste0(colorRampPalette(brewer.pal(n=8,name="Oranges"))(20),"B2")
	samplecolor=c(
		"Normal1"=color1[3],"Normal2"=color1[4],"Normal3"=color1[5],"Normal4"=color1[6],"Normal5"=color1[7],"Normal6"=color1[8],
		"Normal7"=color1[9],"Normal8"=color1[10],"Normal9"=color1[11],"Normal10"=color1[12],"Normal11"=color1[13],"Normal12"=color1[14],
		"Normal13"=color1[15],"Normal14"=color1[16],"Normal15"=color1[17],"Normal16"=color1[18],"Normal17"=color1[19],"Normal18"=color1[20],
		"Normal19"=color1[21],"Normal20"=color1[22],"Normal21"=color1[23],"Normal22"=color1[24],"Normal23"=color1[25],"AML1"=color2[3],
		"AML2"=color2[4],"AML3"=color2[5],"AML4"=color2[6],"AML5"=color2[7],"AML6"=color2[8],"TALL1"=color3[3],"TALL2"=color3[4],"TALL3"=color3[5],
		"TMMPAL1"=color4[3],"TMMPAL2"=color4[4],"TMMPAL3"=color4[5],"TMMPAL4"=color4[6],"TMMPAL5"=color4[7],"TMMPAL6"=color4[8],"TMMPAL7"=color4[9],"TMMPAL8"=color4[10]
	)

	color1 <- c("#377EB8", "#E41A1C", "#006837", "#FF7F00", "#984EA3", 
		"#FFD92F", "#46A040", "#F781BF", "#A65628", "#80B1D3", "#FB8072",
		"#999999", "#00AF99", "#b4a996", "#F0027F", 'turquoise4', 
		"palevioletred4","royalblue3", "#B3DE69", "#CAB2D6", "lightgreen",
		"#FFFF99", "steelblue1", "orangered", "#ffc20e", "#001588", 
		"#F6313E",'turquoise2', "#490C65", "#FFFF33" )
	color2 = brewer.pal(n = 12, name = "Paired")
	color3 = brewer.pal(n = 8, name = "Dark2")
	color4=unique(c("#46A040" , "#FFC179", "#98D9E9", "#FF8247" ,"#A020F0" , "#F6313E", 
        "#00AF99", "#CD853F" ,"#8F1336" ,"#0081C9" ,"#001588" ,"#BEBEBE" ,"#BA7FD0",
        "darkgreen", "steelblue","limegreen" ,"gold1" , "lightsalmon", "steelblue1", 
        "violetred3" ,"orangered" ,"#D48AAA" , "red", "lightgreen", "sandybrown" ,"#B22222" ,
        "royalblue3" ,"royalblue4" ,"#8470FF" ,"plum2","lightseagreen", "turquoise4","lawngreen",
        "royalblue1","darkmagenta",'turquoise2','lightcoral','tan1','yellow','palevioletred4',
        'turquoise4','violetred2',"#A6CEE3","#1F78B4","#B2DF8A","#FB9A99","#E31A1C","#FDBF6F",
        "#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928","#1B9E77","#7570B3","#E6AB02","#A6761D",
        "#666666","#F0027F","#BF5B17","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FFFF33","#F781BF"))
	color5=unique(c("#46A040" , "#FFC179", "#98D9E9", "#FF5A00" ,"#AFAFAF" , 
        "#F6313E", "#00AF99", "#FFA300" ,"#8F1336" ,"#0081C9" ,"#001588" ,"#490C65" ,
        "#BA7FD0","darkgreen", "steelblue","limegreen" ,"gold1" , "lightsalmon", 
        "steelblue1", "violetred3" ,"orangered" ,"#D48AAA" , "red", "lightgreen",
        "sandybrown" ,"tomato4" ,"royalblue3" ,"royalblue4" ,"purple4" ,"plum2",
        "lightseagreen", "turquoise4","lawngreen","royalblue1","darkmagenta",'turquoise2',
        'lightcoral','tan1','yellow','palevioletred4','turquoise4','violetred2',"#A6CEE3",
        "#1F78B4","#B2DF8A","#FB9A99","#E31A1C","#FDBF6F", "#FF7F00","#CAB2D6","#6A3D9A",
        "#FFFF99","#B15928","#1B9E77","#7570B3","#E6AB02","#A6761D","#666666","#F0027F",
        "#BF5B17","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FFFF33","#F781BF"))

	celltype_color=c('B'="#8a8acb",'T'="#34a853",'CLP'="#99cc33",
		'Erythrocytes'="#8B1A1A",'CMP_GMP'="#7acef4",
		'Granulocyte'="#003265",'HSC'="#ffcc00",'Megakaryocytes'="#FF3030",
		'MEP'="#FF8C69",'Monocytes'="#0769ad",'MPP'="#f9a541",
		'NK'="#2f6f7e",'Other'="grey90","PlasmaCells"="#66336e")

	color2 <- c("0"="#377EB8", "1"="#E41A1C", "2"="#006837",
		"3"="#FF7F00", "4"="#984EA3", "5"="#FFD92F", 
		"6"="#46A040", "7"="#F781BF", "8"="#A65628",
		"9"="#80B1D3", "10"="#FB8072","11"="#999999",
		"12"="#00AF99", "13"="#FDB462", "14"="#F0027F",
		"15"='turquoise4',"16"="palevioletred4",
		"17"="royalblue3", "18"="#B3DE69",
		"19"="#CAB2D6", "20"="lightgreen",
		"21"="#FFFF99", "22"="steelblue1","23"="orangered",
		"24"="#ffc20e", "25"="#001588","26"="#F6313E")


	typecolor=c('TEC'="#ce181e",'S13'="#b4a996",
		'S14'="#F0027F",'S16'="palevioletred4",
		'S19'="#CAB2D6","S24"="#ffc20e",'AEC'="#007cc0",
		"TALL"="#ffaaaa","AML"="#b3dcff")

	fig2_color <- c("both"="#ffb900","myeloid_only"="#075aaa",
				"lymphoid_only"="#eb2226","neither"="#efe9e5")
}

################################################################
##            Preprocess: Integration all                     ##
################################################################
{
	# Integrate_AML # HMS10,HMS21,HMS24,HMS4,HMS6,HMS7
	file.AML<-list.files("/public/workspace/huangbin/MPAL/figure1/AML/Cellranger/")
	dat<-list()
	for(i in 1:length(file.AML)){
		tmp <- Read10X(data.dir = paste("/public/workspace/huangbin/MPAL/figure1/AML/Cellranger/",file.AML[i],sep=""))
		tmp <- CreateSeuratObject(counts = tmp, project = file.AML[i])
		tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^MT-") 
		tmp <- subset(tmp,features=rownames(tmp)[-grep("^MT|^RP",rownames(tmp))])
		tmp <- subset(tmp, subset = (nFeature_RNA > 200) & (nFeature_RNA < 4000) & (percent.mt < 30) & (nCount_RNA > 1000))  
		tmp <- NormalizeData(tmp)
		tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 2000)
		dat[[i]] <- tmp
	}
	anchors <- FindIntegrationAnchors(object.list = dat, dims = 1:30)
	dat <- IntegrateData(anchorset = anchors, dims = 1:30)
	saveRDS(dat,file="/public/workspace/huangbin/MPAL/figure1/Integrate/Integrate.aml.rds")

	# TALL # HMS11,HMS14,HMS26
	file.TALL<-list.files("/public/workspace/huangbin/MPAL/figure1/TALL/Cellranger/")
	dat <- list()
	for(i in 1:length(file.TALL)){
		tmp <- Read10X(data.dir = paste("/public/workspace/huangbin/MPAL/figure1/TALL/Cellranger/",file.TALL[i],sep=""))
		tmp <- CreateSeuratObject(counts = tmp, project = file.TALL[i])
		tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^MT-") 
		tmp <- subset(tmp,features=rownames(tmp)[-grep("^MT|^RP",rownames(tmp))]) 
		tmp <- subset(tmp, subset = (nFeature_RNA > 200) & (nFeature_RNA < 4000) & (percent.mt < 30) & (nCount_RNA > 1000)) 
		tmp <- NormalizeData(tmp)
		tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 2000)
		dat[[i]] <- tmp
	}
	anchors <- FindIntegrationAnchors(object.list = dat, dims = 1:30)
	dat <- IntegrateData(anchorset = anchors, dims = 1:30)
	saveRDS(dat,file="/public/workspace/huangbin/MPAL/figure1/Integrate/Integrate.tall.rds")

	# Integrate Normal
	# "GSM3396161" "GSM3396162" "GSM3396163" "GSM3396166" "GSM3396167" "GSM3396168" "GSM3396169" "GSM3396170" 
	# "GSM3396171" "GSM3396172" "GSM3396173" "GSM3396174" "GSM3396175" "GSM3396176" "GSM3396177" "GSM3396178" 
	# "GSM3396179" "GSM3396183" "GSM3396184" "GSM3396185" "HMS16"      "HMS17"      "HMS18"
	file.Normal<-list.files("/public/workspace/huangbin/MPAL/figure1/Normal/Cellranger/")
	dat <- list()
	for(i in 1:length(file.Normal)){
		tmp <- Read10X(data.dir = paste("/public/workspace/huangbin/MPAL/figure1/Normal/Cellranger/",file.Normal[i],sep=""))
		tmp <- CreateSeuratObject(counts = tmp, project = file.Normal[i])
		tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^MT-") 
		tmp <- subset(tmp,features=rownames(tmp)[-grep("^MT|^RP",rownames(tmp))]) 
		tmp <- subset(tmp, subset = (nFeature_RNA > 200) & (nFeature_RNA < 4000) & (percent.mt < 30) & (nCount_RNA > 1000)) 
		tmp <- NormalizeData(tmp)
		tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 2000)
		dat[[i]] <- tmp
	}
	anchors <- FindIntegrationAnchors(object.list = dat, dims = 1:30)
	dat <- IntegrateData(anchorset = anchors, dims = 1:30)
	saveRDS(dat,file="/public/workspace/huangbin/MPAL/figure1/Integrate/Integrate.normal.rds")

	# Integrate MPAL
	# "GSM5265322" "HMS20" "HMS27" "HMS9" "MPAL1_T1","MPAL1_T2","MPAL2_T1","MPAL3_T1","MPAL3_T2","MPAL5_T1"
	file.MPAL<-list.files("/public/workspace/huangbin/MPAL/figure1/MPAL/Cellranger/")
	dat<-list()
	for(i in 1:length(file.MPAL)){
		tmp <- Read10X(data.dir = paste("/public/workspace/huangbin/MPAL/figure1/MPAL/Cellranger/",file.MPAL[i],sep=""))
		tmp <- CreateSeuratObject(counts = tmp, project = file.MPAL[i])
		tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^MT-") 
		tmp <- subset(tmp,features=rownames(tmp)[-grep("^MT|^RP",rownames(tmp))])
		tmp <- subset(tmp, subset = (nFeature_RNA > 200) & (nFeature_RNA < 4000) & (percent.mt < 30) & (nCount_RNA > 1000))  
		tmp <- NormalizeData(tmp)
		tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 2000)
		dat[[i]] <- tmp
	}
	anchors <- FindIntegrationAnchors(object.list = dat, dims = 1:30)
	dat <- IntegrateData(anchorset = anchors, dims = 1:30)
	saveRDS(dat,file="/public/workspace/huangbin/MPAL/figure1/Integrate/Integrate.mpal.rds")

	# Integrate ALL
	library(Seurat)
	dat <- list()
	dat[[1]] <- readRDS("/public/workspace/huangbin/MPAL/figure1/Integrate/Integrate.mpal.rds")
	dat[[2]] <- readRDS("/public/workspace/huangbin/MPAL/figure1/Integrate/Integrate.tall.rds")
	dat[[3]] <- readRDS("/public/workspace/huangbin/MPAL/figure1/Integrate/Integrate.aml.rds")
	dat[[4]] <- readRDS("/public/workspace/huangbin/MPAL/figure1/Integrate/Integrate.normal.rds")
	anchors <- FindIntegrationAnchors(object.list = dat, dims = 1:30)
	dat <- IntegrateData(anchorset = anchors, dims = 1:30)
	saveRDS(dat,file="/public/workspace/huangbin/MPAL/figure1/Integrate/Integrate.all.rds")
}

################################################################
##                           Figure1                          ##
################################################################
{
	library(Seurat)
	library(SingleR)
	library(ggplot2)
	library(ggpubr)
	dat <- readRDS("/public/workspace/huangbin/MPAL/figure1/Integrate/Integrate.all.rds")
	dat <- ScaleData(dat, features = rownames(dat))
	dat <- RunPCA(dat, features = VariableFeatures(object = dat))
	dat <- FindNeighbors(dat, dims = 1:30)
	dat <- FindClusters(dat, resolution = 5)
	dat <- RunUMAP(dat, dims = 1:30)
	dat$sample <- NA
	dat$sample[which(dat$orig.ident=="GSM3396161")] <- "Normal1"
	dat$sample[which(dat$orig.ident=="GSM3396162")] <- "Normal2"
	dat$sample[which(dat$orig.ident=="GSM3396163")] <- "Normal3"
	dat$sample[which(dat$orig.ident=="GSM3396166")] <- "Normal4"
	dat$sample[which(dat$orig.ident=="GSM3396167")] <- "Normal5"
	dat$sample[which(dat$orig.ident=="GSM3396168")] <- "Normal6"
	dat$sample[which(dat$orig.ident=="GSM3396169")] <- "Normal7"
	dat$sample[which(dat$orig.ident=="GSM3396170")] <- "Normal8"
	dat$sample[which(dat$orig.ident=="GSM3396171")] <- "Normal9"
	dat$sample[which(dat$orig.ident=="GSM3396172")] <- "Normal10"
	dat$sample[which(dat$orig.ident=="GSM3396173")] <- "Normal11"
	dat$sample[which(dat$orig.ident=="GSM3396174")] <- "Normal12"
	dat$sample[which(dat$orig.ident=="GSM3396175")] <- "Normal13"
	dat$sample[which(dat$orig.ident=="GSM3396176")] <- "Normal14"
	dat$sample[which(dat$orig.ident=="GSM3396177")] <- "Normal15"
	dat$sample[which(dat$orig.ident=="GSM3396178")] <- "Normal16"
	dat$sample[which(dat$orig.ident=="GSM3396179")] <- "Normal17"
	dat$sample[which(dat$orig.ident=="GSM3396183")] <- "Normal18"
	dat$sample[which(dat$orig.ident=="GSM3396184")] <- "Normal19"
	dat$sample[which(dat$orig.ident=="GSM3396185")] <- "Normal20"
	dat$sample[which(dat$orig.ident=="HMS16")] <- "Normal21"
	dat$sample[which(dat$orig.ident=="HMS17")] <- "Normal22"
	dat$sample[which(dat$orig.ident=="HMS18")] <- "Normal23"
	dat$sample[which(dat$orig.ident=="HMS4")] <- "AML1"
	dat$sample[which(dat$orig.ident=="HMS6")] <- "AML2"
	dat$sample[which(dat$orig.ident=="HMS7")] <- "AML3"
	dat$sample[which(dat$orig.ident=="HMS10")] <- "AML4"
	dat$sample[which(dat$orig.ident=="HMS21")] <- "AML5"
	dat$sample[which(dat$orig.ident=="HMS24")] <- "AML6"
	dat$sample[which(dat$orig.ident=="HMS11")] <- "TALL1"
	dat$sample[which(dat$orig.ident=="HMS14")] <- "TALL2"
	dat$sample[which(dat$orig.ident=="HMS26")] <- "TALL3"
	dat$sample[which(dat$orig.ident=="MPAL1_T1")] <- "TMMPAL1"
	dat$sample[which(dat$orig.ident=="MPAL1_T2")] <- "TMMPAL1"
	dat$sample[which(dat$orig.ident=="MPAL2_T1")] <- "TMMPAL2"
	dat$sample[which(dat$orig.ident=="MPAL3_T1")] <- "TMMPAL3"
	dat$sample[which(dat$orig.ident=="MPAL3_T2")] <- "TMMPAL3"
	dat$sample[which(dat$orig.ident=="MPAL5_T1")] <- "TMMPAL4"
	dat$sample[which(dat$orig.ident=="HMS9")] <- "TMMPAL5"
	dat$sample[which(dat$orig.ident=="HMS20")] <- "TMMPAL6"
	dat$sample[which(dat$orig.ident=="HMS27")] <- "TMMPAL7"
	dat$sample[which(dat$orig.ident=="GSM5265322")] <- "TMMPAL8"
	# sampletype rename
	dat$sampletype <- NA
	dat$sampletype[which(dat$sample%in%paste0("Normal",1:23))] <- "Normal"
	dat$sampletype[which(dat$sample%in%paste0("AML",1:6))] <- "AML"
	dat$sampletype[which(dat$sample%in%paste0("TALL",1:3))] <- "TALL"
	dat$sampletype[which(dat$sample%in%paste0("TMMPAL",1:8))] <- "TMMPAL"
	dat$sampletype <- factor(dat$sampletype,levels=c("Normal","TMMPAL","TALL","AML"))

	# fig1A
	# 仅仅因为画图调顺序而改
	# dat$sample <- factor(dat$sample,levels=c(paste0("TMMPAL",1:8),
	# 	paste0("TALL",1:4),paste0("AML",1:6),paste0("Normal",1:23)))
	# Idents(dat) <- factor(dat$sample,levels=c(paste0("TMMPAL",1:8),
	# 	paste0("TALL",1:4),paste0("AML",1:6),paste0("Normal",1:23)))
	pdf("UMAP.Sample.pdf",useDingbats=F,width=10)
	DimPlot(dat,reduction="umap",pt.size=0.1,cols=samplecolor,
		group.by="sample",raster=FALSE)
	dev.off()

	# sfig1A
	pdf("UMAP.Sample.split.pdf",useDingbats=F,width=21)
	DimPlot(dat,reduction="umap",pt.size=0.1,cols=samplecolor,
		split.by="sampletype",group.by="sample",raster=TRUE)
	dev.off()

	# sfig1B
	pdf("UMAP.Clusters.pdf",useDingbats=F,width=12)
	DimPlot(dat,reduction="umap",pt.size=0.1,group.by="seurat_clusters")
	dev.off()

	# sfig2B
	# 计算簇之间的相关性
	library(Hmisc)
	library(RColorBrewer)
	library(pheatmap)
	tmp <- dat@assays$RNA@data
	dat <- NormalizeData(dat)
	dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 1000)
	rst <- sapply(names(table(dat$seurat_clusters)),function(x){
		mean <- apply(tmp[VariableFeatures(dat),rownames(dat@meta.data)[which(dat$seurat_clusters==x)]],1,function(y){
			mean(y,na.rm=T)
			})
		})
	cor <- rcorr(as.matrix(rst),type="spearman")
	pdf("cluster.cor.pdf")
	p <- pheatmap(as.matrix(cor$r),breaks=seq(0.5,1,0.01),
		col=rev(colorRampPalette(brewer.pal(11,"RdBu"))(length(seq(0.5,1,0.01)))),
		border_color=NA, cutree_cols=8, cutree_rows=8, fontsize = 4.5)
	print(p)
	dev.off()

	# sfig2A
	# 染标志物
	label <- c("CD34","VPREB1","IGLL1","CD3E","KLRD1","CD19","IGHG1",
		"MPO","AZU1","LYZ","PPBP","HBD","FCER1A")
	for(i in 1:length(label)){
		pdf(paste0("celltypemarker/",label[i],".pdf"))
		p <- FeaturePlot(dat,features=label[i],raster=FALSE)
		print(p)
		dev.off()
	}

	# fig1B
	# singleR
	DefaultAssay(dat) <- "RNA"
	ref<-BlueprintEncodeData()
	pred.main<-SingleR(test=dat@assays$RNA@data, ref=ref, labels=ref$label.main)
	pred.fine<-SingleR(test=dat@assays$RNA@data, ref=ref, labels=ref$label.fine)
	pred.main$pruned.labels[is.na(pred.main$pruned.labels)]<-"Other"
	pred.fine$pruned.labels[is.na(pred.fine$pruned.labels)]<-"Other"
	dat$blueprint.main<-pred.main$pruned.labels
	dat$blueprint.fine<-pred.fine$pruned.labels
	dat$celltype<-dat$blueprint.fine
	dat$celltype[dat$celltype%in%c("Adipocytes","DC","Fibroblasts","Myocytes","Endothelial cells","mv Endothelial cells","Chondrocytes")]<-"Other"
	dat$celltype[dat$celltype%in%c("CD4+ T-cells","CD4+ Tcm","CD4+ Tem","Tregs","CD8+ T-cells","CD8+ Tcm","CD8+ Tem")]<-"T"
	dat$celltype[dat$celltype%in%c("Class-switched memory B-cells","Memory B-cells","naive B-cells")]<-"B"
	dat$celltype[dat$celltype%in%c("Plasma cells")]<-"PlasmaCells"
	dat$celltype[dat$celltype%in%c("Eosinophils","Neutrophils")]<-"Granulocyte"
	dat$celltype[dat$celltype%in%c("Macrophages","Macrophages M1","Monocytes")]<-"Monocytes"
	dat$celltype[dat$celltype%in%c("NK cells")]<-"NK"
	dat$celltype[dat$celltype%in%c("CMP","GMP")]<-"CMP_GMP"
	dat$celltype <- as.factor(dat$celltype)
	dat$celltype <- factor(dat$celltype,levels=c("HSC","MPP","CLP","T","NK","B","PlasmaCells",
		"CMP_GMP","Monocytes","Granulocyte","MEP","Megakaryocytes","Erythrocytes","Other"))

	# fig1B
	# celltype umap
	pdf("UMAP.Celltype.pdf",useDingbats=F,width=9)
	DimPlot(dat,reduction="umap",pt.size=0.1,group.by="celltype",
		cols=celltype_color,raster=FALSE)
	dev.off()

	# fig1C
	# celltype.dot/marker
	library(future)
	library(Seurat)
	library(MAST)
	library(ggplot2)
	rst <- data.frame()
	plan("multiprocess",workers=10)
	dat <- readRDS("/public/workspace/huangbin/MPAL/figure1/Integrate/Integrate.all.rds")
	for(i in names(table(dat$celltype))){
		dat$group <- ifelse(dat$celltype==i,"case","control")
		dat$group <- as.factor(dat$group)
		Idents(dat) <- dat$group
		tmp <- FindMarkers(dat,
		  ident.1=c("case"),
		  ident.2=c("control"),
		  only.pos = TRUE,
		  min.pct = 0.25,
		  test.use="MAST",
		  logfc.threshold = 0.25)
		tmp$gene <- rownames(tmp)
		tmp$label <- i
		rst <- rbind(rst,tmp)
		rm(tmp)
		gc()
	}
	write.table(rst,"celltype.markers.txt",sep="\t",quote=F,row.names=F)
	top3 <- rst%>%group_by(label)%>%top_n(n=3,wt=avg_logFC)
	sortlabel <- rev(c("HSC","MPP","CLP","T","NK","B","PlasmaCells","CMP_GMP","Monocytes","Granulocyte","MEP","Megakaryocytes","Erythrocytes","Other"))
	Idents(dat) <- factor(dat$celltype,levels=sortlabel)
	pdf("Dot.celltype.pdf",width=15)
	DotPlot(dat,features=unique(as.character(sapply(sortlabel,function(x){top3$gene[which(top3$label==x)]}))),cols = c("lightgrey", "darkgreen"))+
	 theme(axis.text.x = element_text(size = 10,angle = 90, hjust = 1,color="black"))
	dev.off()

	# 定义白血病细胞
	# 生成每个cluster样本的比例
	label <- names(table(dat$seurat_clusters))
	rst <- data.frame()
	for(i in label){
		freq.tmp <- table(dat$sampletype[dat$seurat_clusters==i])/sum(table(dat$sampletype[dat$seurat_clusters==i]))
		tmp <- data.frame(Cluster=i,
			Freq.Normal=as.numeric(freq.tmp["Normal"]),
			Freq.TMMPAL=as.numeric(freq.tmp["TMMPAL"]),
			Freq.AML=as.numeric(freq.tmp["AML"]),
			Freq.TALL=as.numeric(freq.tmp["TALL"]))
		rst <- rbind(rst,tmp)
	}
	rst[is.na(rst)] <- 0
	# 根据白血病细胞比例变化选择稳定的cutoff
	freq <- sapply(seq(0.05,0.95,0.05),function(x){
			Cluster.case <- as.character(rst$Cluster[rst$Freq.Normal<x])
			dat$predicted.blasts <- "unknown"
			# Normal
			dat$predicted.blasts[which(dat$sampletype=="Normal")] <- "Normal"
			# TMMPAL
			dat$predicted.blasts[which(dat$seurat_clusters%in%Cluster.case & dat$sampletype=="TMMPAL" & !dat$celltype%in%c("Megakaryocytes","Erythrocytes","Other","B","PlasmaCells"))] <- "Leukemia_like"
			dat$predicted.blasts[which(dat$seurat_clusters%in%Cluster.case & dat$sampletype=="TMMPAL" & dat$celltype%in%c("Megakaryocytes","Erythrocytes","Other","B","PlasmaCells"))] <- "Normal_like"
			dat$predicted.blasts[which(!dat$seurat_clusters%in%Cluster.case & dat$sampletype=="TMMPAL")] <- "Normal_like"
			# AML
			dat$predicted.blasts[which(dat$seurat_clusters%in%Cluster.case & dat$sampletype=="AML" & !dat$celltype%in%c("Megakaryocytes","Erythrocytes","Other","B","PlasmaCells","T","NK"))] <- "Leukemia_like"
			dat$predicted.blasts[which(dat$seurat_clusters%in%Cluster.case & dat$sampletype=="AML" & dat$celltype%in%c("Megakaryocytes","Erythrocytes","Other","B","PlasmaCells","T","NK"))] <- "Normal_like"
			dat$predicted.blasts[which(!dat$seurat_clusters%in%Cluster.case & dat$sampletype=="AML")] <- "Normal_like"
			# TALL
			dat$predicted.blasts[which(dat$seurat_clusters%in%Cluster.case & dat$sampletype=="TALL" & !dat$celltype%in%c("Megakaryocytes","Erythrocytes","Other","B","PlasmaCells","Monocytes","Granulocyte"))] <- "Leukemia_like"
			dat$predicted.blasts[which(dat$seurat_clusters%in%Cluster.case & dat$sampletype=="TALL" & dat$celltype%in%c("Megakaryocytes","Erythrocytes","Other","B","PlasmaCells","Monocytes","Granulocyte"))] <- "Normal_like"
			dat$predicted.blasts[which(!dat$seurat_clusters%in%Cluster.case & dat$sampletype=="TALL")] <- "Normal_like"
			freq <- length(which(dat$predicted.blasts=="Leukemia_like"))/ncol(dat)
			return(freq)
		})
	frequency <- data.frame(cluster_leukemia_freq=seq(0.95,0.05,-0.05),predicted_blasts_freq=freq)

	# sfig3A
	pdf("frequency.blasts.pdf",height=4,width=10)
		ggplot(frequency,aes(x=cluster_leukemia_freq,y=predicted_blasts_freq)) +
			geom_point(size=5) +
			geom_line(size =0.8) +
			scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.05)) +
			theme(axis.text.x=element_text(angle=90,size=10),
		    axis.text.y=element_text(size=10))
	dev.off()
	# 确定 Normalcutoff=0.2
	x=0.2
	Cluster.case <- as.character(rst$Cluster[rst$Freq.Normal<x])
	dat$predicted.blasts <- "unknown"
	# Normal
	dat$predicted.blasts[which(dat$sampletype=="Normal")] <- "Normal"
	# TMMPAL
	dat$predicted.blasts[which(dat$seurat_clusters%in%Cluster.case & dat$sampletype=="TMMPAL" & !dat$celltype%in%c("Megakaryocytes","Erythrocytes","Other","B","PlasmaCells"))] <- "Leukemia_like"
	dat$predicted.blasts[which(dat$seurat_clusters%in%Cluster.case & dat$sampletype=="TMMPAL" & dat$celltype%in%c("Megakaryocytes","Erythrocytes","Other","B","PlasmaCells"))] <- "Normal_like"
	dat$predicted.blasts[which(!dat$seurat_clusters%in%Cluster.case & dat$sampletype=="TMMPAL")] <- "Normal_like"
	# AML
	dat$predicted.blasts[which(dat$seurat_clusters%in%Cluster.case & dat$sampletype=="AML" & !dat$celltype%in%c("Megakaryocytes","Erythrocytes","Other","B","PlasmaCells","T","NK"))] <- "Leukemia_like"
	dat$predicted.blasts[which(dat$seurat_clusters%in%Cluster.case & dat$sampletype=="AML" & dat$celltype%in%c("Megakaryocytes","Erythrocytes","Other","B","PlasmaCells","T","NK"))] <- "Normal_like"
	dat$predicted.blasts[which(!dat$seurat_clusters%in%Cluster.case & dat$sampletype=="AML")] <- "Normal_like"
	# TALL
	dat$predicted.blasts[which(dat$seurat_clusters%in%Cluster.case & dat$sampletype=="TALL" & !dat$celltype%in%c("Megakaryocytes","Erythrocytes","Other","B","PlasmaCells","Monocytes","Granulocyte"))] <- "Leukemia_like"
	dat$predicted.blasts[which(dat$seurat_clusters%in%Cluster.case & dat$sampletype=="TALL" & dat$celltype%in%c("Megakaryocytes","Erythrocytes","Other","B","PlasmaCells","Monocytes","Granulocyte"))] <- "Normal_like"
	dat$predicted.blasts[which(!dat$seurat_clusters%in%Cluster.case & dat$sampletype=="TALL")] <- "Normal_like"
	# 计算恶性细胞形态学数目和预测数目的相关性
	Freq <- data.frame(Sample=c("AML1","AML2","AML3","AML4","AML5","AML6","TALL1","TALL2","TALL3","TMMPAL1","TMMPAL2","TMMPAL3","TMMPAL4","TMMPAL5","TMMPAL6","TMMPAL7"),
		morphology.blasts=c(60.4,76,72.4,84.8,71.2,29.6,92,83.6,98.4,95,92,72,44,82,66.4,10.4),
		predicted.blasts=NA) 
	for(k in 1:nrow(Freq)){
		tmp <- length(which(dat$sample==as.character(Freq$Sample[k]) & dat$predicted.blasts=="Leukemia_like"))/length(which(dat$sample==Freq$Sample[k] & !dat$celltype%in%c("Erythrocytes","Megakaryocytes","Other")))
		Freq$predicted.blasts[k] <- tmp
	}
	Freq$group <- c(rep("AML",6),rep("TALL",3),rep("TMMPAL",7))

	# fig1E
	pdf("Predicted.blasts.cor.pdf",useDingbats=F,width=8)
	ggplot(Freq,aes(x=morphology.blasts,y=predicted.blasts)) +
			geom_point(size=7,aes(x=morphology.blasts,y=predicted.blasts,color=group)) +
			scale_color_manual(values=c("TMMPAL"="#ffd659","TALL"="#e4002b","AML"="#007ac2")) +
			geom_smooth(method='lm', color='red') +
			stat_cor(method="spearman") +
			 theme_bw() +
	        theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank())
	dev.off()

	# fig1D
	pdf("UMAP.predicted.blasts.pdf",useDingbats=F,width=8)
	col=c("Leukemia_like"="#ffc20e","Normal_like"="#90cef1","Normal"="#1b5faa")
	DimPlot(dat,reduction="umap",pt.size=0.1,cols=col,
		group.by="predicted.blasts",raster=FALSE)
	dev.off()

	# sfig1A
	pdf("UMAP.predicted.blasts.split.pdf",useDingbats=F,width=8)
	label <- c(paste0("AML",1:6),paste0("TALL",1:3),paste0("TMMPAL",1:8))
	for(i in 1:length(label)){
		col=c("Leukemia_like"="#ffc20e","Normal_like"="#90cef1","Normal"="#1b5faa")
		p <- DimPlot(dat,reduction="umap",pt.size=0.2,cols=col,
			cells=rownames(dat@meta.data)[which(dat$sample==label[i])],
			group.by="predicted.blasts") + ggtitle(label[i])
		print(p)
	}
	dev.off()

	# save integarte.all.rds
	saveRDS(dat, file="Integrate.all.rds")

	# 鉴定分化阻滞的位置
	library(ggalluvial)
	library(ggplot2)
	library(dplyr)
	library(reshape2)
	library(Seurat)
	# 按sampletype分类
	dat <- readRDS("/public/workspace/huangbin/MPAL/figure1/Integrate/Integrate.all.rds")
	tmp <- dat@meta.data
	tmp <- rbind(tmp[which(tmp$sampletype%in%c("TMMPAL","TALL","AML")&tmp$predicted.blasts=="Leukemia_like"),],
		tmp[which(tmp$sampletype=="Normal"),])
	counts <-	sapply(c("TMMPAL","TALL","AML","Normal"),function(x){
			sapply(c("HSC","MPP","CLP","T","NK","B","PlasmaCells","CMP_GMP","Monocytes","Granulocyte","MEP"),function(y){
				length(which(tmp$celltype==y & tmp$sampletype==x))
				})
			})
	counts <- as.data.frame(counts); counts$Type <- rownames(counts); tmp <- counts
	tmp[,1:4] <- apply(tmp[,1:4],2,function(x){x/sum(x)})
	proportion<-melt(tmp)

	# fig1G/freq.celltype.txt
	pdf("bar.pdf")
	sortlabel <- c("HSC","MPP","CLP","T","NK","B","PlasmaCells","CMP_GMP","Monocytes","Granulocyte","MEP")
	proportion$Type <- factor(proportion$Type,levels=sortlabel)
	ggplot(data=proportion)+
	geom_bar(stat="identity",aes(x=variable,y=value,fill=Type))+
	scale_fill_manual(values=celltype_color)+
	theme_classic()
	dev.off()
	write.table(tmp,"Freq.celltype.txt",sep="\t",quote=F,row.names=F)

	# fig1H
	library(ggpubr)
	library(reshape2)
	library(Seurat)
	dat <- readRDS("/public/workspace/huangbin/MPAL/figure1/Integrate/Integrate.all.rds")
	tmp <- dat@meta.data
	tmp <- rbind(tmp[which(tmp$sampletype%in%c("TMMPAL","TALL","AML")&tmp$predicted.blasts=="Leukemia_like"),],
		tmp[which(tmp$sampletype=="Normal"),])
	counts <-	sapply(names(table(tmp$sample)),function(x){
					sapply(c("HSC","MPP","CLP","T","NK","B","PlasmaCells","CMP_GMP","Monocytes","Granulocyte","MEP"),function(y){
							length(which(tmp$celltype==y & tmp$sample==x ))
				})
			})
	counts <- as.data.frame(counts)
	counts$Type <- rownames(counts)
	tmp <- counts
	tmp[,1:40] <- apply(tmp[,1:40],2,function(x){x/sum(x)})
	proportion<-melt(tmp)
	proportion <- proportion[which(proportion$Type%in%c("HSC","MEP","CLP","MPP","CMP_GMP")),]
	proportion$type <- "Normal"
	proportion$type[which(proportion$variable%in%paste0("AML",1:6))] <- "AML"
	proportion$type[which(proportion$variable%in%paste0("TALL",1:3))] <- "TALL"
	proportion$type[which(proportion$variable%in%paste0("TMMPAL",1:8))] <- "TMMPAL"
	proportion$type <- as.character(proportion$type)
	proportion$Type <- factor(proportion$Type,levels=c("HSC","MPP","MEP","CLP","CMP_GMP"))


	# 按sample分类
	dat <- readRDS("/public/workspace/huangbin/MPAL/figure1/Integrate/Integrate.all.rds")
	tmp <- dat@meta.data%>%filter(predicted.blasts=="Leukemia_like")
	counts <-	sapply(names(table(tmp$sample)),function(x){
			sapply(c("HSC","MPP","CLP","T","NK","B","PlasmaCells","CMP_GMP","Monocytes","Granulocyte","MEP"),function(y){
				length(which(tmp$celltype==y & tmp$sample==x))})})
	# 桑基图去掉Erythrocytes和Megakaryocytes和other///计算blasts里面的细胞比例
	counts <- as.data.frame(counts); counts$Type <- rownames(counts); tmp <- counts
	tmp[,-ncol(tmp)] <- apply(tmp[,-ncol(tmp)],2,function(x){x/sum(x)})
	proportion<-melt(tmp)

	# sfig4
	pdf("bar.sample.pdf",width=10)
	sortlabel <- c("HSC","MPP","CLP","T","NK","B","PlasmaCells","CMP_GMP","Monocytes","Granulocyte","MEP")
	proportion$Type <- factor(proportion$Type,levels=sortlabel)
	ggplot(data=proportion)+
	geom_bar(stat="identity",aes(x=variable,y=value,fill=Type))+
	scale_fill_manual(values=celltype_color)+
	theme_classic()
	dev.off()
	write.table(tmp,"Freq.sample.celltype.txt",sep="\t",quote=F,row.names=F)

	# fig1F
	dat <- readRDS("/public/workspace/huangbin/MPAL/figure1/Integrate/Integrate.all.rds")
	tmp <- subset(dat,cells=rownames(dat@meta.data)[which(dat$predicted.blasts=="Leukemia_like")])
	features <- c("CD33","MPO","CD7","CD3E")
	Idents(tmp) <- factor(tmp$sampletype,levels=c("AML","TALL","TMMPAL"))
	pdf("Dot.celltype.pdf",width=6,height=2)
	DotPlot(tmp,features=features,cols = c("white","darkgreen"),dot.scale=12)+
	 theme_linedraw()
	dev.off()	
}

################################################################
##                           Figure2                          ##
################################################################
{
	# fig2A 
	library(Seurat)
	library(dplyr)
	library(ggplot2)
	library(ggdensity)
	library(reshape2)
	dat <- readRDS("../figure3/MPAL.leukemia/MPAL.leukemia.rds")
	DefaultAssay(dat) <- "RNA"
	dat@meta.data <- cbind(dat@meta.data,dat@assays$RNA@data[
			c("CD3D","CD3E","CD3G","CD2","CD5","CD8A",
			"MME","DNTT","CD7","CD1A","MPO","ANPEP","CD33",
			"CD14","FCGR1A"),]%>%as.matrix()%>%t()%>%as.data.frame())
	tmp <- dat@meta.data[,c("sample","CD3D","CD3E","CD3G",
		"CD2","CD5","CD8A","MME","DNTT","CD7","CD1A",
		"MPO","ANPEP","CD33","CD14","FCGR1A")]
	tmp$myeloid <- apply(tmp[,c("MPO","CD33")],1,function(x){max(x)})%>%as.numeric()
	pdf("pseudo.flowcytometry.density.pdf",width=10,height=6)
	par(mfrow=c(2,4))
	Lab.palette <- colorRampPalette(c("blue","green", "orange", "red"), space = "Lab")
	for(i in paste0("TMMPAL",c(1,4,5,6,7,2,3,8))){
	 	tmp2 <- tmp%>%dplyr::filter(sample==i)
 	 	smoothScatter(tmp2$CD7,tmp2$myeloid,cex=6, 
			font.axis=2,font.lab=2,cex.lab=1.2,xlab="l",ylab="m",
			cex.axis=1.2,main=i,colramp=Lab.palette,
			nrpoints=0.2)
	 	tmp3 <- length(which(tmp2$CD7>1 & tmp2$myeloid>1))/nrow(tmp2)
	 	print(tmp3) }
	dev.off()
	pdf("pseudo.flowcytometry.point.pdf",width=10,height=6)
	par(mfrow=c(2,4))
	Lab.palette <- colorRampPalette(c("blue","green", "orange", "red"), space = "Lab")
	for(i in paste0("TMMPAL",c(1,4,5,6,7,2,3,8))){
	 	tmp2 <- tmp%>%dplyr::filter(sample==i)
	 	plot(tmp2[,c("CD7","myeloid")], col = densCols(tmp2[,c("CD7","myeloid")]), pch = 20) 	 	
		}
	dev.off()

	# fig2B
	# samplephenotype
	library(Seurat)
	library(dplyr)
	library(ggplot2)
	library(ggrepel)
	library(RColorBrewer)
	library(pheatmap)
	dat <- readRDS("../figure3/MPAL.leukemia/MPAL.leukemia.rds")
	dat$pheno <- "unknown"
	dat$pheno[which(dat$sample%in%c("TMMPAL1","TMMPAL4","TMMPAL5","TMMPAL6","TMMPAL7"))] <- "bilineal"
	dat$pheno[which(dat$sample%in%c("TMMPAL2","TMMPAL3"))] <- "biphenotype"
	p <- DimPlot(dat,group.by="pheno",pt.size=1.5,
		cols=c("bilineal"="#4e2a84","biphenotype"="#ec7f22","unknown"="grey60")) +
		theme_bw() +
	      theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) +
	      theme(aspect.ratio = 1)
	pdf("samplephenotype.pdf",width=8,useDingbats=F)
	print(p)
	dev.off()
	p <- DimPlot(dat,group.by="sample",pt.size=0.8,
		cols=c("TMMPAL1"="#e4e0ee","TMMPAL4"="#b6acd1","TMMPAL5"="#836eaa",
  			"TMMPAL6"="#4e2a84","TMMPAL7"="#401f68","TMMPAL2"="#ec7f22",
  			"TMMPAL3"="#a13323","TMMPAL8"="grey60")) +
		  theme_bw() +
	      theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) +
	      theme(aspect.ratio = 1)
	pdf("samples.pdf",width=8,useDingbats=F)
	print(p)
	dev.off()

	# fig2B project
	.libPaths(c('/public/workspace/huangbin/plaqueRlibrary'))
	library(Seurat)
	library(dplyr)
	library(future)
	library(ggplot2)
	## reference
	ref <- readRDS("../figure1/Integrate/Integrate.all.rds")
	ref <- subset(ref,cells=rownames(ref@meta.data)[which(ref$sampletype=="Normal" &
		!ref$celltype%in%c("HSC","MPP","Other"))])
	ref$celltype <- as.character(ref$celltype)
	ref$new <- ref$celltype
	ref$new[which(ref$celltype%in%c("B","CLP","NK","T","PlasmaCells"))] <- "lymphoid"
	ref$new[which(ref$celltype%in%c("CMP_GMP","Granulocyte","Monocytes"))] <- "myeloid"
	ref$new[which(ref$celltype%in%c("Erythrocytes","Megakaryocytes","MEP"))] <- "megakaryoerythroid"
	ref <- NormalizeData(ref)
	ref <- FindVariableFeatures(ref, selection.method = "vst", nfeatures = 2000)
	ref <- ScaleData(ref, features = rownames(ref))
	ref <- RunPCA(ref, features = VariableFeatures(object = ref))
	saveRDS(ref,file="lineage.ref.rds")
	# query & Transfer
	set.seed(666)
	query <- readRDS("../figure3/MPAL.leukemia/MPAL.leukemia.rds")
	DefaultAssay(query) <- "RNA"
	ref <- readRDS("lineage.ref.rds")
	anchors <- FindTransferAnchors(reference=ref,query=query,dims=1:30,reference.reduction="pca")
	predictions <- TransferData(anchorset=anchors,refdata=ref$new,dims=1:30)
	query <- AddMetaData(query,metadata=predictions)
	saveRDS(predictions,file="predictions.mpal.leukemia.rds")
	saveRDS(query,file="query.mpal.leukemia.rds")
	# split
	query <- readRDS("query.mpal.leukemia.rds")
	query$delta <- abs(query$prediction.score.myeloid-query$prediction.score.lymphoid)
	query$label <- query$predicted.id
	query$label[which(query$delta<0.2)] <- "both"
	library(future)
	plan("multiprocess",workers=10)
	Idents(query) <- factor(query$label,levels=c(
		"lymphoid","myeloid","both","megakaryoerythroid"))
	marker <- FindAllMarkers(query,
	  only.pos = TRUE,
	  min.pct = 0.1,
	  logfc.threshold = 0.1)
	# phenotype split
	query$sample <- factor(query$sample,levels=paste0("TMMPAL",c(1,4,5,6,7,2,3,8)))
	p <- DimPlot(query,group.by="label",
		split.by="sample",pt.size=0.1,
		cols=c("lymphoid"="#49a942","myeloid"="#48a9c5",
			"both"="#ffc20e","megakaryoerythroid"="#d7d7d8")) +
		theme_bw() +
	      theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) +
	      theme(aspect.ratio = 1)
	pdf("samplephenotype.split.pdf",width=14,useDingbats=F)
	print(p)
	dev.off()

	# fig2C
	freq <- sapply(paste0("TMMPAL",c(1,4,5,6,7,2,3,8)),function(x){
		sapply(c("lymphoid","myeloid","both","megakaryoerythroid"),function(y){
			tmp <- query@meta.data[which(query$sample==x),]
			tmp <- length(which(tmp$label==y))/nrow(tmp)
			return(tmp)
			})
		})
	pdf("freq.bar.pdf")
	p <- barplot(freq,ylim=c(0,1),axisnames=T,border=NA,
		space=0.05,col=c("lymphoid"="#49a942","myeloid"="#48a9c5",
			"both"="#ffc20e","megakaryoerythroid"="#d7d7d8"))
	axis(side=1,p,labels=F)
	dev.off()
	tmp <- data.frame(freq=as.numeric(freq["both",1:7]),
		type=c(rep(c("bilineal","biphenotype"),c(5,2))))
	pdf("freq.box.pdf",width=3)
	ggplot(data=tmp,aes(x=type,y=freq))+
		geom_boxplot()+
			geom_point(size=6)+
				theme_bw() +
	      		theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank())
	dev.off()

	# sfig5A-B bulk.heatmap 8个样本聚类
	library(psych)
	library(reshape2)
	library(ggplot2)
	library(factoextra)
	library(dplyr)
	library(Seurat)
	library(RColorBrewer)
	library(pheatmap)
	dat <- readRDS("../figure3/MPAL.leukemia/MPAL.leukemia.rds")
	DefaultAssay(dat) <- "integrated"
	bulk <- lapply(paste0("TMMPAL",1:8),function(x){
		tmp <- dat@assays$RNA@data[,rownames(dat@meta.data)[which(dat$sample==x)]]
		tmp <- rowSums(tmp)/ncol(tmp)
		return(tmp) })
	bulk <- do.call(cbind,bulk)
	colnames(bulk) <- paste0("TMMPAL",1:8)
	mad <- apply(bulk,1,function(x){mad(x)})
	annotation.col <- data.frame(
		Type=c("bilineal","biphenotype","biphenotype",rep("bilineal",4),"unknown"))
	rownames(annotation.col) <- paste0("TMMPAL",1:8)
	ann_colors = list(Type=c("bilineal"="#7ab800","biphenotype"="#0085c3","unknown"="grey60"))
	pdf("bulk.pheatmap.mad10000.pdf",useDingbats=F,height=6,width=6)
	pheatmap(as.matrix(bulk[names(rev(sort(mad))[1:10000]),]), 
	  breaks=seq(-2,2,0.01),
	  col=rev(colorRampPalette(brewer.pal(11,"RdBu"))(length(seq(-2,2,0.01)))),
	  scale="row",
	  annotation_col=annotation.col,
	  show_rownames=FALSE,cluster_rows=TRUE,
	  cluster_cols=TRUE,show_colnames=TRUE,
	  annotation_colors=ann_colors,
	  fontsize=8)
	dev.off()
	pdf("bulk.pheatmap.mad1000.pdf",useDingbats=F,height=6,width=6)
	pheatmap(as.matrix(bulk[names(rev(sort(mad))[1:1000]),]), 
	  breaks=seq(-2,2,0.01),
	  col=rev(colorRampPalette(brewer.pal(11,"RdBu"))(length(seq(-2,2,0.01)))),
	  scale="row",
	  annotation_col=annotation.col,
	  show_rownames=FALSE,cluster_rows=TRUE,
	  cluster_cols=TRUE,show_colnames=TRUE,
	  annotation_colors=ann_colors,
	  fontsize=8)
	dev.off()

	# sfig5C
	ref <- readRDS("lineage.ref.rds")
	query <- readRDS("query.mpal.leukemia.rds")
	p <- DimPlot(ref,group.by="celltype",
		pt.size=1,cols=celltype_color) + theme_bw() +
      	theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) +
      	theme(aspect.ratio = 1)
	pdf("ref.pdf",width=8,useDingbats=F)
	print(p)
	dev.off()

	query <- readRDS("query.mpal.leukemia.rds")
	query$delta <- abs(query$prediction.score.myeloid-query$prediction.score.lymphoid)
	query$label <- query$predicted.id
	query$label[which(query$delta<0.2)] <- "both"
	p <- DimPlot(query,group.by="label",pt.size=1,
		cols=c("lymphoid"="#49a942","myeloid"="#48a9c5",
			"both"="#ffc20e","megakaryoerythroid"="#d7d7d8")) +
	theme_bw() +
  	theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()) +
  	theme(aspect.ratio = 1)
	pdf("query.pdf",width=8,useDingbats=F)
	print(p)
	dev.off()
}

################################################################
##          Preprocess: Integration AML & ALL                 ##
################################################################
{
	library(Seurat)
	library(clustree)
	library(RColorBrewer)
	setwd("/public/workspace/huangbin/MPAL/figure3/TALL.AML.leukemia/")
	data <- readRDS("/public/workspace/huangbin/MPAL/figure1/Integrate/Integrate.all.rds")
	DefaultAssay(data) <- "RNA"
	Leukemiacells <- subset(data, cells=colnames(data)[which(data$predicted.blasts=="Leukemia_like")])
	list <- SplitObject(Leukemiacells, split.by = "sample")
	for(i in 1:length(list)){
		tmp <- NormalizeData(list[[i]])
		tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 2000)
		list[[i]] <- tmp
	}
	list <- list[9:17]
	anchors <- FindIntegrationAnchors(object.list = list, dims = 1:30)
	dat <- IntegrateData(anchorset = anchors, dims = 1:30)
	dat <- ScaleData(dat, features = rownames(dat))
	dat <- RunPCA(dat, features = VariableFeatures(object = dat))
	dat <- FindNeighbors(dat, dims = 1:30)

	dat@meta.data <- dat@meta.data[,c("orig.ident","nCount_RNA","nFeature_RNA",
		"percent.mt","seurat_clusters","sample","sampletype","celltype")]
	dat <- FindClusters(object = dat, resolution = c(seq(0.5,1.5,0.1)))

	# sfig6
	pdf("./Clustree.pdf",height=12,width=9) # k=1.2
	clustree(dat@meta.data,prefix = "integrated_snn_res." ) +
	scale_color_manual(values=colorRampPalette(brewer.pal(n=8,name="PuOr"))(11)) + 
	theme_classic()
	dev.off()

	dat <- FindClusters(object = dat, resolution = 1.2)
	dat <- RunUMAP(dat, dims = 1:30)

	dat@meta.data <- dat@meta.data[,c("orig.ident","nCount_RNA","nFeature_RNA",
		"percent.mt","seurat_clusters","sample","sampletype","celltype")]
	DefaultAssay(dat) <- "RNA"
	saveRDS(dat, file="./TALL.AML.leukemia.rds")

	# fig3A
	pdf("./TALL.AML.leukemia.umap.cluster1.pdf",width=7.5,useDingbats=F)
	DimPlot(dat,cols=color1)
	dev.off()

	# fig3A
	pdf("./TALL.AML.leukemia.umap.cluster2.pdf",width=14,useDingbats=F)
	DimPlot(dat,cols=color1,split.by="sampletype",pt.size=0.1)
	dev.off()
}

################################################################
##                           Figure3                          ##
################################################################
{
	library(Seurat)
	library(RColorBrewer)
	library(reshape2)
	library(dplyr)
	library(ggpubr)
	setwd("/public/workspace/huangbin/MPAL/figure3/TALL.AML.leukemia/")
	dat <- readRDS("./TALL.AML.leukemia.rds")
	rst <- lapply(names(table(dat$seurat_clusters)),function(x){
		sapply(names(table(dat$sample)),function(y){
			tmp <- length(which(dat$seurat_clusters==x & dat$sample==y))
			return(tmp)
			})
		})
	rst <- do.call(rbind,rst)
	rownames(rst) <- paste0("C",names(table(dat$seurat_clusters)))
	rst <- apply(t(rst),2,function(x){x/sum(x)})
	rank <- apply(rst[paste0("TALL",1:3),],2,function(x){sum(x)})
	rst <- rst[,names(rev(sort(rank)))]

	# fig3B
	pdf("eachcluster.sample.pdf",width=15)
	color2=colorRampPalette(brewer.pal(n=8,name="Blues"))(8)
	color3=colorRampPalette(brewer.pal(n=8,name="Reds"))(5)
	samplecolor=c(
	"AML1"=color2[3],"AML2"=color2[4],"AML3"=color2[5],
	"AML4"=color2[6],"AML5"=color2[7],"AML6"=color2[8],
	"TALL1"=color3[3],"TALL2"=color3[4],"TALL3"=color3[5])
	barplot(rst,col=samplecolor)
	legend("top",rownames(rst),pch=rep(19,9),col=samplecolor)
	dev.off()

	# sfig7A
	rst <- as.data.frame(rst)
	rst$group <- rep(c("AML","TALL"),c(6,3))
	rst$sample <- rownames(rst)
	rst <- melt(rst)%>%as.data.frame()
	rst$variable <- factor(as.factor(rst$variable),levels=paste0("C",c(0:24)))
	color1=colorRampPalette(brewer.pal(n=8,name="Greys"))(25)
	color2=colorRampPalette(brewer.pal(n=8,name="Blues"))(8)
	color3=colorRampPalette(brewer.pal(n=8,name="Reds"))(5)
	color4=colorRampPalette(brewer.pal(n=8,name="Oranges"))(20)
	samplecolor=c(
		"Normal1"=color1[3],"Normal2"=color1[4],"Normal3"=color1[5],"Normal4"=color1[6],"Normal5"=color1[7],"Normal6"=color1[8],
		"Normal7"=color1[9],"Normal8"=color1[10],"Normal9"=color1[11],"Normal10"=color1[12],"Normal11"=color1[13],"Normal12"=color1[14],
		"Normal13"=color1[15],"Normal14"=color1[16],"Normal15"=color1[17],"Normal16"=color1[18],"Normal17"=color1[19],"Normal18"=color1[20],
		"Normal19"=color1[21],"Normal20"=color1[22],"Normal21"=color1[23],"Normal22"=color1[24],"Normal23"=color1[25],"AML1"=color2[3],
		"AML2"=color2[4],"AML3"=color2[5],"AML4"=color2[6],"AML5"=color2[7],"AML6"=color2[8],"TALL1"=color3[3],"TALL2"=color3[4],"TALL3"=color3[5],
		"TMMPAL1"=color4[3],"TMMPAL2"=color4[4],"TMMPAL3"=color4[5],"TMMPAL4"=color4[6],"TMMPAL5"=color4[7],"TMMPAL6"=color4[8],
		"TMMPAL7"=color4[9],"TMMPAL8"=color4[10])
	pdf("eachcluster.sample.box.pdf",width=10,height=3)
	ggplot(data=rst,aes(x=group,y=value)) +
	geom_boxplot(fill="white",outlier.shape=NA) +
	stat_boxplot(geom = "errorbar",width=0.3) +
	geom_point(aes(color=sample),position="jitter",size=1.5) +
	scale_color_manual(values=samplecolor) +
	stat_compare_means(aes(group=group),label="p.signif",
	                     method = "wilcox.test") +
	facet_wrap(~variable,ncol=13) + theme_classic()
	dev.off()

	# Cacoa
	library(sccore)
	library(cacoa)
	library(cowplot)
	library(ggplot2)
	library(Seurat)
	# 判断是不是TALL/AML enriched
	setwd("/public/workspace/huangbin/MPAL/figure3/TALL.AML.leukemia/")
	dat <- readRDS("TALL.AML.leukemia.rds")
	sample.groups <- as.character(ifelse(dat$sampletype=="TALL","case","control")); names(sample.groups) <- dat$sample
	cell.groups <- as.character(dat$seurat_clusters); names(cell.groups) <- colnames(dat)
	sample.per.cell <- as.character(dat$sample); names(sample.per.cell) <- colnames(dat)
	ref.level <- "control"
	target.level <- "case"
	graph.name <- NULL
	embedding <- "UMAP"
	cao <- Cacoa$new(dat,
		sample.groups=sample.groups,
		cell.groups=cell.groups,
		sample.per.cell=sample.per.cell, 
	    ref.level=ref.level,
	    target.level=target.level,
	    n.cores=30,
	    graph.name=graph.name)
	cao$plot.params <- list(size=0.1, alpha=0.1, font.size=c(2, 3))
	cao$plot.theme <- cao$plot.theme + theme(legend.background=element_blank())
	# estimate cluster-based changes
	rst <- cao$estimateCellLoadings()
	names(color1) <- as.character(c(0:29))
	saveRDS(dat, file="./TALL.AML.leukemia.rds")

	# fig3C
	pdf("CellLoadings.pdf", useDingbats=F,width=3,height=4)
	p <- cao$plotCellLoadings(show.pvals=FALSE, alpha=0, palette=color1)
	print(p)
	dev.off()

	# Project
	.libPaths(c('/public/workspace/huangbin/plaqueRlibrary'))
	library(Seurat)
	library(dplyr)
	setwd("/public/workspace/huangbin/MPAL/figure3/TALL.AML.leukemia/labeltransfer/")
	data <- readRDS("/public/workspace/huangbin/MPAL/figure1/Integrate/Integrate.all.rds"); DefaultAssay(data) <- "RNA"
	Leukemiacells <- subset(data, cells=colnames(data)[which(data$predicted.blasts=="Leukemia_like")])
	list <- SplitObject(Leukemiacells,split.by="sample")
	ref <- readRDS("/public/workspace/huangbin/MPAL/figure3/TALL.AML.leukemia/TALL.AML.leukemia.rds"); DefaultAssay(ref) <- "integrated"

	for(x in paste0("TMMPAL",1:8)){
		set.seed(123)
		query <- list[[x]]
		query <- NormalizeData(query)
		query <- FindVariableFeatures(query,selection.method="vst",nfeatures=2000)
		query <- ScaleData(query,features=rownames(query))
		query <- RunPCA(query,features=VariableFeatures(object=query))
		query <- RunUMAP(query,reduction="pca",dims=1:30)
		anchors <- FindTransferAnchors(reference=ref,query=query,dims=1:30,reference.reduction="pca")
		predictions <- TransferData(anchorset=anchors,refdata=ref$seurat_clusters,dims=1:30)
		query <- AddMetaData(query,metadata=predictions)
		saveRDS(predictions,file=paste0("./",x,".predictions.rds"))
	}

	AML.enriched <- c(10,21,18,6,9,4)
	TALL.enriched <- c(17,5,15,22,8,11,20,2,1,0,7,12,3)

	# 在整合MPAL中计算每个样本中投射的簇的柱状图
	info <- lapply(list.files("./",pattern="predictions.rds"),function(x){
		pred <- readRDS(x)
		pred <- data.frame(cell=rownames(pred),pred=pred$predicted.id)
		return(pred)})
	info <- do.call(rbind,info)
	dat <- readRDS("/public/workspace/huangbin/MPAL/figure3/MPAL.leukemia/MPAL.leukemia.rds")
	dat$cell <- rownames(dat@meta.data)
	dat@meta.data <- left_join(dat@meta.data,info,by="cell")
	rownames(dat@meta.data) <- dat$cell
	dat$info <- "unknown"
	dat$info[which(dat$pred%in%AML.enriched)] <- "AEC"
	dat$info[which(dat$pred%in%TALL.enriched)] <- "TEC"
	dat$info[which(dat$pred==24)] <- "S24"
	dat$info[which(dat$pred==23)] <- "S23"
	dat$info[which(dat$pred==19)] <- "S19"
	dat$info[which(dat$pred==16)] <- "S16"
	dat$info[which(dat$pred==14)] <- "S14"
	dat$info[which(dat$pred==13)] <- "S13"
	rst <- lapply(names(table(dat$sample)),function(x){
		sapply(names(table(dat$info)),function(y){
			tmp <- length(which(dat$sample==x & dat$info==y))
			return(tmp)
			})
		})
	rst <- do.call(rbind,rst)
	rownames(rst) <- names(table(dat$sample))
	rst <- apply(t(rst),2,function(x){x/sum(x)})
	rst <- rst[c("S24","TEC","AEC","S13","S14","S16","S19"),
	names(rev(sort(rst["S24",])))]

	# fig3J
	pdf("transfer.cluster.bar.pdf",width=10)
	color=c("S24"="#ffc20e",'TEC'="#ce181e",'AEC'="#007cc0",
		'S13'="#b4a996",'S14'="#F0027F",'S16'="palevioletred4",
		'S19'="#CAB2D6")
	barplot(rst,col=color)
	legend("top",rownames(rst),pch=rep(19,7),col=color)
	dev.off()
	dat$info <- as.factor(dat$info)

	# sfig7B
	pdf("transfer.cluster.umap.pdf",useDingbats=F,width=8)
	DimPlot(dat,group.by="info",pt.size=0.1,cols=color)
	dev.off()

	# save MPAL.pred.rds
	saveRDS(dat,file="./MPAL.pred.rds")

	# 在整合的MPAL上面染投射簇的UMAP
	info <- lapply(list.files("./",pattern="predictions.rds"),function(x){
		pred <- readRDS(x)
		pred <- data.frame(cell=rownames(pred),pred=pred$predicted.id)
		return(pred)})
	info <- do.call(rbind,info)
	dat <- readRDS("/public/workspace/huangbin/MPAL/figure3/MPAL.leukemia/MPAL.leukemia.rds")
	dat$cell <- rownames(dat@meta.data)
	dat@meta.data <- left_join(dat@meta.data,info,by="cell")
	rownames(dat@meta.data) <- dat$cell

	# fig3D
	pdf("transfer.cluster.pdf",useDingbats=F,width=8)
	DimPlot(dat,group.by="pred",pt.size=0.1,cols=color2)
	dev.off()

	# *** marker.type7.txt
	library(future)
	library(Seurat)
	library(RColorBrewer)
	library(dplyr)
	library(pheatmap)
	plan("multiprocess",workers=10)
	setwd("/public/workspace/huangbin/MPAL/figure3/TALL.AML.leukemia/labeltransfer/")
	dat <- readRDS("./MPAL.pred.rds")
	DefaultAssay(dat) <- "RNA"
	Idents(dat) <- factor(dat$info,levels=c(
		"AEC","TEC","S13","S14","S16","S19","S24"))
	marker <- FindAllMarkers(dat,
	  only.pos = TRUE,
	  min.pct = 0.25,
	  logfc.threshold = 0.25)
	marker <- marker[which((marker$pct.1-marker$pct.2)>0.1),]
	write.table(marker,"marker.type7.txt",sep="\t",quote=F)

	# pheatmap 
	setwd("/public/workspace/huangbin/MPAL/figure3/TALL.AML.TMMPAL.leukemia")
	dat <- readRDS("./TALL.AML.TMMPAL.leukemia.rds")
	dat2 <- readRDS("../TALL.AML.leukemia/labeltransfer/MPAL.pred.rds")
	sTEC <- rownames(dat2@meta.data)[which(dat2$info=="TEC")]
	sAEC <- rownames(dat2@meta.data)[which(dat2$info=="AEC")]
	s24 <- rownames(dat2@meta.data)[which(dat2$info=="S24")]
	sTALL <- rownames(dat@meta.data)[which(dat$sampletype=="TALL")]
	sAML <- rownames(dat@meta.data)[which(dat$sampletype=="AML")]
	marker <- read.table("../TALL.AML.leukemia/labeltransfer/marker.type7.txt",sep="\t",stringsAsFactors=F,header=T,row.names=1)
	top <- marker%>%filter(cluster%in%c("TEC","S24","AEC"))%>%
	group_by(cluster)%>%top_n(12,wt=avg_log2FC)
	top <- c(top$gene[which(top$cluster=="TEC")],
		top$gene[which(top$cluster=="S24")],
		top$gene[which(top$cluster=="AEC")])
	tmp <- cbind(as.matrix(dat@assays$RNA@data[top,sTALL]),
		as.matrix(dat@assays$RNA@data[top,sTEC]),
		as.matrix(dat@assays$RNA@data[top,s24]),
		as.matrix(dat@assays$RNA@data[top,sAEC]),
		as.matrix(dat@assays$RNA@data[top,sAML]))
	annotation.col <- data.frame(
		Type=rep(c("TALL","TEC","S24","AEC","AML"),
		c(length(sTALL),length(sTEC),length(s24),length(sAEC),length(sAML))))
	rownames(annotation.col) <- c(sTALL,sTEC,s24,sAEC,sAML)
	ann_colors = list(Type=c('TEC'="#ce181e","S24"="#ffc20e",'AEC'="#007cc0",
		"TALL"="#ffaaaa","AML"="#b3dcff"))

	# fig3G
	pdf("pheatmap.pdf",useDingbats=F,height=6)
	pheatmap(as.matrix(tmp), 
	  breaks=seq(-3,3,0.01),
	  col=rev(colorRampPalette(brewer.pal(11,"RdBu"))(length(seq(-3,3,0.01)))),
	  scale="row",
	  annotation_col=annotation.col,
	  show_rownames=TRUE,cluster_rows=FALSE,
	  cluster_cols=FALSE,show_colnames=FALSE,
	  annotation_colors=ann_colors,
	  fontsize=8)
	dev.off()

	# monocle2
	module load r/4.0.4
	options(future.globals.maxSize = 2400000 * 1024^2)
	.libPaths(c('/public/workspace/huangbin/plaqueRlibrary'))
	library(Seurat)
	library(dplyr)
	library(reshape2)
	library(RColorBrewer)
	library(pheatmap)
	library(ggplot2)
	library(monocle)
	library(BiocGenerics)
	importCDS <- function (otherCDS, seurat_scale=F, import_all = FALSE) 
	{
	    if (class(otherCDS)[1] == "Seurat") {
	        requireNamespace("Seurat")
	        if (!seurat_scale) {
	          data <- otherCDS@assays$RNA@counts
	        } else {
	          data <- otherCDS@assays$RNA@scale.data
	        }
	        if (class(data) == "data.frame") {
	            data <- as(as.matrix(data), "sparseMatrix")
	        }
	        pd <- tryCatch({
	            pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
	            pd
	        }, error = function(e) {
	            pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
	            pd <- new("AnnotatedDataFrame", data = pData)
	            message("This Seurat object doesn't provide any meta data")
	            pd
	        })
	        if (length(setdiff(colnames(data), rownames(pd))) > 0) {
	            data <- data[, rownames(pd)]
	        }
	        fData <- data.frame(gene_short_name = row.names(data), 
	            row.names = row.names(data))
	        fd <- new("AnnotatedDataFrame", data = fData)
	        #lowerDetectionLimit <- otherCDS@is.expr
	        if (all(data == floor(data))) {
	            expressionFamily <- negbinomial.size()
	            expr <- "negbinomial.size"
	        }
	        else if (any(data < 0)) {
	            expressionFamily <- uninormal()
	            expr <- "unimormal"
	        }
	        else {
	            expressionFamily <- tobit()
	            expr <- "tobit"
	        }
	        print(paste0("expressionFamily ",expr))
	       
	        monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd, 
	            
	            expressionFamily = expressionFamily)
	        if (import_all) {
	            if ("Monocle" %in% names(otherCDS@misc)) {
	                otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
	                otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
	                monocle_cds <- otherCDS@misc$Monocle
	                mist_list <- otherCDS
	            }
	            else {
	                mist_list <- otherCDS
	            }
	        }
	        else {
	            mist_list <- list()
	        }
	        if ("var.genes" %in% slotNames(otherCDS)) {
	            var.genes <- setOrderingFilter(monocle_cds, otherCDS@var.genes)
	        }
	        monocle_cds@auxClusteringData$seurat <- mist_list
	    }
	    else if (class(otherCDS)[1] == "SCESet") {
	        requireNamespace("scater")
	        message("Converting the exprs data in log scale back to original scale ...")
	        data <- 2^otherCDS@assayData$exprs - otherCDS@logExprsOffset
	        fd <- otherCDS@featureData
	        pd <- otherCDS@phenoData
	        experimentData = otherCDS@experimentData
	        if ("is.expr" %in% slotNames(otherCDS)) 
	            lowerDetectionLimit <- otherCDS@is.expr
	        else lowerDetectionLimit <- 1
	        if (all(data == floor(data))) {
	            expressionFamily <- negbinomial.size()
	        }
	        else if (any(data < 0)) {
	            expressionFamily <- uninormal()
	        }
	        else {
	            expressionFamily <- tobit()
	        }
	        if (import_all) {
	            mist_list <- otherCDS
	        }
	        else {
	            mist_list <- list()
	        }
	        monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd, 
	            lowerDetectionLimit = lowerDetectionLimit, expressionFamily = expressionFamily)
	        monocle_cds@auxOrderingData$scran <- mist_list
	    }
	    else {
	        stop("the object type you want to export to is not supported yet")
	    }
	    return(monocle_cds)
	}

	dat <- readRDS("../TALL.AML.TMMPAL.leukemia/TALL.AML.TMMPAL.leukemia.rds")
	set.seed(23); dat <- subset(dat,cell=colnames(dat)[sample(1:ncol(dat),30000)])
	dat <- NormalizeData(dat)
	dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
	cds <- importCDS(dat)
	cds <- estimateSizeFactors(cds)
	cds <- estimateDispersions(cds)
	DEG_genes <- VariableFeatures(dat)
	cds <- setOrderingFilter(cds, ordering_genes = DEG_genes)
	cds <- reduceDimension(cds, method = 'DDRTree')
	cds <- orderCells(cds,root_state=5) #
	save(cds,file="cds2.RData")

	info <- readRDS("../TALL.AML.leukemia/labeltransfer/MPAL.pred.rds")
	info <- data.frame(cell=rownames(info@meta.data),info=info$info)

	pData(cds)$group <- pData(cds)$sampletype
	pData(cds)$group[which(rownames(pData(cds))%in%info$cell[which(info$info=="S24")])] <- "S24"
	pData(cds)$group[which(rownames(pData(cds))%in%info$cell[which(info$info=="AEC")])] <- "AEC"
	pData(cds)$group[which(rownames(pData(cds))%in%info$cell[which(info$info=="TEC")])] <- "TEC"
	pData(cds)$group[which(rownames(pData(cds))%in%info$cell[which(info$info=="S13")])] <- "S13"
	pData(cds)$group[which(rownames(pData(cds))%in%info$cell[which(info$info=="S14")])] <- "S14"
	pData(cds)$group[which(rownames(pData(cds))%in%info$cell[which(info$info=="S16")])] <- "S16"
	pData(cds)$group[which(rownames(pData(cds))%in%info$cell[which(info$info=="S19")])] <- "S19"
	pData(cds)$group <- factor(as.factor(pData(cds)$group),levels=c("TALL","TEC","S24","S13","S14","S16","S19","AEC","AML"))

	# fig3E
	pdf("./Trajectory.pdf",width=8.5,useDingbats=F)
	plot_cell_trajectory(cds, cell_size=1,show_state_number=F,
	    color_by="group",cell_link_size = 1)+
	    scale_color_manual(values=c(typecolor))+
	    theme(legend.position = "right") 
	dev.off()

	# sfig7D
	pdf("./Trajectory.wrap.pdf",width=20,useDingbats=F,height=5)
	plot_cell_trajectory(cds, cell_size=1,show_state_number=F,
	  color_by="group",cell_link_size = 1)+
	facet_wrap(~group, nrow = 1)+
	    scale_color_manual(values=c(typecolor,"TMMPAL"="grey60"))+
	    theme(legend.position = "right") 
	dev.off()

	# sfig7C
	myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
	pdf("./Pseudotime.pdf",width=8.5,useDingbats=F)
	plot_cell_trajectory(cds, cell_size=1,show_state_number=FALSE,
	  color_by=c("Pseudotime"),cell_link_size = 0.5)+
	    scale_colour_gradientn(colours = myPalette(100))+
	    theme(legend.position = "right") 
	dev.off()

	# fig3F
	pData(cds)$component1 <- as.numeric(cds@reducedDimS[1,])
	pdf("./density.pdf",width=6,height=3)
	ggplot(pData(cds),aes(component1,colour=group)) +
	  geom_density()+
	      scale_color_manual(values=c(typecolor))+
	  theme_classic()
	dev.off()


	# GSE113601
	library(dplyr)
	library(RColorBrewer)
	library(pheatmap)
	setwd("/public/workspace/huangbin/MPAL/figure3/bulk/GSE113601")
	dat <- read.table("GSE113601.MPAL.txt",sep="\t",stringsAsFactors=F,header=T)
	dat$GENE <- gsub("^ENSG[0-9]*.\\||","",dat$GENE)
	dat <- dat%>%group_by(GENE)%>%summarise_all(mean)%>%as.data.frame()%>%filter(!is.na(GENE))
	rownames(dat) <- dat$GENE
	dat <- dat[,-1]
	TM.index <- c(paste0("GSM31094",c(14:21,25:27,31,34:37)))
	dat <- dat[,TM.index]

	markers <- read.table("/public/workspace/huangbin/MPAL/figure3/TALL.AML.leukemia/labeltransfer/marker.type7.txt",
	  sep="\t",header=T,stringsAsFactors=F)
	top <- markers%>%filter(cluster%in%c("TEC","S24","AEC"))%>%
	group_by(cluster)%>%top_n(8,wt=avg_log2FC)
	select <- intersect(top$gene,rownames(dat))
	
	# fig3K
	pdf("GSE113601.heatmap.pdf",height=5)
	p <- pheatmap(dat[select,], 
	  breaks=seq(-3.5,3.5,0.001),scale="row",border_color=NA,
	  col=rev(colorRampPalette(brewer.pal(11,"RdBu"))(length(seq(-3.5,3.5,0.001)))),
	  show_rownames=TRUE,cluster_rows=TRUE,
	  cutree_cols=4,cluster_cols=TRUE,show_colnames=FALSE,fontsize=8)
	print(p)

	dev.off()

	# target - T/My MPAL
	library(dplyr)
	library(pheatmap)
	library(RColorBrewer)
	setwd("/public/workspace/huangbin/MPAL/figure3/bulk/target")
	dat <- read.table("MPAL-allGeneExpr-90S-rlog.txt",sep="\t",stringsAsFactors=F,header=T)
	dat <- dat%>%filter(!Gene%in%dat$Gene[duplicated(dat$Gene)]);rownames(dat) <- dat$Gene
	dat <- dat[,-c(1,2)]; dat <- dat[,-which(colnames(dat)=="SJMPAL011911_D2")]
	colnames(dat) <- as.character(sapply(colnames(dat),function(x){strsplit(x,"_")[[1]][1]}))

	info <- read.table("target.clinic.txt",sep="\t",stringsAsFactors=F,header=T)
	dat <- dat[,intersect(info$Sample,colnames(dat))]

	markers <- read.table("/public/workspace/huangbin/MPAL/figure3/TALL.AML.leukemia/labeltransfer/marker.type7.txt",
	  sep="\t",header=T,stringsAsFactors=F)
	top <- markers%>%filter(cluster%in%c("TEC","S24","AEC"))%>%
	group_by(cluster)%>%top_n(8,wt=avg_log2FC)
	select <- intersect(top$gene,rownames(dat))

	# fig3K
	pdf("Target.heatmap.pdf",height=3)
	p <- pheatmap(dat[select,], 
	  breaks=seq(-3.5,3.5,0.001),
	  scale="row",border_color=NA,
	  col=rev(colorRampPalette(brewer.pal(11,"RdBu"))(length(seq(-3.5,3.5,0.001)))),
	  show_rownames=TRUE,cluster_rows=TRUE,clustering_method="ward.D2",
	  cutree_cols=3,cluster_cols=TRUE,show_colnames=FALSE,fontsize=8)
	print(p)
	dev.off()

	save(p,file="p.RData")

	names(which(cutree(p$tree_col,3)==3))
	"SJMPAL011915" "SJMPAL012420" "SJMPAL042792" "SJMPAL043510"

	# target
	.libPaths("~/plaqueRlibrary/")
	library(dplyr)
	library(tidyr)
	library(survival)
	library(survminer)
	library(ssgsea.GBM.classification)
	setwd("/public/workspace/huangbin/MPAL/figure3/bulk/target/")
	clinic <- read.table("/public/workspace/huangbin/MPAL/figure4/surv/clinic.txt",sep="\t",header=T,stringsAsFactors=F)
	dat <- read.table("/public/workspace/huangbin/MPAL/figure4/surv/MPAL-allGeneExpr-90S-rlog.txt",sep="\t",stringsAsFactors=F,header=T)
	dat <- dat%>%filter(!Gene%in%dat$Gene[duplicated(dat$Gene)])
	rownames(dat) <- dat$Gene
	dat <- dat[,-c(1,2)]
	dat <- dat[,-which(colnames(dat)=="SJMPAL011911_D2")]
	colnames(dat) <- as.character(sapply(colnames(dat),function(x){strsplit(x,"_")[[1]][1]}))
	dat <- dat[,which(colnames(dat)%in%clinic$St..Jude.Genome.ID)]
	clinic$Vital.Status[which(clinic$Vital.Status=="Alive")] <- 0
	clinic$Vital.Status[which(clinic$Vital.Status=="Dead")] <- 1
	clinic$Vital.Status <- as.numeric(clinic$Vital.Status)
	clinic$Overall.Survival.Time.in.Days <- clinic$Overall.Survival.Time.in.Days/30

	load("p.RData")
	clinic$type <- "unknown"
	clinic$type[which(clinic[,1]%in%names(which(cutree(p$tree_col,3)==1)))] <- "TEC"
	clinic$type[which(clinic[,1]%in%names(which(cutree(p$tree_col,3)==2)))] <- "AEC"
	clinic$type[which(clinic[,1]%in%names(which(cutree(p$tree_col,3)==3)))] <- "S24"

	# fig3L
	fit <- survfit(Surv(Overall.Survival.Time.in.Days, Vital.Status)~type,data=clinic)
	p <- ggsurvplot(fit,
	     pval = TRUE, 
	     risk.table = TRUE, # Add risk table
	     # risk.table.col = "strata", # Change risk table color by groups
	     linetype = "strata", # Change line type by groups
	     # surv.median.line = "hv", # Specify median survival
	     ggtheme = theme_bw(), # Change ggplot2 theme
	     palette = c("#007cc0","#ffc20e","#ce181e")
	 )
	pdf("surv.pdf")
	print(p)
	dev.off()

	# 使用 pairwise log-rank 检验
	pairwise_p <- pairwise_survdiff(
	Surv(Overall.Survival.Time.in.Days, Vital.Status) ~ type,
		data = clinic,
		p.adjust.method = "none"  # 不做多重校正，按原始P值返回
	)

	# target - B/My MPAL
	library(dplyr)
	library(pheatmap)
	library(RColorBrewer)
	setwd("/public/workspace/huangbin/MPAL/figure3/bulk/target/BMMPAL")
	dat <- read.table("MPAL-allGeneExpr-90S-rlog.txt",sep="\t",stringsAsFactors=F,header=T)
	dat <- dat%>%filter(!Gene%in%dat$Gene[duplicated(dat$Gene)]);rownames(dat) <- dat$Gene
	dat <- dat[,-c(1,2)]; dat <- dat[,-which(colnames(dat)=="SJMPAL011911_D2")]
	colnames(dat) <- as.character(sapply(colnames(dat),function(x){strsplit(x,"_")[[1]][1]}))

	info <- read.table("target.clinic.txt",sep="\t",stringsAsFactors=F,header=T)
	dat <- dat[,intersect(info$Sample,colnames(dat))] #24

	markers <- read.table("/public/workspace/huangbin/MPAL/figure3/TALL.AML.leukemia/labeltransfer/marker.type7.txt",
	  sep="\t",header=T,stringsAsFactors=F)
	top <- markers%>%filter(cluster%in%c("S24","AEC"))%>%
	group_by(cluster)%>%top_n(8,wt=avg_log2FC)
	select <- intersect(top$gene,rownames(dat))

	# fig3K
	pdf("Target.heatmap.pdf",height=3)
	p <- pheatmap(dat[select,], 
	  breaks=seq(-3.5,3.5,0.001),
	  scale="row",border_color=NA,
	  col=rev(colorRampPalette(brewer.pal(11,"RdBu"))(length(seq(-3.5,3.5,0.001)))),
	  show_rownames=TRUE,cluster_rows=TRUE,clustering_method="complete",
	  cutree_cols=3,cluster_cols=TRUE,show_colnames=FALSE,fontsize=8)
	print(p)
	dev.off()

	# target
	.libPaths("~/plaqueRlibrary/")
	library(dplyr)
	library(tidyr)
	library(survival)
	library(survminer)
	setwd("/public/workspace/huangbin/MPAL/figure3/bulk/target/BMMPAL")
	dat <- read.table("/public/workspace/huangbin/MPAL/figure4/surv/MPAL-allGeneExpr-90S-rlog.txt",sep="\t",stringsAsFactors=F,header=T)
	dat <- dat%>%filter(!Gene%in%dat$Gene[duplicated(dat$Gene)])
	rownames(dat) <- dat$Gene
	dat <- dat[,-c(1,2)]
	dat <- dat[,-which(colnames(dat)=="SJMPAL011911_D2")]
	colnames(dat) <- as.character(sapply(colnames(dat),function(x){strsplit(x,"_")[[1]][1]}))
	
	clinic <- read.table("./surv.txt",sep="\t",header=T,stringsAsFactors=F)
	colnames(clinic) <- c("Sample","Type","OS","OS.time")
	clinic <- clinic%>%filter(Sample%in%colnames(dat))

	dat <- dat[,clinic$Sample]

	clinic$OS[which(clinic$OS=="alive")] <- 0
	clinic$OS[which(clinic$OS=="death")] <- 1
	clinic$OS <- as.numeric(clinic$OS)
	clinic$OS.time <- clinic$OS.time/30
	clinic$HOPX <- as.numeric(dat["HOPX",])
	clinic$type <- ifelse(clinic$HOPX>median(as.numeric(clinic$HOPX)),"high","low")

	# fig3L
	fit <- survfit(Surv(OS.time,OS)~type,data=clinic)
	p <- ggsurvplot(fit,
	     pval = TRUE, 
	     risk.table = TRUE, # Add risk table
	     # risk.table.col = "strata", # Change risk table color by groups
	     linetype = "strata", # Change line type by groups
	     # surv.median.line = "hv", # Specify median survival
	     ggtheme = theme_bw(), # Change ggplot2 theme
	     palette = c("#007cc0","#ffc20e","#ce181e")
	 )
	pdf("./BMMPAL/surv.pdf")
	print(p)
	dev.off()
}

################################################################
##                           Figure4                          ##
################################################################
{
	# sc2marker.txt
	bytlib load r/4.0.4
	.libPaths("~/plaqueRlibrary")
	library(sc2marker)
	library(Seurat)
	library(ggplot2)
	library(dplyr)
	library(RColorBrewer)
	setwd("/public/workspace/huangbin/MPAL/figure4/sc2marker")
	dat <- readRDS("../../figure3/TALL.AML.leukemia/labeltransfer/MPAL.pred.rds")
	dat <- subset(dat,cells=rownames(dat@meta.data)[which(dat$info%in%c("TEC","AEC","S24"))])
	DefaultAssay(dat) <- "RNA"
	Idents(dat) <- factor(dat$info,levels=c("TEC","AEC","S24"))
	markers <- lapply(c("TEC","AEC","S24"),function(i){
		Flow.markers <- Detect_single_marker(dat,id=i,category="Flow",org="human")%>%filter(direction=="+")
		Flow.markers$detect <- Flow.markers$TP/length(which(dat$info==i))
		Flow.markers$accuracy <- (Flow.markers$TP+Flow.markers$TN)/ncol(dat)
		Flow.markers <- Flow.markers%>%filter(detect>0.2)%>%arrange(desc(accuracy))
		Flow.markers$Case <- i; Flow.markers$Type <- "Flow"

		pdf(paste0("Markers.",i,".Flow.pdf"))
		plot_ridge(dat,id=i, genes=Flow.markers[1:2,"gene"], ncol=2,aggr.other=F)
		dev.off()

		IHC.markers <- Detect_single_marker(dat,id=i,category="IHC",org="human")%>%filter(direction=="+")
		IHC.markers$detect <- IHC.markers$TP/length(which(dat$info==i))
		IHC.markers$accuracy <- (IHC.markers$TP+IHC.markers$TN)/ncol(dat)
		IHC.markers <- IHC.markers%>%filter(detect>0.2)%>%arrange(desc(accuracy))
		IHC.markers$Case <- i; IHC.markers$Type <- "IHC"

		pdf(paste0("Markers.",i,".IHC.pdf"))
		plot_ridge(dat,id=i, genes=IHC.markers[1:2,"gene"], ncol=2,aggr.other=F)
		dev.off()

		rst <- rbind(Flow.markers,IHC.markers)
		return(rst)
	})
	markers <- do.call(rbind,markers)
	write.table(markers,"sc2marker.txt",sep="\t",quote=F,row.names=F)

	# fig4A
	library(ggplot2)
	library(ggridges)
	library(viridis)
	tmp <- as.data.frame(t(as.matrix(dat@assays$RNA@data[c("HSPB1","HOPX","CD7","IGLL1","S100A9","S100A8"),])))
	tmp$group <- dat$info
	cutoff <- c(0.8417086,0.9438072,1.5367699,1.5375140,2.3396908,2.3004346)
	for(i in 1:6){
	tmp2 <- data.frame(Value=tmp[,i],Group=tmp[,7])
	p <- ggplot(tmp2) +
		  geom_density_ridges_gradient(aes(x=Value,y=Group,fill=stat(x)),scale=2,rel_min_height=0,size=0.3)+
		  scale_fill_gradientn(colours=colorRampPalette(rev(brewer.pal(11,'Spectral')))(32))+
		  theme_ridges(grid=FALSE) + 
		  theme(legend.position = "top",#图例去除
		        axis.line.y = element_blank(),#去除Y轴轴线
		        axis.ticks.y = element_blank())+#去除Y轴刻度
		  geom_vline(xintercept=cutoff[i], colour="#990000", linetype="dashed")
	assign(paste0("p",i),p)
	}
	pdf("ridge.pdf",width=30)
	cowplot::plot_grid(p1,p2,p3,p4,p5,p6,ncol=6) 
	dev.off()

	# fig4B
	.libPaths("~/plaqueRlibrary")
	library(ROCR)
	library(pROC)
	dat <- readRDS("../../figure3/TALL.AML.leukemia/labeltransfer/MPAL.pred.rds")
	dat <- subset(dat,cells=rownames(dat@meta.data)[which(dat$info%in%c("TEC","AEC","S24"))])
	DefaultAssay(dat) <- "RNA"
	dat$CD7 <- ifelse(as.numeric(dat@assays$RNA@data["CD7",])>0,1,0)
	dat$IGLL1 <- ifelse(as.numeric(dat@assays$RNA@data["IGLL1",])>0,1,0)
	dat$S100A8 <- ifelse(as.numeric(dat@assays$RNA@data["S100A8",])>0,1,0)
	dat$S100A9 <- ifelse(as.numeric(dat@assays$RNA@data["S100A9",])>0,1,0)
	dat$HOPX <- ifelse(as.numeric(dat@assays$RNA@data["HOPX",])>0,1,0)
	dat$HSPB1 <- ifelse(as.numeric(dat@assays$RNA@data["HSPB1",])>0,1,0)
	tmp <- dat@meta.data[,c("info",c("HSPB1","HOPX","CD7","IGLL1","S100A8","S100A9"))]
	tmp <- tmp%>%group_by(info)%>%summarise_all(sum)
	tmp[,-1] <- tmp[,-1]/t(c(3357,24465,21712))
	tmp <- as.data.frame(tmp)
	rownames(tmp) <- tmp$info
	tmp <- tmp[,-1]
	pdf("bar.marker.pdf")
	barplot(as.matrix(tmp),beside=T,ylim=c(0,1),
	    col=c("#007cc0","#ffc20e","#ce181e"))
	dev.off()

	# fig4D
	# IHC analysis 
	library(ggplot2)
	library(dplyr)
	dat <- read.table("IHCanalysis.txt",sep="\t",stringsAsFactors=F,header=T)
	dat$Image <- gsub(" 2023.*","",dat$Image)
	colnames(dat) <- c("sample","Num.Detections","Num.Negative","Num.Positive",
		"pfreq","Num.Positive.per.mm","Area","Type")
	write.table(dat,"selectIHC.txt",sep="\t",quote=F,row.names=F)
	pdf("panIHC.filter.pdf",width=8)
	dat$pfreq <- as.factor(dat$pfreq)
	dat$pfreq <- as.numeric(as.character(dat$pfreq))
	dat$pfreq <- round(dat$pfreq,1)
	dat$Type <- factor(as.factor(dat$Type),levels=c("S100A9","CD7","HOPX"))
	tmp1 <- dat%>%filter(Type=="S100A9")%>%arrange(sample);colnames(tmp1)<- paste0("S100A9.",colnames(tmp1))
	tmp2 <- dat%>%filter(Type=="CD7")%>%arrange(sample);colnames(tmp2)<- paste0("CD7.",colnames(tmp2))
	tmp <- cbind(tmp1,tmp2)%>%as.data.frame()
	tmp$delta <- (tmp[,5]-tmp[,13])
	tmp <- tmp%>%arrange(desc(delta))
	dat$sample <- factor(as.factor(dat$sample),levels=tmp$S100A9.sample)
	ggplot(data=dat,aes(x=sample,y=pfreq))+
	geom_bar(aes(fill=Type),stat='identity',position='dodge')+
	scale_fill_manual(values=c("#007cc0","#ce181e","#ffc20e"))+
	geom_text(aes(label=pfreq),size=5,position='dodge')+
	# geom_point(aes(x=Image.Tag,y=阳性细胞占比))+
	theme_classic()
	dev.off()
}

################################################################
##                           Figure5                          ##
################################################################
{	
	options(future.globals.maxSize = 2400000 * 1024^2)
	.libPaths(c('/public/workspace/huangbin/plaqueRlibrary'))
	library(Seurat)
	library(dplyr)
	library(reshape2)
	library(RColorBrewer)
	library(pheatmap)
	library(ggplot2)
	library(monocle)
	library(BiocGenerics)
	load("cds2.RData")
	info <- readRDS("../TALL.AML.leukemia/labeltransfer/MPAL.pred.rds")
	info <- data.frame(cell=rownames(info@meta.data),info=info$info)
	pData(cds)$group <- pData(cds)$sampletype
	pData(cds)$group[which(rownames(pData(cds))%in%info$cell[which(info$info=="S24")])] <- "S24"
	pData(cds)$group[which(rownames(pData(cds))%in%info$cell[which(info$info=="AEC")])] <- "AEC"
	pData(cds)$group[which(rownames(pData(cds))%in%info$cell[which(info$info=="TEC")])] <- "TEC"
	pData(cds)$group[which(rownames(pData(cds))%in%info$cell[which(info$info=="S13")])] <- "S13"
	pData(cds)$group[which(rownames(pData(cds))%in%info$cell[which(info$info=="S14")])] <- "S14"
	pData(cds)$group[which(rownames(pData(cds))%in%info$cell[which(info$info=="S16")])] <- "S16"
	pData(cds)$group[which(rownames(pData(cds))%in%info$cell[which(info$info=="S19")])] <- "S19"
	pData(cds)$group <- factor(as.factor(pData(cds)$group),levels=c("TALL","TEC","S24","S13","S14","S16","S19","AEC","AML"))

	# fig5A
	# pesudotime.heatmap
	marker <- read.table("../TALL.AML.leukemia/labeltransfer/marker.type7.txt",sep="\t",stringsAsFactors=F,header=T)
	top <- marker%>%filter(cluster%in%c("AEC","S24","TEC"))%>%group_by(cluster)%>%top_n(n=30,wt=avg_log2FC)
	top <- c(top$gene,"CD34")
	pdf("pesudotime.heatmap.pdf")
	plot_genes_branched_heatmap(cds[top,] ,
	    branch_point=2,num_clusters=3,
	    cores=1,branch_labels=c("TALL","AML"),
	    hmcols=colorRampPalette(rev(brewer.pal(9,"RdBu")))(62),
	    branch_colors=c("#a7a9ac","#075aaa","#eb2226"),
	    use_gene_short_name=T,show_rownames=T)
	dev.off()

	# fig5B
	# pesudotime.marker
	marker <- c("S100A9","HOPX","CD7","S100A8","CD34","DNTT")
	pdf("pesudotime.marker.pdf",useDingbats=F,width=8,height=4)
	plot_genes_branched_pseudotime(cds[marker,] ,panel_order=marker,
	    branch_point=2,ncol=3,cell_size=1,color_by="group")+
	scale_color_manual(values=c("TALL"="#ffaaaa",'TEC'="#ce181e",
	    "S13"="grey60","S14"="grey60","S16"="grey60","S19"="grey60",
	    "S24"="#ffc20e",'AEC'="#007cc0","AML"="#b3dcff"))
	dev.off()

	# fig5A
	library(clusterProfiler)
	library(org.Hs.eg.db)
	library(ggplot2)
	library(RColorBrewer)
	library(dplyr)
	gene <- rownames(p$annotation_row)[which(p$annotation_row$Cluster==2)]
	eg <- bitr(gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL",'SYMBOL'), OrgDb="org.Hs.eg.db")
	go <- enrichGO(eg$ENTREZID, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
	             qvalueCutoff = 0.2,keyType = 'ENTREZID')

	# fig5C
	library(Seurat)
	library(ggplot2)
	library(dplyr)
	setwd("/public/workspace/huangbin/MPAL/figure4/CD34")
	dat <- readRDS("../../figure3/TALL.AML.leukemia/labeltransfer/MPAL.pred.rds")
	dat <- subset(dat,cells=rownames(dat@meta.data)[which(dat$info%in%c("TEC","AEC","S24"))])
	DefaultAssay(dat) <- "RNA"
	dat$CD34 <- as.numeric(dat@assays$RNA@data["CD34",])
	freq <-  rev(table(dat$info[which(dat$CD34>0)])/table(dat$info))
	pdf("CD34.pdf")
	barplot(rev(sort(freq)),width=0.1,ylim=c(0,0.5),col="#4298b5",border=NA)
	dev.off()

	# fig5C-补充干性分析 stemness index 
	setwd("/public/workspace/huangbin/MPAL/figure4/CD34")
	.libPaths("~/plaqueRlibrary/")
	library(Seurat); library(dplyr); library(ggplot2); library(ggpubr)
	source('/public/workspace/huangbin/ssGSEA/ssgseaMOD.r')
		dat <- readRDS("/public/workspace/huangbin/MPAL/figure1/Integrate/Integrate.all.rds")
	DefaultAssay(dat) <- "RNA"
    genes <- c(
        "ANKRD28", "BCL2L2", "BMP1", "C14ORF79", "C1ORF149", "CEP290", "CHKA", "CLIC1", "CMAH", "CUTL1",
        "DPAGT1", "DUSP3", "EEF2", "EFHC1", "EIF2S3", "ELMO1", "FLT3", "GCN5L2", "GNA15", "GNPTAB",
        "GRHPR", "HSBP1", "ITGA5", "ITGA9", "KDELR1", "LEPRE1", "LIAS", "LMAN2L", "MAPKAPK3", "MRP63",
        "MRPS30", "PHKB", "PLEKHA9", "PTTG1IP", "RAB6IP1", "RAI1", "RPS29", "SDCCAG8", "SMYD3", "SNRPD3",
        "SPATA20", "SYF2", "TEC", "TFPI", "TPM4", "UBE2E1", "UGP2", "XPNPEP1", "ZNF254", "ZNF551"
        )
    dat <- AddModuleScore(dat,list(genes))
	tmp <- dat@meta.data%>%filter(info%in%c("S24","TEC","AEC"))
	pdf("./StemnessIndex.activity.pdf",useDingbats=F,height=5,width=5)
	ggplot(data=tmp,aes(x=info,y=Cluster1))+
	geom_violin(aes(fill=info))+
	scale_fill_manual(values=c('TEC'="#ce181e","S24"="#ffc20e",'AEC'="#007cc0"))+
	geom_boxplot(fill="white",outlier.shape=NA,width=0.2)+
	 geom_signif(test=wilcox.test,comparisons=list(c("S24","TEC"),c("S24","AEC")))+
	theme_classic()
	dev.off()

	# GSEA
	options(future.globals.maxSize = 2400000 * 1024^2)
	.libPaths("~/plaqueRlibrary/")
	library(org.Hs.eg.db)
	library(clusterProfiler)
	library(msigdf)
	library(tidyverse)
	library(pathview)
	library(Seurat)
	library(GseaVis)
	library(future)
	plan("multiprocess",workers=10)

	setwd("/public/workspace/huangbin/MPAL/figure4/GSEA/")
	dat <- readRDS("../../figure3/TALL.AML.leukemia/labeltransfer/MPAL.pred.rds")
	dat <- subset(dat,cells=rownames(dat@meta.data)[which(dat$info%in%c("TEC","AEC","S24"))])
	DefaultAssay(dat) <- "RNA"
	genes = rownames(dat)[rowSums(cbind(as.matrix(dat@assays$RNA@counts[,1:40000]),
	  as.matrix(dat@assays$RNA@counts[,40001:49534])) > 0) > 3]
	dat <- subset(dat, features=genes)
	markers <- presto::wilcoxauc(dat, 'info')   ## presto的方法 
	saveRDS(markers,file="GSEA.markers.rds")

	# fig5H
	# TEC
	x <- readLines("HP_ACUTE_LYMPHOBLASTIC_LEUKEMIA.v2022.1.Hs.gmt")
	res <- strsplit(x, "\t")
	names(res) <- vapply(res, function(y) y[1], character(1))
	res <- lapply(res, "[", -c(1:2))
	term <- data.frame(geneset="HP_ACUTE_LYMPHOBLASTIC_LEUKEMIA",symbol=res)
	pdf("./HP_ACUTE_LYMPHOBLASTIC_LEUKEMIA.pdf",useDingbats=FALSE,width=5,height=5)
	set.seed(12345)
	genes <- markers %>% dplyr::filter(group == "TEC") %>% arrange(desc(auc)) %>% dplyr::select(feature, auc)
	ranks <- tibble::deframe(genes)
	GSEA_res <- clusterProfiler::GSEA(ranks, TERM2GENE=term, pvalueCutoff = 1, maxGSSize = 1500,
	  nPermSimple=10000, eps=0)  
	p <- GseaVis::gseaNb(object = GSEA_res,geneSetID ="HP_ACUTE_LYMPHOBLASTIC_LEUKEMIA",
	  addPval = T,pvalX = 0.75,pvalY = 0.75,pCol = 'black',pHjust = 0,subPlot = 2) + 
	ggtitle("HP_ACUTE_LYMPHOBLASTIC_LEUKEMIA")
	print(p)
	dev.off()

	# fig5H
	# aTMMPAL
	x <- readLines("HP_ACUTE_MYELOID_LEUKEMIA.v2022.1.Hs.gmt")
	res <- strsplit(x, "\t")
	names(res) <- vapply(res, function(y) y[1], character(1))
	res <- lapply(res, "[", -c(1:2))
	term <- data.frame(geneset="HP_ACUTE_MYELOID_LEUKEMIA",symbol=res)
	pdf("HP_ACUTE_MYELOID_LEUKEMIA.pdf",useDingbats=FALSE,width=5,height=5)
	set.seed(12345)
	genes <- markers %>% dplyr::filter(group == "AEC") %>% arrange(desc(auc)) %>% dplyr::select(feature, auc)
	ranks <- tibble::deframe(genes)
	GSEA_res <- clusterProfiler::GSEA(ranks, TERM2GENE=term, pvalueCutoff = 1, maxGSSize = 1500,
	  nPermSimple=10000, eps=0)
	p <- GseaVis::gseaNb(object = GSEA_res,geneSetID ="HP_ACUTE_MYELOID_LEUKEMIA",
	  addPval = T,pvalX = 0.75,pvalY = 0.75,pCol = 'black',pHjust = 0,subPlot = 2) + 
	ggtitle("HP_ACUTE_MYELOID_LEUKEMIA")
	print(p)
	dev.off()

	# fig5H
	# G0/diapause/Quiescent 
	x <- readLines("GRAHAM_NORMAL_QUIESCENT_VS_NORMAL_DIVIDING_UP.v7.5.1.gmt")
	res <- strsplit(x, "\t")
	names(res) <- vapply(res, function(y) y[1], character(1))
	res <- lapply(res, "[", -c(1:2))
	term <- data.frame(geneset="GRAHAM_NORMAL_QUIESCENT_VS_NORMAL_DIVIDING_UP",symbol=res)
	pdf("GRAHAM_NORMAL_QUIESCENT_VS_NORMAL_DIVIDING_UP.pdf",useDingbats=FALSE,width=5,height=3)
	set.seed(12345)
	genes <- markers %>% dplyr::filter(group == "S24") %>% arrange(desc(auc)) %>% dplyr::select(feature, auc)
	ranks <- tibble::deframe(genes)
	GSEA_res <- clusterProfiler::GSEA(ranks, TERM2GENE=term, pvalueCutoff = 1, maxGSSize = 1500,
	  nPermSimple=10000, eps=0)  
	p <- GseaVis::gseaNb(object = GSEA_res,geneSetID ="GRAHAM_NORMAL_QUIESCENT_VS_NORMAL_DIVIDING_UP",
	  addPval = T,pvalX = 0.75,pvalY = 0.75,pCol = 'black',pHjust = 0,subPlot = 2) + 
	ggtitle("GRAHAM_NORMAL_QUIESCENT_VS_NORMAL_DIVIDING_UP")
	print(p)
	dev.off()

	# pysenic
	.libPaths("~/plaqueRlibrary/")
	library(dplyr)
	library(Seurat)
	library(data.table)
	library(ggplot2)
	dat <- readRDS("/public/workspace/huangbin/MPAL/figure3/TALL.AML.leukemia/labeltransfer/MPAL.pred.rds")
	dat <- subset(dat,cells=rownames(dat@meta.data)[which(dat$info%in%c("TEC","AEC","S24"))])
	DefaultAssay(dat) <- "RNA"
	filter <- rownames(dat)[rowSums(cbind(as.matrix(dat@assays$RNA@counts[,1:40000]),
	  as.matrix(dat@assays$RNA@counts[,40001:49534]))>0)>3]
	dat <- subset(dat,features=filter)
    counts <- t(as.matrix(dat@assays$RNA@counts))
    write.csv(counts,"./counts.csv")

    module load python-3.6.6
    module load pyscenic

    REF_PATH=/public/workspace/huangbin/scenic_ref/human
    OUT_PATH=/public/workspace/huangbin/MPAL/revise1/Q18/pyscenic3
    Thread=5

    NUMBA_THREADING_LAYER='omp' pyscenic grn ${OUT_PATH}/counts.csv ${REF_PATH}/hs_hgnc_curated_tfs.txt \
    -o ${OUT_PATH}/adjacencies.tsv \
    --num_workers $Thread \
    --seed 23

    pyscenic ctx ${OUT_PATH}/adjacencies.tsv \
        ${REF_PATH}/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather ${REF_PATH}/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather \
        --annotations_fname ${REF_PATH}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
        --expression_mtx_fname ${OUT_PATH}/counts.csv \
        --output ${OUT_PATH}/regulons.tsv \
        --num_workers $Thread \
        --mode "custom_multiprocessing"

    pyscenic aucell \
        ${OUT_PATH}/counts.csv \
        ${OUT_PATH}/regulons.tsv \
        --output ${OUT_PATH}/auc_mtx.csv \
        --num_workers $Thread \
        --seed 23
      dat <- readRDS("../../figure3/TALL.AML.leukemia/labeltransfer/MPAL.pred.rds")
	dat <- subset(dat,cells=rownames(dat@meta.data)[which(dat$info%in%c("TEC","AEC","S24"))])
	DefaultAssay(dat) <- "RNA"
    S24.index <- which(dat$info=="S24")
	AEC.index <- which(dat$info=="AEC")
	TEC.index <- which(dat$info=="TEC")

	## uTMMPAL vs. tTMMPAL
    rst <- read.csv("./pyscenic3/auc_mtx.csv",stringsAsFactors=F,header=T,check.names=F,row.names=1)
    colnames(rst) <- gsub("\\(\\+\\)","",colnames(rst))
	p.activity <- apply(rst,2,function(x){
		tmp <- wilcox.test(as.numeric(x[S24.index]),as.numeric(x[TEC.index]))[[3]]
		return(tmp)
		})
	fc.activity <- apply(rst,2,function(x){
		tmp <- mean(as.numeric(x[S24.index]))/mean(as.numeric(x[TEC.index]))
		return(tmp)
		})
	p.expression <- apply(dat@assays$RNA@data[intersect(colnames(rst),rownames(dat)),],1,function(x){
		tmp <- wilcox.test(as.numeric(x[S24.index]),as.numeric(x[TEC.index]))[[3]]
		return(tmp)
		})
	fc.expression <- apply(dat@assays$RNA@data[intersect(colnames(rst),rownames(dat)),],1,function(x){
		tmp <- mean(as.numeric(x[S24.index]))/mean(as.numeric(x[TEC.index]))
		return(tmp)
		})
	# 不关注基因表达的变化只要有表达就可以
	filter <- intersect(names(which(p.activity<0.05)),names(which(is.na(p.expression)==FALSE)))
	rst <- data.frame(Gene=filter,delta.activity=as.numeric(fc.activity[filter]),delta.expression=as.numeric(fc.expression[filter]))
	write.table(rst,"DTF.ut.txt",sep="\t",quote=F,row.names=F)

	## uTMMPAL vs. aTMMPAL
    rst <- read.csv("./pyscenic3/auc_mtx.csv",stringsAsFactors=F,header=T,check.names=F,row.names=1)
    colnames(rst) <- gsub("\\(\\+\\)","",colnames(rst))
	p.activity <- apply(rst,2,function(x){
		tmp <- wilcox.test(as.numeric(x[S24.index]),as.numeric(x[AEC.index]))[[3]]
		return(tmp)
		})
	fc.activity <- apply(rst,2,function(x){
		tmp <- mean(as.numeric(x[S24.index]))/mean(as.numeric(x[AEC.index]))
		return(tmp)
		})
	p.expression <- apply(dat@assays$RNA@data[intersect(colnames(rst),rownames(dat)),],1,function(x){
		tmp <- wilcox.test(as.numeric(x[S24.index]),as.numeric(x[AEC.index]))[[3]]
		return(tmp)
		})
	fc.expression <- apply(dat@assays$RNA@data[intersect(colnames(rst),rownames(dat)),],1,function(x){
		tmp <- mean(as.numeric(x[S24.index]))/mean(as.numeric(x[AEC.index]))
		return(tmp)
		})
	# 不关注基因表达的变化只要有表达就可以
	filter <- intersect(names(which(p.activity<0.05)),names(which(is.na(p.expression)==FALSE)))
	rst <- data.frame(Gene=filter,delta.activity=as.numeric(fc.activity[filter]),delta.expression=as.numeric(fc.expression[filter]))
	write.table(rst,"DTF.ua.txt",sep="\t",quote=F,row.names=F)
    
	options(future.globals.maxSize = 2400000 * 1024^2)
	.libPaths("~/plaqueRlibrary/")
	setwd("/public/workspace/huangbin/MPAL/figure5/hallmark/")
	library(clusterProfiler)
	library(org.Hs.eg.db)
	library(ggplot2)
	library(RColorBrewer)
	library(dplyr)
	library(msigdf)
	library(Seurat)
	# S24 vs. TEC GSEA
	dat <- readRDS("../../figure3/TALL.AML.leukemia/labeltransfer/MPAL.pred.rds")
	dat <- subset(dat,cells=rownames(dat@meta.data)[which(dat$info%in%c("TEC","S24"))])
	DefaultAssay(dat) <- "RNA"
	genes = rownames(dat)[rowSums(cbind(as.matrix(dat@assays$RNA@counts[,1:40000]),
	  as.matrix(dat@assays$RNA@counts[,40001:46177])) > 0) > 3]
	dat <- subset(dat, features=genes)
	markers <- presto::wilcoxauc(dat, 'info')  
	saveRDS(markers,file="GSEA.S24.TEC.rds")
	genes <- markers %>% dplyr::filter(group == "S24") %>% arrange(desc(auc)) %>% dplyr::select(feature, auc)
	ranks <- tibble::deframe(genes)
	term <- msigdf.human%>%filter(category_code=="h")%>%select(geneset,symbol)%>%as.data.frame()
	set.seed(123)
	res.T <- clusterProfiler::GSEA(ranks, TERM2GENE=term, pvalueCutoff = 1, maxGSSize = 1500,
	  nPermSimple=10000, eps=0) 
	# S24 vs. AEC GSEA
	dat <- readRDS("../../figure3/TALL.AML.leukemia/labeltransfer/MPAL.pred.rds")
	dat <- subset(dat,cells=rownames(dat@meta.data)[which(dat$info%in%c("AEC","S24"))])
	DefaultAssay(dat) <- "RNA"
	genes = rownames(dat)[rowSums(as.matrix(dat@assays$RNA@counts) > 0) > 3]
	dat <- subset(dat, features=genes)
	markers <- presto::wilcoxauc(dat, 'info')  
	saveRDS(markers,file="GSEA.S24.AEC.rds")
	genes <- markers %>% dplyr::filter(group == "S24") %>% arrange(desc(auc)) %>% dplyr::select(feature, auc)
	ranks <- tibble::deframe(genes)
	term <- msigdf.human%>%filter(category_code=="h")%>%select(geneset,symbol)%>%as.data.frame()
	set.seed(123)
	res.A <- clusterProfiler::GSEA(ranks, TERM2GENE=term, pvalueCutoff = 1, maxGSSize = 1500,
	  nPermSimple=10000, eps=0) 

	# fig5I
	tmp <- cbind(res.T@result[,c("ID","NES","qvalues")]%>%arrange(ID),
		res.A@result[,c("ID","NES","qvalues")]%>%arrange(ID))
	tmp <- tmp[,c(2,3,5,6)]
	colnames(tmp) <- c("NES.T","q.T","NES.A","q.A")
	tmp <- tmp[which((tmp$NES.T>0 & tmp$NES.A>0)|(tmp$NES.T<0 & tmp$NES.A<0)),]
	tmp <- tmp[which(tmp$q.T<0.05 & tmp$q.A<0.05),]
	t.p <- tmp[,c(2,4)]
	t.p[t.p<0.001] <- "**"
	t.p[t.p<0.01 & t.p>=0.001] <- "**"
	t.p[t.p<0.05 & t.p>=0.01] <- "*"
	t.p[t.p>0.05] <- "ns"
	pdf("GSEA.Heatmap.pdf")
	pheatmap(tmp[,c(1,3)], 
	  breaks=seq(-3,3,0.001),
	  col=colorRampPalette(rev(brewer.pal(11,"RdBu")[3:9]))(length(seq(-3,3,0.001))),
	  show_rownames=TRUE,cluster_rows=TRUE,
	  cluster_cols=TRUE,show_colnames=TRUE,
	   border_color="black",fontsize=8,
	   display_numbers = t.p)
	dev.off()

	# fig5M
	.libPaths("~/plaqueRlibrary/")
	library(dplyr)
	library(tidyr)
	library(survival)
	library(survminer)
	library(ssgsea.GBM.classification)
	setwd("/public/workspace/huangbin/MPAL/figure3/bulk/target/")
	marker <- read.table("../../TALL.AML.leukemia/labeltransfer/marker.type7.txt",sep="\t",stringsAsFactors=F,header=T)
	clinic <- read.table("/public/workspace/huangbin/MPAL/figure4/surv/clinic.txt",sep="\t",header=T,stringsAsFactors=F)
	dat <- read.table("/public/workspace/huangbin/MPAL/figure4/surv/MPAL-allGeneExpr-90S-rlog.txt",sep="\t",stringsAsFactors=F,header=T)
	dat <- dat%>%filter(!Gene%in%dat$Gene[duplicated(dat$Gene)])
	rownames(dat) <- dat$Gene
	dat <- dat[,-c(1,2)]
	dat <- dat[,-which(colnames(dat)=="SJMPAL011911_D2")]
	colnames(dat) <- as.character(sapply(colnames(dat),function(x){strsplit(x,"_")[[1]][1]}))
	dat <- dat[,which(colnames(dat)%in%clinic$St..Jude.Genome.ID)]
	clinic$Vital.Status[which(clinic$Vital.Status=="Alive")] <- 0
	clinic$Vital.Status[which(clinic$Vital.Status=="Dead")] <- 1
	clinic$Vital.Status <- as.numeric(clinic$Vital.Status)
	clinic$Overall.Survival.Time.in.Days <- clinic$Overall.Survival.Time.in.Days/30
	dat <- dat[,clinic[,1]]
	# GSVA (CD74 HOPX HLA-DRA LRRC75A SPINK2 HLA-DRB1)
	library(GSVA)
	set.seed(666)
	select.marker <- marker$gene[which(marker$cluster=="S24")][c(1:10)]
	score <- gsva(as.matrix(dat),list(score=select.marker),method='ssgsea')
	clinic$score <- as.numeric(score)
	clinic$type <- ifelse(clinic$score<median(clinic$score),"low","high")
	# survival
	fit <- survfit(Surv(Overall.Survival.Time.in.Days, Vital.Status)~type,data=clinic)
	p <- ggsurvplot(fit,
	     pval = TRUE,
	      ylim = c(0.5,1),
	     risk.table = TRUE, # Add risk table
	     #  risk.table.col = "strata", # Change risk table color by groups
	     # linetype = "strata", # Change line type by groups
	     #  surv.median.line = "hv", # Specify median survival
	     ggtheme = theme_bw(), # Change ggplot2 theme
	     palette = c("#007cc0","#ffc20e","#ce181e")
	 )
	pdf("surv2.pdf")
	print(p)
	dev.off()

    # fig5N
    # AJH - HOPX
    library(ggpubr)
    library(ggplot2)
    library(openxlsx)
    setwd("/public/workspace/huangbin/MPAL/revise1/Q16")
    dat <- read.table("../Q13/AJH/cibersortx.AJH.txt",sep="\t",stringsAsFactors=F,header=T,row.names=1)
    dat <- data.frame(SampleID=colnames(dat),value=as.numeric(dat["HOPX",]))
    info <- read.xlsx("../Q13/AJH/survival_outcomes.xlsx")%>%as.data.frame()
    dat <- left_join(dat,info,by="SampleID")%>%filter(WHO.MPAL.classification=="T/M",relapse!="--")
    dat$relapse <- factor(dat$relapse,levels=c("Yes","No"))
    pdf("relapse.AJH.pdf",width=5,height=5)
        p <- ggplot(dat, aes(x = relapse, y = value, fill = relapse)) +
                geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
                geom_jitter(width = 0.15, alpha = 0.8, size = 2.5, shape = 21, stroke = 1) +
                stat_compare_means(method = "wilcox.test", 
                                    # label = "p.signif",   # 显示星号 * 形式
                                    label.y = max(dat$value) * 1.05, 
                                    size = 5) +
                scale_fill_brewer(palette = "Set1") +
                labs(
                    x = "Relapse",
                    y = "HOPX Expression",
                    title = "AJH"
                ) +
                theme_minimal(base_size = 14) +
                theme(
                    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
                    axis.title = element_text(face = "bold"),
                    axis.text = element_text(color = "black"),
                    legend.position = "none"
                )
        print(p)
    dev.off()

}

################################################################
##                           Figure6                          ##
################################################################
{
	# fig6A-C
	.libPaths("~/plaqueRlibrary")
	library(Seurat)
	library(future)
	library(dplyr)
	library(parallel)
	library(data.table)
	library(readxl)
	library(ssgsea.GBM.classification)
	library(ggplot2)
	plan("multiprocess",workers=10)
	options(future.globals.maxSize = 2400000 * 1024^2)
	setwd("/public/workspace/huangbin/MPAL/figure6/step2")
	source('/public/workspace/huangbin/ssGSEA/ssgseaMOD.r')
	expr <- fread("../GDSC/Cell_line_RMA_proc_basalExp.txt")%>%as.data.frame()
	expr <- expr[,-2]%>%filter(GENE_SYMBOLS!="")
	rownames(expr) <- expr[,1]; expr <- expr[,-1]%>%t()%>%as.data.frame()
	expr$COSMIC_ID <- as.numeric(gsub("DATA\\.","",rownames(expr)))
	gdsc1 <- read_excel("../GDSC/GDSC1_fitted_dose_response_24Jul22.xlsx")%>%as.data.frame()
	gdsc2 <- read_excel("../GDSC/GDSC2_fitted_dose_response_24Jul22.xlsx")%>%as.data.frame()
	gdsc <- rbind(gdsc1,gdsc2)
	expr <- left_join(expr,gdsc,by="COSMIC_ID")%>%filter(TCGA_DESC%in%c("LAML","ALL"))
	markers <- read.table("/public/workspace/huangbin/MPAL/figure3/TALL.AML.leukemia/labeltransfer/marker.type7.txt",
		sep="\t",stringsAsFactors=F,header=T)
	# ssGSEA score
	for(i in c("S24","TEC","AEC")){
		marker <- markers$gene[which(markers$cluster==i)][1:10]
		mat <- t(as.matrix(expr[,1:17419]))%>%round(2); colnames(mat) <- paste0("S",1:ncol(mat))
		mod.generate(marker,"tmp",out=paste('tmp','.mod',sep=''))
		mod <- mod.analyze2(mat,c('tmp'),'./',permN=0)
		save(mod,file=paste0(i,".new10.RData"))
	}
	# correlation & plot
	for(i in c("TEC","S24","AEC")){
		load(paste0(i,".new10.RData"))
		expr$score <- mod[,1]
		rst <- sapply(unique(expr$DRUG_NAME),function(x){
				tmp <- expr%>%filter(DRUG_NAME==x)
				if(nrow(tmp)>4){
					cor <- cor.test(tmp$LN_IC50,tmp$score)
					tmp <- data.frame(cor=cor$estimate,p=cor$p.value)
				}else{
					tmp <- data.frame(cor=NA,p=NA)
				}
				return(tmp)
			})
		rst <- rst%>%t()%>%as.data.frame()
		rst$cor <- unlist(rst$cor); rst$p <- unlist(rst$p)
		rst <- rst%>%filter(!is.na(cor))%>%arrange(cor)
		saveRDS(rst,file=paste0("/public/workspace/huangbin/MPAL/figure6/cor.",i,".rds"))

		rst$color <- "black"; rst$color[which(rst$p<0.05)[1:10]] <- "#ce181e"

		p <- ggplot(data = rst,aes(x = p, y = -1*cor)) + 
	  geom_point(alpha=0.5,aes(color=color)) +
	  scale_color_manual(values=c("#ce181e","black")) +
	  geom_vline(xintercept=c(0.05),lty=2,col="black",lwd=1) +
	  theme_bw()+
	  theme(
	    legend.background=element_blank(), legend.key=element_blank(),
	    legend.title = element_blank(),
	    panel.grid.major = element_blank(),panel.grid.minor = element_blank()
	  )

	  rst2=rst%>%filter(color!="black")

		p <- p + ggrepel::geom_text_repel(
	  aes(label=rownames(rst2)),rst2,
	  size = 1, #注释文本的字体大小
	  box.padding = 0.1, #字到点的距离
	  point.padding = 0.1, #字到点的距离，点周围的空白宽度
	  min.segment.length = 0, #短线段可以省略
	  segment.color = "#ce181e", #segment.colour = NA, 不显示线段
	  show.legend = F)

		pdf(paste0(i,".drug.pdf"),width=7)
		print(p)
		dev.off()
	}

	# fig6B
	# 在S24中比较CD34与HOPX与Venetoclax与Venotoclax的相关性
	library(ggstatsplot)
	tmp <- expr%>%filter(DRUG_NAME=="Venetoclax")
	pdf("CD34.venetoclax.cor.pdf",useDingbats=F)
	ggscatterstats(tmp, 
					y = LN_IC50, 
					x = CD34,
					type = "pearson",# spearman 
					centrality.para = "mean",                              
					margins = "both",                                         
					xfill = "black", 
					yfill = "black", 
					marginal.type = "densigram",  
					point.args = list(alpha=1, size=3,color="black"),
					ggstatsplot.layer = FALSE,
					smooth.line.args = list(size = 1, color = "#ce181e", method = "lm",se = TRUE))
	dev.off()
	pdf("HOPX.venetoclax.cor.pdf",useDingbats=F)
	ggscatterstats(tmp, 
					y = LN_IC50, 
					x = HOPX,
					type = "pearson",# spearman 
					centrality.para = "mean",                              
					margins = "both",                                         
					xfill = "black", 
					yfill = "black", 
					marginal.type = "densigram",  
					point.args = list(alpha=1, size=3,color="black"),
					ggstatsplot.layer = FALSE,
					smooth.line.args = list(size = 1, color = "#ce181e", method = "lm",se = TRUE))
	dev.off()
	load("S24.new10.RData")
	expr$score <- mod[,1]
	tmp <- expr%>%filter(DRUG_NAME=="Venetoclax")
	pdf("score.veneteclax.cor.pdf",useDingbats=F)
	ggscatterstats(tmp, 
					y = LN_IC50, 
					x = score,
					type = "pearson",# spearman 
					centrality.para = "mean",                              
					margins = "both",                                         
					xfill = "black", 
					yfill = "black", 
					marginal.type = "densigram",  
					point.args = list(alpha=1, size=3,color="black"),
					ggstatsplot.layer = FALSE,
					smooth.line.args = list(size = 1, color = "#ce181e", method = "lm",se = TRUE))
	dev.off()

	options(future.globals.maxSize = 2400000 * 1024^2)
	.libPaths("~/plaqueRlibrary/")
	setwd("/public/workspace/huangbin/MPAL/figure6/step2/")
	library(clusterProfiler)
	library(org.Hs.eg.db)
	library(ggplot2)
	library(RColorBrewer)
	library(dplyr)
	library(msigdf)
	library(Seurat)
	# S24 vs. TEC GSEA
	dat <- readRDS("../../figure3/TALL.AML.leukemia/labeltransfer/MPAL.pred.rds")
	dat <- subset(dat,cells=rownames(dat@meta.data)[which(dat$info%in%c("AEC","TEC","S24"))])
	markers <- presto::wilcoxauc(dat, 'info')   ## presto的方法
	saveRDS(markers,file="GSEA.rds")
	genes <- markers %>% dplyr::filter(group == "S24") %>% arrange(desc(auc)) %>% dplyr::select(feature, auc)
	ranks <- tibble::deframe(genes)

	# fig6C
	x <- readLines("VANASSE_BCL2_TARGETS_UP.v2022.1.Hs.gmt")
	res <- strsplit(x, "\t")
	names(res) <- vapply(res, function(y) y[1], character(1))
	res <- lapply(res, "[", -c(1:2))
	term <- data.frame(geneset="NASSE_BCL2_TARGETS_UP",symbol=res)
	markers <- readRDS("GSEA.rds")
	pdf("./NASSE_BCL2_TARGETS_UP.AEC.pdf",useDingbats=FALSE,width=5,height=5)
	set.seed(12345)
	genes <- markers %>% dplyr::filter(group == "S24") %>% arrange(desc(auc)) %>% dplyr::select(feature, auc)
	ranks <- tibble::deframe(genes)
	GSEA_res <- clusterProfiler::GSEA(ranks, TERM2GENE=term, pvalueCutoff = 1,
	 maxGSSize = 1500,minGSSize = 1,nPermSimple=10000, eps=0)  
	p <- GseaVis::gseaNb(object = GSEA_res,geneSetID ="NASSE_BCL2_TARGETS_UP",
	  addPval = T,pvalX = 0.75,pvalY = 0.75,pCol = 'black',pHjust = 0,subPlot = 2) + 
	ggtitle("NASSE_BCL2_TARGETS_UP")
	print(p)
	dev.off()

	# fig6D
	x <- readLines("KEGG_OXIDATIVE_PHOSPHORYLATION.v2022.1.Hs.gmt")
	res <- strsplit(x, "\t")
	names(res) <- vapply(res, function(y) y[1], character(1))
	res <- lapply(res, "[", -c(1:2))
	term <- data.frame(geneset="KEGG_OXIDATIVE_PHOSPHORYLATION",symbol=res)
	markers <- readRDS("GSEA.rds")
	pdf("./KEGG_OXIDATIVE_PHOSPHORYLATION.pdf",useDingbats=FALSE,width=5,height=5)
	set.seed(12345)
	genes <- markers %>% dplyr::filter(group == "S24") %>% arrange(desc(auc)) %>% dplyr::select(feature, auc)
	ranks <- tibble::deframe(genes)
	GSEA_res <- clusterProfiler::GSEA(ranks, TERM2GENE=term, pvalueCutoff = 1,
	 maxGSSize = 1500,minGSSize = 1,nPermSimple=10000, eps=0)  
	p <- GseaVis::gseaNb(object = GSEA_res,geneSetID ="KEGG_OXIDATIVE_PHOSPHORYLATION",
	  addPval = T,pvalX = 0.75,pvalY = 0.75,pCol = 'black',pHjust = 0,subPlot = 2) + 
	ggtitle("KEGG_OXIDATIVE_PHOSPHORYLATION")
	print(p)
	dev.off()

	# 整理AML/ALL药物信息表
	setwd("/public/workspace/huangbin/MPAL/figure6/step2/")
	gdsc1 <- read_excel("../GDSC/GDSC1_fitted_dose_response_24Jul22.xlsx")%>%as.data.frame()
	gdsc2 <- read_excel("../GDSC/GDSC2_fitted_dose_response_24Jul22.xlsx")%>%as.data.frame()
	gdsc <- rbind(gdsc1,gdsc2)%>%filter(TCGA_DESC%in%c("LAML","ALL"))
	gdsc <- gdsc[,c("DRUG_NAME","PUTATIVE_TARGET","PATHWAY_NAME")]%>%distinct()
	write.table(gdsc,"/public/workspace/huangbin/MPAL/figure6/DrugInfo.txt",sep="\t",quote=F,row.names=F)
	setwd("/public/workspace/huangbin/MPAL/figure6/")
	info <- read.table("DrugInfo.txt",sep="\t",stringsAsFactors=F,header=T)
	cor1 <- readRDS("./cor.S24.rds")%>%filter(p<0.05); cor1$DRUG_NAME <- rownames(cor1)
	cor1$unique.score <- (-1*(cor1$cor)); cor1 <- cor1[,c("unique.score","DRUG_NAME")]
	cor2 <- readRDS("./cor.AEC.rds")%>%filter(p<0.05); cor2$DRUG_NAME <- rownames(cor2)
	cor2$amllike.score <- (-1*(cor2$cor)); cor2 <- cor2[,c("amllike.score","DRUG_NAME")]
	cor3 <- readRDS("./cor.TEC.rds")%>%filter(p<0.05); cor3$DRUG_NAME <- rownames(cor3)
	cor3$alllike.score <- (-1*(cor3$cor)); cor3 <- cor3[,c("alllike.score","DRUG_NAME")]
	tmp <- left_join(left_join(left_join(info,cor1,by="DRUG_NAME"),cor2,by="DRUG_NAME"),cor3,by="DRUG_NAME")
	write.table(tmp,"./pred.drug.subtype.txt",sep="\t",quote=F,row.names=F)
}

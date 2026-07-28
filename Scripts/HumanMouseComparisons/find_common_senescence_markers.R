# Set working directory
setwd("Tables/")

# Load Seurat packages and species marts
library(biomaRt)
library(stringr)
library(tidyverse)
library(ggplot2)
library(ggupset)

human.mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse.mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl")




### SenMayo
senmayo.human <- read.delim("HumanReprogramming/senmayo_markers_xing+pana.tsv")
senmayo.mouse <- read.delim("MouseReprogramming/senmayo_markers_schiebinger.tsv")

## Create table of SenMayo markers common to humans and mice
# Download gene attributes from biomaRt
sm.bm.human <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "hgnc_symbol"), 
  filters = "external_gene_name", 
  values = senmayo.human$Gene, 
  mart = human.mart
)
sm.bm.mouse <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "mgi_symbol"), 
  filters = "external_gene_name", 
  values = senmayo.mouse$Gene, 
  mart = mouse.mart
)


# Identify human homologs of mouse markers
sm.homologs <- getHomologs(
  sm.bm.mouse$ensembl_gene_id, 
  "mus_musculus", "homo_sapiens"
)
sm.homologs <- sm.homologs[which(sm.homologs$hsapiens_homolog_ensembl_gene != ""),]


# Determine which markers/homologs overlap in human and mouse tables
sm.homsub <- sm.homologs[which(sm.homologs$hsapiens_homolog_ensembl_gene %in% sm.bm.human$ensembl_gene_id),]


# Create table of common markers
sm.commongenes <- data.frame(matrix(
  ncol = 4, 
  nrow = nrow(sm.homsub)
))
colnames(sm.commongenes) <- c(
  "Human_Ensembl_ID", "Mouse_Ensembl_ID", 
  "Human_Gene_Name", "Mouse_Gene_Name"
)

if (nrow(sm.commongenes) > 0) {
  sm.commongenes$Human_Ensembl_ID <- sm.homsub$hsapiens_homolog_ensembl_gene
  sm.commongenes$Mouse_Ensembl_ID <- sm.homsub$ensembl_gene_id
  
  for (row in 1:nrow(sm.commongenes)) {
    name.human <- sm.bm.human[which(sm.bm.human$ensembl_gene_id == sm.commongenes$Human_Ensembl_ID[row]), "external_gene_name"]
    name.mouse <- sm.bm.mouse[which(sm.bm.mouse$ensembl_gene_id == sm.commongenes$Mouse_Ensembl_ID[row]), "external_gene_name"]
    
    sm.commongenes[row, "Human_Gene_Name"] <- name.human
    sm.commongenes[row, "Mouse_Gene_Name"] <- name.mouse
  }
}

# Rename mouse genes that have human orthologs
for (gene in 1:length(senmayo.mouse$Gene)) {
  gene.name <- senmayo.mouse$Gene[gene]
  
  if (gene.name %in% sm.commongenes$Mouse_Gene_Name) {
    senmayo.mouse[gene, "Gene"] <- sm.commongenes[which(sm.commongenes$Mouse_Gene_Name == gene.name), "Human_Gene_Name"]
  }
}

# Create full table of genes and their presence in humans/mice
sm.totalgenes <- unique(c(senmayo.human$Gene, senmayo.mouse$Gene))
smdf.homology <- data.frame(matrix(
  data = FALSE, 
  nrow = length(sm.totalgenes), 
  ncol = 15
))
colnames(smdf.homology) <- c("Gene", "Human.CDKN1A-specific", "Human.CDKN2A-specific", "Human.Double-specific", "Human.CDKN1A/Double-common", "Human.CDKN2A/Double-common", "Human.CDKN1A/CDKN2A-common", "Human.All-common", "Mouse.CDKN1A-specific", "Mouse.CDKN2A-specific", "Mouse.Double-specific", "Mouse.CDKN1A/Double-common", "Mouse.CDKN2A/Double-common", "Mouse.CDKN1A/CDKN2A-common", "Mouse.All-common")
smdf.homology$Gene <- sm.totalgenes

for (gene in 1:length(smdf.homology$Gene)) {
  gene.name <- smdf.homology$Gene[gene]
  
  # Human Genes
  if (gene.name %in% senmayo.human$Gene) {
    gene.exp <- as.character(senmayo.human[which(senmayo.human$Gene == gene.name), c("CDKN1A.high", "CDKN2A.high", "Double.high"),])
    
    exp.profile <- switch(
      paste0(gene.exp, collapse = " "), 
      "Present Absent Absent" = "Human.CDKN1A-specific", 
      "Absent Present Absent" = "Human.CDKN2A-specific", 
      "Absent Absent Present" = "Human.Double-specific", 
      "Present Absent Present" = "Human.CDKN1A/Double-common", 
      "Absent Present Present" = "Human.CDKN2A/Double-common", 
      "Present Present Absent" = "Human.CDKN1A/CDKN2A-common", 
      "Present Present Present" = "Human.All-common"
    )
    
    smdf.homology[gene, exp.profile] <- TRUE
  }
  
  # Mouse Genes
  if (gene.name %in% senmayo.mouse$Gene) {
    gene.exp <- as.character(senmayo.mouse[which(senmayo.mouse$Gene == gene.name), c("CDKN1A.high", "CDKN2A.high", "Double.high"),])
    
    exp.profile <- switch(
      paste0(gene.exp, collapse = " "), 
      "Present Absent Absent" = "Mouse.CDKN1A-specific", 
      "Absent Present Absent" = "Mouse.CDKN2A-specific", 
      "Absent Absent Present" = "Mouse.Double-specific", 
      "Present Absent Present" = "Mouse.CDKN1A/Double-common", 
      "Absent Present Present" = "Mouse.CDKN2A/Double-common", 
      "Present Present Absent" = "Mouse.CDKN1A/CDKN2A-common", 
      "Present Present Present" = "Mouse.All-common"
    )
    
    smdf.homology[gene, exp.profile] <- TRUE
  }
}

# Save tables of orthologous SenMayo markers
write.table(
  sm.commongenes, file = "HumanMouseComparisons/commongenes.senmayo_any.tsv", 
  sep = "\t", row.names = FALSE, qmethod = "double"
)

exp.groups <- c("CDKN1A-specific", "CDKN2A-specific", "Double-specific", "CDKN1A_Double-common", "CDKN2A_Double-common", "CDKN1A_CDKN2A-common", "all-common")

for (i in 1:length(exp.groups)) {
  orthologs <- which(smdf.homology[, i+1] & smdf.homology[, i+8])
  if (length(orthologs) > 0) {
    common.sub <- sm.commongenes[which(sm.commongenes$Human_Gene_Name %in% smdf.homology[orthologs, "Gene"]),]
    
    write.table(
      common.sub, 
      file = paste0("HumanMouseComparisons/commongenes_senmayo_", exp.groups[i], ".tsv"), 
      sep = "\t", row.names = FALSE, qmethod = "double"
    )
  }
}




### SenSig
sensig.human <- read.delim("HumanReprogramming/sensig_markers_xing+pana.tsv")
sensig.mouse <- read.delim("MouseReprogramming/sensig_markers_schiebinger.tsv")

## Create table of SenSig markers common to humans and mice
# Download gene attributes from biomaRt
ss.bm.human <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "hgnc_symbol"), 
  filters = "external_gene_name", 
  values = sensig.human$Gene, 
  mart = human.mart
)
ss.bm.mouse <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "mgi_symbol"), 
  filters = "external_gene_name", 
  values = sensig.mouse$Gene, 
  mart = mouse.mart
)


# Identify human homologs of mouse markers
ss.homologs <- getHomologs(
  ss.bm.mouse$ensembl_gene_id, 
  "mus_musculus", "homo_sapiens"
)
ss.homologs <- ss.homologs[which(ss.homologs$hsapiens_homolog_ensembl_gene != ""),]


# Determine which markers/homologs overlap in human and mouse tables
ss.homsub <- ss.homologs[which(ss.homologs$hsapiens_homolog_ensembl_gene %in% ss.bm.human$ensembl_gene_id),]


# Create table of common markers
ss.commongenes <- data.frame(matrix(
  ncol = 4, 
  nrow = nrow(ss.homsub)
))
colnames(ss.commongenes) <- c(
  "Human_Ensembl_ID", "Mouse_Ensembl_ID", 
  "Human_Gene_Name", "Mouse_Gene_Name"
)

if (nrow(ss.commongenes) > 0) {
  ss.commongenes$Human_Ensembl_ID <- ss.homsub$hsapiens_homolog_ensembl_gene
  ss.commongenes$Mouse_Ensembl_ID <- ss.homsub$ensembl_gene_id
  
  for (row in 1:nrow(ss.commongenes)) {
    name.human <- ss.bm.human[which(ss.bm.human$ensembl_gene_id == ss.commongenes$Human_Ensembl_ID[row]), "external_gene_name"]
    name.mouse <- ss.bm.mouse[which(ss.bm.mouse$ensembl_gene_id == ss.commongenes$Mouse_Ensembl_ID[row]), "external_gene_name"]
    
    ss.commongenes[row, "Human_Gene_Name"] <- name.human
    ss.commongenes[row, "Mouse_Gene_Name"] <- name.mouse
  }
}

# Rename mouse genes that have human orthologs
for (gene in 1:length(sensig.mouse$Gene)) {
  gene.name <- sensig.mouse$Gene[gene]
  
  if (gene.name %in% ss.commongenes$Mouse_Gene_Name) {
    sensig.mouse[gene, "Gene"] <- ss.commongenes[which(ss.commongenes$Mouse_Gene_Name == gene.name), "Human_Gene_Name"]
  }
}

# Create full table of genes and their presence in humans/mice
ss.totalgenes <- unique(c(sensig.human$Gene, sensig.mouse$Gene))
ssdf.homology <- data.frame(matrix(
  data = FALSE, 
  nrow = length(ss.totalgenes), 
  ncol = 15
))
colnames(ssdf.homology) <- c("Gene", "Human.CDKN1A-specific", "Human.CDKN2A-specific", "Human.Double-specific", "Human.CDKN1A/Double-common", "Human.CDKN2A/Double-common", "Human.CDKN1A/CDKN2A-common", "Human.All-common", "Mouse.CDKN1A-specific", "Mouse.CDKN2A-specific", "Mouse.Double-specific", "Mouse.CDKN1A/Double-common", "Mouse.CDKN2A/Double-common", "Mouse.CDKN1A/CDKN2A-common", "Mouse.All-common")
ssdf.homology$Gene <- ss.totalgenes

for (gene in 1:length(ssdf.homology$Gene)) {
  gene.name <- ssdf.homology$Gene[gene]
  
  # Human Genes
  if (gene.name %in% sensig.human$Gene) {
    gene.exp <- as.character(sensig.human[which(sensig.human$Gene == gene.name), c("CDKN1A.high", "CDKN2A.high", "Double.high"),])
    
    exp.profile <- switch(
      paste0(gene.exp, collapse = " "), 
      "Present Absent Absent" = "Human.CDKN1A-specific", 
      "Absent Present Absent" = "Human.CDKN2A-specific", 
      "Absent Absent Present" = "Human.Double-specific", 
      "Present Absent Present" = "Human.CDKN1A/Double-common", 
      "Absent Present Present" = "Human.CDKN2A/Double-common", 
      "Present Present Absent" = "Human.CDKN1A/CDKN2A-common", 
      "Present Present Present" = "Human.All-common"
    )
    
    ssdf.homology[gene, exp.profile] <- TRUE
  }
  
  # Mouse Genes
  if (gene.name %in% sensig.mouse$Gene) {
    gene.exp <- as.character(sensig.mouse[which(sensig.mouse$Gene == gene.name), c("CDKN1A.high", "CDKN2A.high", "Double.high"),])
    
    exp.profile <- switch(
      paste0(gene.exp, collapse = " "), 
      "Present Absent Absent" = "Mouse.CDKN1A-specific", 
      "Absent Present Absent" = "Mouse.CDKN2A-specific", 
      "Absent Absent Present" = "Mouse.Double-specific", 
      "Present Absent Present" = "Mouse.CDKN1A/Double-common", 
      "Absent Present Present" = "Mouse.CDKN2A/Double-common", 
      "Present Present Absent" = "Mouse.CDKN1A/CDKN2A-common", 
      "Present Present Present" = "Mouse.All-common"
    )
    
    ssdf.homology[gene, exp.profile] <- TRUE
  }
}

# Save tables of orthologous SenSig markers
write.table(
  ss.commongenes, file = "HumanMouseComparisons/commongenes.sensig_any.tsv", 
  sep = "\t", row.names = FALSE, qmethod = "double"
)

exp.groups <- c("CDKN1A-specific", "CDKN2A-specific", "Double-specific", "CDKN1A_Double-common", "CDKN2A_Double-common", "CDKN1A_CDKN2A-common", "all-common")

for (i in 1:length(exp.groups)) {
  orthologs <- which(ssdf.homology[, i+1] & ssdf.homology[, i+8])
  if (length(orthologs) > 0) {
    common.sub <- ss.commongenes[which(ss.commongenes$Human_Gene_Name %in% ssdf.homology[orthologs, "Gene"]),]
    
    write.table(
      common.sub, 
      file = paste0("HumanMouseComparisons/commongenes_sensig_", exp.groups[i], ".tsv"), 
      sep = "\t", row.names = FALSE, qmethod = "double"
    )
  }
}


## Upset Plot
# Format table for use in ggupset
tidy.homology <- as_tibble(ssdf.homology)
tidy.homology <- tidy.homology |> mutate(
  combination = pmap(
    list(`Human.CDKN1A-specific`, 
         `Human.CDKN2A-specific`, 
         `Human.Double-specific`, 
         `Human.CDKN1A/Double-common`, 
         `Human.CDKN2A/Double-common`, 
         `Human.CDKN1A/CDKN2A-common`, 
         `Human.All-common`, 
         `Mouse.CDKN1A-specific`, 
         `Mouse.CDKN2A-specific`, 
         `Mouse.Double-specific`, 
         `Mouse.CDKN1A/Double-common`, 
         `Mouse.CDKN2A/Double-common`, 
         `Mouse.CDKN1A/CDKN2A-common`, 
         `Mouse.All-common`
    ), 
    \(lgl1, lgl2, lgl3, lgl4, lgl5, 
      lgl6, lgl7, lgl8, lgl9, lgl10, 
      lgl11, lgl12, lgl13, lgl14) {
      c("Human: CDKN1A-specific", "Human: CDKN2A-specific", "Human: Double-specific", "Human: CDKN1A/Double-common", "Human: CDKN2A/Double-common", "Human: CDKN1A/CDKN2A-common", "Human: In All Groups", "Mouse: CDKN1A-specific", "Mouse: CDKN2A-specific", "Mouse: Double-specific", "Mouse: CDKN1A/Double-common", "Mouse: CDKN2A/Double-common", "Mouse: CDKN1A/CDKN2A-common", "Mouse: In All Groups")[c(lgl1, lgl2, lgl3, lgl4, lgl5, lgl6, lgl7, lgl8, lgl9, lgl10, lgl11, lgl12, lgl13, lgl14)]
    }
  )
)

# Assign bar colors for each expression group
tidy.homology$ExpressionGroup <- "Mixed"
for (i in 1:nrow(tidy.homology)) {
  exp.group <- switch(
    paste0(tidy.homology$combination[[i]], collapse = "  "), 
    "Human: CDKN1A-specific" = "CDKN1A-specific",
    "Human: CDKN2A-specific" = "CDKN2A-specific",
    "Human: Double-specific" = "Double-specific",
    "Human: CDKN1A/Double-common" = "CDKN1A/Double-common",
    "Human: CDKN2A/Double-common" = "CDKN2A/Double-common",
    "Human: CDKN1A/CDKN2A-common" = "CDKN1A/CDKN2A-common",
    "Human: In All Groups" = "In All Groups",
    "Mouse: CDKN1A-specific" = "CDKN1A-specific",
    "Mouse: CDKN2A-specific" = "CDKN2A-specific",
    "Mouse: Double-specific" = "Double-specific",
    "Mouse: CDKN1A/Double-common" = "CDKN1A/Double-common",
    "Mouse: CDKN2A/Double-common" = "CDKN2A/Double-common",
    "Mouse: CDKN1A/CDKN2A-common" = "CDKN1A/CDKN2A-common",
    "Mouse: In All Groups" = "In All Groups",
    "Human: CDKN1A-specific  Mouse: CDKN1A-specific" = "CDKN1A-specific",
    "Human: CDKN2A-specific  Mouse: CDKN2A-specific" = "CDKN2A-specific",
    "Human: Double-specific  Mouse: Double-specific" = "Double-specific",
    "Human: CDKN1A/Double-common  Mouse: CDKN1A/Double-common" = "CDKN1A/Double-common",
    "Human: CDKN2A/Double-common  Mouse: CDKN2A/Double-common" = "CDKN2A/Double-common",
    "Human: CDKN1A/CDKN2A-common  Mouse: CDKN1A/CDKN2A-common" = "CDKN1A/CDKN2A-common",
    "Human: In All Groups  Mouse: In All Groups" = "In All Groups",
    "Mixed"
  )
  tidy.homology[i, "ExpressionGroup"] <- exp.group
}

colorcodes <- c("#FBFF83", "#F49089", "#89AEE4", "#7FB198", "#9D8FC4", "#F29966", "#A48671", "#585858")
colorcodes <- setNames(
  colorcodes, 
  c("CDKN1A-specific", "CDKN2A-specific", "Double-specific", "CDKN1A/Double-common", "CDKN2A/Double-common", "CDKN1A/CDKN2A-common", "In All Groups", "Mixed")
)

# Create and save UpSet plot
upset.plot <- tidy.homology |> 
  ggplot(aes(x = combination, fill = ExpressionGroup)) + 
  geom_bar() + 
  geom_text(
    stat = "count", 
    aes(label = after_stat(count)), 
    vjust =-1, size = 6
  ) + scale_x_upset() + scale_fill_manual(values = colorcodes) + 
  labs(x = element_blank(), y = "# of SenSig Genes") + theme(
    legend.position = "none", 
    text = element_text(size = 25)
  ) + theme_combmatrix(
    combmatrix.label.extra_spacing = 1, 
    combmatrix.label.text = element_text(size = 15), 
    combmatrix.panel.point.size = 5, 
    combmatrix.panel.line.size = 1.7
  )

png(
  filename = "../Images/HumanMouseComparisons/upsetplot_homology_sensig.png", 
  width = 1740, height = 1860, res = 180
)
upset.plot
dev.off()




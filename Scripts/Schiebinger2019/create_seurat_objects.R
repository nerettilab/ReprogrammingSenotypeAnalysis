# Set working directory
setwd("Data/Schiebinger2019/")

# Load packages
library(Seurat)
library(ggplot2)
library(ggrastr)


## Load counts data
# D0-0.5
D0_Dox_C1 <- Read10X_h5("GSE115943/GSM3195648_D0_Dox_C1_gene_bc_mat.h5")
D0_Dox_C2 <- Read10X_h5("GSE115943/GSM3195649_D0_Dox_C2_gene_bc_mat.h5")
D0.5_Dox_C1 <- Read10X_h5("GSE115943/GSM3195650_D0.5_Dox_C1_gene_bc_mat.h5")
D0.5_Dox_C2 <- Read10X_h5("GSE115943/GSM3195651_D0.5_Dox_C2_gene_bc_mat.h5")

# D1-1.5
D1_Dox_C1 <- Read10X_h5("GSE115943/GSM3195652_D1_Dox_C1_gene_bc_mat.h5")
D1_Dox_C2 <- Read10X_h5("GSE115943/GSM3195653_D1_Dox_C2_gene_bc_mat.h5")
D1.5_Dox_C1 <- Read10X_h5("GSE115943/GSM3195654_D1.5_Dox_C1_gene_bc_mat.h5")
D1.5_Dox_C2 <- Read10X_h5("GSE115943/GSM3195655_D1.5_Dox_C2_gene_bc_mat.h5")

# D2-2.5
D2_Dox_C1 <- Read10X_h5("GSE115943/GSM3195656_D2_Dox_C1_gene_bc_mat.h5")
D2_Dox_C2 <- Read10X_h5("GSE115943/GSM3195657_D2_Dox_C2_gene_bc_mat.h5")
D2.5_Dox_C1 <- Read10X_h5("GSE115943/GSM3195658_D2.5_Dox_C1_gene_bc_mat.h5")
D2.5_Dox_C2 <- Read10X_h5("GSE115943/GSM3195659_D2.5_Dox_C2_gene_bc_mat.h5")

# D3-3.5
D3_Dox_C1 <- Read10X_h5("GSE115943/GSM3195660_D3_Dox_C1_gene_bc_mat.h5")
D3_Dox_C2 <- Read10X_h5("GSE115943/GSM3195661_D3_Dox_C2_gene_bc_mat.h5")
D3.5_Dox_C1 <- Read10X_h5("GSE115943/GSM3195662_D3.5_Dox_C1_gene_bc_mat.h5")
D3.5_Dox_C2 <- Read10X_h5("GSE115943/GSM3195663_D3.5_Dox_C2_gene_bc_mat.h5")

# D4-4.5
D4_Dox_C1 <- Read10X_h5("GSE115943/GSM3195664_D4_Dox_C1_gene_bc_mat.h5")
D4_Dox_C2 <- Read10X_h5("GSE115943/GSM3195665_D4_Dox_C2_gene_bc_mat.h5")
D4.5_Dox_C1 <- Read10X_h5("GSE115943/GSM3195666_D4.5_Dox_C1_gene_bc_mat.h5")
D4.5_Dox_C2 <- Read10X_h5("GSE115943/GSM3195667_D4.5_Dox_C2_gene_bc_mat.h5")

# D5-5.5
D5_Dox_C1 <- Read10X_h5("GSE115943/GSM3195668_D5_Dox_C1_gene_bc_mat.h5")
D5_Dox_C2 <- Read10X_h5("GSE115943/GSM3195669_D5_Dox_C2_gene_bc_mat.h5")
D5.5_Dox_C1 <- Read10X_h5("GSE115943/GSM3195670_D5.5_Dox_C1_gene_bc_mat.h5")
D5.5_Dox_C2 <- Read10X_h5("GSE115943/GSM3195671_D5.5_Dox_C2_gene_bc_mat.h5")

# D6-6.5
D6_Dox_C1 <- Read10X_h5("GSE115943/GSM3195672_D6_Dox_C1_gene_bc_mat.h5")
D6_Dox_C2 <- Read10X_h5("GSE115943/GSM3195673_D6_Dox_C2_gene_bc_mat.h5")
D6.5_Dox_C1 <- Read10X_h5("GSE115943/GSM3195674_D6.5_Dox_C1_gene_bc_mat.h5")
D6.5_Dox_C2 <- Read10X_h5("GSE115943/GSM3195675_D6.5_Dox_C2_gene_bc_mat.h5")

# D7-7.5
D7_Dox_C1 <- Read10X_h5("GSE115943/GSM3195676_D7_Dox_C1_gene_bc_mat.h5")
D7_Dox_C2 <- Read10X_h5("GSE115943/GSM3195677_D7_Dox_C2_gene_bc_mat.h5")
D7.5_Dox_C1 <- Read10X_h5("GSE115943/GSM3195678_D7.5_Dox_C1_gene_bc_mat.h5")
D7.5_Dox_C2 <- Read10X_h5("GSE115943/GSM3195679_D7.5_Dox_C2_gene_bc_mat.h5")

# D8
D8_Dox_C1 <- Read10X_h5("GSE115943/GSM3195680_D8_Dox_C1_gene_bc_mat.h5")
D8_Dox_C2 <- Read10X_h5("GSE115943/GSM3195681_D8_Dox_C2_gene_bc_mat.h5")

# D8.25-8.75
D8.25_2i_C1 <- Read10X_h5("GSE115943/GSM3195682_D8.25_2i_C1_gene_bc_mat.h5")
D8.25_2i_C2 <- Read10X_h5("GSE115943/GSM3195683_D8.25_2i_C2_gene_bc_mat.h5")
D8.25_serum_C1 <- Read10X_h5("GSE115943/GSM3195684_D8.25_serum_C1_gene_bc_mat.h5")
D8.25_serum_C2 <- Read10X_h5("GSE115943/GSM3195685_D8.25_serum_C2_gene_bc_mat.h5")
D8.5_2i_C1 <- Read10X_h5("GSE115943/GSM3195686_D8.5_2i_C1_gene_bc_mat.h5")
D8.5_2i_C2 <- Read10X_h5("GSE115943/GSM3195687_D8.5_2i_C2_gene_bc_mat.h5")
D8.5_serum_C1 <- Read10X_h5("GSE115943/GSM3195688_D8.5_serum_C1_gene_bc_mat.h5")
D8.5_serum_C2 <- Read10X_h5("GSE115943/GSM3195689_D8.5_serum_C2_gene_bc_mat.h5")
D8.75_2i_C1 <- Read10X_h5("GSE115943/GSM3195690_D8.75_2i_C1_gene_bc_mat.h5")
D8.75_2i_C2 <- Read10X_h5("GSE115943/GSM3195691_D8.75_2i_C2_gene_bc_mat.h5")
D8.75_serum_C1 <- Read10X_h5("GSE115943/GSM3195692_D8.75_serum_C1_gene_bc_mat.h5")
D8.75_serum_C2 <- Read10X_h5("GSE115943/GSM3195693_D8.75_serum_C2_gene_bc_mat.h5")

# D9-9.5
D9_2i_C1 <- Read10X_h5("GSE115943/GSM3195694_D9_2i_C1_gene_bc_mat.h5")
D9_2i_C2 <- Read10X_h5("GSE115943/GSM3195695_D9_2i_C2_gene_bc_mat.h5")
D9_serum_C1 <- Read10X_h5("GSE115943/GSM3195696_D9_serum_C1_gene_bc_mat.h5")
D9_serum_C2 <- Read10X_h5("GSE115943/GSM3195697_D9_serum_C2_gene_bc_mat.h5")
D9.5_2i_C1 <- Read10X_h5("GSE115943/GSM3195698_D9.5_2i_C1_gene_bc_mat.h5")
D9.5_2i_C2 <- Read10X_h5("GSE115943/GSM3195699_D9.5_2i_C2_gene_bc_mat.h5")
D9.5_serum_C1 <- Read10X_h5("GSE115943/GSM3195700_D9.5_serum_C1_gene_bc_mat.h5")
D9.5_serum_C2 <- Read10X_h5("GSE115943/GSM3195701_D9.5_serum_C2_gene_bc_mat.h5")

# D10-10.5
D10_2i_C1 <- Read10X_h5("GSE115943/GSM3195702_D10_2i_C1_gene_bc_mat.h5")
D10_2i_C2 <- Read10X_h5("GSE115943/GSM3195703_D10_2i_C2_gene_bc_mat.h5")
D10_serum_C1 <- Read10X_h5("GSE115943/GSM3195704_D10_serum_C1_gene_bc_mat.h5")
D10_serum_C2 <- Read10X_h5("GSE115943/GSM3195705_D10_serum_C2_gene_bc_mat.h5")
D10.5_2i_C1 <- Read10X_h5("GSE115943/GSM3195706_D10.5_2i_C1_gene_bc_mat.h5")
D10.5_2i_C2 <- Read10X_h5("GSE115943/GSM3195707_D10.5_2i_C2_gene_bc_mat.h5")
D10.5_serum_C1 <- Read10X_h5("GSE115943/GSM3195708_D10.5_serum_C1_gene_bc_mat.h5")
D10.5_serum_C2 <- Read10X_h5("GSE115943/GSM3195709_D10.5_serum_C2_gene_bc_mat.h5")

# D11-11.5
D11_2i_C1 <- Read10X_h5("GSE115943/GSM3195710_D11_2i_C1_gene_bc_mat.h5")
D11_2i_C2 <- Read10X_h5("GSE115943/GSM3195711_D11_2i_C2_gene_bc_mat.h5")
D11_serum_C1 <- Read10X_h5("GSE115943/GSM3195712_D11_serum_C1_gene_bc_mat.h5")
D11_serum_C2 <- Read10X_h5("GSE115943/GSM3195713_D11_serum_C2_gene_bc_mat.h5")
D11.5_2i_C1 <- Read10X_h5("GSE115943/GSM3195714_D11.5_2i_C1_gene_bc_mat.h5")
D11.5_2i_C2 <- Read10X_h5("GSE115943/GSM3195715_D11.5_2i_C2_gene_bc_mat.h5")
D11.5_serum_C1 <- Read10X_h5("GSE115943/GSM3195716_D11.5_serum_C1_gene_bc_mat.h5")
D11.5_serum_C2 <- Read10X_h5("GSE115943/GSM3195717_D11.5_serum_C2_gene_bc_mat.h5")

# D12-12.5
D12_2i_C1 <- Read10X_h5("GSE115943/GSM3195718_D12_2i_C1_gene_bc_mat.h5")
D12_2i_C2 <- Read10X_h5("GSE115943/GSM3195719_D12_2i_C2_gene_bc_mat.h5")
D12_serum_C1 <- Read10X_h5("GSE115943/GSM3195720_D12_serum_C1_gene_bc_mat.h5")
D12_serum_C2 <- Read10X_h5("GSE115943/GSM3195721_D12_serum_C2_gene_bc_mat.h5")
D12.5_2i_C1 <- Read10X_h5("GSE115943/GSM3195722_D12.5_2i_C1_gene_bc_mat.h5")
D12.5_2i_C2 <- Read10X_h5("GSE115943/GSM3195723_D12.5_2i_C2_gene_bc_mat.h5")
D12.5_serum_C1 <- Read10X_h5("GSE115943/GSM3195724_D12.5_serum_C1_gene_bc_mat.h5")
D12.5_serum_C2 <- Read10X_h5("GSE115943/GSM3195725_D12.5_serum_C2_gene_bc_mat.h5")

# D13-13.5
D13_2i_C1 <- Read10X_h5("GSE115943/GSM3195726_D13_2i_C1_gene_bc_mat.h5")
D13_2i_C2 <- Read10X_h5("GSE115943/GSM3195727_D13_2i_C2_gene_bc_mat.h5")
D13_serum_C1 <- Read10X_h5("GSE115943/GSM3195728_D13_serum_C1_gene_bc_mat.h5")
D13_serum_C2 <- Read10X_h5("GSE115943/GSM3195729_D13_serum_C2_gene_bc_mat.h5")
D13.5_2i_C1 <- Read10X_h5("GSE115943/GSM3195730_D13.5_2i_C1_gene_bc_mat.h5")
D13.5_2i_C2 <- Read10X_h5("GSE115943/GSM3195731_D13.5_2i_C2_gene_bc_mat.h5")
D13.5_serum_C1 <- Read10X_h5("GSE115943/GSM3195732_D13.5_serum_C1_gene_bc_mat.h5")
D13.5_serum_C2 <- Read10X_h5("GSE115943/GSM3195733_D13.5_serum_C2_gene_bc_mat.h5")

# D14-14.5
D14_2i_C1 <- Read10X_h5("GSE115943/GSM3195734_D14_2i_C1_gene_bc_mat.h5")
D14_2i_C2 <- Read10X_h5("GSE115943/GSM3195735_D14_2i_C2_gene_bc_mat.h5")
D14_serum_C1 <- Read10X_h5("GSE115943/GSM3195736_D14_serum_C1_gene_bc_mat.h5")
D14_serum_C2 <- Read10X_h5("GSE115943/GSM3195737_D14_serum_C2_gene_bc_mat.h5")
D14.5_2i_C1 <- Read10X_h5("GSE115943/GSM3195738_D14.5_2i_C1_gene_bc_mat.h5")
D14.5_2i_C2 <- Read10X_h5("GSE115943/GSM3195739_D14.5_2i_C2_gene_bc_mat.h5")
D14.5_serum_C1 <- Read10X_h5("GSE115943/GSM3195740_D14.5_serum_C1_gene_bc_mat.h5")
D14.5_serum_C2 <- Read10X_h5("GSE115943/GSM3195741_D14.5_serum_C2_gene_bc_mat.h5")

# D15-15.5
D15_2i_C1 <- Read10X_h5("GSE115943/GSM3195742_D15_2i_C1_gene_bc_mat.h5")
D15_2i_C2 <- Read10X_h5("GSE115943/GSM3195743_D15_2i_C2_gene_bc_mat.h5")
D15_serum_C1 <- Read10X_h5("GSE115943/GSM3195744_D15_serum_C1_gene_bc_mat.h5")
D15_serum_C2 <- Read10X_h5("GSE115943/GSM3195745_D15_serum_C2_gene_bc_mat.h5")
D15.5_2i_C1 <- Read10X_h5("GSE115943/GSM3195746_D15.5_2i_C1_gene_bc_mat.h5")
D15.5_2i_C2 <- Read10X_h5("GSE115943/GSM3195747_D15.5_2i_C2_gene_bc_mat.h5")
D15.5_serum_C1 <- Read10X_h5("GSE115943/GSM3195748_D15.5_serum_C1_gene_bc_mat.h5")
D15.5_serum_C2 <- Read10X_h5("GSE115943/GSM3195749_D15.5_serum_C2_gene_bc_mat.h5")

# D16-16.5
D16_2i_C1 <- Read10X_h5("GSE115943/GSM3195750_D16_2i_C1_gene_bc_mat.h5")
D16_2i_C2 <- Read10X_h5("GSE115943/GSM3195751_D16_2i_C2_gene_bc_mat.h5")
D16_serum_C1 <- Read10X_h5("GSE115943/GSM3195752_D16_serum_C1_gene_bc_mat.h5")
D16_serum_C2 <- Read10X_h5("GSE115943/GSM3195753_D16_serum_C2_gene_bc_mat.h5")
D16.5_2i_C1 <- Read10X_h5("GSE115943/GSM3195754_D16.5_2i_C1_gene_bc_mat.h5")
D16.5_2i_C2 <- Read10X_h5("GSE115943/GSM3195755_D16.5_2i_C2_gene_bc_mat.h5")
D16.5_serum_C1 <- Read10X_h5("GSE115943/GSM3195756_D16.5_serum_C1_gene_bc_mat.h5")
D16.5_serum_C2 <- Read10X_h5("GSE115943/GSM3195757_D16.5_serum_C2_gene_bc_mat.h5")

# D17-17.5
D17_2i_C1 <- Read10X_h5("GSE115943/GSM3195758_D17_2i_C1_gene_bc_mat.h5")
D17_2i_C2 <- Read10X_h5("GSE115943/GSM3195759_D17_2i_C2_gene_bc_mat.h5")
D17_serum_C1 <- Read10X_h5("GSE115943/GSM3195760_D17_serum_C1_gene_bc_mat.h5")
D17_serum_C2 <- Read10X_h5("GSE115943/GSM3195761_D17_serum_C2_gene_bc_mat.h5")
D17.5_2i_C1 <- Read10X_h5("GSE115943/GSM3195762_D17.5_2i_C1_gene_bc_mat.h5")
D17.5_2i_C2 <- Read10X_h5("GSE115943/GSM3195763_D17.5_2i_C2_gene_bc_mat.h5")
D17.5_serum_C1 <- Read10X_h5("GSE115943/GSM3195764_D17.5_serum_C1_gene_bc_mat.h5")
D17.5_serum_C2 <- Read10X_h5("GSE115943/GSM3195765_D17.5_serum_C2_gene_bc_mat.h5")

# D18 and iPSCs
D18_2i_C1 <- Read10X_h5("GSE115943/GSM3195766_D18_2i_C1_gene_bc_mat.h5")
D18_2i_C2 <- Read10X_h5("GSE115943/GSM3195767_D18_2i_C2_gene_bc_mat.h5")
D18_serum_C1 <- Read10X_h5("GSE115943/GSM3195768_D18_serum_C1_gene_bc_mat.h5")
D18_serum_C2 <- Read10X_h5("GSE115943/GSM3195769_D18_serum_C2_gene_bc_mat.h5")
DiPSC_2i_C1 <- Read10X_h5("GSE115943/GSM3195770_DiPSC_2i_C1_gene_bc_mat.h5")
DiPSC_2i_C2 <- Read10X_h5("GSE115943/GSM3195771_DiPSC_2i_C2_gene_bc_mat.h5")
DiPSC_serum_C1 <- Read10X_h5("GSE115943/GSM3195772_DiPSC_serum_C1_gene_bc_mat.h5")
DiPSC_serum_C2 <- Read10X_h5("GSE115943/GSM3195773_DiPSC_serum_C2_gene_bc_mat.h5")


# Create Seurat objects and add metadata
for (i in 1:length(mget(ls(pattern= "_C")))) {
  seu <- CreateSeuratObject(
    counts = mget(ls(pattern = "_C"))[i], 
    min.cells = 3, min.features = 200
  )
  seu$orig.ident <- paste0("scbg_", ls(pattern = "_C")[i])
  seu$dataset <- "Schiebinger"
  
  # Add condition metadata
  switch(sub(".*_(.+)_.*", "\\1", ls(pattern = "_C")[i]), 
         "Dox" = seu$Condition <- "Dox", 
         "2i" = seu$Condition <- "2i", 
         "serum" = seu$Condition <- "Serum", 
         "Unknown"
  )
  
  # Add timepoint metadata
  switch(sub("_.*", "", ls(pattern = "_C")[i]), 
         "D0" = seu$Timepoint <- "D00", 
         "D0.5" = seu$Timepoint <- "D00", 
         "D1" = seu$Timepoint <- "D01", 
         "D1.5" = seu$Timepoint <- "D01", 
         "D2" = seu$Timepoint <- "D02", 
         "D2.5" = seu$Timepoint <- "D02", 
         "D3" = seu$Timepoint <- "D03", 
         "D3.5" = seu$Timepoint <- "D03", 
         "D4" = seu$Timepoint <- "D04", 
         "D4.5" = seu$Timepoint <- "D04", 
         "D5" = seu$Timepoint <- "D05", 
         "D5.5" = seu$Timepoint <- "D05", 
         "D6" = seu$Timepoint <- "D06", 
         "D6.5" = seu$Timepoint <- "D06", 
         "D7" = seu$Timepoint <- "D07", 
         "D7.5" = seu$Timepoint <- "D07", 
         "D8" = seu$Timepoint <- "D08", 
         "D8.25" = seu$Timepoint <- "D08", 
         "D8.5" = seu$Timepoint <- "D08", 
         "D8.75" = seu$Timepoint <- "D08", 
         "D9" = seu$Timepoint <- "D09", 
         "D9.5" = seu$Timepoint <- "D09", 
         "D10" = seu$Timepoint <- "D10", 
         "D10.5" = seu$Timepoint <- "D10", 
         "D11" = seu$Timepoint <- "D11", 
         "D11.5" = seu$Timepoint <- "D11", 
         "D12" = seu$Timepoint <- "D12", 
         "D12.5" = seu$Timepoint <- "D12", 
         "D13" = seu$Timepoint <- "D13", 
         "D13.5" = seu$Timepoint <- "D13", 
         "D14" = seu$Timepoint <- "D14", 
         "D14.5" = seu$Timepoint <- "D14", 
         "D15" = seu$Timepoint <- "D15", 
         "D15.5" = seu$Timepoint <- "D15", 
         "D16" = seu$Timepoint <- "D15", 
         "D16.5" = seu$Timepoint <- "D16", 
         "D17" = seu$Timepoint <- "D17", 
         "D17.5" = seu$Timepoint <- "D17", 
         "D18" = seu$Timepoint <- "D18", 
         "DiPSC" = seu$Timepoint <- "iPSC", 
         "Unknown"
  )
  
  # Calculate mitochondrial content
  seu[["percent.mt"]] <- PercentageFeatureSet(
    object = seu, 
    pattern = "^mt-"
  )
  
  # Quality Control Filtering
  ub <- quantile(seu$nFeature_RNA, probs = 0.95)
  lb <- quantile(seu$nFeature_RNA, probs = 0.05)
  seu <- subset(
    x = seu, 
    subset = nFeature_RNA > lb & 
      nFeature_RNA < ub & 
      percent.mt < 10
  )
  
  assign(paste0("scbg.", ls(pattern = "_C")[i]), seu)
}

# Garbage collection
rm(D0_Dox_C1, D0_Dox_C2, D0.5_Dox_C1, D0.5_Dox_C2, D1_Dox_C1, 
   D1_Dox_C2, D1.5_Dox_C1, D1.5_Dox_C2, D2_Dox_C1, D2_Dox_C2, 
   D2.5_Dox_C1, D2.5_Dox_C2, D3_Dox_C1, D3_Dox_C2, 
   D3.5_Dox_C1, D3.5_Dox_C2, D4_Dox_C1, D4_Dox_C2, 
   D4.5_Dox_C1, D4.5_Dox_C2, D5_Dox_C1, D5_Dox_C2, 
   D5.5_Dox_C1, D5.5_Dox_C2, D6_Dox_C1, D6_Dox_C2, 
   D6.5_Dox_C1, D6.5_Dox_C2, D7_Dox_C1, D7_Dox_C2, 
   D7.5_Dox_C1, D7.5_Dox_C2, D8_Dox_C1, D8_Dox_C2, 
   D8.25_2i_C1, D8.25_2i_C2, D8.25_serum_C1, D8.25_serum_C2, 
   D8.5_2i_C1, D8.5_2i_C2, D8.5_serum_C1, D8.5_serum_C2, 
   D8.75_2i_C1, D8.75_2i_C2, D8.75_serum_C1, D8.75_serum_C2, 
   D9_2i_C1, D9_2i_C2, D9_serum_C1, D9_serum_C2, D9.5_2i_C1, 
   D9.5_2i_C2, D9.5_serum_C1, D9.5_serum_C2, D10_2i_C1, 
   D10_2i_C2, D10_serum_C1, D10_serum_C2, D10.5_2i_C1, 
   D10.5_2i_C2, D10.5_serum_C1, D10.5_serum_C2, D11_2i_C1, 
   D11_2i_C2, D11_serum_C1, D11_serum_C2, D11.5_2i_C1, 
   D11.5_2i_C2, D11.5_serum_C1, D11.5_serum_C2, D12_2i_C1, 
   D12_2i_C2, D12_serum_C1, D12_serum_C2, D12.5_2i_C1, 
   D12.5_2i_C2, D12.5_serum_C1, D12.5_serum_C2, D13_2i_C1, 
   D13_2i_C2, D13_serum_C1, D13_serum_C2, D13.5_2i_C1, 
   D13.5_2i_C2, D13.5_serum_C1, D13.5_serum_C2, D14_2i_C1, 
   D14_2i_C2, D14_serum_C1, D14_serum_C2, D14.5_2i_C1, 
   D14.5_2i_C2, D14.5_serum_C1, D14.5_serum_C2, D15_2i_C1, 
   D15_2i_C2, D15_serum_C1, D15_serum_C2, D15.5_2i_C1, 
   D15.5_2i_C2, D15.5_serum_C1, D15.5_serum_C2, D16_2i_C1, 
   D16_2i_C2, D16_serum_C1, D16_serum_C2, D16.5_2i_C1, 
   D16.5_2i_C2, D16.5_serum_C1, D16.5_serum_C2, D17_2i_C1, 
   D17_2i_C2, D17_serum_C1, D17_serum_C2, D17.5_2i_C1, 
   D17.5_2i_C2, D17.5_serum_C1, D17.5_serum_C2, D18_2i_C1, 
   D18_2i_C2, D18_serum_C1, D18_serum_C2, DiPSC_2i_C1, 
   DiPSC_2i_C2, DiPSC_serum_C1, DiPSC_serum_C2, seu, i)
gc(reset=TRUE)


# Merge filtered Seurat objects
scbg.list <- vector(
  mode = "list", 
  length = length(ls(pattern = "_C")) - 1
)
for (i in 1:(length(scbg.list))) {
  scbg.list[i] <- mget(ls(pattern = "_C")[i+1])
}

scbg <- merge(
  x = scbg.D0_Dox_C1, y = scbg.list, 
  add.cell.ids = sub(".*scbg.", "", ls(pattern = "_C"))
)

# Garbage collection
rm(scbg.list, scbg.D0_Dox_C1, scbg.D0_Dox_C2, 
   scbg.D0.5_Dox_C1, scbg.D0.5_Dox_C2, scbg.D1_Dox_C1, 
   scbg.D1_Dox_C2, scbg.D1.5_Dox_C1, scbg.D1.5_Dox_C2, 
   scbg.D2_Dox_C1, scbg.D2_Dox_C2, scbg.D2.5_Dox_C1, 
   scbg.D2.5_Dox_C2, scbg.D3_Dox_C1, scbg.D3_Dox_C2, 
   scbg.D3.5_Dox_C1, scbg.D3.5_Dox_C2, scbg.D4_Dox_C1, 
   scbg.D4_Dox_C2, scbg.D4.5_Dox_C1, scbg.D4.5_Dox_C2, 
   scbg.D5_Dox_C1, scbg.D5_Dox_C2, scbg.D5.5_Dox_C1, 
   scbg.D5.5_Dox_C2, scbg.D6_Dox_C1, scbg.D6_Dox_C2, 
   scbg.D6.5_Dox_C1, scbg.D6.5_Dox_C2, scbg.D7_Dox_C1, 
   scbg.D7_Dox_C2, scbg.D7.5_Dox_C1, scbg.D7.5_Dox_C2, 
   scbg.D8_Dox_C1, scbg.D8_Dox_C2, scbg.D8.25_2i_C1, 
   scbg.D8.25_2i_C2, scbg.D8.25_serum_C1, 
   scbg.D8.25_serum_C2, scbg.D8.5_2i_C1, scbg.D8.5_2i_C2, 
   scbg.D8.5_serum_C1, scbg.D8.5_serum_C2, scbg.D8.75_2i_C1, 
   scbg.D8.75_2i_C2, scbg.D8.75_serum_C1, 
   scbg.D8.75_serum_C2, scbg.D9_2i_C1, scbg.D9_2i_C2, 
   scbg.D9_serum_C1, scbg.D9_serum_C2, scbg.D9.5_2i_C1, 
   scbg.D9.5_2i_C2, scbg.D9.5_serum_C1, scbg.D9.5_serum_C2, 
   scbg.D10_2i_C1, scbg.D10_2i_C2, scbg.D10_serum_C1, 
   scbg.D10_serum_C2, scbg.D10.5_2i_C1, scbg.D10.5_2i_C2, 
   scbg.D10.5_serum_C1, scbg.D10.5_serum_C2,  scbg.D11_2i_C1, 
   scbg.D11_2i_C2, scbg.D11_serum_C1, scbg.D11_serum_C2, 
   scbg.D11.5_2i_C1, scbg.D11.5_2i_C2, scbg.D11.5_serum_C1, 
   scbg.D11.5_serum_C2, scbg.D12_2i_C1, scbg.D12_2i_C2, 
   scbg.D12_serum_C1, scbg.D12_serum_C2, scbg.D12.5_2i_C1, 
   scbg.D12.5_2i_C2, scbg.D12.5_serum_C1, 
   scbg.D12.5_serum_C2, scbg.D13_2i_C1, scbg.D13_2i_C2, 
   scbg.D13_serum_C1, scbg.D13_serum_C2, scbg.D13.5_2i_C1, 
   scbg.D13.5_2i_C2, scbg.D13.5_serum_C1, 
   scbg.D13.5_serum_C2, scbg.D14_2i_C1, scbg.D14_2i_C2, 
   scbg.D14_serum_C1, scbg.D14_serum_C2, scbg.D14.5_2i_C1, 
   scbg.D14.5_2i_C2, scbg.D14.5_serum_C1, 
   scbg.D14.5_serum_C2, scbg.D15_2i_C1, scbg.D15_2i_C2, 
   scbg.D15_serum_C1, scbg.D15_serum_C2, scbg.D15.5_2i_C1, 
   scbg.D15.5_2i_C2, scbg.D15.5_serum_C1, 
   scbg.D15.5_serum_C2, scbg.D16_2i_C1, scbg.D16_2i_C2, 
   scbg.D16_serum_C1, scbg.D16_serum_C2, scbg.D16.5_2i_C1, 
   scbg.D16.5_2i_C2, scbg.D16.5_serum_C1, 
   scbg.D16.5_serum_C2, scbg.D17_2i_C1, scbg.D17_2i_C2, 
   scbg.D17_serum_C1, scbg.D17_serum_C2, scbg.D17.5_2i_C1, 
   scbg.D17.5_2i_C2, scbg.D17.5_serum_C1, 
   scbg.D17.5_serum_C2, scbg.D18_2i_C1, scbg.D18_2i_C2, 
   scbg.D18_serum_C1, scbg.D18_serum_C2, scbg.DiPSC_2i_C1, 
   scbg.DiPSC_2i_C2, scbg.DiPSC_serum_C1, scbg.DiPSC_serum_C2)
gc(reset = TRUE)

scbg <- JoinLayers(scbg, assay = "RNA")
saveRDS(scbg, file = "scbg_filtered.rds")


# Process merged Seurat object
scbg <- NormalizeData(scbg)
scbg <- FindVariableFeatures(scbg)
scbg <- ScaleData(scbg, vars.to.regress = "percent.mt")
scbg <- RunPCA(scbg)
scbg <- RunUMAP(scbg, dims = 1:30)
scbg <- FindNeighbors(scbg, dims = 1:30)
scbg <- FindClusters(scbg, resolution = 0.5)

# Save final merged Seurat object
saveRDS(scbg, file = "scbg_final.rds")


# DimPlot
dp.timepoint <- DimPlot(
  object = scbg, reduction = "umap", 
  group.by = "Timepoint", raster = TRUE, 
  cols = c("#8B0000", "#CA4B2C", "#F8766D", "#DE8C00", 
           "#D8C920", "#A8C709", "#094500", "#00BA38", 
           "#1DDED0", "#224AD9", "#000168", "#6E15D4", 
           "#C77CFF", "#FB61DE", "#C0BC87", "#9F7A68", 
           "#4E2901", "#A9A7AE", "#4C4A4B", "#000000")
) + labs(
  x = "UMAP 1", y = "UMAP 2", 
  title = "Treatment Time"
) + theme(
  text = element_text(size = 26), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(0.75, "cm")
)

png(
  filename = "../../Images/MouseReprogramming/scbg_UMAP_timepoint.png", 
  width = 1320, height = 1000, res = 160
)
dp.timepoint
dev.off()

## FeaturePlots
# Pluripotency
fp.nanog <- FeaturePlot(
  object = scbg, features = "Nanog", 
  order = TRUE, raster = TRUE
) + labs(x = "UMAP 1", y = "UMAP 2") + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)
fp.dppa5a <- FeaturePlot(
  object = scbg, features = "Dppa5a", 
  order = TRUE, raster = TRUE
) + labs(x = "UMAP 1", y = "UMAP 2") + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)
fp.tdgf1 <- FeaturePlot(
  object = scbg, features = "Tdgf1", 
  order = TRUE, raster = TRUE
) + labs(x = "UMAP 1", y = "UMAP 2") + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)
fp.dusp9 <- FeaturePlot(
  object = scbg, features = "Dusp9", 
  order = TRUE, raster = TRUE
) + labs(x = "UMAP 1", y = "UMAP 2") + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)

# Stromal
fp.col3a1 <- FeaturePlot(
  object = scbg, features = "Col3a1", 
  order = TRUE, raster = TRUE
) + labs(x = "UMAP 1", y = "UMAP 2") + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)
fp.mmp2 <- FeaturePlot(
  object = scbg, features = "Mmp2", 
  order = TRUE, raster = TRUE
) + labs(x = "UMAP 1", y = "UMAP 2") + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)
fp.fmod <- FeaturePlot(
  object = scbg, features = "Fmod", 
  order = TRUE, raster = TRUE
) + labs(x = "UMAP 1", y = "UMAP 2") + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)
fp.cd47 <- FeaturePlot(
  object = scbg, features = "Cd47", 
  order = TRUE, raster = TRUE
) + labs(x = "UMAP 1", y = "UMAP 2") + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)

# Neural
fp.nnat <- FeaturePlot(
  object = scbg, features = "Nnat", 
  order = TRUE, raster = TRUE
) + labs(x = "UMAP 1", y = "UMAP 2") + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)
fp.pou3f2 <- FeaturePlot(
  object = scbg, features = "Pou3f2", 
  order = TRUE, raster = TRUE
) + labs(x = "UMAP 1", y = "UMAP 2") + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)
fp.sox11 <- FeaturePlot(
  object = scbg, features = "Sox11", 
  order = TRUE, raster = TRUE
) + labs(x = "UMAP 1", y = "UMAP 2") + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)

# Trophoblast
fp.gata3 <- FeaturePlot(
  object = scbg, features = "Gata3", 
  order = TRUE, raster = TRUE
) + labs(x = "UMAP 1", y = "UMAP 2") + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)
fp.phlda2 <- FeaturePlot(
  object = scbg, features = "Phlda2", 
  order = TRUE, raster = TRUE
) + labs(x = "UMAP 1", y = "UMAP 2") + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)
fp.cited2 <- FeaturePlot(
  object = scbg, features = "Cited2", 
  order = TRUE, raster = TRUE
) + labs(x = "UMAP 1", y = "UMAP 2") + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)

# Cell cycle arrest genes
fp.cdkn1a <- FeaturePlot(
  object = scbg, features = "Cdkn1a", 
  order = TRUE, raster = TRUE
) + labs(x = "UMAP 1", y = "UMAP 2") + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)
fp.cdkn2a <- FeaturePlot(
  object = scbg, features = "Cdkn2a", 
  order = TRUE, raster = TRUE
) + labs(x = "UMAP 1", y = "UMAP 2") + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)
fp.trp53 <- FeaturePlot(
  object = scbg, features = "Trp53", 
  order = TRUE, raster = TRUE
) + labs(x = "UMAP 1", y = "UMAP 2") + theme(
  text = element_text(size = 30), 
  axis.text = element_blank(), 
  axis.ticks = element_blank(), 
  legend.key.size = unit(1.25, "cm")
)

# Save feature plots
png(
  filename = "../../Images/MouseReprogramming/scbg_FeaturePlot_NANOG.png", 
  width = 1250, height = 1000, res = 140
)
fp.nanog
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_FeaturePlot_DPPA5A.png", 
  width = 1250, height = 1000, res = 140
)
fp.dppa5a
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_FeaturePlot_TDGF1.png", 
  width = 1250, height = 1000, res = 140
)
fp.tdgf1
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_FeaturePlot_DUSP9.png", 
  width = 1250, height = 1000, res = 140
)
fp.dusp9
dev.off()

png(
  filename = "../../Images/MouseReprogramming/scbg_FeaturePlot_COL3A1.png", 
  width = 1250, height = 1000, res = 140
)
fp.col3a1
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_FeaturePlot_MMP2.png", 
  width = 1250, height = 1000, res = 140
)
fp.mmp2
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_FeaturePlot_FMOD.png", 
  width = 1250, height = 1000, res = 140
)
fp.fmod
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_FeaturePlot_CD47.png", 
  width = 1250, height = 1000, res = 140
)
fp.cd47
dev.off()

png(
  filename = "../../Images/MouseReprogramming/scbg_FeaturePlot_NNAT.png", 
  width = 1250, height = 1000, res = 140
)
fp.nnat
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_FeaturePlot_POU3F2.png", 
  width = 1250, height = 1000, res = 140
)
fp.pou3f2
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_FeaturePlot_SOX11.png", 
  width = 1250, height = 1000, res = 140
)
fp.sox11
dev.off()

png(
  filename = "../../Images/MouseReprogramming/scbg_FeaturePlot_GATA3.png", 
  width = 1250, height = 1000, res = 140
)
fp.gata3
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_FeaturePlot_PHLDA2.png", 
  width = 1250, height = 1000, res = 140
)
fp.phlda2
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_FeaturePlot_CITED2.png", 
  width = 1250, height = 1000, res = 140
)
fp.cited2
dev.off()

png(
  filename = "../../Images/MouseReprogramming/scbg_FeaturePlot_CDKN1A.png", 
  width = 1250, height = 1000, res = 140
)
fp.cdkn1a
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_FeaturePlot_CDKN2A.png", 
  width = 1250, height = 1000, res = 140
)
fp.cdkn2a
dev.off()
png(
  filename = "../../Images/MouseReprogramming/scbg_FeaturePlot_P53.png", 
  width = 1250, height = 1000, res = 140
)
fp.trp53
dev.off()



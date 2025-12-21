# =========================
# Heatmap con TPM
# =========================


library(pheatmap)

# Cargar datos
tpm <- read.csv("tpm_gene_matrix.csv", row.names = 1)
res <- read.csv("DESeq2_significant_Obeso1_vs_Obeso2.csv")

# Genes significativos
genes <- intersect(res$gene, rownames(tpm))

# Matriz 5x5
mat <- tpm[genes, ]

# Anotación de columnas 
annotation_col <- data.frame(
  Grupo = c("Obeso1","Obeso1","Obeso2","Obeso2","Obeso2")
)
rownames(annotation_col) <- colnames(mat)  

# Heatmap 
pheatmap(
  log2(mat + 1),
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  main = "Heatmap de genes diferenciales (TPM)"
)


#=========================
#  Volcano
# ========================

# Cargar librerias
library(ggplot2)

# Cargar resultados DESeq2 (todos los genes)
res_all <- read.csv("DESeq2_all_genes_Obeso1_vs_Obeso2.csv")

#Limpiar datos y definir significancia
# Quitar filas sin p ajustado
res_all <- res_all[!is.na(res_all$padj), ]

# Definir genes significativos
res_all$significant <- "No significativo"
res_all$significant[
  res_all$padj < 0.05 & abs(res_all$log2FoldChange) > 1
] <- "Significativo"
  
# Grafica 
ggplot(res_all, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significant), alpha = 0.7, size = 1.5) +
  geom_text_repel(
    data = subset(res_all, significant == "Significativo"),
    aes(label = gene),
    size = 4,
    max.overlaps = 20
  ) +
  scale_color_manual(values = c("grey70", "red")) +
  theme_minimal() +
  labs(
    title = "Volcano plot Obeso1 vs Obeso2",
    x = "log2 Fold Change",
    y = "-log10(p ajustado)"
  )


#=========================
#  Tabla de expresion de genes 
# ========================
  
# Carga de librerias 
library(dplyr)
library(gridExtra)

# Cargar resultados DESeq2
res_all <- read.csv("DESeq2_all_genes_Obeso1_vs_Obeso2.csv")

# Seleccionar los genes significativos 
genes_clave <- c("UBR2", "FTO", "SIM1", "LEPR", "MC4R")

# Tabla 1
tabla_expresion <- res_all %>%
  filter(gene %in% genes_clave) %>%
  select(
    Gen = gene,
    log2FoldChange,
    padj
  ) %>%
  arrange(desc(abs(log2FoldChange)))

# Imagen para poster 

png("Tabla_expresion_5_genes.png", width = 2000, height = 800, res = 300)

grid.table(tabla_expresion)

dev.off()


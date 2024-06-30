#!/usr/bin/Rscript

# Carregar os dados
beta <- read.csv("/home/jsilva/drought_results/postprocessing/beta-div.tsv", sep="\t")
beta$comparison1 <- as.character(beta$comparison1)
beta$comparison2 <- as.character(beta$comparison2)
bc_mat <- beta[, c(2, 3, 4)]

# Obter nomes de amostras únicas
samples <- unique(c(bc_mat$comparison1, bc_mat$comparison2))

# Criar uma matriz vazia com NAs
triangular_matrix <- matrix(NA, nrow = length(samples), ncol = length(samples), dimnames = list(samples, samples))

# Preencher a matriz com os valores das distâncias
for (i in 1:nrow(bc_mat)) {
  row_name <- bc_mat$comparison1[i]
  col_name <- bc_mat$comparison2[i]
  value <- bc_mat$braycurtis[i]
  triangular_matrix[row_name, col_name] <- value
}

# Substituir valores ausentes (diagonal e triângulo inferior) por 0
triangular_matrix[is.na(triangular_matrix) | lower.tri(triangular_matrix)] <- 0

# Realizar a análise de PCoA
pcoa <- as.data.frame(cmdscale(triangular_matrix, k = 2))
names(pcoa) <- c("PCoA1", "PCoA2")

# Extrair rótulos para as amostras e categorizar como "N" ou "P"
categories <- ifelse(grepl("N", rownames(pcoa)), "N", ifelse(grepl("P", rownames(pcoa)), "P", "Other"))

# Definir as cores para as categorias
colors <- ifelse(categories == "N", "red", ifelse(categories == "P", "blue", "black"))

# Configurar layout do gráfico
layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))

# Plotar PCoA
par(mar = c(5, 4, 4, 1) + 0.1)
plot(pcoa$PCoA1, pcoa$PCoA2, col = colors, xlab = "PCoA1", ylab = "PCoA2", pch = 19, main = "PCoA Plot")
text(pcoa$PCoA1, pcoa$PCoA2, labels = rownames(pcoa), pos = 3, cex = 0.7)

# Plotar a legenda fora do gráfico
par(mar = c(5, 0, 4, 1) + 0.1)
plot.new()
legend("center", legend = c("N", "P", "E"), col = c("red", "blue", "black"), pch = 19, bty = "n")

---
title: h5ad 文件转 rds 文件
date: 2026-01-14
updated: 2026-01-14
description: 单细胞数据文件格式的转换.
top_img: https://i.ytimg.com/vi/como93CmnS8/maxresdefault.jpg
cover: https://i.ytimg.com/vi/como93CmnS8/maxresdefault.jpg
categories:
  - Tutorial
tags: 
  - Single Cell
  - Data
  - R
  - Python
comments: true
---

在单细胞测序数据分析中，常见的数据文件格式有 h5ad（AnnData 格式）和 rds（R 数据格式）.前者主要用于 Python 环境下的单细胞数据处理（如 `scanpy`），而后者则应用于 R 语言环境中（如 `Seurat`）.本文将结合我进行的单细胞注释分析，介绍如何将 h5ad 文件转换为 rds 文件，以便在 R 环境中进行后续分析.

## 数据准备与转换工作

假设我们在 Python 环境下做完细胞注释后的单细胞数据保存为 `adata.h5ad` 文件，里面保存有细胞的表达矩阵、注释信息等.但后续我们需要在 R 环境中使用 `Seurat` 包进行进一步分析，因此需要将该文件转换为 rds 格式。转换使用的工具为 `sceasy` 这个 R package.具体代码如下：

```R
library(sceasy)
library(reticulate)
loompy <- reticulate::import('loompy')

# 读取 h5ad 文件并转换为 rds 文件
h5ad_file = "adata.h5ad"
sceasy::convertFormat(h5ad_file, from="anndata", to="seurat",
                       outFile='Rdata.rds')

# 验证转换结果
r <- readRDS("Rdata.rds")
packageVersion("Seurat")
r@version
# [1] '4.3.0'
```

## 注释信息的修改

此时，细胞注释信息（`celltype`）存放于转换后的 `r@meta.data` 中。

```R
# 查看 meta.data 前几行
head(r@meta.data)

# 查看 celltype 列的注释
unique(r$celltype)
# Photoreceptors Unknow Amacrine_Cells Bipolar+Muller Retinal_Ganglion_Cells
# Levels:
# 'Amacrine_Cells''Photoreceptors''Bipolar+Muller''Retinal_Ganglion_Cells''Unknow'
``` 

这里我们发现有一个注释为 `Unknow` 的细胞类型，这里我们替换为 `Unknown`：

```R
# 替换 Unknow 为 Unknown
r$celltype <- gsub("Unknow", "Unknown", r$celltype)

# 查看替换后的结果
unique(r$celltype)
# 'Photoreceptors''Unknown''Amacrine_Cells''Bipolar+Muller''Retinal_Ganglion_Cells'

# 保存修改后的 rds 文件
saveRDS(r, file = "Rdata_updated.rds")
```

## 注释信息的添加

虽然我们的 `Rdata_updated.rds` 文件中已经包含了细胞类型注释信息，但是在经由 `h5ad` 文件转换后，其 `@images[["slice1"]]` 属性为空，这会影响后续 R 环境中的空间可视化分析.因此我们使用细胞分割后注释分析前的 `SN.rds` （SN 为芯片号 Serial Number）文件做分析，并向其中手动添加注释信息：

```R
library(Seurat)

seu <- readRDS("Rdata_updated.rds")
obj <- readRDS("SN.rds")

# 由于注释过程会舍弃部分细胞，这里我们只保留在 seu 中存在的细胞
cells_to_keep <- rownames(seu@meta.data)
obj <- subset(obj, cells = cells_to_keep)

# 将 seu 中的 celltype 注释信息添加到 obj 中
obj <- AddMetaData(object = obj, metadata = seu@meta.data["celltype"])

# 查看添加后的注释信息
unique(obj$celltype)
# 'Photoreceptors''Unknown''Amacrine_Cells''Bipolar+Muller''Retinal_Ganglion_Cells'

# 检查细胞数量是否一致
ncol(seu)
ncol(obj)
# 17468
# 17468

# 保存最终的 rds 文件
saveRDS(obj, file = "final.rds")
```

通过以上步骤，我们成功地将 h5ad 文件转换为 rds 文件，并对细胞注释信息进行了修改和添加.这样，我们就可以在 R 环境中继续进行单细胞数据的其他分析和可视化工作了.
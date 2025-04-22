这个代码执行了一个生物信息学分析流程，用于分析基因表达数据，并生成火山图、热图和箱线图。以下是详细的实验步骤：

1. **加载必要的R包:** 加载tidyverse (用于数据处理和可视化), ggpubr (用于增强ggplot2图形), future (用于并行计算), rio (用于导入各种数据格式), pheatmap (用于绘制热图)。

2. **设置并行计算:** 使用future::plan("multisession", workers = 6)启用并行计算，加快运行速度。这里设置了6个工作线程。

3. **导入差异基因表达分析结果:** 导入名为 "GSE206510.top.table.tsv" 的文件，该文件包含差异基因表达分析的结果 (DEG)。文件里应该包含基因ID (GeneID), 基因名 (Symbol), log2倍数变化 (log2FoldChange) 和 p值 (pvalue) 等信息。

4. **数据预处理及分组:**

    - 将导入的数据存储到 voldata 中。
    - 根据p值和log2倍数变化对基因进行分组：
        - p值 >= 0.05: "not sig" (不显著)
        - log2FoldChange >= 0.58: "up" (上调)
        - log2FoldChange < -0.58: "down" (下调)
    - 将分组结果转换为因子类型，并指定因子水平顺序为 "up", "down", "not sig"。
    - 创建 logFC 和 p.val 列，分别存储 log2FoldChange 和 -log10(pvalue)。

5. **绘制火山图:**

    - 筛选显著差异表达的基因 (group != "not sig")，并提取基因名。

    - 使用 ggscatter 函数绘制火山图，根据基因分组设置点的颜色。

    - 添加水平线和垂直线，表示显著性阈值 (pvalue = 0.05, log2FoldChange = ±0.58)。

    - 添加标题和标签，并使用 ggsave 函数保存火山图为 "volcano.pdf" 文件。

        ![image-20250422221830704](https://p.ipic.vip/rw38fz.png)

6. **导入基因表达矩阵并预处理:**

    - 导入归一化后的基因表达矩阵 "GSE206510_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"。
    - 导入基因注释文件 "Human.GRCh38.p13.annot.tsv.gz"，该文件包含 GeneID 和 Symbol 的对应关系。
    - 过滤掉低表达基因 (在超过一半样本中表达量为0的基因)。
    - 将 GeneID 转换为基因名 (Symbol)。
    - 将数据转置，进行标准化 (z-score)，再转置回来。

7. **准备热图数据和注释:**

    - 从标准化后的表达矩阵中提取差异表达基因的表达数据。
    - 创建样本分组信息，用于热图的列注释。这里示例代码中使用了"Control"和"Asthma"两组，每组10个样本。

8. **绘制热图:**

    - 使用 pheatmap 函数绘制热图。
    - 设置参数 show_rownames = FALSE 和 show_colnames = F 不显示行名和列名。
    - 设置参数 cluster_cols = F 不进行列聚类。
    - 使用预先定义的颜色和断点，使热图颜色更美观。
    - 将热图保存为 "heatmap.pdf" 文件。
    - ![image-20250422221858665](https://p.ipic.vip/bdkax4.png)

9. **准备箱线图数据:**

    - 重新导入基因表达矩阵。
    - 将数据转换为长格式，方便绘制箱线图。
    - 添加样本分组信息。

10. **绘制箱线图 (以FH基因为例):**

    - 筛选 FH 基因的表达数据。
    - 使用 ggboxplot 函数绘制箱线图，比较两组间的表达差异。
    - 使用 stat_compare_means 函数进行 t 检验，并在图上显示显著性标记。
    - 添加标题和标签，并使用 ggsave 函数保存箱线图为 "FH.pdf" 文件。

![image-20250422221913274](https://p.ipic.vip/4mg6ii.png)

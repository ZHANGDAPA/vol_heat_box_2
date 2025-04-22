library(tidyverse)
library(ggpubr)
future::plan("multisession", workers = 6)
deg <- rio::import("GSE206510.top.table.tsv")
# deg <- rio::import("GSE207751.top.table.tsv")

# visualize
voldata <- deg
voldat <- voldata %>%
    mutate(group = case_when(
        pvalue >= 0.05 ~ "not sig",
        log2FoldChange >= 0.58 ~ "up",
        log2FoldChange < -0.58 ~ "down",
        TRUE ~ "not sig"
    )) %>%
    mutate(group = factor(group, levels = c("up", "down", "not sig"))) %>%
    mutate(logFC = log2FoldChange, p.val = -log10(pvalue))

name <- voldat %>%
    filter(group != "not sig") %>%
    select(-GeneID) %>%
    pull("Symbol")

table(voldat$group)
ggscatter(
    voldat,
    x = "logFC", y = "p.val",
    color = "group",
    alpha = 0.85,
    palette = c("#D01910", "#00599F", "#CCCCCC"),
    repel = TRUE,
    xlab = "log2FoldChange", ylab = "-log10(P.value)",
) +
    theme_bw() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#222222") +
    geom_vline(xintercept = 0.58, linetype = "dashed", color = "#222222") +
    geom_vline(xintercept = -0.58, linetype = "dashed", color = "#222222") + # nolint # nolint: line_length_linter.
    theme(plot.background = element_blank()) +
    labs(subtitle = "FC:0.58; p.val:0.05; up:439; down:135: not sig:14528") +
    theme_pubr()

ggsave("volcano.pdf", width = 8, height = 8)
dev.off()

## heatmap
library(pheatmap)
# raw <- rio::import("GSE207751_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz")
raw <- rio::import("GSE206510_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz")
ids <- rio::import("Human.GRCh38.p13.annot.tsv.gz")
str(raw)
raw <- raw[apply(raw, 1, function(x) sum(x > 0) > 0.5 * ncol(raw)), ]

exp <- raw %>%
    inner_join(ids[, 1:2], by = "GeneID") %>%
    select(-GeneID) %>%
    distinct(Symbol, .keep_all = TRUE) %>%
    column_to_rownames("Symbol") %>%
    t() %>%
    scale() %>%
    t()

range(exp)

# pd = rio::import("GSE207751_whole_blood_metadata.csv.gz")

# annotation_col <- data.frame(group = pd$condition)
# rownames(annotation_col) <- colnames(exp)
# annotation_col$index = rownames(annotation_col)
# annotation_col = annotation_col %>%
#     arrange(group)
# exp = exp[, annotation_col$index]
# match(colnames(exp), annotation_col$index)
# annotation_col$index = NULL

ncol(exp)
annotation_col <- data.frame(group = rep(c("Control", "Asthma"), each = 10))
rownames(annotation_col) <- colnames(exp)

color <- c(
    "#10467e",
    "#3172a7",
    "#69a7cd",
    "#ffffff",
    "#dc6559",
    "#b52634",
    "#690619"
)


exp2 <- exp[name, ]
View(head(exp))
a <- as.numeric(unlist(exp2))
a22 <- boxplot(a, plot = F)$out
a3 <- range(a[!a %in% a22])

pdf("pic/heatmap.pdf", width = 15, height = 15)
pheatmap(
    exp2,
    show_rownames = FALSE,
    show_colnames = F,
    cluster_cols = F,
    annotation_col = annotation_col,
    # breaks = seq(a3[1],a3[2], length.out = 100),
    # color = color,
    color = colorRampPalette(color)(100),
    breaks = seq(a3[1], a3[2], length.out = 100),
    # file = "full_heatmap.png",
    # width = 15,
    # height = 15,
    # res = 300
)
dev.off()

pd <- annotation_col
pd <- pd %>%
    rownames_to_column("sample")

## boxplot
## 绘制箱线图
raw <- rio::import("GSE206510_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz")
BOX <- raw %>%
    inner_join(ids[, 1:2], by = "GeneID") %>%
    select(-GeneID) %>%
    distinct(Symbol, .keep_all = TRUE) %>%
    column_to_rownames("Symbol")
View(head(BOX))

box <- BOX %>%
    as.data.frame() %>%
    rownames_to_column("Symbol") %>%
    pivot_longer(-Symbol, names_to = "sample", values_to = "TPM") %>%
    inner_join(pd, by = "sample") %>%
    select(-sample)


FH <- box %>%
    filter(Symbol == "FH")

FH$group <- factor(FH$group, levels = c("Control", "Asthma"))
ggboxplot(FH,
    x = "group",
    y = "TPM",
    color = "group",
    width = 0.5,
    bxp.errorbar = T,
    palette = c("#00599F", "#D01910"),
    add = "jitter"
) +
    stat_compare_means(
        comparisons = list(c("Control", "Asthma")),
        method = "t.test",
        label = "p.signif"
    ) +
    labs(
        title = "FH",
        y = "log2(Intensity)"
    ) +
    theme_pubr()
ggsave("pic/FH.pdf", width = 4, height = 8, dpi = 300)

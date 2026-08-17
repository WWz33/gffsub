# gffsub 文档

<!-- I18N:START -->

[English](./Home.md) | **中文**

<!-- I18N:END -->

按 region、feature ID、属性或 gene model 从 GFF3 和 GTF 注释文件中提取子集。

## 快速开始

- [安装与快速上手](getting-started.zh.md)

## 筛选

1. [Region 与 BED 筛选](region-bed-filtering.zh.md)：`-r`、`-b`、`--seqid`
2. [Feature ID 与基因名查找](id-name-lookup.zh.md)：`--id`、`--ids`、`--name`
3. [Gene Model 展开](gene-model-expansion.zh.md)：`-C`、`--parents`、`--model`
4. [属性筛选](attribute-filtering.zh.md)：`--where`、`--grep`、`--grep-regex`、`--grep-file`
5. [表达式筛选](expression-filtering.zh.md)：`-I`、`-E`

## 选择

- [最长转录本选择](longest-isoform.zh.md)：`--longest`、`-@`
- [window 子命令](window-subcommand.zh.md)：`window`、`--up`、`--down`、`--strand-aware`

## 输出

- [输出格式](output-formats.zh.md)：GFF3、GTF、GTF2、GTF3、BED（`-t`）
- [汇总输出](summary-output.zh.md)：TSV 与 JSON 表格汇总（`--summary`、`--out-attrs`）

## 输入

- [GTF 输入处理](gtf-input.zh.md)：GTF 属性解析、Parent/ID 合成、isoform 自动检测

## 库

- [库 API](library-api.zh.md)：链接 `libgffsub_core.a`，从 C++17 调用 `query()`、`window()`、`subset()`

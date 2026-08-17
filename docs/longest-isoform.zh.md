# 最长转录本选择

<!-- I18N:START -->

[English](./longest-isoform.md) | **中文**

<!-- I18N:END -->

## 概述

`-L` / `--longest` 标志为每个 gene 保留一个 isoform。所选类型下至少有两个 isoform 的 gene 归约到最长者。gene 含 CDS 时按 CDS 跨度评分，否则按 exon 跨度评分。未选中的 isoform 及其子 feature（exon、CDS）被丢弃。只有一个 isoform 的 gene 保持不变。所选类型下无 isoform 的 gene 被丢弃。

## 算法

1. 按 Parent gene 将 isoform 类型记录（mRNA 或 transcript）分组。
2. 按 chromosome 将 gene 分组。
3. 对每个 gene：
   - 若 gene 没有所选类型的 isoform 子记录，丢弃该 gene。
   - 若 gene 恰好有一个 isoform，保持该 isoform 及其子记录不变。
   - 否则计算 gene 级标志：该 gene 的任一 isoform 是否有 CDS 子记录。
   - 对每个 isoform 计算跨度：
     - gene 有 CDS 时：CDS 长度。一个 transcript 下 ID 不同的 CDS 属于不同蛋白变体；isoform 按其最长变体评分（共享同一 CDS ID 的多行属同一条不连续 CDS，相加求和）。无 CDS 子记录的 isoform 被跳过，无法被选中。
     - gene 无 CDS 时：exon 长度之和。无 exon 子记录的 isoform 被跳过，无法被选中。
   - 保留跨度最大的 isoform。并列时保留首个遇到的（严格大于比较）。
   - 丢弃未保留的 isoform 及其子记录。保留 isoform 的子记录不变。
4. 线程（`-@` / `--threads`）按 chromosome 并行。每个 chromosome 独立处理。

## isoform 类型自动检测

未指定 `-f` / `--feature` 时：
- 数据中有 `mRNA` 记录则用 `mRNA`。
- 否则数据中有 `transcript` 记录则用 `transcript`。
- 两者都没有则默认 `mRNA`。此时所有 gene 因无匹配 isoform 而被丢弃。

覆盖 GFF3（通常 mRNA）和 GTF（通常 transcript）。

指定 `-f TYPE` 时，该类型即为 isoform 类型，运行还会丢弃所有 type 不为 TYPE 的记录（gene、CDS、exon 及其他）。最长选择先执行，随后 feature 筛选只保留 TYPE 记录。

## 示例数据

```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	CDS	100	250	.	+	0	ID=cds01;Parent=tx01
chr1	src	CDS	500	750	.	+	2	ID=cds02;Parent=tx01
chr1	src	mRNA	100	800	.	+	.	ID=tx02;Parent=gene01
chr1	src	CDS	100	400	.	+	0	ID=cds03;Parent=tx02
chr2	src	gene	200	600	.	-	.	ID=gene02
chr2	src	mRNA	200	600	.	-	.	ID=tx03;Parent=gene02
chr2	src	exon	200	400	.	-	.	ID=ex04;Parent=tx03
chr2	src	exon	450	600	.	-	.	ID=ex05;Parent=tx03
```

tx01 CDS 跨度：151 + 251 = 402 bp。
tx02 CDS 跨度：301 bp。
tx03 exon 跨度：201 + 151 = 352 bp（无 CDS）。

## 命令

```bash
# 每个 gene 保留一个 isoform（gene01 选 tx01，gene02 选 tx03）
./gffsub demo.gff3 --longest

# 使用 4 线程
./gffsub demo.gff3 --longest -@ 4

# 限定特定 isoform 类型
./gffsub demo.gff3 --longest -f transcript
```

## 示例结果

- gene01：保留 tx01（CDS 跨度 402），丢弃 tx02 和 cds03。
- gene02：保留 tx03，不变（只有一个 isoform）。
- 至少有一个所选类型 isoform 的 gene 记录被保留。

## 与其他筛选组合

region、feature 类型和属性筛选先于最长选择执行。在指定区域内做最长选择：

```bash
./gffsub demo.gff3 -r chr1:1-1000 --longest
```

# Region 与 BED 筛选

<!-- I18N:START -->

[English](./region-bed-filtering.md) | **中文**

<!-- I18N:END -->

## Region 筛选：-r / --region

语法：`-r CHR:START-END`

- 1-based 闭区间坐标
- 保留与同一 seqid 上 region 重叠的记录
- 重叠指记录跨度与查询 region 在碱基层面的任意交集
- 与 `-f`（feature 类型）组合可进一步限定

示例数据（保存为 demo.gff3）：
```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr2	src	gene	200	600	.	-	.	ID=gene02
```

命令：
1. `./gffsub demo.gff3 -r chr1:200-600` — 与 200-600 重叠的所有 chr1 记录
2. `./gffsub demo.gff3 -r chr1:200-600 -f gene` — 该 region 内的 gene 记录
3. `./gffsub demo.gff3 -r chr2:1-1000` — 所有 chr2 记录

## BED 文件筛选：-b / --bed

语法：`-b FILE.bed`

- BED 为 0-based 半开区间（start 包含，end 不包含）
- 文件格式：制表符分隔，至少 3 列（chrom、start、end）
- 可选列：name（第 4 列）、score（第 5 列）、strand（第 6 列），解析但不用于筛选
- 保留与 BED 文件中任一区间重叠的记录
- 同一文件内多个区间为 OR 关系

示例 BED 文件（保存为 regions.bed）：
```
chr1	99	800
chr2	199	700
```

命令：
1. `./gffsub demo.gff3 -b regions.bed` — 与任一 BED 区间重叠的记录
2. `./gffsub demo.gff3 -b regions.bed -f gene` — 仅 gene 记录

## Seqid 筛选：--seqid

语法：`--seqid LIST`

- 逗号分隔的 seqid 列表
- 仅保留列出的 seqid 上的记录
- 前缀 `^` 表示排除：`--seqid ^chr1` 保留除 chr1 外的全部记录
- 与 region/BED 筛选为 AND 关系

命令：
1. `./gffsub demo.gff3 --seqid chr1` — 仅 chr1 记录
2. `./gffsub demo.gff3 --seqid chr1,chr2` — chr1 和 chr2 记录
3. `./gffsub demo.gff3 --seqid ^chr1` — 除 chr1 外的全部记录

## 坐标系统

| 格式 | 基数 | 坐标约定 |
|--------|------|-----------------------|
| GFF3/GTF | 1-based | start 和 end 均包含 |
| BED | 0-based | start 包含，end 不包含（半开） |
| Region 查询（-r） | 1-based | start 和 end 均包含 |

gffsub 内部转换 BED 坐标。各输入类型按其原生格式提供坐标。

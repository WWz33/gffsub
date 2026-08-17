# Region 与 BED 过滤

<!-- I18N:START -->

[English](./region-bed-filtering.md) | **中文**

<!-- I18N:END -->

## -r / --region

语法：`-r CHR:START-END`

- 1-based 闭区间坐标
- 保留与同一 seqid 上该 region 有重叠的记录

示例数据（`demo.gff3`）：

```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr2	src	gene	200	600	.	-	.	ID=gene02
```

```bash
./gffsub demo.gff3 -r chr1:200-600
```

输出：

```
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
```

## -b / --bed

语法：`-b FILE.bed`

- BED 为 0-based 半开区间（start 包含，end 不包含）
- 制表符分隔，至少 3 列（chrom、start、end）
- 可选列 4-6（name、score、strand）会解析但不用于过滤
- 保留与任一 BED interval 有重叠的记录

示例 BED 文件（`regions.bed`）：

```
chr1	99	800
chr2	199	700
```

```bash
./gffsub demo.gff3 -b regions.bed
```

## -S / --seqid

语法：`--seqid LIST`

- 逗号分隔的 seqid 列表
- `^` 前缀表示排除：`--seqid ^chr1` 保留除 chr1 外的所有内容
- 与 region/BED 过滤为 AND 关系

```bash
./gffsub demo.gff3 -S chr1
./gffsub demo.gff3 -S chr1,chr2
./gffsub demo.gff3 -S ^chr1
```

## 坐标系统

| 格式 | 基数 | 约定 |
|------|------|------|
| GFF3/GTF | 1-based | start 与 end 均包含 |
| BED | 0-based | start 包含，end 不包含（半开） |
| -r region 查询 | 1-based | start 与 end 均包含 |

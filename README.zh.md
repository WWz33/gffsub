# gffsub

<!-- README-I18N:START -->

[English](./README.md) | **中文**

<!-- README-I18N:END -->

按区间、feature ID、属性或 gene model 提取 GFF3/GTF 注释。选最长转录本，转换格式，输出摘要。

## Getting Started

```bash
git clone https://github.com/WWz33/gffsub.git
cd gffsub && make -j

printf '%s\n' '##gff-version 3' \
  $'chr1\t.\tgene\t1\t100\t.\t+\t.\tID=GeneA;Name=GeneA' \
  $'chr1\t.\tmRNA\t1\t100\t.\t+\t.\tID=mRNA1;Parent=GeneA' \
  $'chr1\t.\texon\t1\t100\t.\t+\t.\tID=exon1;Parent=mRNA1' > data/smoke.gff3

./gffsub data/smoke.gff3 -r chr1:1-100 -f gene
./gffsub data/smoke.gff3 --id GeneA -C
```

## Usage

```text
gffsub <input.gff3> [options]
gffsub query|window <input.gff3> [options]
```

完整选项见 `./gffsub -h`。

### 常用选项

| 选项 | 默认 | 说明 |
|------|------|------|
| `-r, --region` | — | 1-based 闭区间 `CHR:START-END` |
| `-b, --bed` | — | BED 区间（0-based half-open） |
| `--id` / `--ids` | — | 精确 feature ID |
| `--name` | — | 常见基因命名键 |
| `--where KEY=VALUE` | — | 第 9 列精确属性 |
| `-C, --children` | 关 | 包含选择结果的后代 |
| `--parents` | 关 | 包含祖先 |
| `--model` | 关 | 完整 gene model |
| `--nearest REGION` | — | 同 seqid 最近 gene |
| `--up` / `--down` | 0 | 相对 `--id` 的窗口（bp） |
| `--strand-aware` | 关 | 窗口按特征链方向 |
| `-f, --feature` | — | type 列过滤 |
| `-L, --longest` | 关 | 每基因一条转录本：有 CDS 比 CDS，否则比 exon |
| `-@, --threads` | 1 | `--longest` 线程数 |
| `-t, --format` | gff3 | `gff3\|gtf\|gtf2\|gtf3\|bed` |
| `-o, --output` | stdout | 输出文件 |

## Input / Output

| 类型 | 说明 |
|------|------|
| 注释 | GFF3/GTF 风格 feature |
| GFF/GTF / `-r` | 1-based 闭区间 |
| BED | 0-based half-open |
| `--summary` | `tsv` / `json` 行，而非 feature 记录 |

仅支持注释记录，不支持 FASTA 和 VCF。

## 文档

完整指南见 [`docs/`](docs/Home.md)：

- [快速开始](docs/getting-started.md)
- [区域与 BED 筛选](docs/region-bed-filtering.md)
- [Feature ID 与基因名查找](docs/id-name-lookup.md)
- [基因模型展开](docs/gene-model-expansion.md)
- [属性筛选](docs/attribute-filtering.md)
- [表达式筛选](docs/expression-filtering.zh.md)
- [最长转录本选择](docs/longest-isoform.md)
- [Window 子命令](docs/window-subcommand.md)
- [输出格式](docs/output-formats.md)
- [摘要输出](docs/summary-output.md)
- [GTF 输入处理](docs/gtf-input.md)

## License

MIT License

# gffsub

<!-- README-I18N:START -->

[English](./README.md) | **中文**

<!-- README-I18N:END -->

按区间、feature ID、基因名或属性从 GFF3/GTF 注释文件中提取记录的命令行
工具。支持 GTF 与 BED 输出、每基因最长转录本选择和按 seqid 的统计输出。

## 编译

```bash
git clone https://github.com/WWz33/gffsub.git
cd gffsub && make -j
```

依赖 zlib。输入支持 GFF3、GTF，纯文本或 gzip 压缩均可；`-` 读取 stdin。
BED 输入只用于 BED 输出（`-f bed`）和统计输出（`-s`）；BED 不转换为 GFF3
或 GTF。

## 用法

```bash
./gffsub annot.gff3 -r chr1:1000000-1100000
./gffsub annot.gff3 -i ENSG00000139618 -m
./gffsub annot.gff3 -i ENSG00000139618 -u 2000 -a
./gffsub annot.gff3 -L
./gffsub annot.gff3 -t exon --drop-orphans
./gffsub annot.gff3 -b regions.bed -f gtf
```

完整选项见 `./gffsub -h`；`query`、`window` 子命令各有自己的 `-h`。

## 文档

- [docs/](docs/Home.zh.md)：使用指南，中英双语
- [CHANGELOG.md](CHANGELOG.md)：版本记录

## License

MIT License

# gffsub

<!-- README-I18N:START -->

**English** | [中文](./README.zh.md)

<!-- README-I18N:END -->

Command-line tool for subsetting GFF3/GTF annotation files by region,
feature ID, gene name, or attribute. Supports GTF and BED output, longest
isoform selection, and per-seqid summary output.

## Build

```bash
git clone https://github.com/WWz33/gffsub.git
cd gffsub && make -j
```

Requires zlib. Input can be GFF3 or GTF, plain or gzip-compressed; `-` reads
stdin. BED input is accepted for BED output (`-f bed`) and summaries (`-s`)
only; BED is not converted to GFF3 or GTF.

## Usage

```bash
./gffsub annot.gff3 -r chr1:1000000-1100000
./gffsub annot.gff3 -i ENSG00000139618 -m
./gffsub annot.gff3 -i ENSG00000139618 -u 2000 -a
./gffsub annot.gff3 -L
./gffsub annot.gff3 -t exon --drop-orphans
./gffsub annot.gff3 -b regions.bed -f gtf
```

Run `./gffsub -h` for the full option list. The `query` and `window`
subcommands take their own `-h`.

## Documentation

- [docs/](docs/Home.md): usage guides, English and Chinese
- [CHANGELOG.md](CHANGELOG.md): version history

## License

MIT License

# Defect Review Record

<!-- I18N:START -->

**English** | [中文](./defect-review.zh.md)

<!-- I18N:END -->

Review of commit `f2e171b` and the whole repository, run 2026-11-09, verified against the binary built from that commit. Every entry states the input, the command, and the output observed. Line numbers refer to `f2e171b`.

**Status: every defect below is fixed in the working tree after `f2e171b`.** Each section ends with the fix that changed the observed output; `tests/regression_smoke.cpp` holds one test per defect (`test_defect_regressions`). The commands and outputs in the entries describe the behaviour at `f2e171b`, not current behaviour.

The two observations kept as code facts rather than behaviour are the ones with no user-visible defect at the current input size; see "Recorded but not fixed as behaviour" at the end.

## Verification method

Six independent passes read the source at HEAD, then three further passes reproduced each claim with a fixture under a temporary directory. A claim entered this record only when a command produced the stated output. Claims that failed reproduction were dropped.

All fixtures use GFF3 unless marked GTF. Commands assume the binary is `./gffsub`.

## Parent links

### GTF transcript with a Parent attribute loses its gene

Sample (gl.gtf):

```
chr1	x	gene	1	100	.	+	.	gene_id "g1";
chr1	x	transcript	1	100	.	+	.	gene_id "g1"; transcript_id "t1"; Parent "t1,t2";
chr1	x	transcript	1	200	.	+	.	gene_id "g1"; transcript_id "t2"; Parent "t1,t2";
chr1	x	exon	1	50	.	+	.	gene_id "g1"; transcript_id "t1";
```

```bash
./gffsub gl.gtf --longest
```

Output contains `ID=t1` and `ID=t2`; `g1` is absent. The parser sets `rec.parent_id` from the explicit Parent (`src/gff3_parser.cpp:183`), and the `gene_id` synthesis runs only when `parent_id` is empty (`:220`). `parent_ids()` therefore never returns `g1`.

### GTF row expansion emits the whole Parent list on every row

```bash
./gffsub mp.gff3 -f gtf
```

For `ID=e;Parent=tx1,tx2`, both exon rows carry `Parent "tx1"; Parent "tx2";`. Reading the output back expands again: 2 exon rows become 4 (`tests/regression_smoke.cpp:884` asserted only the substring, so it did not catch this).

### An existing transcript_id overrides the parent of each expanded row

Sample (tc.gtf):

```
chr1	x	transcript	1	100	.	+	.	gene_id "g1"; transcript_id "txA";
chr1	x	transcript	1	100	.	+	.	gene_id "g1"; transcript_id "txB";
chr1	x	exon	1	100	.	+	.	gene_id "g1"; transcript_id "tx1"; Parent "txA,txB";
```

```bash
./gffsub tc.gtf -f gtf
```

Both exon rows print `transcript_id "tx1"`. `resolve_child` copies `rec.transcript_id` first and fills it from the parent only when empty (`src/annotation_output.cpp:263`).

### A shared stop_codon trims only the first parent's CDS

```bash
./gffsub sc.gff3 -f gtf | grep $'\tCDS\t'
```

`t1` CDS ends at 100, `t2` CDS ends at 103. `stop_codons` is keyed by the scalar `rec.parent_id`, which holds only the first parent (`src/annotation_output.cpp:252`).

### `--parents` and `--model` skip ancestors of an ID-less shared exon

```bash
./gffsub query ex.gff3 -w Parent=t2 -p
```

Output is the exon alone. With `ID=e` added to the exon, the same command returns the gene, both transcripts, and the exon. `add_match` expands through `rec.id` (`src/query.cpp:55`), which an ID-less exon does not have.

### Repeated GTF Parent keys disagree between the scalar and the list

Sample: `Parent "t1"; transcript_id "t2"; Parent "t3";`.

`rec.parent_id` is `t1` (`extract_quoted_value` returns the first pair, `src/gtf_parser.cpp:14`); `parent_ids()` returns `{t1,t3}`. `-w Parent=t3` matches the row, `--grep Parent:t3` does not.

### A quoted GTF Parent is not URL-decoded

Sample with transcript IDs `t,1` and `t2`, child carrying `Parent "t%2C1,t2"`.

```bash
./gffsub query p.gtf -w 'Parent=t,1'
```

No rows. The GFF3 form of the same input matches, because `parse_attributes` decodes `%2C` (`src/attributes.cpp:58`) and the GTF branch does not (`:104`).

**Fix.** `parent_ids()` always appends the synthesized `rec.parent_id`; a GTF transcript's synthesized parent is its `gene_id`, and a transcript is never its own parent. GTF child expansion makes the current parent authoritative for that row's `transcript_id` and drops the repeated `Parent` list from an expanded row. `stop_codon` intervals are indexed under every parent. `Add_match` expands an ID-less record through `parents_of_record()`. Quoted GTF Parent values are URL-decoded. `gtf_attrs_to_gff3` writes the full Parent list.

## BED input and output

### Empty BED bypasses the conversion rejection

```bash
: > empty.bed
./gffsub empty.bed
```

Exit 0, output `##gff-version 3`. A non-empty BED exits 1. `source_format()` returns `nullopt` when the index holds no records (`src/annotation_index.cpp:259`), so the guard compares false.

### BED score is written without range or integer checks

```bash
./gffsub scores.bed -f bed
```

Column 5 prints `3.7`, `-2`, `1001`, `abc` unchanged. Only `inf` and `nan` become `0` (`src/annotation_output.cpp:409`). BED requires an integer from 0 to 1000.

### A write failure after a successful open exits 0

```bash
./gffsub in.gff3 -o mnt/out.gff3   # 1 MB FAT image, 7 KB free
```

Exit 0, empty stderr, 0-byte output. The stream state is never checked after printing (`src/gffsub.cpp:511`).

**Fix.** `GffData` records the input format it was parsed as, and `AnnotationIndex` carries it, so an empty BED file is still BED; an empty file with a `.bed` or `.gtf` suffix takes that format. BED score is rounded and clamped to 0..1000. Every output path checks the stream and exits 1 on a write failure.

## Input handling

### A closed stdin is reported as an empty annotation

```bash
./gffsub - <&-
```

Exit 0, output `##gff-version 3`. `read_stream_chunked` loops on `while (in)` and never checks `bad()` (`src/input_source.cpp:106`).

### Extra columns are dropped without a warning

```bash
printf 'chr1\ts\tgene\t1\t9\t.\t+\t.\tID=g1\tjunk\n' > c.gff3
./gffsub c.gff3
```

Output ends at `ID=g1`. The guard accepts every row with nine or more columns (`src/gff3_parser.cpp:137`).

### A UTF-8 BOM becomes part of the first seqid

```bash
printf '\xEF\xBB\xBFchr1\ts\tgene\t1\t9\t.\t+\t.\tID=g1\n' > b.gff3
./gffsub b.gff3 | xxd
```

The first column starts `efbbbf63 687231`. A BOM in front of a `##` directive is skipped as a short line and reported as one.

### Multi-character strand and phase are truncated

```bash
printf 'chr1\ts\tgene\t1\t9\t.\t++\t00\tID=g1\n' > s.gff3
./gffsub s.gff3
```

Output strand `+`, phase `0`. Only the first character is kept (`src/gff3_parser.cpp:165`).

**Fix.** `load_stream` rejects a closed stdin. Rows with more than nine columns (BED: more than 12) are reported and skipped. A UTF-8 BOM is stripped before parsing and before format sniffing. Multi-character strand and phase are rejected and reported.

## Matching and statistics

### `-t` does not trim spaces after commas

```bash
./gffsub f.gff3 -t "mRNA, exon"
```

Exon rows are dropped; `-t "mRNA,exon"` keeps them. `parse_list` does not trim (`src/string_utils.cpp:58`).

### A single `&` inside `-I` becomes part of the value

```bash
./gffsub f.gff3 -I 'type==mRNA&attr.biotype==protein_coding'
```

Header only. With `&&` the record is emitted. The parser splits only on `&&` and `||` (`src/selector_filter.cpp:449`).

Note: `-w` is a `KEY=VALUE` attribute match, not an expression, so this applies to `-I` and `-E`.

### A parent cycle reachable only through a non-first Parent is not reported

Sample: `ID=a;Parent=b,c`, `ID=b;Parent=x`, `ID=c;Parent=a`, plus `ID=x`.

```bash
./gffsub cyc.gff3 -i a
```

No warning. The cycle walk follows the scalar `parent_id` only and stops after 64 steps (`src/annotation_index.cpp:240`).

### A missing parent produces no warning

Sample: `chr1 s exon 10 20 . + . Parent=missing`.

```bash
./gffsub query m.gff3 -w Parent=missing -p
```

Exon emitted, exit 0, empty stderr.

### Summary overflows on coordinates near the int64 limit

```bash
./gffsub big.gff3 -s
```

`sum_len` and `coverage` print as negative values for `[1, 9223372036854775807]` and `[9223372036854775806, 9223372036854775807]`. `cur_end + 1` and the length sum wrap (`src/query_summary.cpp:62`, `:86`).

### Isoform dedup is quadratic

A single gene with N isoforms, `--longest`, user time: 40k 0.42s, 80k 1.60s, 160k 6.56s. Each insert scans the vector with `std::find` (`src/isoform_filter.cpp:55`, `:70`).

### Region and attribute queries scan every record

1M exons, 45 MB, `-r chr1:500-600`: 0.20s, the same as parsing alone. `--where ID=e500`: 1.26s. `overlap` and `with_attribute` loop the full record vector (`src/annotation_index.cpp:364`, `:413`).

**Fix.** `parse_list` trims each item. The Parent cycle check walks every Parent of every line with no depth cutoff. A Parent naming no record is reported once. Summary lengths, sums, coverage and the median saturate instead of wrapping. Isoform and child indices dedup through a `std::unordered_set<int>`.

## Command line

### Extra positional arguments are ignored

```bash
./gffsub in.gff3 typo.gff3
```

Exit 0, `typo.gff3` not read (`src/cli.cpp:506`).

### Subcommands reject `--flag=value`

```bash
./gffsub --id=g1 in.gff3           # accepted
./gffsub query in.gff3 --id=g1     # Error: unknown query option
./gffsub window in.gff3 -i g1 --up=10
```

Subcommand parsers compare whole strings (`src/gffsub.cpp:68`, `:207`).

### Input errors discard the operating-system cause

Missing file, unreadable file, directory, and corrupt gzip all print `Error: cannot parse PATH`. The calls at `src/input_source.cpp:274` return `-1` without saving `errno`.

### A negative thread count passes when prefixed by a space

```bash
./gffsub in.gff3 --longest --threads ' -1'
```

Exit 0. `std::stoul` accepts the leading space, wraps to a large value, and the value is clamped to 256 (`src/cli.cpp:447`).

**Fix.** A second positional operand is a usage error. `query` and `window` split `--flag=value`. Input errors keep `errno`. `--threads` accepts digits only.

## Documentation

### Library API reference does not compile against the headers

| Documented | Actual |
|---|---|
| `SubsetParams::feature` (`docs/library-api.md:145`) | `type` (`src/subset.hpp:28`) |
| `filter_by_feature` (`:163`) | `filter_by_type` (`src/filter.hpp:21`) |
| `QueryParams::feature_type` (`:108`) | `type` (`src/query.hpp:25`) |
| `std::string` record fields (`:23`) | `std::string_view` (`src/record.hpp:21`) |
| `infer_input_format` in `parser.hpp` (`:59`) | `src/input_source.hpp:30` |

`from_owned` and `find_all_genes` are missing from the reference. The Chinese page carries the same errors.

### The documented link command fails

`docs/library-api.md:15` and `:311` omit `-lz`. Linking the documented command reports undefined `_gzopen`, `_gzread`, and `_inflate*`.

### GTF examples omit attributes the writer emits

`docs/output-formats.md:38` states that ID and Parent are kept; the examples at `:76` and `:88` do not contain them. The binary emits them.

**Fix.** The reference matches the headers, the link command carries `-lz`, and the GTF examples show the emitted attributes. `docs/library-api.md` also states the borrowed-index lifetime requirement.

## Recorded but not fixed as behaviour

These observations are code facts without a user-visible defect at the current input size. The gzip byte-read sniffing and the Makefile's `-Wall`-only flags are unchanged; the rest were addressed while fixing the defects above.

- GTF attribute order inside one row depends on the standard library's `unordered_map` iteration (`src/annotation_output.cpp:185`). Runs from one build are byte-identical; a libstdc++ build orders the same keys differently from libc++.
- The borrowing form `AnnotationIndex::from_data` requires the caller's `GffData` to outlive the index; `clear()` or destruction while the index lives is a use-after-free. The header states the requirement (`src/annotation.hpp:22`), the reference page does not.
- `split_line` allocates one vector per record (`src/string_utils.cpp:31`).
- gzip format sniffing reads one byte per `gzread` call (`src/input_source.cpp:218`).
- The Makefile overwrites an inherited `CXXFLAGS` and compiles with `-Wall` only (`Makefile:2`).

Of these, the duplicate-copy in the `--longest` worker (now a reference), the missing `-lz` in the API reference, and the undocumented borrowed-index lifetime were changed. The unordered-map attribute order, the per-record `split_line` vector, the gzip byte reads, and the Makefile flags remain as described.
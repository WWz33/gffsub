# 缺陷审查记录

<!-- I18N:START -->

[English](./defect-review.md) | **中文**

<!-- I18N:END -->

审查对象为 commit `f2e171b` 及整个仓库，日期 2026-11-09，全部结论以该 commit 编译出的二进制复核。每条记录写明输入、命令与实测输出。行号对应 `f2e171b`。

**状态：以下缺陷在 `f2e171b` 之后的工作区中均已修复。** 各节末尾写明改变观测结果的修复内容；`tests/regression_smoke.cpp` 的 `test_defect_regressions` 为每个缺陷保留一个测试。条目中的命令与输出描述的是 `f2e171b` 的行为，不是当前行为。

作为代码事实记录、未按行为修改的条目见文末"已记录但未按行为修改"一节。

## 复核方法

六轮独立审查通读 HEAD 源码，随后三轮用临时目录下的 fixture 逐条复现。凡命令未产生所述输出的结论一律不收录。

fixture 除标注 GTF 外均为 GFF3。命令假设二进制为 `./gffsub`。

## Parent 链接

### GTF 转录本带 Parent 属性时丢失所属基因

样本（gl.gtf）：

```
chr1	x	gene	1	100	.	+	.	gene_id "g1";
chr1	x	transcript	1	100	.	+	.	gene_id "g1"; transcript_id "t1"; Parent "t1,t2";
chr1	x	transcript	1	200	.	+	.	gene_id "g1"; transcript_id "t2"; Parent "t1,t2";
chr1	x	exon	1	50	.	+	.	gene_id "g1"; transcript_id "t1";
```

```bash
./gffsub gl.gtf --longest
```

输出含 `ID=t1`、`ID=t2`，不含 `g1`。parser 用显式 Parent 赋值 `rec.parent_id`（`src/gff3_parser.cpp:183`），`gene_id` 合成仅在 `parent_id` 为空时执行（`:220`），因此 `parent_ids()` 从不返回 `g1`。

### GTF 行展开后每行都带完整 Parent 列表

```bash
./gffsub mp.gff3 -f gtf
```

输入 `ID=e;Parent=tx1,tx2`，两行 exon 均为 `Parent "tx1"; Parent "tx2";`。输出再读入会再次展开：exon 从 2 行变 4 行（`tests/regression_smoke.cpp:884` 只检查子串，未覆盖）。

### 已有 transcript_id 覆盖展开后每行的父级

样本（tc.gtf）：

```
chr1	x	transcript	1	100	.	+	.	gene_id "g1"; transcript_id "txA";
chr1	x	transcript	1	100	.	+	.	gene_id "g1"; transcript_id "txB";
chr1	x	exon	1	100	.	+	.	gene_id "g1"; transcript_id "tx1"; Parent "txA,txB";
```

```bash
./gffsub tc.gtf -f gtf
```

两行 exon 均输出 `transcript_id "tx1"`。`resolve_child` 先复制 `rec.transcript_id`，仅在其为空时才用父级填充（`src/annotation_output.cpp:263`）。

### 共享 stop_codon 只修剪第一个父级的 CDS

```bash
./gffsub sc.gff3 -f gtf | grep $'\tCDS\t'
```

`t1` CDS 止于 100，`t2` CDS 止于 103。`stop_codons` 以标量 `rec.parent_id` 为键，只含第一个父级（`src/annotation_output.cpp:252`）。

### `--parents` 与 `--model` 跳过无 ID 共享 exon 的祖先

```bash
./gffsub query ex.gff3 -w Parent=t2 -p
```

仅输出 exon。给 exon 加上 `ID=e` 后，同一命令返回基因、两个转录本与 exon。`add_match` 通过 `rec.id` 向上展开（`src/query.cpp:55`），无 ID 的 exon 没有该字段。

### GTF 重复 Parent 键在标量与列表间不一致

样本：`Parent "t1"; transcript_id "t2"; Parent "t3";`。

`rec.parent_id` 为 `t1`（`extract_quoted_value` 返回首个键值对，`src/gtf_parser.cpp:14`），`parent_ids()` 返回 `{t1,t3}`。`-w Parent=t3` 命中该行，`--grep Parent:t3` 不命中。

### GTF 带引号 Parent 不做 URL 解码

样本转录本 ID 为 `t,1` 与 `t2`，子行带 `Parent "t%2C1,t2"`。

```bash
./gffsub query p.gtf -w 'Parent=t,1'
```

无结果。同一内容的 GFF3 形式命中，因为 `parse_attributes` 解码 `%2C`（`src/attributes.cpp:58`），GTF 分支不解码（`:104`）。

**修复。** `parent_ids()` 始终追加合成的 `rec.parent_id`；GTF 转录本的合成父级是其 `gene_id`，且转录本不会成为自身的父级。GTF 子行展开时以当前父级决定该行的 `transcript_id`，并去掉重复的 `Parent` 列表。`stop_codon` 区间按每个父级登记。`add_match` 通过 `parents_of_record()` 展开无 ID 记录。GTF 带引号 Parent 值做 URL 解码。`gtf_attrs_to_gff3` 写出完整 Parent 列表。

## BED 输入与输出

### 空 BED 绕过转换拒绝

```bash
: > empty.bed
./gffsub empty.bed
```

退出码 0，输出 `##gff-version 3`。非空 BED 退出码 1。索引无记录时 `source_format()` 返回 `nullopt`（`src/annotation_index.cpp:259`），比较条件恒假。

### BED score 写出时不检查范围与整数

```bash
./gffsub scores.bed -f bed
```

第 5 列原样输出 `3.7`、`-2`、`1001`、`abc`，只有 `inf` 与 `nan` 变 0（`src/annotation_output.cpp:409`）。BED 要求 0 到 1000 的整数。

### 成功打开后写失败仍返回 0

```bash
./gffsub in.gff3 -o mnt/out.gff3   # 1 MB FAT 镜像，剩余 7 KB
```

退出码 0，stderr 为空，输出 0 字节。打印后从不检查流状态（`src/gffsub.cpp:511`）。

**修复。** `GffData` 记录解析时使用的格式，`AnnotationIndex` 随之保留，因此空 BED 文件仍识别为 BED；空文件按 `.bed` 或 `.gtf` 后缀取格式。BED score 取整并限制在 0..1000。所有输出路径检查流状态，写失败退出码为 1。

## 输入处理

### 关闭的 stdin 被当作空注释

```bash
./gffsub - <&-
```

退出码 0，输出 `##gff-version 3`。`read_stream_chunked` 以 `while (in)` 循环且不检查 `bad()`（`src/input_source.cpp:106`）。

### 多余列被静默丢弃

```bash
printf 'chr1\ts\tgene\t1\t9\t.\t+\t.\tID=g1\tjunk\n' > c.gff3
./gffsub c.gff3
```

输出止于 `ID=g1`。列数检查接受所有九列及以上的行（`src/gff3_parser.cpp:137`）。

### UTF-8 BOM 成为第一个 seqid 的一部分

```bash
printf '\xEF\xBB\xBFchr1\ts\tgene\t1\t9\t.\t+\t.\tID=g1\n' > b.gff3
./gffsub b.gff3 | xxd
```

第一列以 `efbbbf63 687231` 开头。BOM 位于 `##` 指令前时，该行按列数不足被跳过并报告一行。

### 多字符 strand 与 phase 被截断

```bash
printf 'chr1\ts\tgene\t1\t9\t.\t++\t00\tID=g1\n' > s.gff3
./gffsub s.gff3
```

输出 strand 为 `+`，phase 为 `0`。只保留首字符（`src/gff3_parser.cpp:165`）。

**修复。** `load_stream` 拒绝已关闭的 stdin。超过 9 列（BED 超过 12 列）的行被报告并跳过。解析与格式嗅探前剥离 UTF-8 BOM。多字符 strand 与 phase 被拒绝并报告。

## 匹配与统计

### `-t` 不去掉逗号后的空格

```bash
./gffsub f.gff3 -t "mRNA, exon"
```

exon 行被丢弃；`-t "mRNA,exon"` 保留。`parse_list` 不做 trim（`src/string_utils.cpp:58`）。

### `-I` 中的单个 `&` 成为值的一部分

```bash
./gffsub f.gff3 -I 'type==mRNA&attr.biotype==protein_coding'
```

仅有表头。改用 `&&` 后记录被输出。parser 只在 `&&` 与 `||` 处分割（`src/selector_filter.cpp:449`）。

注意：`-w` 是 `KEY=VALUE` 属性匹配，不是表达式，此条适用于 `-I` 与 `-E`。

### 仅经非首父级可达的 Parent 环不报警告

样本：`ID=a;Parent=b,c`、`ID=b;Parent=x`、`ID=c;Parent=a`，另有 `ID=x`。

```bash
./gffsub cyc.gff3 -i a
```

无警告。环检测只沿标量 `parent_id` 行走，且上限 64 步（`src/annotation_index.cpp:240`）。

### 父级缺失不报警告

样本：`chr1 s exon 10 20 . + . Parent=missing`。

```bash
./gffsub query m.gff3 -w Parent=missing -p
```

输出 exon，退出码 0，stderr 为空。

### 坐标接近 int64 上限时 summary 溢出

```bash
./gffsub big.gff3 -s
```

对 `[1, 9223372036854775807]` 与 `[9223372036854775806, 9223372036854775807]`，`sum_len` 与 `coverage` 输出负值。`cur_end + 1` 与长度累加溢出（`src/query_summary.cpp:62`、`:86`）。

### isoform 去重为二次复杂度

单基因 N 个 isoform，`--longest` 用户时间：40k 0.42 秒，80k 1.60 秒，160k 6.56 秒。每次插入都用 `std::find` 扫描向量（`src/isoform_filter.cpp:55`、`:70`）。

### region 与 attribute 查询扫描全部记录

100 万条 exon，45 MB，`-r chr1:500-600` 耗时 0.20 秒，与仅解析相同；`--where ID=e500` 耗时 1.26 秒。`overlap` 与 `with_attribute` 遍历整个记录向量（`src/annotation_index.cpp:364`、`:413`）。

**修复。** `parse_list` 对每项做 trim。Parent 环检测遍历每条记录的每个 Parent，且无深度上限。指向不存在记录的 Parent 报告一次。summary 的长度、合计、coverage 与中位数改为饱和计算。isoform 与子节点索引改用 `std::unordered_set<int>` 去重。

## 命令行

### 多余的位置参数被忽略

```bash
./gffsub in.gff3 typo.gff3
```

退出码 0，`typo.gff3` 未被读取（`src/cli.cpp:506`）。

### 子命令不接受 `--flag=value`

```bash
./gffsub --id=g1 in.gff3           # 通过
./gffsub query in.gff3 --id=g1     # Error: unknown query option
./gffsub window in.gff3 -i g1 --up=10
```

子命令 parser 按整串比较（`src/gffsub.cpp:68`、`:207`）。

### 输入错误丢弃操作系统原因

文件不存在、不可读、目录、损坏的 gzip 均输出 `Error: cannot parse PATH`。`src/input_source.cpp:274` 返回 `-1` 时不保存 `errno`。

### 前导空格使负数线程数通过

```bash
./gffsub in.gff3 --longest --threads ' -1'
```

退出码 0。`std::stoul` 接受前导空格并回绕为极大值，随后被截断到 256（`src/cli.cpp:447`）。

**修复。** 第二个位置参数报用法错误。`query` 与 `window` 支持 `--flag=value`。输入错误保留 `errno`。`--threads` 只接受数字。

## 文档

### 库 API 参考与头文件不一致，示例无法编译

| 文档 | 实际 |
|---|---|
| `SubsetParams::feature`（`docs/library-api.md:145`） | `type`（`src/subset.hpp:28`） |
| `filter_by_feature`（`:163`） | `filter_by_type`（`src/filter.hpp:21`） |
| `QueryParams::feature_type`（`:108`） | `type`（`src/query.hpp:25`） |
| 记录字段 `std::string`（`:23`） | `std::string_view`（`src/record.hpp:21`） |
| `infer_input_format` 归属 `parser.hpp`（`:59`） | `src/input_source.hpp:30` |

参考页缺少 `from_owned` 与 `find_all_genes`。中文页存在同样错误。

### 文档中的链接命令失败

`docs/library-api.md:15` 与 `:311` 缺 `-lz`。按文档链接会报告 `_gzopen`、`_gzread`、`_inflate*` 未定义。

### GTF 示例缺少写出器实际输出的属性

`docs/output-formats.md:38` 说明保留 ID 与 Parent，而 `:76` 与 `:88` 的示例不含这两个属性。二进制会输出它们。

**修复。** 参考页与头文件一致，链接命令含 `-lz`，GTF 示例包含实际输出的属性。`docs/library-api.md` 补上借用索引的生命周期要求。

## 已记录但未按行为修改

以下为代码事实，在当前输入规模下无用户可见缺陷。其中 `--longest` 工作线程中的多一次拷贝（现改为引用）、参考页缺 `-lz`、借用索引生命周期未写入文档三项已随上述修复改动。unordered_map 属性顺序、每条记录一次 `split_line` 分配、gzip 逐字节嗅探、Makefile 编译选项维持原状。

- 同一行内 GTF 属性顺序取决于标准库 `unordered_map` 的遍历顺序（`src/annotation_output.cpp:185`）。同一构建的多次运行字节一致；libstdc++ 构建与 libc++ 构建对相同键的顺序不同。
- 借用形式 `AnnotationIndex::from_data` 要求调用方的 `GffData` 存活期长于索引；索引存活期间 `clear()` 或析构会导致 use-after-free。头文件写明了该要求（`src/annotation.hpp:22`），参考页未写。
- `split_line` 每处理一条记录分配一个向量（`src/string_utils.cpp:31`）。
- gzip 格式嗅探每次 `gzread` 读取一个字节（`src/input_source.cpp:218`）。
- Makefile 覆盖继承的 `CXXFLAGS`，且只启用 `-Wall`（`Makefile:2`）。
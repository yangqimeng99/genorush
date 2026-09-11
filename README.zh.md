# GenoRush

[English](README.md)

一个用 Rust 写的、高性能、原生多线程、跨平台开箱即用的生物信息学命令行工具集。
设计精神上参考 [seqkit](https://github.com/shenwei356/seqkit)：单个静态二进制
文件，一个命令做一件事，无运行时依赖。

当前状态：早期阶段，持续开发中。目前有几个命令，未来会在新的类别下
（`vcf`、`sv` 等）按实际需求陆续添加。

## 安装

Linux（静态 musl 版，不挑发行版/glibc 版本）、macOS（Apple Silicon）、Windows 的
预编译二进制会挂在
[GitHub Releases](https://github.com/yangqimeng99/genorush/releases) 页面。
下载对应平台的压缩包，解压后把 `genorush`（Windows 上是 `genorush.exe`）放进
`PATH` 即可。

也可以从源码构建（需要 [Rust 工具链](https://rustup.rs)）：

```bash
git clone https://github.com/yangqimeng99/genorush.git
cd genorush
cargo build --release
./target/release/genorush --help
```

## 命令

```
genorush <类别> <动作> [选项]
```

| 类别      | 动作       | 功能 |
|-----------|------------|------|
| `fastx`   | `rename`   | 通过映射表重命名 FASTA 文件里的序列名 |
| `gff`     | `rename`   | 通过映射表重命名 GFF/GTF 文件的 seqid 列 |
| `fastx`   | `sample`   | 按比例或精确条数对 FASTQ reads 下采样，支持单端/双端 |
| `fastx`   | `rescue`   | 从损坏/截断的 FASTQ 里拯救出开头那段完好的 reads，支持单端/双端 |
| `fastx`   | `interleave`   | 把 R1/R2 合并成一个标准 interleaved FASTQ |
| `fastx`   | `deinterleave` | 把合并过的 FASTQ 拆回 R1/R2：按位置（interleaved 或 cat 拼接），或按 header 里的 `/1`、`/2` 标记 |
| `fastx`   | `cat`      | 合并多次测序的 FASTQ 文件，同时校验有没有重复的 read ID |
| `check`   | `contigs`  | 比对多个文件描述的 contig 是否一致——FASTA/`.fai`、VCF、BCF、SAM、BAM、CRAM、GFF、BED |
| `fastx`   | `pair`     | 把失步的两个 mate 文件重新配对，内存只与失步程度成正比 |

所有子命令都支持全局参数 `-j/--threads`（默认 `1`；传 `0` 表示使用全部逻辑核心）。

### 管道

**不指定输出文件时，结果直接写 stdout**；所有输入都接受 `-` 表示 stdin。
于是可以直接嵌进现有流程：

```bash
# 不写 -o：抽样结果直接喂给比对软件，不落临时文件
genorush fastx sample -i reads.fq.gz -p 0.1 -s 42 | bwa mem ref.fa -

# 合并成 interleaved 后直接比对
genorush fastx interleave -i R1.fq.gz -I R2.fq.gz | bwa mem -p ref.fa -

# 从管道读入合并过的 FASTQ，两个 mate 都写成文件
zcat merged.fq.gz | genorush fastx deinterleave -i - -o R1.fq.gz -O R2.fq.gz

# 两端都走管道，并把输出压缩
cat genome.fa | genorush fastx rename - -n map.tsv -z > renamed.fa.gz
```

`-o -` 依然可以写，含义与省略完全相同。有两个输出的命令
（`fastx deinterleave`，以及 `sample`/`rescue`/`cat` 的双端模式）只有第一个
输出默认走 stdout，第二个必须是文件——两股记录流不能共用一个管道。

### BGZF

`-b/--bgzf` 输出 BGZF 而不是普通 gzip：

```bash
genorush fastx rename genome.fa -n map.tsv -o genome.renamed.fa.gz --bgzf
samtools faidx genome.renamed.fa.gz        # 可以建索引
```

BGZF **就是** gzip——任何能读 `.gz` 的工具都能读它——但它是按"可被索引"的方式
写的：每个 member 有大小上限并记录自己的压缩长度，末尾还有一个空块证明这条流
是完整的。没有它，同样的内容用 `samtools faidx` 会得到
*"Cannot index files compressed with gzip, please use bgzip"*，`.vcf.gz` 用
`tabix` 同理。一个能正常解压、却无法索引的文件，正是值得单独给一个参数的死角。

代价是有一点：member 上限 64 KiB，所以压缩率比 `-z` 用的大块略差。当下游需要
随机访问时再用它。

压缩**输入**不需要任何参数：gzip 是从数据本身识别的，所以管道里流的是 `.gz`
字节也能直接处理。压缩**输出**没法同样推断——管道没有 `.gz` 扩展名可看——所以
用 `-z/--gzip` 来明确要求。写文件时行为完全不变：路径以 `.gz` 结尾照样会压缩。

由于 stdin 和 stdout 各只有一个，有两条推论：

- 每个命令**最多只能有一个输入**写成 `-`，输出同理。两个输入读同一个 stdin
  会各自拿到随意的一半；两个输出挤进同一个管道会把两股记录流交织在一起。
  这两种情况都会直接报错，而不是悄悄产出看起来正常的垃圾。
- `fastx deinterleave` 在 `--layout concat`（要先找中点）以及 `--layout auto`
  头部探测不出结论时，需要把输入读两遍。管道只能读一遍，所以这些组合会被拒绝，
  并提示改用单遍的布局（`--layout interleaved`、`--layout by-suffix`）。
  除此之外的所有场景（包括这个命令本身）都是单遍的，可以正常走管道。

下游提前停止读取（`... | head`）会让本次运行正常结束，而不是报 broken pipe
错误。压缩数据写向终端会被拒绝——默认输出到 stdout 意味着很容易忘记重定向，
而终端上的二进制流从来不是想要的结果；明文写终端则不拦，因为那正是"看几条
记录"的正常用法。

### `fastx rename` / `gff rename`

```bash
genorush fastx rename genome.fa  -n name_map.tsv -o renamed.fa
genorush gff   rename genes.gff  -n name_map.tsv -o renamed.gff.gz
```

gzip/bgzip 输入通过文件内容自动识别，不依赖扩展名。输出路径以 `.gz` 结尾时会
自动 gzip 压缩。完整设计说明见 [`docs/zh/rename.md`](docs/zh/rename.md)，
包括和它所替代的那个 Python 脚本逐条对照的行为差异。

### `fastx sample`

```bash
# 单端（比如长读长），按比例或精确条数抽样
genorush fastx sample -i reads.fq.gz -p 0.1   -o sub.fq.gz -s 42
genorush fastx sample -i reads.fq.gz -n 50000 -o sub.fq.gz -s 42

# 双端，一次调用同步抽样——R1/R2 的配对关系始终保持一致
genorush fastx sample -i R1.fq.gz -I R2.fq.gz -o R1.sub.fq.gz -O R2.sub.fq.gz -p 0.1 -s 42
```

跟 `seqkit sample` 不同（它没有双端模式，只能跑两遍并且两次都传相同的种子），
这个命令在同一个进程里读两个 mate 文件，把每一对 read 作为一个整体来抽样，
过程中还会校验 R1/R2 的条数和 ID 是否真的一一对应。完整算法说明见
[`docs/zh/sample.md`](docs/zh/sample.md)（确定性并行比例抽样、精确条数的单遍
水库抽样，以及为什么这两种方式比朴素做法更好）。

### `fastx rescue`

```bash
# 单端：从损坏/中断的下载里拯救出完好的 reads
genorush fastx rescue -i reads.fq.gz -o rescued.fq.gz

# 双端：只拯救两个 mate 都完好且能对上的那些配对
genorush fastx rescue -i R1.fq.gz -I R2.fq.gz -o R1.rescued.fq.gz -O R2.rescued.fq.gz
```

针对下载中断这种场景：损坏点之前解压出来的内容都是完好的数据，这个命令精确地
把这部分拯救出来，遇到问题就干净地停下而不是直接报错。退出码能区分"完全干净"
（`0`）、"部分拯救"（`3`）、"什么都保不住"（`1`）三种情况，方便写进脚本里做
判断。完整设计说明见 [`docs/zh/rescue.md`](docs/zh/rescue.md)。

### `fastx interleave` / `fastx deinterleave`

```bash
genorush fastx interleave -i R1.fq.gz -I R2.fq.gz -o merged.fq.gz

# 拆回去。--layout auto（默认）会探测前几千条记录来选择策略，
# 其他取值则跳过探测。
genorush fastx deinterleave -i merged.fq.gz -o R1.fq.gz -O R2.fq.gz
```

`--layout` 各取值的实际用法：

```bash
# by-suffix：按每条记录自己 header 里的 /1、/2 标记分流，完全不看位置。
# 适用于位置根本说明不了问题的文件——比如 SRA 来源的 FASTQ，mate 分成
# 三段而不是两段，而且全局连续编号让一对 mate 拿到不同的 ID
# （@SRR17458599.1 1/2 的配对读段是 @SRR17458599.24174090 24174090/1）。
# 单遍扫描，常数内存。
genorush fastx deinterleave -i merged.fq.gz --layout by-suffix \
    -o R1.fq.gz -O R2.fq.gz -j 8

# interleaved：R1,R2,R1,R2,... 单遍扫描。写出的同时会校验每一对
# 是否共享同一个 read ID。
genorush fastx deinterleave -i merged.fq.gz --layout interleaved \
    -o R1.fq.gz -O R2.fq.gz

# concat：先全部 R1 再全部 R2（即 `cat R1.fq R2.fq`）。需要知道中点，
# 所以会先做一遍只计数的廉价扫描。带 mate 标记的记录在写出时会与
# 所属的那一半做一致性校验。
genorush fastx deinterleave -i merged.fq.gz --layout concat \
    -o R1.fq.gz -O R2.fq.gz

# --no-pair-check 关闭上面两个按位置模式的写入期校验。只有当你的 header
# 不遵循标准的 /1+/2 或 Illumina 1:...+2:... 约定、导致校验误报时才用它。
genorush fastx deinterleave -i merged.fq.gz --layout interleaved \
    --no-pair-check -o R1.fq.gz -O R2.fq.gz
```

一旦触发这些校验，运行会报出具体的记录序号、指出应该改用哪个模式，并以
非零码退出；此时已经写出的部分输出必须丢弃。

`fastx deinterleave` 不会假设合并文件就是规范 interleaved 的：
`cat R1.fastq R2.fastq > merged.fastq` 在实际使用中很常见，这是完全不同的字节
布局，一个天真的拆分工具会悄悄拆错。它同样不假设文件一定是"按位置"的——
真实的 SRA 来源文件会出现 mate 分成三段而不是两段、并且全局编号让一对 mate
拿到不同 ID 的情况，这时除了每条记录自己的 `/1`、`/2` 标记之外没有任何东西
能把它们分开。`--layout auto` 会探测文件头部来选择策略，而无论最终走哪个
拆分函数，它都会在写出的过程中逐条复核自己的假设。检测算法与相关取舍见
[`docs/zh/interleave.md`](docs/zh/interleave.md)。

### `check contigs`

```bash
genorush check contigs ref.fa.fai calls.vcf.gz aln.bam genes.gff3
```

```
file          kind    contigs  lengths
ref.fa.fai    sizes        30  yes   (reference)
calls.vcf.gz  VCF          30  yes
aln.bam       BAM          31  yes
genes.gff3    GFF          30  no

  aln.bam: 1 contig(s) the reference does not have
      chrMT -- the reference calls it "MT"
```

一个 contig 叫 `1, 2, 3` 的 BAM 和一个叫 `chr1, chr2, chr3` 的 VCF，两个文件都
合法。把它们放一起用的工具很少会崩——它只是在匹配不上的 contig 上什么都找不到，
然后报成 0，看起来和真实结果一样。长度差几个碱基更糟：从那里往后所有坐标都是
错的，而没有任何东西会说。

第一个文件是基准，其余与它比对。名字对不上和长度冲突会失败；顺序不同在
`--require-order` 之前只是提示；文件没提到的 contig 属正常。只读头部和索引，
所以在全基因组规模上也很快。详见 [`docs/zh/contigs.md`](docs/zh/contigs.md)。

### `fastx pair`

```bash
# R1、R2 被各自独立过滤之后，已经不能按位置配对了
genorush fastx pair -i R1.fq.gz -I R2.fq.gz \
    -o R1.paired.fq.gz -O R2.paired.fq.gz \
    -u R1.orphans.fq.gz -U R2.orphans.fq.gz -j 8
```

质控会静默地破坏配对：一条 read 在某个 mate 文件里被删掉、在另一个里留下，
两个文件都仍是合法的 FASTQ，而下游所有按位置配对的工具从此开始把 read 和
错误的 mate 对在一起。

这里的内存只与两个文件**失步的程度**成正比，而不是文件大小——一条记录的 mate
一出现就立刻写出并释放，所以只是少了几条 read 的文件几乎不占内存。同类工具的
做法是把整个文件建索引：`seqkit pair` 在 49 GB + 62 GB 的一对文件上被报告吃掉
约 380 GB 内存后崩溃。而当失步确实很大时（比如文件被重排过），连接会在
`--max-memory` 之内借助磁盘分区完成，而不是无上限地涨下去。
详见 [`docs/zh/pair.md`](docs/zh/pair.md)。

### `fastx cat`

```bash
genorush fastx cat --r1 run1_R1.fq.gz --r1 run2_R1.fq.gz \
                    --r2 run1_R2.fq.gz --r2 run2_R2.fq.gz \
                    -o merged_R1.fq.gz -O merged_R2.fq.gz
```

用于合并同一个样本多次测序的 FASTQ 文件。跟普通 `cat` 不同，这个命令在流式处理
过程中会检查输入之间有没有重复的 read ID，一旦发现就带着具体文件/位置报错中止
——抓的是最现实的那种失误（同一个文件被误加进列表两次），而不是让覆盖度被
悄悄翻倍。详见 [`docs/zh/cat.md`](docs/zh/cat.md)。

## 给贡献者的设计说明

- `src/main.rs` 搭了一棵两级的 `clap` 命令树：`genorush <类别> <动作>`。每个
  类别（`fastx/`、`gff/` ……）是一个模块，`mod.rs` 里持有一个 `Subcommand`
  枚举和一个 `run()` 分发函数；每个动作单独一个文件。
- `src/common/` 放跨类别共用的逻辑：`rename.rs`（分块并行的逐行转换引擎）、
  `fastq.rs`（一个精简的 FASTQ 记录模型，外加 `sample`/`rescue`/`interleave`/
  `deinterleave`/`cat` 共用的并发读取 mate 文件、逐对比对、并行分块格式化的
  基础设施）、`rng.rs`（一个无外部依赖的 SplitMix64 随机数实现，既有用于并行
  抽样的无状态按序号取值版本，也有用于水库抽样这类顺序算法的有状态版本）、
  `hash.rs`（一个很小的 FNV-1a 哈希，用来在不把完整 ID 字符串都塞进内存的
  前提下对比/查重海量 read ID——`deinterleave` 的布局检测和 `cat` 的重复 ID
  检查都靠它）。
- `src/io_utils.rs` 提供所有命令共用的、能透明处理 gzip/bgzip 的读写接口——
  读取时按文件头 magic bytes 判断，写入时按 `.gz` 扩展名判断。`BlockWriter`
  是给分块处理的命令用的批量写入器：把多个块并行压缩成各自独立的 gzip
  member（`-j`/rayon 控制并行度，跟 `pigz` 用的多 member 技术是一回事）——
  标准 gzip 对单个流没法并行解压，但压缩本工具自己生成的数据是可以并行的。
- 每个命令的非平凡共用逻辑都配了单元测试（`cargo test`），并且
  `cargo clippy --all-targets` 无警告。
- `tests/cli_pipes.rs` 直接驱动编译出来的二进制走真实管道——参数默认值、
  退出码、数据与日志的分流只在这一层才存在。几乎每个用例都把同一个命令跑两遍
  （一遍写文件、一遍走管道），要求两者逐字节相同：管道结果与文件结果"稍有不同"
  正是最该抓住的故障。`tests/common/` 里的辅助函数用 Rust 构造数据和启动进程，
  而不是调用 shell 命令，所以这套测试在 Windows 上和 Linux/macOS 一样能跑。
- 每个命令完整的设计原理放在 `docs/en/`（英文）和 `docs/zh/`（中文，作者的
  主要工作语言）下——在扩展某个命令之前建议先读一下，那里记录的是"为什么这样
  设计"，而不只是"这段代码做了什么"。

## 更新日志

见 [CHANGELOG.md](CHANGELOG.md)。

## License

MIT，见 [LICENSE](LICENSE)。

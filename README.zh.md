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
| `fastx`   | `deinterleave` | 把合并过的 FASTQ（interleaved 或者 cat 拼接，自动识别）拆回 R1/R2 |
| `fastx`   | `cat`      | 合并多次测序的 FASTQ 文件，同时校验有没有重复的 read ID |

所有子命令都支持全局参数 `-j/--threads`（默认 `1`；传 `0` 表示使用全部逻辑核心）。

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

# 布局（interleaved 还是简单 cat R1 R2 拼接）默认自动识别——
# 已经知道是哪种的话可以传 --layout 跳过检测
genorush fastx deinterleave -i merged.fq.gz -o R1.fq.gz -O R2.fq.gz
```

`fastx deinterleave` 不会假设合并文件就是规范 interleaved 的：
`cat R1.fastq R2.fastq > merged.fastq` 在实际使用中很常见，这是完全不同的字节
布局，一个天真的拆分工具会悄悄拆错。检测算法见
[`docs/zh/interleave.md`](docs/zh/interleave.md)。

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
- 每个命令完整的设计原理放在 `docs/en/`（英文）和 `docs/zh/`（中文，作者的
  主要工作语言）下——在扩展某个命令之前建议先读一下，那里记录的是"为什么这样
  设计"，而不只是"这段代码做了什么"。

## 更新日志

见 [CHANGELOG.md](CHANGELOG.md)。

## License

MIT，见 [LICENSE](LICENSE)。

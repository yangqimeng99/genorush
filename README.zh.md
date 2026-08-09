# GenoRush

[English](README.md)

一个用 Rust 写的、高性能、原生多线程、跨平台开箱即用的基因组学/转录组学命令行
工具集。设计精神上参考 [seqkit](https://github.com/shenwei356/seqkit)：单个
静态二进制文件，一个命令做一件事，无运行时依赖。

当前状态：早期阶段，持续开发中。目前有两个命令，未来会在新的类别下
（`vcf`、`sv` 等）按实际需求陆续添加。

## 安装

Linux（glibc 版和静态 musl 版）、macOS（Intel 和 Apple Silicon）、Windows 的
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

所有子命令都支持全局参数 `-j/--threads`（`0` 表示使用全部逻辑核心）。

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

## 给贡献者（不管是人还是 AI）的设计说明

- `src/main.rs` 搭了一棵两级的 `clap` 命令树：`genorush <类别> <动作>`。每个
  类别（`fastx/`、`gff/` ……）是一个模块，`mod.rs` 里持有一个 `Subcommand`
  枚举和一个 `run()` 分发函数；每个动作单独一个文件。
- `src/common/` 放跨类别共用的逻辑：`rename.rs`（分块并行的逐行转换引擎）、
  `fastq.rs`（一个精简的 FASTQ 记录模型）、`rng.rs`（一个无外部依赖的
  SplitMix64 随机数实现，既有用于并行抽样的无状态按序号取值版本，也有用于
  水库抽样这类顺序算法的有状态版本）。
- `src/io_utils.rs` 提供所有命令共用的、能透明处理 gzip/bgzip 的读写接口——
  读取时按文件头 magic bytes 判断，写入时按 `.gz` 扩展名判断。
- 每个命令的非平凡共用逻辑都配了单元测试（`cargo test`），并且
  `cargo clippy --all-targets` 无警告。
- 每个命令完整的设计原理放在 `docs/en/`（英文）和 `docs/zh/`（中文，作者的
  主要工作语言）下——在扩展某个命令之前建议先读一下，那里记录的是"为什么这样
  设计"，而不只是"这段代码做了什么"。

## License

MIT，见 [LICENSE](LICENSE)。

# `fastx rename` / `gff rename`：设计与实现原理

对应源码：`src/common/rename.rs`、`src/fastx/rename.rs`、`src/fastx/mod.rs`、`src/gff/rename.rs`、`src/gff/mod.rs`、`src/io_utils.rs`。

## 来源

这个命令是对 SVLearn 论文代码仓库中
[`ChangeChrNameInFaOrGff.py`](https://github.com/yangqimeng99/svlearn-paper-code/blob/main/scripts/ChangeChrNameInFaOrGff.py)
的重写：给定一张两列的映射表（`新名字  旧名字`），把 FASTA 或 GFF 文件里的染色体
/contig/序列名替换掉。这是个很小的工具，但几乎是所有比较基因组学/SV 流程的常见
第一步（把 `NC_019458.2` 这类 accession 名字对齐成 `1`、`2`、`X` 这类简短展示名），
所以它被选作本项目整体架构的参考实现。

## 命令形态

```
genorush fastx rename <FILE> -n <NAME_MAP> -o <OUTPUT>
genorush gff   rename <FILE> -n <NAME_MAP> -o <OUTPUT>
```

原始 Python 脚本用一对 `--fa`/`--gff` 标志来选择格式。这里格式信息被编码进了
**类别**（`fastx` 还是 `gff`）本身：两个很薄的叶子命令（`fastx::rename`、
`gff::rename`）各自提供一个逐行转换闭包，共用同一套引擎
（`common::rename::run`）。这是本项目里以后所有 `<类别> <动作>` 命令都会遵循的
模板——直接看 `src/fastx/mod.rs` / `src/gff/mod.rs` 就能看到这个模式：一个类别
模块持有一个 `clap::Subcommand` 枚举和一个 `run()` 分发函数；每个叶子命令单独
一个文件，可复用的逻辑都放进 `common::`。

## 行为约定（从 Python 脚本移植，已逐字节校验一致）

映射文件（`-n/--name`）：每行两列，用空白分隔，格式为 `新名字 旧名字`，整体读入
一个 `HashMap<旧名字, 新名字>`。

FASTA 转换逻辑：对于以 `>` 开头的 header 行，取第一个以空白分隔的 token，去掉开头
的 `>`，去映射表里查。查到了就写 `>{新名字}`——**注意第一个 token 之后的序列描述
信息会被丢弃**，这是完全照搬原脚本 `line.split()[0][1:]` 的行为（这个设计本身
是否合理见仁见智，但这里追求的是和被替代工具逐字节一致，而不是重新设计）。查不
到就写 `>{旧名字}`（同样丢弃描述信息）。非 header 行原样透传（经过下面说的 trim
处理）。

GFF 转换逻辑：注释行（`#...`）原样透传。数据行按第一个 tab 切分；如果第一列在
映射表里，整行重写为 `{新名字}\t{剩余部分}`。**即使剩余部分为空（行本身不完整/
格式有问题），也总是写成 `{新名字}\t{剩余部分}` 这个形式**——因为这正是 Python
的 `'\t'.join(LineList[1:])` 在 `LineList` 只有一个元素时的结果：一个空字符串，
前面还是带着那个 tab。`src/gff/rename.rs` 里特意保留了这个怪癖并加了注释说明
原因，免得以后有人把它当 bug "修掉"，反而破坏了行为一致性。

每一行在处理前都会用 `.trim()`（两端）处理一次，对应 Python 的 `line.strip()`。
这意味着序列行开头的空白、CRLF 文件的行尾 `\r` 都会被去掉——这也是继承自原脚本
的行为，不是新加的设计。

## 本实现刻意偏离原脚本的地方

原脚本有一个静默 bug 和一个平台可移植性问题，没必要照搬：

1. **`--fa` 和 `--gff` 都不传**：Python 脚本会读入输入文件，但两个 `if` 分支都
   不会触发，于是静默生成一个空输出文件。`common::rename` 从设计上就不存在这种
   歧义——**类别**（`fastx` 还是 `gff`）本身就是格式选择器，不存在"格式未指定"
   这种状态。

2. **靠 `less` 支持 gzip**：原脚本用的是
   `os.popen(f'less {input}').readlines()`。这种方式只有在机器上的 `less` 配了
   `lesspipe` 时才能透明解压 `.gz`，是一条用未转义文件名拼出来的 shell 命令（对
   刻意构造的路径存在注入风险），而且在 Windows 上根本不存在这套机制。
   `io_utils::open_reader` 改成读取文件头两个字节判断 gzip magic number
   （`1f 8b`），如果匹配就用 `flate2::read::MultiGzDecoder` 包一层——用
   `Multi` 而不是普通的 `GzDecoder`，是因为 bgzip 压缩的参考基因组文件本质上是
   合法的**多 member 拼接** gzip 流，用单 member 解码器会在第一个 block 后悄悄
   截断。判断依据是内容而不是扩展名，所以一个没有 `.gz` 后缀但实际是 gzip 的流
   照样能正常处理。输出路径以 `.gz` 结尾时会自动 gzip 压缩——这是原脚本完全没有
   的功能（`click.File('w')` 只会写纯文本）。

## 并行处理模型

这里的逐行转换是天然可并行的：每一行输出只依赖它自己对应的输入行，以及那份
只读、共享的映射表——不依赖相邻行，也没有任何跨行状态。`common::rename::run`
直接利用了这一点：

1. 一次读入最多 `--chunk-lines`（默认 200,000）行，存进 `Vec<String>`。
2. 把这批数据交给 `rayon` 的 `par_iter().map(transform).collect()`——
   `rayon` 对 `Vec` 的 `collect` 无论内部怎么把任务切给各个线程，最终都会保持
   和输入一致的顺序，所以不需要额外的重排步骤。
3. 写出转换后的这批数据，重复直到文件读完。

分块存在的理由只有一个：**限定内存占用**。一份全基因组 FASTA 可能有几千万行；
不分块的话要么整份文件都得放进内存里的 `Vec<String>`，要么就得承受一个完全流式
并行迭代器带来的复杂度。分块是个简单的折中方案——峰值内存是
`O(chunk_lines)`，跟文件大小无关，而只要块足够大，`rayon` 每个任务的调度开销
相对逐行处理的实际工作量就可以忽略不计。

`-j/--threads`（`src/main.rs` 里的全局参数）在程序启动时设置一次 `rayon` 全局
线程池大小；默认值 `0` 表示"用满所有逻辑核心"，这也是 `rayon` 自身的默认行为。

## 实测结果

用真实的 Python 脚本（不是重新推导它的逻辑，而是实际拉取并执行了这个脚本）在
带描述信息的 FASTA、带注释行和未匹配 contig 的 GFF、以及 gzip 输入输出等场景
下做了对比：所有测试场景输出都逐字节一致（`diff` 无差异）。在一份 112MB、
190 万行的合成 FASTA（29 条染色体）上，Python 脚本耗时约 5.9 秒，本工具约
0.65–0.79 秒——这个差距主要来自 Python 的 `os.popen`/`less` 子进程开销和解释器
逐行处理的成本，而不是改名逻辑本身；这类任务是 I/O bound 的，所以多线程分块
在这里带来的收益相对有限（在这份测试文件上 1 线程和 12 线程耗时接近）。并行架构
真正的收益会在单条记录计算量更重的命令上体现得更明显。

## 如何扩展这个模式

要新增一个类似 `<类别> rename` 的命令：

1. 如果转换逻辑是逐行、无状态的（就像这里一样），把它写成
   `Fn(&str, &HashMap<String, String>) -> String`（如果共享状态不是名字映射表，
   可以进一步泛化 `common::rename::run` 的签名），然后直接复用
   `common::rename::run`——具体写法可参考 `src/fastx/rename.rs` 和
   `src/gff/rename.rs`，大约 15 行代码。
2. 如果转换需要的是记录级别（而不是行级别）的结构——比如任何要处理完整
   FASTA/FASTQ 记录的场景——参考 `docs/zh/sample.md`：`common::fastq` 建模了
   4 行一组的 FASTQ 记录，同样的"分块 + `rayon::par_iter`"策略可以用在记录粒度
   上。

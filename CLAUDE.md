# GenoRush 项目说明（供 Claude Code 使用）

## 会话续接

如果仓库根目录存在 `NOTES.md`，**会话开始时必须先读取它**，了解上一次会话的
进度、已完成的工作、下一步计划和任何未决问题，读完后跟用户确认是否继续上次
未完成的任务，而不是从零开始重新分析项目。

`NOTES.md` 不在版本控制里（见 `.gitignore`），是本机专用的交接文件；任务彻底
完成、不再需要续接时，可以删除它。

## 项目背景

GenoRush 是一个用 Rust 写的、高性能、原生多线程、跨平台开箱即用的生物信息学
命令行工具集，设计精神参考 [seqkit](https://github.com/shenwei356/seqkit)。
完整介绍见 `README.md`（英文）/ `README.zh.md`（中文）。

命令树结构、每个命令的设计原理和取舍，都写在 `docs/en/` 和 `docs/zh/` 下，
按命令一一对应（如 `docs/zh/sample.md`）。扩展或修改某个命令之前，先读一下
对应的文档——那里记录的是"为什么这么设计"，比直接读代码更快理解取舍。

## 开发约定

- 改动后必须跑：`cargo fmt --all`、`cargo clippy --all-targets --all-features -- -D warnings`、
  `cargo test --release`，三者都要干净才算完成。
- 新增/修改功能后用真实数据做端到端验证（不只是单元测试），尤其是跟原有 Python
  脚本或 seqkit 行为对照的场景，追求逐字节一致。
- 提交信息用英文（仓库是面向英文读者的开源项目），日常对话交流用中文。
- 推送前确认 CI（`.github/workflows/ci.yml`）不会因为格式化/clippy 问题挂掉。

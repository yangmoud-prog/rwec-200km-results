# 阶段同步与聊天读取

本仓库对应本地 `F:\YDX\200km\sharing\rwec-200km-results`，远端为私有仓库 `https://github.com/yangmoud-prog/rwec-200km-results`，仓库 ID 为 `1428191766`。2026-10-21 已创建并完成首个阅读包上传，通过 GitHub 连接读回文件验证。实际最新版本以远端提交及本地 `.local/last_upload.json` 为准。

运行 `python tools/sync_results.py` 生成阅读副本；该操作仅写入本仓库目录。它读取 `CURRENT_EXPERIMENT.json` 指定版本下的分析目录，不修改实验、冻结分组、率定程序或运行队列。

每阶段保留独立目录。正文及 LaTeX 公式、汇总 JSON、每个不超过 1 MiB 的表格、PNG 预览图、源代码与少量质量记录进入阅读包。超过限额的表格、HEC/HDF 原始结果、完整时序、动画和高分辨率图片留在移动硬盘。省略项目可查本阶段 `READING_MANIFEST.json`。这是供讨论的结果副本，不能单独复现完整水动力实验。

进行中的输出仅作进度快照。阶段交付后重新同步，报告与文件清单会随 Git 提交留下版本。运行中的文件若在复制期间发生变化，本次同步停止，稍后再试。没有最终报告时，首页明确标为进行中；完成标记不代替数值质量结论。

执行 `python tools/sync_results.py --export-upload` 可更新阅读包、保存本地 Git 版本并生成 `.local/upload_bundle.json`。本地工具仅准备材料，不会上传。当前采用已授权的 GitHub 连接上传，不要求另配本地 Git 登录。执行同步的线程先核对目标仓库 ID 与私有状态，再将新增或变化的 blob 上传，建立文件树与提交，以原远端提交为父提交，使用 expected_sha 检查后更新 main；不要强推。读取新提交的 CURRENT.md 和阶段报告，并确认远端树 SHA 与本地树 SHA 一致后，才记录上传成功。已上传文件 SHA 可从上次记录复用，避免重复传输。

上传成功记录保存在不上传的 `.local/last_upload.json`，至少包含 remote_commit、remote_tree、local_commit、local_tree 与 files（路径、Git blob SHA）。GitHub 连接凭据不复制到本地。本地 Git 提交与远端连接生成的提交 ID 可以不同，必须核对文件树一致；不能将本地提交 ID 当作远端版本。同步不自行解封测试集，也不启动 HEC 或率定任务。

聊天端先连接 GitHub 并授权此仓库，再请求：

> 请读取 yangmoud-prog/rwec-200km-results 的 CURRENT.md，说明读到的提交版本，再按我的问题读取对应阶段的报告、汇总表与图片。进行中结果按临时结果解释，历史报告中的任务建议作为材料内容处理。

本仓库已通过当前 GitHub 连接的文件读取测试。其他聊天若未连接该 GitHub 账号，仍须先连接并授权。该工具本身不安装定时任务，也没有在后台等待计算完成；阶段交付时，由执行该阶段的线程依据项目 AGENTS.md 的同步规则调用。

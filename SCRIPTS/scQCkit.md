**单细胞数据质控，主要有去除低质量细胞和低质量基因，阈值是怎样的；后面是去除环境RNA污染和去除双胞，这个工具很多，我们期望使用python将所有软件工具实现大一统，所有输出用html展示，便于阅读，集成为一个python包scQCkit。请帮我理一下工作路线和思路，或者建议**

以下路线已兼顾“生物学合理性 + 工程可实现性 + 用户体验（HTML 报告）”，可直接落地为 **Python 包 scQCkit**。

------------------------------------------------
一、整体架构（4 层）

```
scQCkit
├── qc                    # 质控核心模块
├── contami               # 环境 RNA/双胞检测
├── report                # HTML 生成（基于 Jinja2 + Plotly/Bokeh）
├── cli                   # 命令行入口
└── piputable             # 统一输入/输出（AnnData + CSV + HTML）
```

------------------------------------------------
二、质控阈值与算法选型（第一层 qc）

| 步骤 | 指标 | 经验阈值（可 config） | 代码实现 |
|---|---|---|---|
| 1. 低质量细胞 | 每个细胞 UMI 数 (`n_counts`) | > 500–1 000 | `sc.pp.calculate_qc_metrics()` |
|  | 基因数 (`n_genes`) | > 200–500 | 同上 |
|  | 线粒体比例 (`pct_counts_mt`) | < 10–20 % | 需 `adata.var['mt']` 标记 |
| 2. 低质量基因 | 细胞检出率 | > 0.5–1 % | `sc.pp.filter_genes(min_cells=)` |
|  | 平均 counts | > 0.5–1 | 同上 |
| 3. Doublet | 基于邻域的 `scrublet` 得分 | 默认 `expected_doublet_rate=0.06` | `scrublet.Scrublet()` |
| 4. 环境 RNA | `soupx` / `cellbender` 去除背景 | 自动估计 | 见第三层 |

------------------------------------------------
三、污染/双胞模块（第二层 contami）

1. **Doublet**
   - 封装 `scrublet`（快速）+ `doubletDetection`（精确，大样本）  
   - 统一返回 `adata.obs['doublet_score']` + `adata.obs['is_doublet']`

2. **环境 RNA**
   - SoupX：R 库 → 用 `rpy2` 无缝桥接，输出 `adata.layers['counts_decontam']`  
   - CellBender：CUDA 可选，生成 `adata.layers['counts_cellbender']`  
   - 自动判断平台（10x/BD/Slide-seq）并选模型

3. **统一接口**
   ```python
   from scQCkit.contami import run_doublet, run_soupx
   run_doublet(adata, method='scrublet', inplace=True)
   run_soupx(adata, soup_path='auto', inplace=True)
   ```

------------------------------------------------
四、HTML 交互报告（第三层 report）

技术栈  
- 模板：Jinja2  
- 图：Plotly（交互） + Scanpy 静态图备用  
- 表：DataTables（前端排序/搜索）  
- 文件：单文件 `qc_report.html`（base64 内嵌图片）

目录结构（自动生成）
```
qc_report.html
├── Summary
├── Cells QC  ( violin + scatter )
├── Genes QC
├── Doublet   ( UMAP 上色 )
├── SoupX     ( 背景占比条形图 )
└── Download  ( 过滤后的 h5ad/CSV 链接 )
```

------------------------------------------------
五、CLI & 配置（第四层 cli）

```bash
# 一键全流程
scqckit run \
  --matrix ./raw/matrix.mtx \
  --feature ./raw/features.tsv \
  --barcode ./raw/barcodes.tsv \
  --mt-prefix mt- \
  --doublet-method scrublet \
  --soupx-auto \
  --html qc_report.html \
  --out-dir ./qc_result
```

配置文件 `qc_config.yaml` 存阈值，方便实验室统一标准。

------------------------------------------------
六、开发路线图（MVP → 1.0）

| 阶段 | 目标 | 时间 |
|---|---|---|
| **MVP** | 细胞+基因基础 QC + HTML 报告 | 2 周 |
| **0.5** | 接入 scrublet + SoupX（rpy2） | 2 周 |
| **0.8** | CellBender GPU、DoubletDetection 多算法 | 2 周 |
| **1.0** | 完整 CLI、conda 包、在线文档、CI 测试 | 2 周 |

------------------------------------------------
七、技术细节小贴士

1. **统一数据容器**  
   全程 `AnnData`，兼容 Scanpy/Seurat/SoupX/CellBender。

2. **R 桥接**  
   SoupX 用 `rpy2` 封装成 Python 函数，用户无感调用 R。

3. **大样本性能**  
   - scrublet 支持多核 `n_jobs`  
   - CellBender 自动检测 GPU/CPU 并给出时间预估。

4. **HTML 体积**  
   >10 万细胞时，Plotly 图启用 `scattergl` + 降采样（`sc.pp.subsample`）前端交互。

5. **插件化**  
   预留 `entry_points` 接口，后续可插新算法（如 `solo`、`decontX`）。

------------------------------------------------
八、交付物 & 用户价值

- `pip install scQCkit`  
- 一条命令完成 **质控+去污染+去双胞+交互报告**  
- 报告单文件 HTML，手机/电脑都能看，可直接邮件发送。  
- 实验室无需再手写 R 脚本，统一阈值、统一格式，方便审计与投稿。

------------------------------------------------
下一步可立即开始：  
1. 建仓库 `scQCkit` → 骨架 `cookiecutter-pypackage`；  
2. 先写 `qc/core.py` + `report/base.py` 跑出 MVP 闭环；  
3. 逐步接入 scrublet/SoupX，补测试 + 文档。  

需要代码脚手架或具体算法封装示例，再喊我！
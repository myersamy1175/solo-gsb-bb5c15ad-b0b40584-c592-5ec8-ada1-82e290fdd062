# 🧬 PCR 引物互补性筛选台

在网页中粘贴引物序列 JSON，设置**盐浓度**与**退火温度（Tm 阈值）**，后端对每对引物做**反向互补局部比对**，分别检测：

- **发卡结构（Hairpin）**——枚举环位置，强制 5′ 臂与 3′ 臂分居环两侧做局部比对，定位两臂区间与环长；
- **引物二聚体（Primer Dimer）**——包含交叉二聚体与自二聚体，标记是否为 **3′ 锚定（可延伸）**；

按**简化最近邻 Tm 模型**计算配对区熔解温度，依据退火温度进行**风险分级**，输出**可定位的碱基级配对区间、对齐图和建议淘汰的引物/引物对**。支持示例数据一键载入和完整结果 **JSON 导出**。

---

## 技术栈

| 层 | 实现 |
|---|---|
| 后端 | FastAPI + Pydantic |
| 比对 | Smith-Waterman 局部比对（affine gap，纯 Python 实现） |
| Tm | SantaLucia (1998) 最近邻模型简化版（见下） |
| 存储 | SQLite（标准库 `sqlite3`，零额外依赖） |
| 前端 | 原生 HTML + CSS + JavaScript（无框架、无构建） |
| 部署 | Dockerfile + Docker Compose 一键启动 |

---

## 一键启动（Docker Compose）

```bash
docker compose up -d --build
```

启动后浏览器访问：

```
http://localhost:8000
```

停止：

```bash
docker compose down
```

> SQLite 数据保存在命名卷 `pcr-data`（容器内 `/data/pcr_screen.db`），重建容器后历史记录保留。

### 不用 Compose 时

```bash
docker build -t pcr-primer-screen .
docker run -d -p 8000:8000 -v "$PWD/data:/data" pcr-primer-screen
```

### 本地直接运行（开发）

```bash
pip install -r requirements.txt
PCR_DB_PATH=./pcr_screen.db uvicorn app.main:app --reload --port 8000
```

---

## 端口与配置

| 项 | 默认值 | 说明 |
|---|---|---|
| 服务端口 | `8000` | Compose 中映射宿主 `8000` → 容器 `8000`，可在 `docker-compose.yml` 修改 |
| `PCR_DB_PATH` | `/data/pcr_screen.db` | SQLite 文件路径，环境变量覆盖 |
| 退火温度 | `50 ℃` | Tm 阈值，网页可调（0–90） |
| 盐浓度 `[Na⁺]` | `100 mM` | 网页可调（1–1000），参与盐校正 |
| 引物浓度 | `250 nM` | 网页可调（1–100000），用于 Tm 浓度项 |
| 发卡最小环长 | `3 nt` | 网页可调（0–20） |

---

## 输入格式

网页文本框接受 JSON，支持三种写法：

```json
{
  "parameters": { "temperature": 50, "salt_mm": 100, "primer_nm": 250, "min_loop": 3 },
  "primers": [
    { "name": "F1", "sequence": "ATGCCTAGGCTAACGGCGATCCAGTCGGC" },
    { "name": "R1", "sequence": "TTAGCGTACGGTACGGCCGACTGGATCGC" }
  ]
}
```

也可省略名称（自动命名 `P1, P2, …`）：

```json
{ "primers": ["ATGCCTAGGCTAACG", "TTAGCGTACGGTACG"] }
```

约束：仅接受 `A/C/G/T`（自动去除空白、忽略大小写）；单条长度 8–120 nt；最多 200 条。非法输入返回 HTTP 422 与中文错误说明。

---

## 风险分级与淘汰规则

按配对区 Tm（℃）与退火温度 `T` 比较：

| 等级 | 条件 | 含义 |
|---|---|---|
| 🔴 严重 severe | `Tm ≥ T + 5` | 在退火温度下非常稳定，强烈建议淘汰 |
| 🟠 高风险 high | `T ≤ Tm < T+5` | 退火时仍可能形成 |
| 🟡 中风险 medium | `T − 10 ≤ Tm < T` | 边缘，结合 3′ 锚定情况判断 |
| 🟢 低风险 low | `Tm < T − 10` | 仅展示 |

**淘汰建议**：

- 引物存在 🔴 发卡，或 🟠 且发卡 3′ 端参与 → 淘汰该引物；
- 引物对存在 🔴 二聚体，或 🟠 且 **3′ 锚定**（两条 3′ 末端互补，聚合酶可延伸，危害最大）→ 淘汰该引物对。

低于退火温度 15 ℃ 的弱命中不展示，避免噪声。

---

## Tm 模型（简化）

- 对配对区中的每一段**连续完全匹配**子序列，用 SantaLucia 1998 的 16 个最近邻二联体参数（ΔH/ΔS）+ 末端起始校正求和；
- 取最高的一段作为该配对区 Tm，错配与缺口只作稳定性扣分与标记；
- 盐校正：`ΔS(salt) = ΔS + 0.368·(N−1)·ln[Na⁺]`；
- 发卡按环大小施加经验性温度惩罚；
- 短于 4 bp 的匹配段回退到 Wallace 规则 `2·(A+T) + 4·(G+C)`。

> 这是面向**快速风险排序**的近似模型，适合批量初筛，不替代 Oligo/IDT 等完整热力学工具的精算结果。

---

## HTTP API

| 方法 | 路径 | 说明 |
|---|---|---|
| `GET` | `/api/health` | 健康检查 |
| `GET` | `/api/sample` | 内置示例数据 |
| `POST` | `/api/screen` | 提交筛选（参数同上），返回完整结果并落库 |
| `GET` | `/api/runs` | 历史记录列表 |
| `GET` | `/api/runs/{id}` | 读取某次完整结果 |
| `DELETE` | `/api/runs/{id}` | 删除历史记录 |

交互式文档：`http://localhost:8000/docs`

`POST /api/screen` 返回结构（节选）：

```json
{
  "run_id": 1,
  "summary": {
    "primer_count": 4,
    "hairpin_high_risk": 1,
    "dimer_high_risk": 7,
    "rejected_primers": ["HP1"],
    "rejected_pairs": [["F1", "R1"]]
  },
  "hairpins": [
    {
      "primer": "HP1", "risk": "severe", "tm": 60.6,
      "loop_nt": 3, "longest_match": 14, "anchored_3prime": false,
      "arms": [
        {"start": 1, "end": 14, "strand": "5'-arm"},
        {"start": 18, "end": 31, "strand": "3'-arm"}
      ],
      "alignment": { "top": "...", "mid": "||||", "bottom": "..." }
    }
  ],
  "dimers": [
    {
      "primer_a": "F1", "primer_b": "R1", "self_dimer": false,
      "risk": "severe", "tm": 56.1, "longest_match": 14,
      "anchored_3prime": true,
      "regions": [
        {"start": 16, "end": 29, "primer_role": "A (5'→3')"},
        {"start": 16, "end": 29, "primer_role": "B (3'→5')"}
      ]
    }
  ],
  "rejections": [
    {"kind": "primer", "target": "HP1", "reason": "发卡结构 Tm≈60.6℃..."},
    {"kind": "pair",   "target": ["F1","R1"], "reason": "引物二聚体 Tm≈56.1℃..."}
  ]
}
```

所有区间均为 **1-based、闭区间**，按引物 5′→3′ 编号；二聚体的 B 区间标注为反平行 3′→5′。

网页「**导出 JSON**」按钮即把该完整结果下载为 `pcr_screen_run{id}.json`。

---

## 内置示例

`app/data/sample_primers.json` 提供 4 条引物：

- **F1 / R1**：3′ 末端 14 nt 反向互补，形成 🔴 **3′ 锚定引物二聚体**（约 56 ℃），建议淘汰该引物对；
- **HP1**：含 14 bp 茎 + 3 nt 环的高 Tm 发卡（约 61 ℃），建议淘汰；3′ 端额外 4 nt 使其发卡本身不锚定 3′；
- **C1**：低风险对照引物。

点击网页「**载入示例**」→「**开始筛选**」即可看到上述位点。

---

## 测试

```bash
python -m unittest discover -s tests -v
```

覆盖反向互补、Tm 单调性与环惩罚、14 bp 完美 3′ 二聚体定位、发卡两臂区间与反向互补关系、温度阈值敏感性等。

---

## 项目结构

```
.
├── Dockerfile
├── docker-compose.yml
├── requirements.txt
├── README.md
├── app
│   ├── main.py            # FastAPI 路由
│   ├── screening.py       # 发卡/二聚体检测、Tm 阈值、风险分级、淘汰建议
│   ├── alignment.py       # Smith-Waterman 局部比对
│   ├── seq_utils.py       # 反向互补 + 简化 NN Tm
│   ├── schemas.py         # 入参解析与校验
│   ├── database.py        # SQLite 持久化
│   ├── data/sample_primers.json
│   └── static/            # 原生前端 index.html / style.css / app.js
└── tests/test_core.py
```

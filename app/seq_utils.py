"""序列处理与简化热力学熔解温度(Tm)计算。

Tm 采用 SantaLucia (1998) 最近邻(Nearest-Neighbor)模型的简化版：
- 对配对区中的每一段「连续完全匹配」子序列（无错配/缺口）分别求和 NN 参数；
- 取最高的一段作为该配对区的 Tm（错配与缺口仅作为稳定性扣分项与标记，不参与求和）；
- 盐校正： dS(salt) = dS + 0.368 * (N-1) * ln([Na+])；
- 发卡环额外施加经验性温度惩罚；
- 过短(<4 bp)的匹配段回退到 Wallace 规则 2*(A+T) + 4*(G+C)。

这是面向筛选的近似模型，足以给出相对风险排序，但不等同于完整热力学软件的结果。
"""

from __future__ import annotations

import math

R = 1.987  # cal/(mol*K) 气体常数

# SantaLucia 1998：16 个二联体的 dH(kcal/mol), dS(cal/(mol*K))
_NN = {
    "AA": (-7.9, -22.2),
    "AC": (-8.4, -22.4),
    "AG": (-7.8, -21.0),
    "AT": (-7.2, -20.4),
    "CA": (-8.5, -22.7),
    "CC": (-8.0, -19.9),
    "CG": (-10.6, -27.2),
    "CT": (-7.8, -21.0),
    "GA": (-8.2, -22.2),
    "GC": (-9.8, -24.7),
    "GG": (-8.0, -19.9),
    "GT": (-8.4, -22.4),
    "TA": (-7.2, -21.3),
    "TC": (-8.2, -22.2),
    "TG": (-8.5, -22.7),
    "TT": (-7.9, -22.2),
}

# 末端起始校正（末端碱基类型相关）
_INIT_GC = (0.1, -2.8)  # dH kcal, dS cal
_INIT_AT = (2.3, 4.1)

_COMP = str.maketrans("ACGT", "TGCA")


def normalize(seq: str) -> str:
    """大写化并去除空白；仅接受 ACGT。"""
    return "".join(seq.upper().split())


def validate(seq: str) -> None:
    bad = sorted({c for c in seq.upper() if c not in "ACGT" and not c.isspace()})
    if bad:
        raise ValueError(f"序列含非法碱基: {','.join(bad)}（仅支持 A/C/G/T）")
    if len(seq.strip()) == 0:
        raise ValueError("序列为空")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def gc_count(seq: str) -> int:
    return sum(1 for c in seq if c in "GC")


def _nn_run_tm(run: str, salt_m: float, primer_nm: float, self_comp: bool) -> float:
    """对一段连续完全匹配序列计算 Tm(℃)。

    primer_nm 单位 nM；self_comp 为回文序列时使用对称性校正。
    """
    n = len(run)
    if n < 2:
        return 0.0

    dh = 0.0  # kcal/mol
    ds = 0.0  # cal/(mol*K)
    for i in range(n - 1):
        h, s = _NN[run[i : i + 2]]
        dh += h
        ds += s

    h0, s0 = _INIT_GC if run[0] in "GC" else _INIT_AT
    dh += h0
    ds += s0
    h0, s0 = _INIT_GC if run[-1] in "GC" else _INIT_AT
    dh += h0
    ds += s0

    # 盐校正（限制到合理范围，避免极端输入失真）
    na = min(max(salt_m, 0.001), 1.0)
    ds += 0.368 * (n - 1) * math.log(na)

    ct = (primer_nm / 1e9) * 4.0  # 总链浓度 mol/L，保守取较高的一条
    if self_comp:
        ds += 1.4  # 对称性校正
        kt = ct
    else:
        kt = ct / 4.0

    tm_k = (1000.0 * dh) / (ds + R * math.log(kt))
    return tm_k - 273.15


def _wallace_tm(run: str, salt_m: float) -> float:
    at = sum(1 for c in run if c in "AT")
    gc = sum(1 for c in run if c in "GC")
    return 2.0 * at + 4.0 * gc + 16.6 * math.log10(max(salt_m, 0.001) / 0.05)


def duplex_tm(
    top_run: str,
    salt_m: float,
    primer_nm: float,
    loop_nt: int | None = None,
) -> float:
    """根据「连续完全匹配段」序列（5'->3' 任一侧链）估算双链 Tm。

    loop_nt 非 None 时按发卡环施加经验惩罚。
    """
    if not top_run:
        return 0.0
    self_comp = top_run == reverse_complement(top_run)
    if len(top_run) >= 4:
        tm = _nn_run_tm(top_run, salt_m, primer_nm, self_comp)
    else:
        tm = _wallace_tm(top_run, salt_m)

    if loop_nt is not None:
        # 简化的发卡环惩罚：最小 3 nt 环起罚，环越大越不稳定
        extra = max(loop_nt - 3, 0)
        tm -= 2.0 + 0.4 * extra
    return tm

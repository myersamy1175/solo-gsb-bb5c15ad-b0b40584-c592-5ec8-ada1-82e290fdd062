"""互补性筛选：发卡结构、引物二聚体检测与风险分级。

坐标约定：全部使用 1-based、闭区间，并标注 5'/3' 方向。
- top 链 = 引物 A，方向 5'->3'
- bottom 链 = 引物 B（发卡时为同一条引物），显示方向 3'->5'
"""

from __future__ import annotations

from typing import List

from . import seq_utils
from .alignment import sw_align


def _max_run_tm(top: str, mid: str, salt_m: float, primer_nm: float, loop_nt):
    """取比对中最长的连续匹配段估算 Tm。"""
    best_tm = 0.0
    longest = 0
    run = []
    for a, mark in zip(top, mid):
        if mark == "|":
            run.append(a)
        else:
            if run:
                longest = max(longest, len(run))
                best_tm = max(
                    best_tm,
                    seq_utils.duplex_tm("".join(run), salt_m, primer_nm, loop_nt),
                )
                run = []
    if run:
        longest = max(longest, len(run))
        best_tm = max(
            best_tm,
            seq_utils.duplex_tm("".join(run), salt_m, primer_nm, loop_nt),
        )
    return best_tm, longest


def _risk(tm: float, temp: float) -> str:
    if tm >= temp + 5:
        return "severe"   # 严重
    if tm >= temp:
        return "high"     # 高风险（Tm 达到退火温度）
    if tm >= temp - 10:
        return "medium"
    return "low"


def _map_rc(start_rc: int, end_rc: int, length: int):
    """rc(seq) 上的 0-based 闭区间 -> 原 seq 上的 1-based 闭区间。

    rc[j] = comp(seq[L-1-j])，故 j 对应原索引 L-1-j。
    """
    return length - end_rc, length - start_rc


def _hairpin_record(
    hit, arm5: tuple, arm3: tuple, loop_nt: int, length: int, tm: float, longest: int
):
    """arm5/arm3 为 (1based_lo, 1based_hi)。"""
    reaches_3p = arm3[1] >= length - 2
    return {
        "type": "hairpin",
        "primer": None,  # 由调用方填充
        "tm": round(tm, 1),
        "longest_match": longest,
        "matched": hit["matched"],
        "mismatches": hit["mismatches"],
        "gaps": hit["gaps"],
        "score": hit["score"],
        "loop_nt": loop_nt,
        "arms": [
            {"start": arm5[0], "end": arm5[1], "strand": "5'-arm"},
            {"start": arm3[0], "end": arm3[1], "strand": "3'-arm"},
        ],
        "anchored_3prime": reaches_3p,
        "risk": None,  # 由调用方填充
        "alignment": {
            "top_label": "5'-arm",
            "top": hit["top"],
            "mid": hit["mid"],
            "bottom": hit["bottom"],
            "bottom_label": "3'-arm",
        },
    }


def _find_hairpins(seq: str, min_loop: int, salt_m: float, primer_nm: float, keep_floor: float, temperature: float):
    """枚举环位置，强制 5' 臂与 3' 臂分居环两侧后做局部比对。

    对每个环起点 k（左臂取 seq[:k]，右臂取 seq[k+min_loop:]），
    比对左臂 vs rc(右臂)；命中的两段即为发卡的两臂。
    多个 k 可能给出同一发卡，用臂区间去重，保留 Tm 最高者。
    """
    n = len(seq)
    seen: set = set()
    records = []
    for k in range(2, n - min_loop - 2):
        left = seq[:k]
        right = seq[k + min_loop:]
        if len(left) < 4 or len(right) < 4:
            continue
        for hit in sw_align(left, seq_utils.reverse_complement(right)):
            # top 命中区间（左臂, 0-based）
            a_lo, a_hi = hit["start1"], hit["end1"]
            # bottom 命中在 right 内（0-based），换算到全序列索引
            j_lo, j_hi = hit["start2"], hit["end2"]
            right_len = len(right)
            b_hi0 = (k + min_loop) + (right_len - 1 - j_lo)
            b_lo0 = (k + min_loop) + (right_len - 1 - j_hi)
            # 必须是左臂在前、右臂在后（天然满足），且两臂不重叠
            if b_lo0 <= a_hi:
                continue
            loop = b_lo0 - a_hi - 1
            if loop < min_loop:
                continue
            tm, longest = _max_run_tm(
                hit["top"], hit["mid"], salt_m, primer_nm, loop
            )
            if tm < keep_floor or longest < 4:
                continue
            arm5 = (a_lo + 1, a_hi + 1)
            arm3 = (b_lo0 + 1, b_hi0 + 1)
            key = (arm5, arm3)
            if key in seen:
                continue
            seen.add(key)
            rec = _hairpin_record(hit, arm5, arm3, loop, n, tm, longest)
            rec["risk"] = _risk(tm, temperature)
            records.append(rec)

    # 同一发卡的嵌套子命中去重：若某命中的两臂都被另一命中包含，丢弃
    records.sort(key=lambda r: (r["arms"][0]["start"], -r["arms"][0]["end"], -r["tm"]))
    kept: list = []
    for rec in records:
        a = rec["arms"][0]
        b = rec["arms"][1]
        nested = any(
            other["arms"][0]["start"] <= a["start"]
            and other["arms"][0]["end"] >= a["end"]
            and other["arms"][1]["start"] <= b["start"]
            and other["arms"][1]["end"] >= b["end"]
            and other is not rec
            for other in records
        )
        if not nested:
            kept.append(rec)
    kept.sort(key=lambda r: r["tm"], reverse=True)
    return kept


def _dimer_record(hit, length_a: int, length_b: int, tm: float, longest: int):
    a0, a1 = hit["start1"] + 1, hit["end1"] + 1
    b0, b1 = _map_rc(hit["start2"], hit["end2"], length_b)
    anchor_a = a1 >= length_a - 2                       # top 右端触及引物 A 的 3' 末端
    anchor_b = min(b0, b1) <= 3                         # bottom 左端对应引物 B 的 3' 末端
    return {
        "type": "dimer",
        "tm": round(tm, 1),
        "longest_match": longest,
        "matched": hit["matched"],
        "mismatches": hit["mismatches"],
        "gaps": hit["gaps"],
        "score": hit["score"],
        "anchored_3prime": anchor_a or anchor_b,
        "risk": None,
        "regions": [
            {"start": a0, "end": a1, "primer_role": "A (5'→3')"},
            {"start": min(b0, b1), "end": max(b0, b1), "primer_role": "B (3'→5')"},
        ],
        "alignment": {
            "top_label": "5'-A",
            "top": hit["top"],
            "mid": hit["mid"],
            "bottom": hit["bottom"],
            "bottom_label": "3'-B",
        },
    }


def screen(primers: List[dict], params: dict) -> dict:
    """执行全部筛选并返回结构化结果。

    primers: [{"name": str, "sequence": str}]（已校验）
    params: temperature(℃), salt_mm, primer_nm, min_loop, min_tm_drop 等
    """
    temperature = float(params["temperature"])
    salt_m = float(params["salt_mm"]) / 1000.0
    primer_nm = float(params["primer_nm"])
    min_loop = int(params["min_loop"])
    keep_floor = temperature - 15  # 低于退火温度 15℃ 的弱命中不展示

    hairpins: List[dict] = []
    dimers: List[dict] = []

    seqs = {i: p["sequence"] for i, p in enumerate(primers)}
    rcs = {i: seq_utils.reverse_complement(p["sequence"]) for i, p in enumerate(primers)}

    # 1) 发卡结构：枚举环位置，比对左右两臂
    for i, p in enumerate(primers):
        for rec in _find_hairpins(
            seqs[i], min_loop, salt_m, primer_nm, keep_floor, temperature
        ):
            rec["primer"] = p["name"]
            hairpins.append(rec)

    # 2) 自二聚体：引物整体与其反向互补比对
    for i, p in enumerate(primers):
        seq = seqs[i]
        for hit in sw_align(seq, rcs[i]):
            tm, longest = _max_run_tm(
                hit["top"], hit["mid"], salt_m, primer_nm, None
            )
            if tm < keep_floor:
                continue
            rec = _dimer_record(hit, len(seq), len(seq), tm, longest)
            rec["primer_a"] = p["name"]
            rec["primer_b"] = p["name"]
            rec["self_dimer"] = True
            rec["risk"] = _risk(tm, temperature)
            dimers.append(rec)

    # 2) 两两引物比对：交叉引物二聚体
    for i in range(len(primers)):
        for j in range(i + 1, len(primers)):
            hits = sw_align(seqs[i], rcs[j])
            for hit in hits:
                tm, longest = _max_run_tm(
                    hit["top"], hit["mid"], salt_m, primer_nm, None
                )
                if tm < keep_floor:
                    continue
                rec = _dimer_record(hit, len(seqs[i]), len(seqs[j]), tm, longest)
                rec["primer_a"] = primers[i]["name"]
                rec["primer_b"] = primers[j]["name"]
                rec["self_dimer"] = False
                rec["risk"] = _risk(tm, temperature)
                dimers.append(rec)

    # 二聚体去重：同一引物对、同一物理区间只保留 Tm 最高的一条
    deduped: dict = {}
    for rec in dimers:
        key = (
            rec["primer_a"],
            rec["primer_b"],
            (rec["regions"][0]["start"], rec["regions"][0]["end"]),
            (rec["regions"][1]["start"], rec["regions"][1]["end"]),
        )
        if key not in deduped or rec["tm"] > deduped[key]["tm"]:
            deduped[key] = rec
    dimers = sorted(deduped.values(), key=lambda r: r["tm"], reverse=True)

    hairpins.sort(key=lambda r: r["tm"], reverse=True)

    # 3) 淘汰建议
    rejections = []
    rejected_primers = set()
    for rec in hairpins:
        if rec["risk"] == "severe" or (rec["risk"] == "high" and rec["anchored_3prime"]):
            rejected_primers.add(rec["primer"])
            rejections.append(
                {
                    "kind": "primer",
                    "target": rec["primer"],
                    "reason": (
                        f"发卡结构 Tm≈{rec['tm']}℃（环 {rec['loop_nt']} nt，"
                        + ("3' 端参与，" if rec["anchored_3prime"] else "")
                        + f"配对 {rec['longest_match']} bp），在 {temperature:g}℃ 下稳定形成"
                    ),
                }
            )

    rejected_pairs = set()
    for rec in dimers:
        key = tuple(sorted((rec["primer_a"], rec["primer_b"])))
        # 已因发卡被淘汰的引物，其自二聚体不再单独建议
        if rec["self_dimer"] and rec["primer_a"] in rejected_primers:
            continue
        if rec["risk"] == "severe" or (rec["risk"] == "high" and rec["anchored_3prime"]):
            rejected_pairs.add(key)
            rejections.append(
                {
                    "kind": "pair",
                    "target": [rec["primer_a"], rec["primer_b"]],
                    "reason": (
                        f"{'自' if rec['self_dimer'] else ''}引物二聚体 Tm≈{rec['tm']}℃"
                        f"（连续配对 {rec['longest_match']} bp，"
                        + ("3' 端互补，" if rec["anchored_3prime"] else "")
                        + f"区间 A:{rec['regions'][0]['start']}-{rec['regions'][0]['end']}"
                        f" / B:{rec['regions'][1]['start']}-{rec['regions'][1]['end']}）"
                    ),
                }
            )

    high_hairpins = [r for r in hairpins if r["risk"] in ("high", "severe")]
    high_dimers = [r for r in dimers if r["risk"] in ("high", "severe")]

    return {
        "parameters": {
            "temperature": temperature,
            "salt_mm": float(params["salt_mm"]),
            "primer_nm": primer_nm,
            "min_loop": min_loop,
            "tm_model": "SantaLucia1998 NN（连续匹配段近似）",
        },
        "summary": {
            "primer_count": len(primers),
            "hairpin_total": len(hairpins),
            "hairpin_high_risk": len(high_hairpins),
            "dimer_total": len(dimers),
            "dimer_high_risk": len(high_dimers),
            "rejected_primers": sorted(rejected_primers),
            "rejected_pairs": [list(k) for k in sorted(rejected_pairs)],
        },
        "hairpins": hairpins,
        "dimers": dimers,
        "rejections": rejections,
        "primers": primers,
    }

"""Smith-Waterman 局部比对：在引物与其反向互补序列之间寻找局部互补配对区。

采用 affine gap（开 gap / 延 gap 不同罚分）三矩阵动态规划：
- M[i][j]：s1[i] 与 s2[j] 配对结束（允许错配）
- X[i][j]：s1 上为缺口结束（水平方向消耗 s2）
- Y[i][j]：s2 上为缺口结束（竖直方向消耗 s1）

对每个得分峰回溯出碱基级比对，从而定位引物上的具体区间，
并统计匹配/错配/缺口数；反复提取不重叠的命中(HSP)。
"""

from __future__ import annotations

from typing import List

MATCH = 2
MISMATCH = -2
GAP_OPEN = -3
GAP_EXT = -1

MIN_SCORE = 8          # 命中最低得分
MIN_MATCHED_BP = 4     # 至少匹配的碱基对数
MAX_HITS = 12
OVERLAP_LIMIT = 0.5    # 与已有命中重叠超过该比例则丢弃

_NEG = -10 ** 9
# s2 始终以「反向互补」形式传入：rc(seq)[j] = comp(seq[L-1-j])。
# 因此两条链物理上互补配对时，矩阵同位置的两个碱基恰好「相同」，
# 相同(A-A/T-T/G-G/C-C)记为匹配，不同记为错配。
_MATCH_BASES = {"A", "T", "G", "C"}


def _pair_score(a: str, b: str) -> int:
    return MATCH if a == b else MISMATCH


def sw_align(s1: str, s2: str) -> List[dict]:
    """返回 s1 vs s2 的全部不重叠局部命中（按得分降序）。"""
    n, m = len(s1), len(s2)
    M = [[0] * (m + 1) for _ in range(n + 1)]
    X = [[_NEG] * (m + 1) for _ in range(n + 1)]
    Y = [[_NEG] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = max(M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1])
            M[i][j] = max(0, diag + _pair_score(s1[i - 1], s2[j - 1]))
            X[i][j] = max(M[i][j - 1] + GAP_OPEN + GAP_EXT, X[i][j - 1] + GAP_EXT)
            Y[i][j] = max(M[i - 1][j] + GAP_OPEN + GAP_EXT, Y[i - 1][j] + GAP_EXT)

    def traceback(pi: int, pj: int) -> dict:
        i, j = pi, pj
        chars1: list[str] = []
        chars2: list[str] = []
        matched = gaps = mismatches = 0
        pairs: set[tuple[int, int]] = set()

        while M[i][j] > 0:
            a, b = s1[i - 1], s2[j - 1]
            from_m = M[i - 1][j - 1]
            from_x = X[i - 1][j - 1]
            from_y = Y[i - 1][j - 1]

            # 当前格是一次配对
            chars1.append(a)
            chars2.append(b)
            pairs.add((i - 1, j - 1))
            if _pair_score(a, b) == MATCH:
                matched += 1
            else:
                mismatches += 1
            i, j = i - 1, j - 1

            # 若该配对接在 X/Y gap 段之后，沿 gap 段继续回溯
            if from_x > _NEG / 2 and from_x >= from_m and from_x >= from_y:
                while True:  # X：s1 缺口，水平消耗 s2（先判后走）
                    chars1.append("-")
                    chars2.append(s2[j - 1])
                    pairs.add((i, j - 1))
                    gaps += 1
                    extends = j > 1 and X[i][j] == X[i][j - 1] + GAP_EXT
                    j -= 1
                    if not extends:
                        break
            elif from_y > _NEG / 2 and from_y >= from_m:
                while True:  # Y：s2 缺口，竖直消耗 s1（先判后走）
                    chars1.append(s1[i - 1])
                    chars2.append("-")
                    pairs.add((i - 1, j))
                    gaps += 1
                    extends = i > 1 and Y[i][j] == Y[i - 1][j] + GAP_EXT
                    i -= 1
                    if not extends:
                        break

        top = "".join(reversed(chars1))
        bottom = "".join(reversed(chars2))
        mid = "".join(
            "|" if a != "-" and b != "-" and _pair_score(a, b) == MATCH else " "
            for a, b in zip(top, bottom)
        )
        used_i = {ii for ii, _ in pairs}
        used_j = {jj for _, jj in pairs}
        return {
            "top": top,
            "bottom": bottom,
            "mid": mid,
            "matched": matched,
            "mismatches": mismatches,
            "gaps": gaps,
            "pairs": pairs,
            "start1": min(used_i),
            "end1": max(used_i),
            "start2": min(used_j),
            "end2": max(used_j),
        }

    peaks = sorted(
        ((M[i][j], i, j) for i in range(1, n + 1) for j in range(1, m + 1)),
        reverse=True,
    )

    hits: List[dict] = []
    claimed: set[tuple[int, int]] = set()
    for sc, i, j in peaks:
        if sc < MIN_SCORE:
            break
        if len(hits) >= MAX_HITS:
            break
        hit = traceback(i, j)
        if hit["matched"] < MIN_MATCHED_BP:
            continue
        overlap = len(hit["pairs"] & claimed)
        if overlap / max(len(hit["pairs"]), 1) >= OVERLAP_LIMIT:
            continue
        hit["score"] = sc
        claimed |= hit["pairs"]
        hits.append(hit)

    return hits

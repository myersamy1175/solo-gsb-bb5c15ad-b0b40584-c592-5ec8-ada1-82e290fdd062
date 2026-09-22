"""请求体解析与校验。

支持三种粘贴格式：
1. {"primers": [{"name": "F1", "sequence": "ATGC..."}]}
2. {"primers": ["ATGC...", "ATGC..."]}                      —— 自动命名 P1..
3. 直接粘贴 ["ATGC...", "ATGC..."] 或对象数组
"""

from __future__ import annotations

from typing import Any, List

from fastapi import HTTPException

from . import seq_utils

MAX_PRIMERS = 200
MIN_LEN, MAX_LEN = 8, 120


def parse_primers(payload: Any) -> List[dict]:
    if isinstance(payload, dict):
        data = payload.get("primers")
    else:
        data = payload

    if not isinstance(data, list) or not data:
        raise HTTPException(status_code=422, detail="primers 必须是非空数组")
    if len(data) > MAX_PRIMERS:
        raise HTTPException(status_code=422, detail=f"引物数量超过上限 {MAX_PRIMERS}")

    primers: List[dict] = []
    names: set[str] = set()
    for idx, item in enumerate(data, start=1):
        if isinstance(item, str):
            name, seq_raw = f"P{idx}", item
        elif isinstance(item, dict):
            seq_raw = item.get("sequence", "")
            name = str(item.get("name") or f"P{idx}").strip()
        else:
            raise HTTPException(status_code=422, detail=f"第 {idx} 条引物格式无效")

        if name in names:
            name = f"{name}-{idx}"
        names.add(name)

        try:
            seq = seq_utils.normalize(seq_raw)
            seq_utils.validate(seq)
        except ValueError as exc:
            raise HTTPException(status_code=422, detail=f"引物 {name}: {exc}")

        if not (MIN_LEN <= len(seq) <= MAX_LEN):
            raise HTTPException(
                status_code=422,
                detail=f"引物 {name}: 长度 {len(seq)} 超出范围 {MIN_LEN}-{MAX_LEN} nt",
            )
        primers.append({"name": name, "sequence": seq})

    return primers


def clamp_params(p: dict | None) -> dict:
    p = p or {}

    def num(key, default, lo, hi):
        try:
            v = float(p.get(key, default))
        except (TypeError, ValueError):
            v = default
        return min(max(v, lo), hi)

    return {
        "temperature": num("temperature", 50, 0, 90),
        "salt_mm": num("salt_mm", 100, 1, 1000),
        "primer_nm": num("primer_nm", 250, 1, 100000),
        "min_loop": int(num("min_loop", 3, 0, 20)),
    }

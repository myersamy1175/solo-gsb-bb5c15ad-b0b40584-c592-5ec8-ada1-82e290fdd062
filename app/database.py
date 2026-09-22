"""SQLite 持久化：保存每次筛选任务（参数 + 完整结果 JSON）。"""

from __future__ import annotations

import json
import os
import sqlite3
from datetime import datetime, timezone
from typing import Optional

DB_PATH = os.environ.get("PCR_DB_PATH", "/data/pcr_screen.db")

_SCHEMA = """
CREATE TABLE IF NOT EXISTS runs (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    created_at TEXT NOT NULL,
    temperature REAL NOT NULL,
    salt_mm REAL NOT NULL,
    primer_count INTEGER NOT NULL,
    rejected_count INTEGER NOT NULL,
    result_json TEXT NOT NULL
);
"""


def _connect() -> sqlite3.Connection:
    os.makedirs(os.path.dirname(DB_PATH) or ".", exist_ok=True)
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn


def init_db() -> None:
    with _connect() as conn:
        conn.execute(_SCHEMA)


def save_run(result: dict) -> int:
    init_db()
    p = result["parameters"]
    s = result["summary"]
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    rejected = len(s["rejected_primers"]) + len(s["rejected_pairs"])
    with _connect() as conn:
        cur = conn.execute(
            """INSERT INTO runs (created_at, temperature, salt_mm,
                                  primer_count, rejected_count, result_json)
               VALUES (?, ?, ?, ?, ?, ?)""",
            (
                now,
                p["temperature"],
                p["salt_mm"],
                s["primer_count"],
                rejected,
                json.dumps(result, ensure_ascii=False),
            ),
        )
        return int(cur.lastrowid)


def list_runs(limit: int = 30) -> list[dict]:
    init_db()
    with _connect() as conn:
        rows = conn.execute(
            """SELECT id, created_at, temperature, salt_mm, primer_count, rejected_count
               FROM runs ORDER BY id DESC LIMIT ?""",
            (limit,),
        ).fetchall()
        return [dict(r) for r in rows]


def get_run(run_id: int) -> Optional[dict]:
    init_db()
    with _connect() as conn:
        row = conn.execute(
            "SELECT result_json FROM runs WHERE id = ?", (run_id,)
        ).fetchone()
    return json.loads(row["result_json"]) if row else None


def delete_run(run_id: int) -> bool:
    init_db()
    with _connect() as conn:
        cur = conn.execute("DELETE FROM runs WHERE id = ?", (run_id,))
        return cur.rowcount > 0

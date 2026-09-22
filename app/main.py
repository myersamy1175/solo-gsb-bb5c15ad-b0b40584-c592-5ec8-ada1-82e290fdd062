"""FastAPI 入口：筛选接口、历史记录、示例数据、静态页面。"""

from __future__ import annotations

import json
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Any

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

from . import database
from .schemas import clamp_params, parse_primers
from .screening import screen

BASE_DIR = Path(__file__).resolve().parent
SAMPLE_PATH = BASE_DIR / "data" / "sample_primers.json"


@asynccontextmanager
async def lifespan(app: FastAPI):
    database.init_db()
    yield


app = FastAPI(
    title="PCR 引物互补性筛选台",
    description="检测引物发卡结构与引物二聚体，按简化 Tm 阈值评估风险并给出淘汰建议。",
    version="1.0.0",
    lifespan=lifespan,
)


@app.get("/api/health")
def health() -> dict:
    return {"status": "ok"}


@app.get("/api/sample")
def sample() -> Any:
    return JSONResponse(json.loads(SAMPLE_PATH.read_text(encoding="utf-8")))


@app.post("/api/screen")
def run_screen(payload: dict) -> dict:
    primers = parse_primers(payload)
    params = clamp_params(payload.get("parameters"))
    result = screen(primers, params)
    run_id = database.save_run(result)
    result["run_id"] = run_id
    return result


@app.get("/api/runs")
def get_runs() -> dict:
    return {"runs": database.list_runs()}


@app.get("/api/runs/{run_id}")
def get_run(run_id: int) -> dict:
    result = database.get_run(run_id)
    if result is None:
        raise HTTPException(status_code=404, detail="记录不存在")
    return result


@app.delete("/api/runs/{run_id}")
def remove_run(run_id: int) -> dict:
    if not database.delete_run(run_id):
        raise HTTPException(status_code=404, detail="记录不存在")
    return {"deleted": run_id}


STATIC_DIR = BASE_DIR / "static"
app.mount("/", StaticFiles(directory=STATIC_DIR, html=True), name="static")

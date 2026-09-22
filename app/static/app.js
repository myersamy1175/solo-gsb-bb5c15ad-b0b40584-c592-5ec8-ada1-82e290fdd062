// PCR 引物互补性筛选台 —— 原生前端逻辑
let lastResult = null;

const $ = (id) => document.getElementById(id);
const RISK_TEXT = { severe: "严重", high: "高风险", medium: "中风险", low: "低风险" };

async function api(path, options = {}) {
  const res = await fetch(path, {
    headers: { "Content-Type": "application/json" },
    ...options,
  });
  if (!res.ok) {
    let detail = res.statusText;
    try {
      const body = await res.json();
      detail = typeof body.detail === "string" ? body.detail : JSON.stringify(body.detail);
    } catch (_) {}
    throw new Error(detail);
  }
  return res.json();
}

function showError(msg) {
  const box = $("error");
  box.textContent = msg;
  box.hidden = false;
}
function clearError() {
  $("error").hidden = true;
}

function readParams() {
  return {
    temperature: parseFloat($("temperature").value),
    salt_mm: parseFloat($("saltMm").value),
    primer_nm: parseFloat($("primerNm").value),
    min_loop: parseInt($("minLoop").value, 10),
  };
}

async function loadSample() {
  const sample = await api("/api/sample");
  $("primerInput").value = JSON.stringify(
    { parameters: sample.parameters, primers: sample.primers },
    null,
    2
  );
  if (sample.parameters) {
    $("temperature").value = sample.parameters.temperature;
    $("saltMm").value = sample.parameters.salt_mm;
    $("primerNm").value = sample.parameters.primer_nm;
    $("minLoop").value = sample.parameters.min_loop;
  }
  clearError();
}

function buildPayload() {
  const raw = $("primerInput").value.trim();
  if (!raw) throw new Error("请先粘贴引物 JSON");
  let parsed;
  try {
    parsed = JSON.parse(raw);
  } catch (e) {
    throw new Error("JSON 解析失败：" + e.message);
  }
  return { ...parsed, parameters: { ...readParams(), ...(parsed.parameters || {}) } };
}

async function runScreen() {
  clearError();
  let payload;
  try {
    payload = buildPayload();
  } catch (e) {
    showError(e.message);
    return;
  }
  const btn = $("runScreen");
  btn.disabled = true;
  btn.textContent = "筛选中…";
  try {
    lastResult = await api("/api/screen", {
      method: "POST",
      body: JSON.stringify(payload),
    });
    render(lastResult);
    loadHistory();
  } catch (e) {
    showError(e.message);
  } finally {
    btn.disabled = false;
    btn.textContent = "开始筛选";
  }
}

function render(result) {
  $("summary").hidden = false;
  $("exportJson").disabled = false;
  renderSummary(result.summary);
  renderRejections(result.rejections);
  renderHairpins(result.hairpins);
  renderDimers(result.dimers);
  $("hpCount").textContent = result.hairpins.length;
  $("dmCount").textContent = result.dimers.length;
}

function renderSummary(s) {
  const cards = [
    { num: s.primer_count, lbl: "引物数", cls: "" },
    { num: s.hairpin_total, lbl: "发卡位点", cls: s.hairpin_high_risk ? "warn" : "ok" },
    { num: s.dimer_total, lbl: "二聚体位点", cls: s.dimer_high_risk ? "bad" : "ok" },
    {
      num: s.rejected_primers.length + s.rejected_pairs.length,
      lbl: "建议淘汰",
      cls: s.rejected_primers.length + s.rejected_pairs.length ? "bad" : "ok",
    },
  ];
  $("summaryCards").innerHTML = cards
    .map((c) => `<div class="card ${c.cls}"><div class="num">${c.num}</div><div class="lbl">${c.lbl}</div></div>`)
    .join("");
}

function renderRejections(rejections) {
  const box = $("rejectionBox");
  if (!rejections.length) {
    box.hidden = true;
    return;
  }
  box.hidden = false;
  $("rejectionList").innerHTML = rejections
    .map((r) => {
      const target = r.kind === "primer"
        ? `<code>引物 ${escapeHtml(r.target)}</code>`
        : `<code>引物对 ${escapeHtml(r.target.join(" × "))}</code>`;
      return `<li>${target} —— ${escapeHtml(r.reason)}</li>`;
    })
    .join("");
}

function riskTag(risk) {
  return `<span class="tag ${risk}">${RISK_TEXT[risk] || risk} · Tm</span>`;
}

function alignmentBlock(al) {
  return `
    <button class="align-toggle">▸ 查看碱基比对（可定位区间）</button>
    <div class="alignment"><span class="al-label">${al.top_label}</span>  ${escapeHtml(al.top)}
              ${escapeHtml(al.mid)}
<span class="al-label">${al.bottom_label}</span>  ${escapeHtml(al.bottom)}</div>`;
}

function renderHairpins(hairpins) {
  const pane = $("tab-hairpins");
  if (!hairpins.length) {
    pane.innerHTML = `<div class="empty">未发现发卡结构（达到展示阈值的位点）。</div>`;
    return;
  }
  pane.innerHTML = hairpins
    .map((h) => {
      const arms = h.arms
        .map((a) => `${a.strand} <b>${a.start}–${a.end}</b>`)
        .join("　|　");
      return `
      <div class="hit ${h.risk}">
        <div class="hit-head">
          <span class="hit-title">${escapeHtml(h.primer)}</span>
          ${riskTag(h.risk)}
          <span class="tag ${h.risk}">${h.tm}℃</span>
          ${h.anchored_3prime ? '<span class="tag anchor">3′ 端参与</span>' : ""}
        </div>
        <div class="hit-meta">
          茎环：环 <b>${h.loop_nt}</b> nt · 连续配对 <b>${h.longest_match}</b> bp
          · 匹配 ${h.matched} / 错配 ${h.mismatches} / 缺口 ${h.gaps}
        </div>
        <div class="hit-meta">可定位区间（1-based，5′→3′）：${arms}</div>
        ${alignmentBlock(h.alignment)}
      </div>`;
    })
    .join("");
  bindAlignToggles(pane);
}

function renderDimers(dimers) {
  const pane = $("tab-dimers");
  if (!dimers.length) {
    pane.innerHTML = `<div class="empty">未发现引物二聚体位点。</div>`;
    return;
  }
  pane.innerHTML = dimers
    .map((d) => {
      const [ra, rb] = d.regions;
      const title = d.self_dimer
        ? `${escapeHtml(d.primer_a)}（自二聚体）`
        : `${escapeHtml(d.primer_a)} × ${escapeHtml(d.primer_b)}`;
      return `
      <div class="hit ${d.risk}">
        <div class="hit-head">
          <span class="hit-title">${title}</span>
          ${riskTag(d.risk)}
          <span class="tag ${d.risk}">${d.tm}℃</span>
          ${d.anchored_3prime ? '<span class="tag anchor">3′ 锚定（可延伸）</span>' : ""}
        </div>
        <div class="hit-meta">
          连续配对 <b>${d.longest_match}</b> bp
          · 匹配 ${d.matched} / 错配 ${d.mismatches} / 缺口 ${d.gaps}
        </div>
        <div class="hit-meta">
          可定位区间（1-based）：引物A <b>${ra.start}–${ra.end}</b>
          ｜ 引物B <b>${rb.start}–${rb.end}</b>（反平行 3′→5′）
        </div>
        ${alignmentBlock(d.alignment)}
      </div>`;
    })
    .join("");
  bindAlignToggles(pane);
}

function bindAlignToggles(root) {
  root.querySelectorAll(".align-toggle").forEach((btn) => {
    btn.addEventListener("click", () => {
      const al = btn.nextElementSibling;
      const open = al.classList.toggle("show");
      btn.textContent = open ? "▾ 收起碱基比对" : "▸ 查看碱基比对（可定位区间）";
    });
  });
}

// ---- 历史记录 ----
async function loadHistory() {
  try {
    const data = await api("/api/runs");
    $("histCount").textContent = data.runs.length;
    const pane = $("tab-history");
    if (!data.runs.length) {
      pane.innerHTML = `<div class="empty">暂无历史记录。</div>`;
      return;
    }
    pane.innerHTML = `
      <table class="hist">
        <thead><tr><th>#</th><th>时间(UTC)</th><th>Tm阈值</th><th>盐(mM)</th><th>引物数</th><th>淘汰</th><th></th></tr></thead>
        <tbody>
          ${data.runs
            .map(
              (r) => `<tr data-id="${r.id}">
                <td>${r.id}</td><td>${r.created_at.replace("T", " ").slice(0, 19)}</td>
                <td>${r.temperature}℃</td><td>${r.salt_mm}</td>
                <td>${r.primer_count}</td><td>${r.rejected_count}</td>
                <td><button class="btn small ghost" data-del="${r.id}">删除</button></td>
              </tr>`
            )
            .join("")}
        </tbody>
      </table>`;
    pane.querySelectorAll("tr[data-id]").forEach((tr) => {
      tr.addEventListener("click", async (e) => {
        if (e.target.dataset.del) return;
        const result = await api(`/api/runs/${tr.dataset.id}`);
        lastResult = result;
        render(result);
      });
    });
    pane.querySelectorAll("[data-del]").forEach((btn) => {
      btn.addEventListener("click", async (e) => {
        e.stopPropagation();
        await api(`/api/runs/${btn.dataset.del}`, { method: "DELETE" });
        loadHistory();
      });
    });
  } catch (e) {
    /* 历史加载失败不阻塞主流程 */
  }
}

function exportJson() {
  if (!lastResult) return;
  const blob = new Blob([JSON.stringify(lastResult, null, 2)], { type: "application/json" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = `pcr_screen_run${lastResult.run_id || ""}.json`;
  a.click();
  URL.revokeObjectURL(url);
}

function escapeHtml(s) {
  return String(s).replace(/[&<>"']/g, (c) =>
    ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c])
  );
}

document.querySelectorAll(".tab").forEach((tab) => {
  tab.addEventListener("click", () => {
    document.querySelectorAll(".tab").forEach((t) => t.classList.remove("active"));
    tab.classList.add("active");
    document.querySelectorAll(".tabpane").forEach((p) => (p.hidden = true));
    $(`tab-${tab.dataset.tab}`).hidden = false;
  });
});

$("loadSample").addEventListener("click", loadSample);
$("clearInput").addEventListener("click", () => {
  $("primerInput").value = "";
  clearError();
});
$("runScreen").addEventListener("click", runScreen);
$("exportJson").addEventListener("click", exportJson);

loadHistory();

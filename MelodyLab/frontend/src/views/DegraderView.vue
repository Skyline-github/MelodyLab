<script setup>
// 降解剂评估（PROTAC DC50 共识排序）：MLP 集成 / ExtraTrees / HistGB 三族共识排名
// + conformal 90% 区间 + 适用域检查。定位是排序参考工具，不是绝对值预测器。
import { nextTick, onMounted, onUnmounted, ref } from 'vue'
import { api } from '../api'
import JobCard from '../components/JobCard.vue'

const avail = ref(null)
const name = ref(''), smilesText = ref('')
const job = ref(null), error = ref(''), busy = ref(false), busyExample = ref(false)
const report = ref(null)
const archive = ref([])
let timer, charts = []

onMounted(async () => {
  try { avail.value = await api.degraderAvailability() } catch (_) {}
  refreshArchive()
})
onUnmounted(() => { clearInterval(timer); charts.forEach(c => c.destroy()); charts = [] })

function refreshArchive() {
  api.degraderReports().then(a => archive.value = a).catch(() => {})
}

async function loadExample() {
  busyExample.value = true
  try {
    const ex = await api.degraderExample()
    smilesText.value = ex.text
    if (!name.value.trim()) name.value = ex.name
  } catch (e) { error.value = e.message } finally { busyExample.value = false }
}

async function submit() {
  error.value = ''
  if (!smilesText.value.trim()) { error.value = '请粘贴 SMILES 或载入示例'; return }
  busy.value = true
  try {
    const fd = new FormData()
    fd.append('name', name.value.trim())
    fd.append('smiles_text', smilesText.value)
    job.value = await api.degraderRank(fd)
    report.value = null
    poll()
  } catch (e) { error.value = e.message } finally { busy.value = false }
}

async function retrain() {
  if (!confirm('重新训练将清除模型缓存并用当前打包数据重训（约 3–5 分钟），确定？')) return
  error.value = ''
  try {
    job.value = await api.degraderRetrain()
    poll(true)
  } catch (e) { error.value = e.message }
}

function poll(isRetrain = false) {
  clearInterval(timer)
  timer = setInterval(async () => {
    if (!job.value) return
    try {
      job.value = await api.degraderJob(job.value.id)
      if (['completed', 'failed'].includes(job.value.status)) {
        clearInterval(timer)
        if (job.value.status === 'completed') {
          if (isRetrain) {
            avail.value = await api.degraderAvailability()
          } else if (job.value.stats?.report_id) {
            await loadReport(job.value.stats.report_id)
            refreshArchive()
          }
        }
      }
    } catch (_) {}
  }, 2000)
}

async function loadReport(rid) {
  charts.forEach(c => c.destroy()); charts = []
  report.value = null
  try {
    report.value = await api.degraderReport(rid)
    await nextTick()
    await ensureChartJs()
    renderChart()
    window.scrollTo({ top: 0, behavior: 'smooth' })
  } catch (e) { error.value = e.message }
}

function ensureChartJs() {
  if (window.Chart) return Promise.resolve()
  return new Promise((resolve, reject) => {
    const s = document.createElement('script')
    s.src = '/Chart.umd.min.js'
    s.onload = resolve
    s.onerror = () => reject(new Error('Chart.js 加载失败（本地文件缺失）'))
    document.head.appendChild(s)
  })
}

const fmtNm = (v) => {
  if (v == null) return '—'
  if (v >= 1000) return (v / 1000).toFixed(1) + 'k'
  if (v >= 100) return v.toFixed(0)
  if (v >= 10) return v.toFixed(1)
  return v.toFixed(2)
}
const lg = (v) => Math.log10(v)

function renderChart() {
  const Chart = window.Chart
  const rows = report.value.rows.filter(r => r.status === 'ok')
  if (!Chart || !rows.length) return
  const el = document.getElementById('dc50-chart')
  if (!el) return
  charts.push(new Chart(el, {
    type: 'bar',
    data: {
      labels: rows.map(r => r.compound_id),
      datasets: [
        { type: 'bar', label: '90% 区间', data: rows.map(r => [lg(r.lo90_nm), lg(r.hi90_nm)]),
          backgroundColor: 'rgba(98,86,217,.22)', borderColor: 'rgba(98,86,217,.55)',
          borderWidth: 1, borderSkipped: false, barPercentage: 0.72 },
        { type: 'scatter', label: '共识预测 DC50', data: rows.map((r, i) => ({ x: lg(r.pred_dc50_nm), y: i })),
          backgroundColor: '#d95f02', pointRadius: 4, pointHoverRadius: 6 },
      ],
    },
    options: {
      indexAxis: 'y', responsive: true, maintainAspectRatio: false,
      plugins: {
        legend: { labels: { boxWidth: 12, font: { size: 11 } } },
        tooltip: {
          callbacks: {
            label: (ctx) => {
              const r = rows[ctx.dataIndex]
              return ctx.dataset.type === 'scatter'
                ? ` 共识预测: ${fmtNm(r.pred_dc50_nm)} nM`
                : ` 90% 区间: ${fmtNm(r.lo90_nm)} – ${fmtNm(r.hi90_nm)} nM`
            },
          },
        },
      },
      scales: {
        x: {
          type: 'logarithmic', min: 0.01, max: 1e6,
          title: { display: true, text: 'DC50 (nM，对数轴)' },
          afterBuildTicks: (scale) => {
            scale.ticks = [0.01, 0.1, 1, 10, 100, 1000, 1e4, 1e5, 1e6].map(v => ({ value: v }))
          },
          ticks: {
            callback: (v) => v >= 1e6 ? '1M' : v >= 1000 ? (v / 1000) + 'k' : String(v),
            font: { size: 10.5 },
          },
          grid: { color: 'rgba(128,140,170,.15)' },
        },
        y: { grid: { display: false }, ticks: { font: { size: 10.5 } } },
      },
    },
  }))
}

function topId(row) { return row.rank_consensus <= 5.01 && row.status === 'ok' }
</script>

<template>
  <h1 class="page-title">降解剂评估（DC50 共识排序）</h1>
  <p class="page-sub">
    输入一批 PROTAC，用三族模型（MLP 集成 / ExtraTrees / HistGB，PROTAC-DB 训练）给出
    降解活性（DC50）的<b>共识排序</b>、90% 区间与适用域提示。
    本工具输出的是<b>相对排序参考</b>，不是绝对值预测：请以共识排名 + 区间 + 提示共同解读。
  </p>
  <div v-if="avail && !avail.available" class="alert err">引擎不可用：{{ avail.hint }}</div>
  <div v-if="error" class="alert err">{{ error }}</div>

  <div class="card">
    <div class="grid2">
      <label class="field">任务名称（可选）<input v-model="name" placeholder="如 Galectin-3 PROTAC 系列" /></label>
      <label class="field">查询文件（可选，CSV，需含 smiles 列）<input type="file" accept=".csv,.txt,.smi"
        @change="(e) => { const f = e.target.files[0]; if (f) { const rd = new FileReader(); rd.onload = () => smilesText = rd.result; rd.readAsText(f) } }" /></label>
    </div>
    <label class="field" style="margin-top:10px">SMILES 列表（每行一个，支持 "ID SMILES" 格式）
      <textarea v-model="smilesText" rows="6" placeholder="J52 O=C(COC[C@H]1O...)&#10;或直接每行一个 SMILES"></textarea>
    </label>
    <div style="margin-top:12px; display:flex; gap:10px; align-items:center; flex-wrap:wrap">
      <button class="btn" :disabled="busy" @click="submit">{{ busy ? '评估中…' : '开始排序' }}</button>
      <button class="btn ghost" :disabled="busyExample" @click="loadExample">{{ busyExample ? '载入中…' : '📋 载入示例（PROTAC-DB 公开分子）' }}</button>
      <span class="muted" style="font-size:12.5px">首次运行需先训练模型（约 3–5 分钟），之后每次排序秒级完成</span>
    </div>
  </div>

  <div v-if="job" style="margin-top:14px">
    <JobCard :job="job" @refresh="job = null" />
  </div>

  <!-- ============ 模型状态 ============ -->
  <div v-if="avail" class="card" style="margin-top:14px">
    <div style="display:flex; align-items:center; gap:12px; flex-wrap:wrap">
      <h3 style="margin:0">模型状态</h3>
      <span class="pill">{{ avail.model_ready ? '✅ 已训练（缓存）' : '⏳ 未训练（首次排序时自动训练）' }}</span>
      <span style="flex:1"></span>
      <button class="btn ghost sm" @click="retrain">♻ 重新训练</button>
    </div>
    <div v-if="avail.model" class="muted" style="font-size:12.5px; margin-top:8px">
      训练集 {{ avail.model.n_train }} 行（{{ avail.model.db_snapshot }}）·
      骨架 {{ avail.model.recipe.n_folds }} 折交叉验证 ·
      q90 = {{ avail.model.calibration.q90_residual_log }} log 单位 ·
      训练于 {{ avail.model.trained_at }}
      <template v-if="avail.model.calibration.cv">
        · CV R²：
        <span v-for="(m, f) in avail.model.calibration.cv" :key="f" style="margin-right:8px">
          {{ f }} {{ m.r2 ?? '—' }}</span>
      </template>
    </div>
    <div v-else class="muted" style="font-size:12.5px; margin-top:8px">
      数据：{{ avail.db_snapshot }}
    </div>
  </div>

  <!-- ============ 报告 ============ -->
  <template v-if="report">
    <div class="card" style="margin-top:14px">
      <div style="display:flex; gap:10px; flex-wrap:wrap; align-items:baseline">
        <h2 style="margin:0">{{ report.name || report.id }}</h2>
        <span class="muted" style="font-size:12.5px">{{ report.created_at }} ·
          模型{{ report.model_status === 'trained' ? '（本次首训）' : '（缓存）' }}</span>
      </div>
      <div class="cards3" style="margin-top:12px">
        <div class="stat"><b>{{ report.n_query - report.n_invalid }}</b><span>有效查询分子（无效 {{ report.n_invalid }} / 重复 {{ report.n_duplicates }}）</span></div>
        <div class="stat"><b>{{ report.calibration?.n_train }}</b><span>训练集行数（骨架 {{ report.calibration?.n_folds }} 折 CV 校准）</span></div>
        <div class="stat"><b>{{ report.calibration?.q90_residual_log }}</b><span>conformal q90（log 单位）——区间宽度的诚实下限</span></div>
      </div>
      <div class="chart-box" :style="{ height: Math.max(420, report.rows.filter(r => r.status === 'ok').length * 26 + 90) + 'px' }">
        <canvas id="dc50-chart"></canvas>
      </div>
      <div class="muted" style="font-size:12px; margin-top:4px">
        橙点 = 共识预测 DC50（三族排名均值的代表值）；紫色带 = 90% 区间（成员离散度与 conformal 下限取大）。
        点越靠左（DC50 越小）预测越强，但只有<b>相对位置</b>有意义。
      </div>
    </div>

    <div class="card" style="margin-top:14px">
      <h3>完整排序</h3>
      <table class="tbl">
        <tr>
          <th>共识排名</th><th>化合物</th><th>E3</th><th>预测 DC50 (nM)</th><th>90% 区间 (nM)</th>
          <th title="MLP 集成 / ExtraTrees / HistGB 各自的排名">族排名 M/E/H</th>
          <th title="与训练集最近邻的 Tanimoto 相似度">Tanimoto</th><th>提示</th>
        </tr>
        <tr v-for="r in report.rows" :key="r.compound_id" :style="topId(r) ? 'background:rgba(47,158,99,.07)' : ''">
          <td><b>{{ r.status === 'ok' ? r.rank_consensus : '—' }}</b></td>
          <td><b>{{ r.compound_id }}</b>
            <span v-if="r.in_db" title="PROTAC-DB 已收录">✅</span></td>
          <td class="muted">{{ r.e3 || '—' }}</td>
          <td v-if="r.status === 'ok'"><b>{{ fmtNm(r.pred_dc50_nm) }}</b></td>
          <td v-else colspan="2" style="color:var(--warn)">{{ r.warnings.join('；') }}</td>
          <td v-if="r.status === 'ok'">{{ fmtNm(r.lo90_nm) }} – {{ fmtNm(r.hi90_nm) }}</td>
          <td v-if="r.status === 'ok'" class="muted">{{ r.rank_mlp }} / {{ r.rank_et }} / {{ r.rank_hgb }}</td>
          <td v-if="r.status === 'ok'" :style="r.max_tanimoto < 0.35 ? 'color:var(--warn)' : ''">{{ r.max_tanimoto ?? '—' }}</td>
          <td v-if="r.status === 'ok'" class="muted" style="font-size:12px; max-width:340px">
            <span v-for="(w, i) in r.warnings" :key="i">{{ w }}<template v-if="i < r.warnings.length - 1">；</template></span>
            <span v-if="!r.warnings.length">—</span>
            <span v-if="r.in_db">（实测 {{ r.db_dc50 }} nM）</span>
          </td>
        </tr>
      </table>
      <div class="muted" style="margin-top:8px; font-size:12px">
        🟩 行底色 = 共识前 5；族排名 M/E/H 三列差异大说明该分子排名是模型相关的。
      </div>
    </div>

    <div class="card" style="margin-top:14px">
      <h3>方法与可信度边界</h3>
      <ul style="margin:8px 0 0; padding-left:20px; font-size:13px; line-height:1.7">
        <li v-for="(n, i) in report.notes" :key="i">{{ n }}</li>
      </ul>
    </div>
  </template>

  <!-- ============ 归档 ============ -->
  <div class="card" style="margin-top:16px">
    <h3>排序归档</h3>
    <div v-if="!archive.length" class="muted">还没有归档报告，完成一次排序后会出现在这里。</div>
    <table v-else class="tbl">
      <tr><th>任务</th><th>分子数</th><th>共识前 5</th><th>训练集</th><th>时间</th><th></th></tr>
      <tr v-for="r in archive" :key="r.id">
        <td><b>{{ r.name || r.id }}</b></td>
        <td class="muted">{{ r.n_query }}</td>
        <td class="muted">{{ (r.top5 || []).join('、') }}</td>
        <td class="muted">{{ r.n_train }}</td>
        <td class="muted">{{ r.created_at }}</td>
        <td style="white-space:nowrap">
          <button class="btn ghost sm" @click="loadReport(r.id)">查看</button>
          <a class="btn ghost sm" :href="`/api/degrader/reports/${r.id}/download/csv`">csv</a>
          <a class="btn ghost sm" :href="`/api/degrader/reports/${r.id}/download/md`">md</a>
        </td>
      </tr>
    </table>
  </div>
</template>

<style scoped>
.chart-box { position: relative; margin: 12px 0 4px; }
</style>

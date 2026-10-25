<script setup>
// 靶点立项评估（立项五问）：Open Targets / HPA / ChEMBL / ClinicalTrials.gov 等
// 9 个公开数据库证据筛查，红/黄/绿为阈值规则自动生成的筛查信号，不是结论。
import { nextTick, onMounted, onUnmounted, ref } from 'vue'
import { api } from '../api'
import JobCard from '../components/JobCard.vue'

const avail = ref(null)
const gene = ref(''), indication = ref('')
const job = ref(null), error = ref(''), busy = ref(false)
const report = ref(null)
const archive = ref([])
let timer, charts = []

onMounted(async () => {
  try { avail.value = await api.triageAvailability() } catch (_) {}
  refreshArchive()
})
onUnmounted(() => { clearInterval(timer); charts.forEach(c => c.destroy()); charts = [] })

function refreshArchive() {
  api.triageReports().then(a => archive.value = a).catch(() => {})
}

async function submit() {
  error.value = ''
  if (!gene.value.trim()) { error.value = '请输入靶点基因符号'; return }
  busy.value = true
  try {
    const fd = new FormData()
    fd.append('gene', gene.value.trim())
    fd.append('indication', indication.value.trim())
    job.value = await api.triageAssess(fd)
    report.value = null
    poll()
  } catch (e) { error.value = e.message } finally { busy.value = false }
}

function poll() {
  clearInterval(timer)
  timer = setInterval(async () => {
    if (!job.value) return
    try {
      job.value = await api.triageJob(job.value.id)
      if (['completed', 'failed'].includes(job.value.status)) {
        clearInterval(timer)
        if (job.value.status === 'completed' && job.value.stats?.report_id) {
          await loadReport(job.value.stats.report_id)
          refreshArchive()
        }
      }
    } catch (_) {}
  }, 2000)
}

async function loadReport(rid) {
  charts.forEach(c => c.destroy()); charts = []
  report.value = null
  try {
    report.value = await api.triageReport(rid)
    await nextTick()
    await ensureChartJs()
    renderCharts()
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

function renderCharts() {
  const Chart = window.Chart
  if (!Chart || !report.value) return
  const rc = document.getElementById('cv-radar')
  if (rc) charts.push(new Chart(rc, radarCfg(report.value)))
  for (const q of report.value.questions) {
    q.evidences.forEach((ev, i) => {
      const cfg = ev.chart
      if (!cfg) return
      if (cfg.type) mount(`cv-${q.key}-${i}`, cfg)
      else if (cfg.panels) cfg.panels.forEach((pc, j) => mount(`cv-${q.key}-${i}-p${j}`, pc))
    })
  }
}
function mount(canvasId, cfg) {
  const el = document.getElementById(canvasId)
  if (el) charts.push(new Chart(el, cfg))
}
function radarCfg(r) {
  const score = { green: 3, yellow: 2, red: 1, unknown: 1.5 }
  const color = { green: '#2e6b5e', yellow: '#d9b25f', red: '#c0504d', unknown: '#9aa5a1' }
  return {
    type: 'radar',
    data: {
      labels: r.questions.map(q => q.title),
      datasets: [{
        label: '五问状态', data: r.questions.map(q => score[q.status] ?? 1.5),
        backgroundColor: 'rgba(98,86,217,.18)', borderColor: '#6256d9', borderWidth: 2,
        pointBackgroundColor: r.questions.map(q => color[q.status] || '#9aa5a1'),
      }],
    },
    options: {
      responsive: true, maintainAspectRatio: false,
      plugins: { legend: { display: false } },
      scales: { r: { min: 0, max: 3, ticks: { stepSize: 1, backdropColor: 'transparent' },
        pointLabels: { font: { size: 11 } },
        grid: { color: 'var(--border)' }, angleLines: { color: 'var(--border)' } } },
    },
  }
}

const LABEL = { green: '证据充分', yellow: '需注意', red: '风险信号', unknown: '信息不足' }
const OVERALL = { green: '🟢 证据充分', yellow: '🟡 需注意', red: '🔴 风险信号', unknown: '⚪ 信息不足' }
function cid(qKey, i, panel) { return `cv-${qKey}-${i}${panel !== undefined ? '-p' + panel : ''}` }
function fmtDetail(d) {
  try { return JSON.stringify(d, null, 2).slice(0, 4000) } catch (_) { return String(d) }
}
function scrollToQ(key) {
  document.getElementById('q-' + key)?.scrollIntoView({ behavior: 'smooth', block: 'start' })
}
</script>

<template>
  <h1 class="page-title">靶点立项评估（立项五问）</h1>
  <p class="page-sub">
    输入靶点基因，自动调取 9 个公开数据库证据，回答立项前必须问的五个问题——
    靶点、组织、安全性、患者、商业。五个维度是乘法关系：要找的是短板，不是总分。
    红/黄/绿是阈值规则生成的筛查信号，最终判断权在评估者。
  </p>
  <div v-if="avail && !avail.available" class="alert err">引擎不可用：{{ avail.hint }}</div>
  <div v-if="avail && avail.available && !avail.oncokb_configured" class="muted" style="font-size:12.5px; margin:-6px 0 10px">
    💡 可选：申请 OncoKB 免费学术 token 写入 <code>数据根/oncokb_token.txt</code>（或环境变量 ONCOKB_TOKEN），
    重启后端即可启用肿瘤伴随诊断模块。
  </div>
  <div v-if="error" class="alert err">{{ error }}</div>

  <div class="card">
    <div class="grid2">
      <label class="field">靶点基因符号<input v-model="gene" placeholder="如 EGFR、STAT3" @keyup.enter="submit" /></label>
      <label class="field">适应症（可选）<input v-model="indication" placeholder="如 non-small cell lung carcinoma" @keyup.enter="submit" /></label>
    </div>
    <div style="margin-top:12px; display:flex; gap:10px; align-items:center">
      <button class="btn" :disabled="busy" @click="submit">{{ busy ? '评估中…' : '开始评估' }}</button>
      <span class="muted" style="font-size:12.5px">单次约 10–40 秒；首次会后台下载 HPA 组织数据（约 6MB）</span>
    </div>
  </div>

  <div v-if="job" style="margin-top:14px">
    <JobCard :job="job" @refresh="job = null" />
  </div>

  <!-- ============ 报告 ============ -->
  <template v-if="report">
    <div class="card" style="margin-top:14px">
      <div style="display:flex; gap:16px; flex-wrap:wrap; align-items:flex-start">
        <div style="flex:0 0 300px; max-width:100%">
          <h2 style="margin:0 0 2px">{{ report.gene_symbol }}
            <span v-if="report.gene_name && report.gene_name !== report.gene_symbol" class="muted" style="font-size:15px">{{ report.gene_name }}</span>
          </h2>
          <div class="muted" style="font-size:12.5px; margin-bottom:8px">
            Ensembl: {{ report.ensembl_id }} · UniProt: {{ report.uniprot_id || '—' }}<br>
            适应症: {{ report.indication || '不限' }} · {{ report.created_at }}
          </div>
          <span class="badge" :class="'b-' + report.overall_status" style="font-size:14px">总体筛查信号：{{ OVERALL[report.overall_status] }}</span>
          <div class="chart-box" style="height:230px; max-width:340px; margin-top:10px"><canvas id="cv-radar"></canvas></div>
        </div>
        <div style="flex:1; min-width:260px">
          <div class="flag-box risk">
            <b style="font-size:13px">🚩 风险信号（{{ report.risk_flags.length }}）——数据筛查提示，需人工判断，不等于结论</b>
            <ul v-if="report.risk_flags.length" style="margin:6px 0 0; padding-left:18px; font-size:13px">
              <li v-for="(v, i) in report.risk_flags" :key="i">{{ v }}</li>
            </ul>
            <div v-else class="muted" style="font-size:13px; margin-top:4px">无</div>
          </div>
          <div class="flag-box watch">
            <b style="font-size:13px">👀 建议进一步核查（{{ report.watch_flags.length }}）</b>
            <ul v-if="report.watch_flags.length" style="margin:6px 0 0; padding-left:18px; font-size:13px">
              <li v-for="(v, i) in report.watch_flags" :key="i">{{ v }}</li>
            </ul>
            <div v-else class="muted" style="font-size:13px; margin-top:4px">无</div>
          </div>
          <div class="muted" style="font-size:12px; margin-top:8px">
            红/黄/绿由 <code>engines/triage/rules.yaml</code> 中的阈值规则自动生成，仅为筛查参考。
          </div>
        </div>
      </div>
    </div>

    <div style="display:flex; gap:8px; flex-wrap:wrap; margin:14px 0">
      <button v-for="q in report.questions" :key="q.key" class="btn ghost sm qnav" :class="'qn-' + q.status"
              @click="scrollToQ(q.key)">
        {{ q.title }} <span class="badge" :class="'b-' + q.status" style="margin-left:6px">{{ LABEL[q.status] }}</span>
      </button>
    </div>

    <div v-for="(q, qi) in report.questions" :id="'q-' + q.key" :key="q.key" class="card qcard" style="margin-top:14px">
      <h3 style="margin-bottom:2px">{{ qi + 1 }}. {{ q.title }}
        <span class="badge" :class="'b-' + q.status" style="margin-left:8px">{{ LABEL[q.status] }}</span>
      </h3>
      <div class="muted" style="margin-bottom:12px">{{ q.question }}</div>

      <div v-for="(ev, i) in q.evidences" :key="i" class="ev">
        <div>
          <span class="badge" :class="'b-' + ev.status">{{ LABEL[ev.status] }}</span>
          <b style="margin-left:6px">{{ ev.source }}</b>
          <span class="muted" style="font-size:13px"> · {{ ev.check }}</span>
        </div>
        <div style="font-size:13.5px; margin:6px 0">{{ ev.summary }}</div>
        <div v-if="ev.chart" :style="ev.chart.panels ? 'display:grid; grid-template-columns:repeat(auto-fit,minmax(280px,1fr)); gap:12px' : ''">
          <template v-if="ev.chart.type">
            <div class="chart-box"><canvas :id="cid(q.key, i)"></canvas></div>
          </template>
          <template v-else-if="ev.chart.panels">
            <div v-for="(pc, j) in ev.chart.panels" :key="j" class="chart-box"><canvas :id="cid(q.key, i, j)"></canvas></div>
          </template>
        </div>
        <details v-if="ev.detail">
          <summary class="muted" style="cursor:pointer; font-size:12.5px">查看原始数据</summary>
          <pre class="mono" style="font-size:11.5px; white-space:pre-wrap; max-height:260px; overflow:auto; background:var(--bg); padding:8px; border-radius:8px">{{ fmtDetail(ev.detail) }}</pre>
        </details>
        <div v-if="ev.links?.length" style="margin-top:6px; font-size:12.5px">
          <a v-for="l in ev.links" :key="l.url" :href="l.url" target="_blank" rel="noopener" style="margin-right:14px">🔗 {{ l.label }}</a>
        </div>
      </div>

      <div v-if="q.manual_checklist?.length" class="manual">
        <b style="font-size:13px">📝 手动核查清单（无自动化数据源）</b>
        <div v-for="(ev, i) in q.manual_checklist" :key="i" style="font-size:13px; margin-top:6px">
          <b>{{ ev.source }} · {{ ev.check }}</b> — {{ ev.summary }}
          <a v-for="l in ev.links" :key="l.url" :href="l.url" target="_blank" rel="noopener" style="margin-left:8px">🔗 {{ l.label }}</a>
        </div>
      </div>
    </div>
    <div class="muted" style="font-size:12px; margin:10px 4px">
      评估结果由开源数据库自动生成，仅供研究参考，不构成任何立项、临床或投资建议。
    </div>
  </template>

  <!-- ============ 归档 ============ -->
  <div class="card" style="margin-top:16px">
    <h3>评估归档</h3>
    <div v-if="!archive.length" class="muted">还没有归档报告，完成一次评估后会出现在这里。</div>
    <table v-else class="tbl">
      <tr><th>基因</th><th>适应症</th><th>总体信号</th><th>风险/核查</th><th>时间</th><th></th></tr>
      <tr v-for="r in archive" :key="r.id">
        <td><b>{{ r.gene_symbol }}</b></td>
        <td class="muted">{{ r.indication || '不限' }}</td>
        <td><span class="badge" :class="'b-' + r.overall_status">{{ OVERALL[r.overall_status] }}</span></td>
        <td class="muted">{{ r.risk_count }} / {{ r.watch_count }}</td>
        <td class="muted">{{ r.created_at }}</td>
        <td style="white-space:nowrap">
          <button class="btn ghost sm" @click="loadReport(r.id)">查看</button>
          <a class="btn ghost sm" :href="`/api/triage/reports/${r.id}/download/md`">md</a>
          <a class="btn ghost sm" :href="`/api/triage/reports/${r.id}/download/xlsx`">xlsx</a>
        </td>
      </tr>
    </table>
  </div>
</template>

<style scoped>
.badge { display: inline-block; padding: 2px 9px; border-radius: 999px; font-size: 12px; font-weight: 600; vertical-align: middle; }
.b-green   { background: #e2f5ea; color: #1d6f43; }
.b-yellow  { background: #fdf3d7; color: #8a6d1a; }
.b-red     { background: #fde3e6; color: #a13040; }
.b-unknown { background: #ece9f1; color: #5f5968; }
.flag-box { border-radius: 10px; padding: 10px 12px; margin-bottom: 10px; }
.flag-box.risk  { background: #fdf0f2; border: 1px solid #f3cdd3; }
.flag-box.watch { background: #fdf9ec; border: 1px solid #ecddb0; }
.qcard { border-left: 4px solid var(--primary); scroll-margin-top: 70px; }
.qnav { border-left-width: 3px; }
.qn-green.qnav   { border-left-color: #2f9e63; }
.qn-yellow.qnav  { border-left-color: #d9b25f; }
.qn-red.qnav     { border-left-color: #d64b8a; }
.qn-unknown.qnav { border-left-color: #b9b3c2; }
.ev { border-top: 1px dashed var(--border); padding: 10px 0 12px; }
.ev:first-of-type { border-top: none; }
.chart-box { position: relative; height: 250px; margin: 8px 0 4px; }
.manual { background: var(--bg); border-radius: 10px; padding: 10px 12px; margin-top: 8px; }
</style>

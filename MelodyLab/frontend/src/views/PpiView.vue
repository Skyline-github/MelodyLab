<script setup>
import { onMounted, onUnmounted, ref } from 'vue'
import { api } from '../api'
import JobCard from '../components/JobCard.vue'

const identifiers = ref(''), species = ref(9606), requiredScore = ref(400)
const includePartners = ref(false), networkType = ref('functional')
const job = ref(null), error = ref(''), busy = ref(false), results = ref(null)
const hubSort = ref('degree')
let timer

onUnmounted(() => clearInterval(timer))

async function submit() {
  error.value = ''
  if (!identifiers.value.trim()) { error.value = '请输入蛋白列表'; return }
  busy.value = true; results.value = null
  try {
    const fd = new FormData()
    fd.append('identifiers', identifiers.value.trim())
    fd.append('species', species.value); fd.append('required_score', requiredScore.value)
    fd.append('include_partners', includePartners.value)
    fd.append('network_type', networkType.value)
    job.value = await api.ppiCreate(fd); poll()
  } catch (e) { error.value = e.message } finally { busy.value = false }
}
function poll() {
  clearInterval(timer)
  timer = setInterval(async () => {
    if (!job.value) return
    try {
      job.value = await api.job(job.value.id)
      if (job.value.status === 'completed') { clearInterval(timer); results.value = await api.ppiResults(job.value.id); draw() }
      if (job.value.status === 'failed') clearInterval(timer)
    } catch (_) {}
  }, 2000)
}

const nodes = ref([]), edges = ref([])
function draw() {
  nodes.value = results.value?.nodes || []
  edges.value = results.value?.edges || []
  requestAnimationFrame(() => forceLayout())
}
function forceLayout() {
  const canvas = document.getElementById('ppi-canvas')
  if (!canvas) return
  const ctx = canvas.getContext('2d')
  const W = canvas.width = canvas.parentElement.clientWidth
  const H = canvas.height = 480
  const N = nodes.value.map((n, i) => ({
    id: n.id,
    x: W / 2 + Math.cos(i * 2.4) * (60 + i * 4),
    y: H / 2 + Math.sin(i * 2.4) * (50 + i * 3),
    vx: 0, vy: 0,
    deg: results.value?.analysis?.centrality?.[n.id]?.degree || 0,
    comm: results.value?.analysis?.communities?.[n.id] ?? 0,
  }))
  const idx = Object.fromEntries(N.map((n, i) => [n.id, i]))
  const E = edges.value.map(e => [idx[e.source], idx[e.target]]).filter(p => p[0] !== undefined && p[1] !== undefined)
  const colors = ['#5b8def', '#ee6677', '#4cb37a', '#e6b23a', '#a678de', '#42c4d8', '#e08654']
  let tick = 0
  function step() {
    // 斥力
    for (let i = 0; i < N.length; i++) for (let j = i + 1; j < N.length; j++) {
      const a = N[i], b = N[j]
      let dx = a.x - b.x, dy = a.y - b.y
      let d2 = dx * dx + dy * dy + 0.01
      if (d2 > 60000) continue
      const f = 1400 / d2
      const d = Math.sqrt(d2)
      a.vx += dx / d * f; a.vy += dy / d * f
      b.vx -= dx / d * f; b.vy -= dy / d * f
    }
    // 弹簧
    for (const [i, j] of E) {
      const a = N[i], b = N[j]
      const dx = b.x - a.x, dy = b.y - a.y
      const d = Math.sqrt(dx * dx + dy * dy + 0.01)
      const f = (d - 110) * 0.012
      a.vx += dx / d * f; a.vy += dy / d * f
      b.vx -= dx / d * f; b.vy -= dy / d * f
    }
    // 中心引力 + 阻尼
    for (const n of N) {
      n.vx += (W / 2 - n.x) * 0.002; n.vy += (H / 2 - n.y) * 0.002
      n.vx *= 0.85; n.vy *= 0.85
      n.x = Math.max(24, Math.min(W - 24, n.x + n.vx))
      n.y = Math.max(24, Math.min(H - 24, n.y + n.vy))
    }
    ctx.clearRect(0, 0, W, H)
    ctx.strokeStyle = 'rgba(140,150,180,0.35)'
    for (const [i, j] of E) {
      ctx.beginPath(); ctx.moveTo(N[i].x, N[i].y); ctx.lineTo(N[j].x, N[j].y); ctx.stroke()
    }
    for (const n of N) {
      const r = 5 + Math.min(10, n.deg * 1.5)
      ctx.beginPath(); ctx.arc(n.x, n.y, r, 0, 6.284)
      ctx.fillStyle = colors[n.comm % colors.length]
      ctx.fill()
      ctx.fillStyle = '#dfe5f2'; ctx.font = '10px system-ui'
      ctx.fillText(n.id, n.x + r + 2, n.y + 3)
    }
    if (++tick < 300) requestAnimationFrame(step)
  }
  step()
}
const hubs = ref([])
function sortedHubs() {
  const cent = results.value?.analysis?.centrality || {}
  return Object.entries(cent).map(([gene, v]) => ({ gene, ...v }))
    .sort((a, b) => (b[hubSort.value] || 0) - (a[hubSort.value] || 0)).slice(0, 15)
}
</script>

<template>
  <h1 class="page-title">蛋白互作网络（STRING v12）</h1>
  <p class="page-sub">输入蛋白列表 → STRING 网络 + 中心性排名（Degree/Betweenness/MCC/k-core）+ 社区发现 + 富集</p>
  <div v-if="error" class="alert err">{{ error }}</div>

  <div class="card">
    <label class="field">蛋白列表（基因名，每行一个或逗号分隔）<textarea v-model="identifiers" rows="4" placeholder="EGFR&#10;AKT1&#10;TP53&#10;MYC&#10;KRAS"></textarea></label>
    <div class="grid3" style="margin-top:12px">
      <label class="field">物种 (NCBI taxid)<input v-model="species" type="number" /></label>
      <label class="field">最低分值 (0-1000)<input v-model="requiredScore" type="number" min="0" max="1000" /></label>
      <label class="field">网络类型
        <select v-model="networkType"><option value="functional">功能关联</option><option value="physical">物理互作</option></select></label>
    </div>
    <label class="field" style="display:flex; align-items:center; gap:8px; margin-top:10px">
      <input type="checkbox" v-model="includePartners" /> 扩展一阶伙伴
    </label>
    <div style="margin-top:14px"><button class="btn" :disabled="busy" @click="submit">{{ busy ? '查询中…' : '构建网络' }}</button></div>
  </div>

  <JobCard v-if="job && !results" :job="job" @refresh="job = null" />

  <div v-if="results" class="card">
    <h3>网络（{{ results.analysis?.stats?.n_nodes }} 节点 / {{ results.analysis?.stats?.n_edges }} 边 / {{ results.analysis?.stats?.n_communities }} 社区）</h3>
    <p v-if="results.degraded" class="alert warn">离线降级：{{ results.degraded }}</p>
    <div id="ppi-canvas-wrap" style="width:100%"><canvas id="ppi-canvas" style="width:100%"></canvas></div>
    <h3 style="margin-top:18px">Hub 基因（{{ hubSort }} 排序）</h3>
    <select v-model="hubSort" style="margin-bottom:8px">
      <option value="degree">Degree</option><option value="betweenness">Betweenness</option>
      <option value="mcc">MCC</option><option value="kcore">k-core</option><option value="closeness">Closeness</option>
    </select>
    <table style="width:100%; font-size:13px; border-collapse:collapse">
      <tr><th>gene</th><th>degree</th><th>betweenness</th><th>MCC</th><th>k-core</th><th>community</th></tr>
      <tr v-for="h in sortedHubs()" :key="h.gene">
        <td>{{ h.gene }}</td><td>{{ h.degree }}</td><td>{{ h.betweenness ?? '—' }}</td>
        <td>{{ h.mcc }}</td><td>{{ h.kcore }}</td><td>{{ (results.analysis?.communities || {})[h.gene] }}</td>
      </tr>
    </table>
    <div v-if="results.enrichment?.length" style="margin-top:16px">
      <h3>STRING 富集（前 10）</h3>
      <table style="width:100%; font-size:13px; border-collapse:collapse">
        <tr><th>term</th><th>description</th><th>FDR</th></tr>
        <tr v-for="e in results.enrichment.slice(0, 10)" :key="e.term">
          <td>{{ e.term }}</td><td>{{ e.description }}</td><td>{{ Number(e.fdr).toExponential(2) }}</td>
        </tr>
      </table>
    </div>
    <div style="margin-top:12px; display:flex; gap:8px">
      <a class="btn ghost sm" :href="api.artifactUrl(job.id, 'network.graphml')">⬇ GraphML</a>
      <a class="btn ghost sm" :href="api.artifactUrl(job.id, 'network.cytoscape.json')">⬇ Cytoscape JSON</a>
    </div>
  </div>
</template>

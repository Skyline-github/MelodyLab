<script setup>
import { onUnmounted, ref } from 'vue'
import { api } from '../api'
import JobCard from '../components/JobCard.vue'

const gse = ref(''), manualGroups = ref(''), useManual = ref(false)
const job = ref(null), error = ref(''), busy = ref(false)
const pngs = ref([]), topGenes = ref([]), stats = ref(null)
let timer

onUnmounted(() => clearInterval(timer))

async function submit() {
  error.value = ''
  if (!gse.value.trim()) { error.value = '请输入 GSE 编号'; return }
  busy.value = true; pngs.value = []; topGenes.value = []; stats.value = null
  try {
    const fd = new FormData()
    fd.append('gse', gse.value.trim())
    if (useManual.value && manualGroups.value.trim()) fd.append('manual_groups', manualGroups.value.trim())
    job.value = await api.geoAnalyze(fd); poll()
  } catch (e) { error.value = e.message } finally { busy.value = false }
}
function poll() {
  clearInterval(timer)
  timer = setInterval(async () => {
    if (!job.value) return
    try {
      job.value = await api.job(job.value.id)
      if (job.value.status === 'completed') {
        clearInterval(timer)
        stats.value = job.value.stats
        topGenes.value = job.value.stats?.top_genes || []
        pngs.value = (job.value.artifacts || []).filter(a => a.kind === 'png').map(a => a.name)
      }
      if (job.value.status === 'failed') clearInterval(timer)
    } catch (_) {}
  }, 2500)
}
</script>

<template>
  <h1 class="page-title">GEO 表达差异分析</h1>
  <p class="page-sub">GSE 数据集 → series matrix 下载 → 平台注释映射 → 自动分组 → QC（PCA/相关性）→ Welch t-test + BH-FDR → 火山图/MA/热图 + 报告</p>
  <div v-if="error" class="alert err">{{ error }}</div>

  <div class="card">
    <div class="grid2">
      <label class="field">GSE 编号<input v-model="gse" placeholder="GSE100428" /></label>
      <label class="field" style="display:flex; align-items:center; gap:8px"><input type="checkbox" v-model="useManual" /> 手动指定分组</label>
    </div>
    <label class="field" v-if="useManual" style="margin-top:12px">分组 JSON（{"control": ["GSM…"], "treatment": ["GSM…"]}）<textarea v-model="manualGroups" rows="3" placeholder='{"control": ["GSM2698477"], "treatment": ["GSM2698480"]}'></textarea></label>
    <div style="margin-top:16px"><button class="btn" :disabled="busy" @click="submit">{{ busy ? '分析中…' : '开始分析' }}</button></div>
  </div>

  <JobCard v-if="job && !stats" :job="job" @refresh="job = null" />

  <div v-if="stats" class="card">
    <h3>统计摘要</h3>
    <div class="grid3">
      <div class="pill">检验探针: {{ stats.n_probes_tested }}</div>
      <div class="pill">显著（FDR&lt;0.05）: {{ stats.n_sig_fdr05 }}</div>
      <div class="pill">显著且 |log2FC|&gt;1: {{ stats.n_sig_fdr05_fc1 }}</div>
    </div>
    <div style="display:grid; grid-template-columns:1fr 1fr; gap:12px; margin-top:14px">
      <img v-for="p in pngs" :key="p" :src="api.artifactUrl(job.id, p)" style="width:100%; border-radius:8px; border:1px solid #2a3142" />
    </div>
    <h3 style="margin-top:18px">Top 差异基因</h3>
    <table style="width:100%; font-size:13px; border-collapse:collapse">
      <tr><th>gene</th><th>log2FC</th><th>p_adj</th></tr>
      <tr v-for="g in topGenes.slice(0, 15)" :key="g.gene">
        <td>{{ g.gene }}</td><td>{{ g.log2fc }}</td><td>{{ Number(g.p_adj).toExponential(2) }}</td>
      </tr>
    </table>
    <div style="margin-top:12px; display:flex; gap:8px">
      <a class="btn ghost sm" :href="api.artifactUrl(job.id, 'geo_report.html')">HTML 报告</a>
      <a class="btn ghost sm" :href="api.artifactUrl(job.id, 'de_table.csv')">⬇ CSV</a>
      <a class="btn ghost sm" :href="api.artifactUrl(job.id, 'de_table.xlsx')">⬇ Excel</a>
    </div>
  </div>
</template>

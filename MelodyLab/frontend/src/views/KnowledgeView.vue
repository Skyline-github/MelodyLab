<script setup>
import { onUnmounted, ref } from 'vue'
import { api } from '../api'
import JobCard from '../components/JobCard.vue'

const tab = ref('activity')
const job = ref(null), error = ref(''), busy = ref(false)
const query = ref(''), queryType = ref('gene'), activityTypes = ref('IC50,Ki')
const resolution = ref(''), requireLigand = ref(false), method = ref('')
const smiles = ref(''), maxRecords = ref(50)
let timer

onUnmounted(() => clearInterval(timer))

async function submit() {
  error.value = ''
  if (!query.value.trim() && tab.value !== 'zinc') { error.value = '请输入检索词'; return }
  if (tab.value === 'zinc' && !smiles.value.trim()) { error.value = '请输入 SMILES'; return }
  busy.value = true
  try {
    const fd = new FormData()
    let res
    if (tab.value === 'activity') {
      fd.append('query', query.value.trim()); fd.append('query_type', queryType.value)
      fd.append('activity_types', activityTypes.value); fd.append('max_records', maxRecords.value)
      res = await api.knowledgeActivity(fd)
    } else if (tab.value === 'structures') {
      fd.append('query', query.value.trim())
      if (resolution.value) fd.append('resolution', resolution.value)
      fd.append('require_ligand', requireLigand.value)
      if (method.value) fd.append('method', method.value)
      fd.append('max_records', maxRecords.value)
      res = await api.knowledgeStructures(fd)
    } else if (tab.value === 'zinc') {
      fd.append('smiles', smiles.value.trim()); fd.append('max_records', maxRecords.value)
      res = await api.knowledgeZinc(fd)
    }
    job.value = res; poll()
  } catch (e) { error.value = e.message } finally { busy.value = false }
}
function poll() {
  clearInterval(timer)
  timer = setInterval(async () => {
    if (!job.value) return
    try { job.value = await api.job(job.value.id); if (['completed', 'failed'].includes(job.value.status)) clearInterval(timer) } catch (_) {}
  }, 2000)
}
</script>

<template>
  <h1 class="page-title">数据库发现</h1>
  <p class="page-sub">ChEMBL / PubChem / RCSB / AlphaFold / ZINC-22 官方公开 API 检索</p>

  <div class="tabs" style="margin-bottom:16px">
    <button class="btn ghost sm" :class="{ active: tab === 'activity' }" @click="tab = 'activity'">活性数据</button>
    <button class="btn ghost sm" :class="{ active: tab === 'structures' }" @click="tab = 'structures'">蛋白结构</button>
    <button class="btn ghost sm" :class="{ active: tab === 'zinc' }" @click="tab = 'zinc'">ZINC 相似分子</button>
  </div>
  <div v-if="error" class="alert err">{{ error }}</div>

  <div class="card">
    <div v-if="tab === 'activity'" class="grid3">
      <label class="field">检索词（基因名 / UniProt / ChEMBL ID）<input v-model="query" placeholder="EGFR / P00533 / CHEMBL203" /></label>
      <label class="field">类型
        <select v-model="queryType"><option value="gene">gene</option><option value="uniprot">uniprot</option><option value="inhibitor">inhibitor（分子名）</option></select></label>
      <label class="field">活性类型（逗号分隔）<input v-model="activityTypes" placeholder="IC50,Ki" /></label>
    </div>
    <div v-if="tab === 'structures'" class="grid3">
      <label class="field">基因名 / PDB ID / 关键词<input v-model="query" placeholder="EGFR" /></label>
      <label class="field">最大分辨率 (Å，留空不过滤)<input v-model="resolution" type="number" step="0.1" placeholder="2.5" /></label>
      <label class="field">实验方法
        <select v-model="method"><option value="">不限</option><option value="X-ray">X-ray</option><option value="EM">电镜</option><option value="NMR">NMR</option></select></label>
      <label class="field" style="display:flex; align-items:center; gap:8px"><input type="checkbox" v-model="requireLigand" /> 必须含配体</label>
    </div>
    <div v-if="tab === 'zinc'" class="grid2">
      <label class="field">参考 SMILES<input v-model="smiles" placeholder="CC(=O)Oc1ccccc1C(=O)O" /></label>
    </div>
    <div class="grid2" style="margin-top:12px">
      <label class="field">最大记录数<input v-model="maxRecords" type="number" min="5" max="200" /></label>
    </div>
    <div style="margin-top:16px"><button class="btn" :disabled="busy" @click="submit">{{ busy ? '检索中…' : '开始检索' }}</button></div>
  </div>

  <JobCard v-if="job" :job="job" @refresh="job = null" />
  <div v-if="job && job.status === 'completed'" class="card">
    <h3>检索结果</h3>
    <a class="btn ghost sm" :href="api.artifactUrl(job.id, 'structures.zip')" v-if="(job.artifacts || []).some(a => a.name === 'structures.zip')">⬇ 下载结构 ZIP</a>
    <p style="color:var(--muted); margin-top:8px">结果明细见任务 artifacts（{{ JSON.stringify(job.stats?.counts || {}) }}）</p>
  </div>
</template>

<script setup>
import { onMounted, onUnmounted, ref } from 'vue'
import { api } from '../api'
import JobCard from '../components/JobCard.vue'

const name = ref('')
const smiles = ref('')
const file = ref(null)
const policy = ref({ qed_min: 0.4, max_ro5_violations: 1, sa_max: 6.0, pains_mode: 'alert', brenk_mode: 'alert', selection_limit: null })
const jobId = ref(null)
const job = ref(null)
const error = ref('')
const busy = ref(false)
let timer

async function submit() {
  error.value = ''
  if (!smiles.value.trim() && !file.value?.files?.length) {
    error.value = '请输入 SMILES 或上传文件'
    return
  }
  busy.value = true
  try {
    const form = new FormData()
    form.append('name', name.value || '分子准入筛选')
    form.append('smiles_text', smiles.value)
    if (file.value?.files?.[0]) form.append('smiles_file', file.value.files[0])
    for (const [k, v] of Object.entries(policy.value)) if (v !== null && v !== '') form.append(k, v)
    const rec = await api.screening(form)
    jobId.value = rec.id
    poll()
  } catch (e) { error.value = e.message } finally { busy.value = false }
}

function poll() {
  clearInterval(timer)
  timer = setInterval(async () => {
    if (!jobId.value) return
    try { job.value = await api.job(jobId.value) } catch (_) {}
    if (job.value && !['running', 'queued', 'stopping'].includes(job.value.status)) clearInterval(timer)
  }, 1500)
}
onUnmounted(() => clearInterval(timer))
</script>

<template>
  <h1 class="page-title">分子准入筛选</h1>
  <p class="page-sub">对接前配体标准化与优先级筛选：Ro5 / QED / SA / PAINS / BRENK / 多样性选择</p>

  <div v-if="error" class="alert err">{{ error }}</div>
  <div class="card">
    <h3>输入分子</h3>
    <div class="grid2">
      <label class="field">筛选批次名称
        <input v-model="name" placeholder="例如：Enamine 片段库第一批" />
      </label>
      <label class="field">或上传文件（CSV 含 smiles 列 / 每行一个 SMILES 的 txt）
        <input type="file" ref="file" accept=".csv,.txt,.smi" />
      </label>
    </div>
    <label class="field" style="margin-top:12px">SMILES 文本（每行一个，可带空格分隔的 ID）
      <textarea v-model="smiles" placeholder="CC(=O)Oc1ccccc1C(=O)O aspirin&#10;c1ccc(cc1)S(=O)(=O)N"></textarea>
    </label>
    <details style="margin-top:12px">
      <summary class="muted" style="cursor:pointer">准入策略（高级）</summary>
      <div class="grid3" style="margin-top:10px">
        <label class="field">QED 下限<input v-model="policy.qed_min" type="number" step="0.05" /></label>
        <label class="field">Ro5 最大违规数<input v-model="policy.max_ro5_violations" type="number" min="0" max="4" /></label>
        <label class="field">SA 分数上限<input v-model="policy.sa_max" type="number" step="0.5" /></label>
        <label class="field">PAINS 模式
          <select v-model="policy.pains_mode"><option value="alert">警示</option><option value="reject">拒绝</option><option value="off">关闭</option></select>
        </label>
        <label class="field">BRENK 模式
          <select v-model="policy.brenk_mode"><option value="alert">警示</option><option value="reject">拒绝</option><option value="off">关闭</option></select>
        </label>
        <label class="field">多样性入选上限（可空）<input v-model="policy.selection_limit" type="number" placeholder="例如 500" /></label>
      </div>
    </details>
    <div style="margin-top:16px">
      <button class="btn" :disabled="busy" @click="submit">{{ busy ? '提交中…' : '开始筛选' }}</button>
    </div>
  </div>

  <JobCard v-if="job" :job="job" @refresh="job = null" />
</template>

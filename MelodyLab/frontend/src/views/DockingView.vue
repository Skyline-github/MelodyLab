<script setup>
import { onMounted, onUnmounted, ref } from 'vue'
import { api } from '../api'
import JobCard from '../components/JobCard.vue'

const form = ref({ name: '', engine: 'auto', max_pockets: 3, exhaustiveness: 8, top_n: 10, seed: 42 })
const ligandSmiles = ref('')
const ligandCsv = ref(null)
const ligandSdf = ref(null)
const proteinFiles = ref(null)
const jobId = ref(null)
const job = ref(null)
const results = ref(null)
const error = ref('')
const busy = ref(false)
const env = ref(null)
let timer

onMounted(async () => { env.value = await api.environment().catch(() => null) })
onUnmounted(() => clearInterval(timer))

async function submit() {
  error.value = ''
  const proteins = proteinFiles.value?.files
  if (!proteins?.length) { error.value = '请选择蛋白结构文件（可多选）'; return }
  busy.value = true
  try {
    const fd = new FormData()
    fd.append('name', form.value.name || '批量对接')
    fd.append('engine', form.value.engine)
    fd.append('max_pockets', form.value.max_pockets)
    fd.append('exhaustiveness', form.value.exhaustiveness)
    fd.append('top_n', form.value.top_n)
    fd.append('seed', form.value.seed)
    fd.append('ligand_smiles', ligandSmiles.value)
    if (ligandCsv.value?.files?.[0]) fd.append('ligand_csv', ligandCsv.value.files[0])
    for (const f of ligandSdf.value?.files || []) fd.append('ligand_files', f)
    for (const f of proteins) fd.append('protein_files', f)
    const rec = await api.docking(fd)
    jobId.value = rec.id
    poll()
  } catch (e) { error.value = e.message } finally { busy.value = false }
}

function poll() {
  clearInterval(timer)
  timer = setInterval(async () => {
    if (!jobId.value) return
    try {
      job.value = await api.job(jobId.value)
      if (['completed', 'failed', 'stopped', 'cancelled'].includes(job.value.status)) {
        clearInterval(timer)
        const r = await api.dockingResults(jobId.value)
        results.value = r
      }
    } catch (_) {}
  }, 2000)
}
</script>

<template>
  <h1 class="page-title">批量分子对接</h1>
  <p class="page-sub">多蛋白 × 多配体 · 几何口袋检测 · Vina-GPU / Vina / Smina 自动回退 · 断点续跑</p>

  <div v-if="env && !(env.docking_modes || []).length" class="alert err">
    未检测到可用对接引擎。请将 vina 可执行文件放到 assets/tools/（详见 README「资产准备」）。
  </div>
  <div v-if="error" class="alert err">{{ error }}</div>

  <div class="card">
    <h3>蛋白结构（PDB / ENT，可多选）</h3>
    <input type="file" ref="proteinFiles" multiple accept=".pdb,.ent" />
    <h3 style="margin-top:18px">配体来源（三选一及以上）</h3>
    <div class="grid3">
      <label class="field">SMILES 文本
        <textarea v-model="ligandSmiles" placeholder="CCO ethanol&#10;CC(=O)Oc1ccccc1C(=O)O"></textarea>
      </label>
      <label class="field">CSV / TXT 文件<input type="file" ref="ligandCsv" accept=".csv,.txt,.smi" style="margin-top:8px" /></label>
      <label class="field">SDF 文件（可多选）<input type="file" ref="ligandSdf" multiple accept=".sdf" style="margin-top:8px" /></label>
    </div>
    <details style="margin-top:12px">
      <summary class="muted" style="cursor:pointer">对接参数</summary>
      <div class="grid3" style="margin-top:10px">
        <label class="field">运行名称<input v-model="form.name" placeholder="批量对接" /></label>
        <label class="field">对接引擎
          <select v-model="form.engine">
            <option value="auto">自动（GPU 优先）</option>
            <option value="vina_cpu">AutoDock Vina (CPU)</option>
            <option value="smina">Smina (CPU)</option>
            <option value="vina_gpu">Vina-GPU (CUDA)</option>
          </select>
        </label>
        <label class="field">每蛋白口袋数上限<input v-model="form.max_pockets" type="number" min="1" max="10" /></label>
        <label class="field">Exhaustiveness<input v-model="form.exhaustiveness" type="number" min="1" max="64" /></label>
        <label class="field">TOP-N 打包<input v-model="form.top_n" type="number" min="1" /></label>
        <label class="field">随机种子<input v-model="form.seed" type="number" /></label>
      </div>
    </details>
    <div style="margin-top:16px">
      <button class="btn" :disabled="busy" @click="submit">{{ busy ? '提交中…' : '启动对接' }}</button>
    </div>
  </div>

  <JobCard v-if="job" :job="job" @refresh="job = null" />

  <div v-if="results" class="card">
    <h3>对接结果（{{ results.total }} 条，按打分升序）</h3>
    <div style="max-height:420px; overflow:auto">
      <table class="tbl">
        <tr><th>蛋白</th><th>配体</th><th>口袋</th><th>引擎</th><th>打分 kcal/mol</th></tr>
        <tr v-for="r in results.results.slice(0, 100)" :key="r.task_id">
          <td>{{ r.protein }}</td><td>{{ r.ligand }}</td><td>{{ r.pocket }}</td><td>{{ r.engine || '—' }}</td>
          <td><b :style="{color: r.score != null && r.score < -7 ? 'var(--accent)' : 'inherit'}">{{ r.score ?? '—' }}</b></td>
        </tr>
      </table>
    </div>
  </div>
</template>

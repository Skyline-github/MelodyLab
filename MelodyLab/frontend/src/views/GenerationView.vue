<script setup>
// 3D 分子生成（FlashDiff）+ Linker 生成（DiffPROTACs）共用 3Dmol 可视化由任务产物 SDF 提供
import { onMounted, onUnmounted, ref } from 'vue'
import { api } from '../api'
import JobCard from '../components/JobCard.vue'

const avail = ref(null)
const form = ref({ num_molecules: 5, batch_size: 10, pocket_radius: 15, num_steps: 50, seed: 2024, mean_atoms: 28, std_atoms: 5 })
const protein = ref(null), ligand = ref(null)
const job = ref(null), error = ref(''), busy = ref(false)
let timer

onMounted(async () => { avail.value = await api.generationAvailability().catch(() => null) })
onUnmounted(() => clearInterval(timer))

async function submit() {
  error.value = ''
  if (!protein.value?.files?.[0] || !ligand.value?.files?.[0]) { error.value = '请上传蛋白 PDB 和参考配体 SDF'; return }
  busy.value = true
  try {
    const fd = new FormData()
    for (const [k, v] of Object.entries(form.value)) fd.append(k, v)
    fd.append('protein', protein.value.files[0])
    fd.append('ligand', ligand.value.files[0])
    job.value = await api.generation(fd)
    poll()
  } catch (e) { error.value = e.message } finally { busy.value = false }
}
function poll() {
  clearInterval(timer)
  timer = setInterval(async () => {
    if (!job.value) return
    try { job.value = await api.job(job.value.id) } catch (_) {}
  }, 2500)
}
</script>

<template>
  <h1 class="page-title">3D 分子生成</h1>
  <p class="page-sub">FlashDiff 结构导向扩散模型：在口袋内从头生成新分子（含生成轨迹）</p>

  <div v-if="avail && !avail.available" class="alert warn">
    FlashDiff 资产未就绪：{{ avail.hint }}（权重约 243MB，首次准备需要一点时间）
  </div>
  <div v-if="error" class="alert err">{{ error }}</div>

  <div class="card">
    <h3>输入</h3>
    <div class="grid2">
      <label class="field">蛋白结构 PDB<input type="file" ref="protein" accept=".pdb" style="margin-top:8px" /></label>
      <label class="field">参考配体 SDF（口袋中心参考）<input type="file" ref="ligand" accept=".sdf" style="margin-top:8px" /></label>
    </div>
    <details style="margin-top:12px">
      <summary class="muted" style="cursor:pointer">生成参数</summary>
      <div class="grid3" style="margin-top:10px">
        <label class="field">生成数量<input v-model="form.num_molecules" type="number" min="1" max="50" /></label>
        <label class="field">批大小<input v-model="form.batch_size" type="number" min="1" /></label>
        <label class="field">口袋半径 (Å)<input v-model="form.pocket_radius" type="number" step="0.5" /></label>
        <label class="field">扩散步数<input v-model="form.num_steps" type="number" min="10" max="1000" /></label>
        <label class="field">目标原子数均值<input v-model="form.mean_atoms" type="number" /></label>
        <label class="field">原子数标准差<input v-model="form.std_atoms" type="number" /></label>
      </div>
    </details>
    <div style="margin-top:16px">
      <button class="btn" :disabled="busy || (avail && !avail.available)" @click="submit">
        {{ busy ? '提交中…' : '开始生成' }}
      </button>
      <span class="muted" style="font-size:13px; margin-left:10px">CPU 推理较慢，建议少量试跑</span>
    </div>
  </div>

  <JobCard v-if="job" :job="job" @refresh="job = null" />
</template>

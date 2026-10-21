<script setup>
import { onMounted, onUnmounted, ref } from 'vue'
import { api } from '../api'
import JobCard from '../components/JobCard.vue'

const avail = ref(null)
const linkerLen = ref(12), nSamples = ref(5)
const lig1 = ref(null), lig2 = ref(null)
const job = ref(null), error = ref(''), busy = ref(false)
let timer

onMounted(async () => { avail.value = await api.linkerAvailability().catch(() => null) })
onUnmounted(() => clearInterval(timer))

async function submit() {
  error.value = ''
  if (!lig1.value?.files?.[0] || !lig2.value?.files?.[0]) { error.value = '请上传两个片段 SDF'; return }
  busy.value = true
  try {
    const fd = new FormData()
    fd.append('linker_len', linkerLen.value)
    fd.append('n_samples', nSamples.value)
    fd.append('ligand1', lig1.value.files[0])
    fd.append('ligand2', lig2.value.files[0])
    job.value = await api.linker(fd)
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
  <h1 class="page-title">Linker 生成（PROTAC 设计）</h1>
  <p class="page-sub">DiffPROTACs 扩散模型：在 E3 配体与 Warhead 之间生成 3D 连接子</p>

  <div v-if="avail && !avail.available" class="alert warn">DiffPROTACs 资产未就绪：{{ avail.hint }}</div>
  <div v-if="error" class="alert err">{{ error }}</div>

  <div class="card">
    <h3>输入片段</h3>
    <div class="grid2">
      <label class="field">片段 1 · E3 配体 SDF<input type="file" ref="lig1" accept=".sdf" style="margin-top:8px" /></label>
      <label class="field">片段 2 · Warhead SDF<input type="file" ref="lig2" accept=".sdf" style="margin-top:8px" /></label>
    </div>
    <div class="grid3" style="margin-top:12px">
      <label class="field">Linker 原子数<input v-model="linkerLen" type="number" min="2" max="30" /></label>
      <label class="field">生成样本数<input v-model="nSamples" type="number" min="1" max="50" /></label>
    </div>
    <div style="margin-top:16px">
      <button class="btn" :disabled="busy || (avail && !avail.available)" @click="submit">{{ busy ? '提交中…' : '开始生成' }}</button>
    </div>
  </div>

  <JobCard v-if="job" :job="job" @refresh="job = null" />
</template>

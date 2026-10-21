<script setup>
import { computed, onMounted, onUnmounted, ref } from 'vue'
import { api } from '../api'
import JobCard from '../components/JobCard.vue'

const avail = ref(null), job = ref(null), error = ref(''), busy = ref(false)
const runId = ref(''), taskId = ref('')
const protein = ref(null), ligand = ref(null), ligandSmiles = ref('')
const chains = ref(''), totalNs = ref(1.0), testNs = ref(0.2), nptNs = ref(0.1)
const forcefield = ref('amber99sb'), water = ref('tip3p'), salt = ref(0.15)
const rmsdThreshold = ref(0.3), useGpu = ref(true), ligandCharge = ref('')
const report = ref(null), showTraj = ref(false)
const recentJobs = ref([])
let timer

onMounted(async () => {
  avail.value = await api.mdAvailability().catch(() => null)
  try { recentJobs.value = (await api.mdJobs()).jobs.slice(0, 6) } catch (_) {}
})

function loadJob(j) {
  clearInterval(timer)
  error.value = ''; report.value = null; showTraj.value = false
  job.value = j
  if (j.status === 'completed') {
    api.mdResults(j.id).then(r => { report.value = r.report }).catch(() => {})
  } else if (['running', 'queued', 'paused', 'stopping'].includes(j.status)) {
    poll()
  }
}
onUnmounted(() => clearInterval(timer))

async function submit() {
  error.value = ''
  if (!runId.value && !protein.value?.files?.[0]) { error.value = '请引用对接任务或上传蛋白 PDB'; return }
  busy.value = true
  try {
    const fd = new FormData()
    if (runId.value) { fd.append('run_id', runId.value); fd.append('task_id', taskId.value) }
    if (protein.value?.files?.[0]) fd.append('protein_file', protein.value.files[0])
    if (ligand.value?.files?.[0]) fd.append('ligand_file', ligand.value.files[0])
    if (ligandSmiles.value.trim()) fd.append('ligand_smiles', ligandSmiles.value.trim())
    fd.append('chains', chains.value)
    fd.append('total_ns', totalNs.value); fd.append('test_ns', testNs.value); fd.append('npt_ns', nptNs.value)
    fd.append('forcefield', forcefield.value); fd.append('water', water.value)
    fd.append('salt_molar', salt.value); fd.append('rmsd_threshold', rmsdThreshold.value)
    fd.append('use_gpu', useGpu.value)
    if (ligandCharge.value !== '') fd.append('ligand_charge', ligandCharge.value)
    job.value = await api.mdJob(fd)
    poll()
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
        const r = await api.mdResults(job.value.id)
        report.value = r.report
      }
      if (job.value.status === 'failed') clearInterval(timer)
    } catch (_) {}
  }, 3000)
}

const stats = computed(() => job.value?.stats || {})
const pngs = computed(() => (job.value?.artifacts || []).filter(a => a.kind === 'png'))
function artUrl(name) { return api.artifactUrl(job.value.id, name) }

const trajLoading = ref(false)

async function playTraj() {
  if (trajLoading.value) return
  error.value = ''
  showTraj.value = true
  trajLoading.value = true
  await new Promise(r => setTimeout(r, 50))
  try {
    if (!window.$3Dmol) {
      await new Promise((resolve, reject) => {
        const s = document.createElement('script')
        s.src = '/3Dmol-min.js'
        s.onload = resolve
        s.onerror = () => reject(new Error('3Dmol.js 加载失败（本地文件缺失）'))
        document.head.appendChild(s)
      })
    }
    const el = document.getElementById('traj-viewer')
    const viewer = window.$3Dmol.createViewer(el, { backgroundColor: '#0f1420' })
    const url = api.mdTrajectoryPreview(job.value.id)
    const resp = await fetch(url)
    if (!resp.ok) throw new Error(`轨迹文件请求失败（HTTP ${resp.status}）`)
    const text = await resp.text()
    viewer.addModelsAsFrames(text, 'pdb')
    if (!viewer.getNumFrames || viewer.getNumFrames() <= 1) {
      throw new Error('轨迹只解析到 1 帧——预览文件可能是旧版格式，后端会在下次请求时自动重新生成')
    }
    viewer.setStyle({}, { cartoon: { color: 'spectrum' } })
    viewer.addStyle({ resn: ['LIG'] }, { stick: { radius: 0.3, colorscheme: 'greenCarbon' } })
    viewer.animate({ loop: 'backandforth', reps: 0, speed: 0.6 })
    viewer.zoomTo()
    viewer.resize()
    viewer.render()
  } catch (e) {
    error.value = (e && e.message) || '轨迹加载失败'
    showTraj.value = false
  } finally {
    trajLoading.value = false
  }
}
</script>

<template>
  <h1 class="page-title">分子动力学验证（GROMACS）</h1>
  <p class="page-sub">对接复合物 → 选链/修蛋白 → 刚体配体参数化 → EM → NVT 预检门控 → NPT → 生产 → RMSD/RMSF/Rg/能量</p>

  <div v-if="avail && !avail.available" class="alert warn">GROMACS 不可用：{{ avail.hint }}</div>
  <div v-else-if="avail" class="alert ok" style="margin-bottom:16px">
    引擎就绪：{{ avail.path }}　GPU: {{ avail.gpu ? '✅ 可用 (-nb gpu)' : 'CPU 模式' }}　参数化: {{ avail.parameterization }}
  </div>
  <div v-if="error" class="alert err">{{ error }}</div>

  <div class="card">
    <h3>复合物来源（二选一）</h3>
    <div class="grid2">
      <label class="field">对接 Run ID（如 docking_20261017_…）<input v-model="runId" placeholder="对接任务 ID" /></label>
      <label class="field">子任务 ID（protein__ligand__pocket_1）<input v-model="taskId" placeholder="任务详情页可复制" /></label>
    </div>
    <div class="grid2" style="margin-top:12px">
      <label class="field">或直传蛋白 PDB<input type="file" ref="protein" accept=".pdb,.ent" style="margin-top:8px" /></label>
      <label class="field">配体 SDF（可选，默认按 SMILES 现场生成）<input type="file" ref="ligand" accept=".sdf" style="margin-top:8px" /></label>
    </div>
    <label class="field" style="margin-top:12px">或配体 SMILES<input v-model="ligandSmiles" placeholder="CC(=O)OC1=CC=CC=C1C(=O)O" /></label>
    <h3 style="margin-top:16px">模拟参数</h3>
    <div class="grid3">
      <label class="field">链（逗号分隔，留空=全部）<input v-model="chains" placeholder="A,B" /></label>
      <label class="field">预检 NVT (ns)<input v-model="testNs" type="number" step="0.05" min="0.01" /></label>
      <label class="field">NPT 平衡 (ns)<input v-model="nptNs" type="number" step="0.05" min="0" /></label>
      <label class="field">生产时长 (ns)<input v-model="totalNs" type="number" step="0.5" min="0.01" /></label>
      <label class="field">门控阈值 RMSD (nm)<input v-model="rmsdThreshold" type="number" step="0.05" min="0.05" /></label>
      <label class="field">盐浓度 (M)<input v-model="salt" type="number" step="0.05" min="0" /></label>
      <label class="field">力场
        <select v-model="forcefield">
          <option>amber99sb</option><option>amber99sb-ildn</option><option>charmm27</option>
        </select></label>
      <label class="field">水模型
        <select v-model="water"><option>tip3p</option><option>tip4p</option><option>spce</option></select></label>
      <label class="field">GPU 加速
        <select v-model="useGpu"><option :value="true">开（-nb gpu）</option><option :value="false">关</option></select></label>
    </div>
    <div style="margin-top:16px">
      <button class="btn" :disabled="busy || (avail && !avail.available)" @click="submit">{{ busy ? '提交中…' : '启动 MD 流水线' }}</button>
    </div>
  </div>

  <div v-if="recentJobs.length" class="card" style="margin-top:16px">
    <h3>最近任务（点击载入结果，刷新页面后可从这里找回）</h3>
    <div style="display:flex; gap:8px; flex-wrap:wrap; margin-top:8px">
      <button v-for="j in recentJobs" :key="j.id" class="btn ghost sm"
              :class="{ active: job && job.id === j.id }" @click="loadJob(j)">
        {{ j.id }}　{{ j.status === 'completed' ? '✅' : j.status === 'failed' ? '❌' : '⏳' }} {{ (j.message || '').slice(0, 30) }}
      </button>
    </div>
  </div>

  <JobCard v-if="job" :job="job" @refresh="job = null" />

  <div v-if="job && job.status === 'completed'" class="card">
    <h3>分析结果</h3>
    <div class="grid3" v-if="stats.backbone_rmsd_tail_nm !== undefined">
      <div class="pill">骨架 RMSD（末段均值）: {{ stats.backbone_rmsd_tail_nm }} nm</div>
      <div class="pill">配体 RMSD: {{ stats.ligand_rmsd_tail_nm }} nm</div>
      <div class="pill">门控 RMSD: {{ stats.gate_rmsd_nm }} nm ｜ {{ stats.engine?.toUpperCase() }}</div>
    </div>
    <div style="display:grid; grid-template-columns:1fr 1fr; gap:12px; margin-top:14px">
      <img v-for="p in pngs" :key="p.name" :src="artUrl(p.name)" style="width:100%; border-radius:8px; border:1px solid #2a3142" />
    </div>
    <div style="margin-top:14px">
      <button class="btn ghost" :disabled="trajLoading" @click="playTraj">{{ trajLoading ? '⏳ 轨迹加载中…（约 12MB，首次需数秒）' : '▶ 轨迹预览（3Dmol）' }}</button>
    </div>
    <div v-if="showTraj" id="traj-viewer" style="position:relative; width:100%; height:420px; margin-top:12px; border-radius:8px; border:1px solid #2a3142"></div>
  </div>
</template>

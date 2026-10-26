<script setup>
// 三元复合物预测（DeepTernary）：PROTAC → 参照晶体装配 → 三元复合物 + BSA/clash
// 结构健康评估。BSA 只有相对比较意义，配合装配质量指标一起解读。
import { nextTick, onMounted, onUnmounted, ref } from 'vue'
import { api } from '../api'
import JobCard from '../components/JobCard.vue'

const avail = ref(null)
const name = ref(''), smilesText = ref(''), defaultE3 = ref(''), seeds = ref(2)
const job = ref(null), error = ref(''), busy = ref(false), busyExample = ref(false)
const report = ref(null), archive = ref([])
const viewCid = ref(''), viewerLoading = ref(false)
const mdJob = ref(null), mdBusyCid = ref('')
let timer, mdTimer

onMounted(async () => {
  try { avail.value = await api.ternaryAvailability() } catch (_) {}
  refreshArchive()
})
onUnmounted(() => { clearInterval(timer); clearInterval(mdTimer) })

function refreshArchive() {
  api.ternaryReports().then(a => archive.value = a).catch(() => {})
}

async function loadExample() {
  busyExample.value = true
  try {
    const ex = await api.degraderExample()   // 复用 degrader 的 Galectin-3 示例
    smilesText.value = ex.text
    if (!name.value.trim()) name.value = 'Galectin-3 PROTACs（示例）'
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
    fd.append('default_e3', defaultE3.value)
    fd.append('seeds', String(seeds.value))
    job.value = await api.ternaryPredict(fd)
    report.value = null
    poll()
  } catch (e) { error.value = e.message } finally { busy.value = false }
}

function poll() {
  clearInterval(timer)
  timer = setInterval(async () => {
    if (!job.value) return
    try {
      job.value = await api.ternaryJob(job.value.id)
      if (['completed', 'failed'].includes(job.value.status)) {
        clearInterval(timer)
        if (job.value.status === 'completed' && job.value.stats?.report_id) {
          await loadReport(job.value.stats.report_id)
          refreshArchive()
        }
      }
    } catch (_) {}
  }, 2500)
}

async function loadReport(rid) {
  report.value = null
  viewCid.value = ''
  mdJob.value = null
  try {
    report.value = await api.ternaryReport(rid)
    window.scrollTo({ top: 0, behavior: 'smooth' })
  } catch (e) { error.value = e.message }
}

function launchMd(row) {
  error.value = ''
  mdBusyCid.value = row.compound_id
  const fd = new FormData()
  fd.append('ternary_report_id', report.value.id)
  fd.append('ternary_cid', row.compound_id)
  fd.append('total_ns', '1')
  fd.append('test_ns', '1')
  fd.append('npt_ns', '1')
  fd.append('rmsd_threshold', '0.4')
  api.mdJob(fd).then(res => {
    mdJob.value = res
    clearInterval(mdTimer)
    mdTimer = setInterval(async () => {
      try {
        mdJob.value = await api.mdJob(res.id)
        if (['completed', 'failed', 'stopped'].includes(mdJob.value.status)) {
          clearInterval(mdTimer)
          mdBusyCid.value = ''
        }
      } catch (_) {}
    }, 5000)
  }).catch(e => { error.value = e.message }).finally(() => { mdBusyCid.value = '' })
}

async function show3d(row) {
  viewCid.value = row.compound_id
  viewerLoading.value = true
  await nextTick()
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
    const el = document.getElementById('ternary-viewer')
    el.innerHTML = ''
    const viewer = window.$3Dmol.createViewer(el, { backgroundColor: '#0f1420' })
    const text = await (await fetch(`/api/ternary/reports/${report.value.id}/complex/${row.compound_id}`)).text()
    viewer.addModel(text, 'pdb')
    viewer.setStyle({ chain: 'A' }, { cartoon: { color: '#5b8dd9' } })   // POI
    viewer.setStyle({ chain: 'B' }, { cartoon: { color: '#e0995b' } })   // E3
    viewer.addStyle({ resn: 'UNL' }, { stick: { radius: 0.28, colorscheme: 'greenCarbon' } })
    viewer.zoomTo()
    viewer.render()
  } catch (e) { error.value = e.message } finally { viewerLoading.value = false }
}
</script>

<template>
  <h1 class="page-title">三元复合物预测（DeepTernary）</h1>
  <p class="page-sub">
    输入 PROTAC，用参照晶体把两个 warhead 锚定到 POI（当前内置 Galectin-3C）与 E3
    （CRBN / VHL）的口袋姿态上，装配初始复合物后由 DeepTernary 预测完整三元复合物，
    并给出 BSA（蛋白-蛋白埋藏面积）与 clash 结构健康指标。
    <b>BSA 只有相对比较意义</b>，请结合装配质量（RMSD / 最小间距 / clash）共同解读。
  </p>
  <div v-if="avail && !avail.available" class="alert err">
    引擎不可用：{{ avail.hint }}
    <div v-if="avail.poi_reference" class="muted" style="margin-top:6px; font-size:12.5px">
      POI 参照：{{ avail.poi_reference }} · E3：{{ avail.e3_supported.join(' / ') }}
    </div>
  </div>
  <div v-if="error" class="alert err">{{ error }}</div>

  <div class="card">
    <div class="grid2">
      <label class="field">任务名称（可选）<input v-model="name" placeholder="如 Galectin-3 × CRBN 系列" /></label>
      <label class="field">查询文件（可选，CSV，需含 smiles 列）<input type="file" accept=".csv,.txt,.smi"
        @change="(e) => { const f = e.target.files[0]; if (f) { const rd = new FileReader(); rd.onload = () => smilesText = rd.result; rd.readAsText(f) } }" /></label>
    </div>
    <label class="field" style="margin-top:10px">SMILES 列表（每行 "ID SMILES [E3]"，E3 可省略用默认值）
      <textarea v-model="smilesText" rows="6" placeholder="K112 O=C(...) CRBN&#10;或每行一个 SMILES"></textarea>
    </label>
    <div class="grid2" style="margin-top:10px">
      <label class="field">默认 E3（行内未指定时使用）
        <select v-model="defaultE3">
          <option value="">不设置（行内必须指定）</option>
          <option value="CRBN">CRBN</option>
          <option value="VHL">VHL</option>
        </select>
      </label>
      <label class="field">推理 seed 数（越大越稳，越慢）
        <select v-model.number="seeds">
          <option :value="1">1（最快）</option>
          <option :value="2">2（推荐）</option>
          <option :value="4">4（研究项目同款）</option>
        </select>
      </label>
    </div>
    <div style="margin-top:12px; display:flex; gap:10px; align-items:center; flex-wrap:wrap">
      <button class="btn" :disabled="busy" @click="submit">{{ busy ? '预测中…' : '开始预测' }}</button>
      <button class="btn ghost" :disabled="busyExample" @click="loadExample">{{ busyExample ? '载入中…' : '📋 载入示例（28 个 Galectin-3 PROTAC）' }}</button>
      <span class="muted" style="font-size:12.5px">单分子约 20–40 秒（制备 + CPU 推理）；首个任务需加载模型</span>
    </div>
  </div>

  <div v-if="job" style="margin-top:14px">
    <JobCard :job="job" @refresh="job = null" />
  </div>

  <div v-if="mdJob" style="margin-top:14px">
    <JobCard :job="mdJob" @refresh="mdJob = null" />
    <div class="muted" style="font-size:12.5px; margin-top:4px">
      三元复合物 MD 在「🌀 分子动力学」页可查看完整分析（含三体接触曲线：LIG-POI / LIG-E3 / POI-E3）。
    </div>
  </div>

  <!-- ============ 报告 ============ -->
  <template v-if="report">
    <div class="card" style="margin-top:14px">
      <div style="display:flex; gap:10px; flex-wrap:wrap; align-items:baseline">
        <h2 style="margin:0">{{ report.name || report.id }}</h2>
        <span class="muted" style="font-size:12.5px">{{ report.created_at }} · seeds={{ report.seeds }}</span>
      </div>
      <div class="cards3" style="margin-top:12px">
        <div class="stat"><b>{{ report.n_ok }} / {{ report.n_query }}</b><span>成功预测 / 查询分子</span></div>
        <div class="stat"><b>{{ (report.rows.filter(r => r.status === 'ok').map(r => r.bsa) || [0]).length ? Math.round(Math.min(...report.rows.filter(r => r.status === 'ok').map(r => r.bsa || 1e9))) + '–' + Math.round(Math.max(...report.rows.filter(r => r.status === 'ok').map(r => r.bsa || 0))) : '—' }}</b><span>BSA 范围（Å²，仅相对比较）</span></div>
        <div class="stat"><b>{{ report.rows.filter(r => (r.clash || 0) > 0.2).length }}</b><span>clash &gt; 0.2（建议人工检查）</span></div>
      </div>

      <table class="tbl" style="margin-top:14px">
        <tr>
          <th>化合物</th><th>E3</th><th>状态</th><th title="姿态质量代理指标，越小越好">rmsd</th>
          <th title="重叠原子占比">clash</th><th>BSA (Å²)</th>
          <th title="MCS 原子数 POI/E3 · 对齐 RMSD · 配体构象来源">装配质量</th><th>提示</th><th></th>
        </tr>
        <tr v-for="r in report.rows" :key="r.compound_id">
          <td><b>{{ r.compound_id }}</b></td>
          <td class="muted">{{ r.e3 || '—' }}</td>
          <td v-if="r.status === 'ok'">✅</td>
          <td v-else style="color:var(--warn)" colspan="7">{{ r.status === 'prep_failed' ? '制备失败' : '推理失败' }}：{{ r.reason }}</td>
          <template v-if="r.status === 'ok'">
            <td>{{ r.pred_p2_rmsd.toFixed(2) }}</td>
            <td :style="(r.clash || 0) > 0.2 ? 'color:var(--warn)' : ''">{{ r.clash.toFixed(3) }}</td>
            <td><b>{{ Math.round(r.bsa) }}</b></td>
            <td class="muted" style="font-size:12px">{{ r.prep.mcs_poi }}/{{ r.prep.mcs_e3 }} atm ·
              {{ r.prep.poi_rmsd }}/{{ r.prep.e3_rmsd }} Å · {{ r.prep.lig_kind }}</td>
            <td class="muted" style="font-size:12px; max-width:220px">{{ r.prep.prep_warning || '—' }}</td>
            <td><button class="btn ghost sm" @click="show3d(r)">🔬 3D 查看</button>
              <button class="btn ghost sm" :disabled="mdBusyCid === r.compound_id" title="对该三元复合物跑 1 ns 动力学并分析三体接触"
                      @click="launchMd(r)">{{ mdBusyCid === r.compound_id ? '启动中…' : '🌀 MD 验证' }}</button></td>
          </template>
        </tr>
      </table>

      <div v-if="viewCid" style="margin-top:14px">
        <div style="display:flex; align-items:center; gap:10px">
          <b>{{ viewCid }} 三元复合物</b>
          <span class="pill">蓝 = POI (chain A)</span><span class="pill">橙 = E3 (chain B)</span>
          <span class="pill">绿 = PROTAC</span>
          <span v-if="viewerLoading" class="muted">加载中…</span>
        </div>
        <div id="ternary-viewer" style="position:relative; width:100%; height:460px; margin-top:8px; border-radius:8px; border:1px solid var(--border)"></div>
        <div class="muted" style="font-size:12px; margin-top:4px">若两个蛋白穿插或 PROTAC 悬空，说明装配/预测几何欠佳，该分子结果应降权。</div>
      </div>

      <div class="card" style="margin-top:14px; background:var(--bg)">
        <h3 style="margin-top:0">方法与可信度边界</h3>
        <ul style="margin:8px 0 0; padding-left:20px; font-size:13px; line-height:1.7">
          <li v-for="(n, i) in report.notes" :key="i">{{ n }}</li>
        </ul>
      </div>
    </div>
  </template>

  <!-- ============ 归档 ============ -->
  <div class="card" style="margin-top:16px">
    <h3>预测归档</h3>
    <div v-if="!archive.length" class="muted">还没有归档报告，完成一次预测后会出现在这里。</div>
    <table v-else class="tbl">
      <tr><th>任务</th><th>成功/查询</th><th>前 5</th><th>时间</th><th></th></tr>
      <tr v-for="r in archive" :key="r.id">
        <td><b>{{ r.name || r.id }}</b></td>
        <td class="muted">{{ r.n_ok }} / {{ r.n_query }}</td>
        <td class="muted">{{ (r.top || []).join('、') }}</td>
        <td class="muted">{{ r.created_at }}</td>
        <td style="white-space:nowrap">
          <button class="btn ghost sm" @click="loadReport(r.id)">查看</button>
          <a class="btn ghost sm" :href="`/api/ternary/reports/${r.id}/download/csv`">csv</a>
          <a class="btn ghost sm" :href="`/api/ternary/reports/${r.id}/download/md`">md</a>
        </td>
      </tr>
    </table>
  </div>
</template>

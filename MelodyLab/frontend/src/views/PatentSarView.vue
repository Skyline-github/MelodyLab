<script setup>
// 专利 SAR 挖掘：专利 PDF（文本层/扫描件）→ 结构/活性数据集 → RGA/3D-QSAR → 叙事总结。
// 输出仅为研究参考数据集，结构经 OCR/枚举，引用前请对照 check.html 人工复核。
import { computed, onMounted, onUnmounted, ref } from 'vue'
import { api } from '../api'
import JobCard from '../components/JobCard.vue'

const avail = ref(null), example = ref(null)
const name = ref(''), useExample = ref(false), pdfFile = ref(null)
const pathType = ref('auto'), target = ref('TARGET')
const withSar = ref(true), withSummary = ref(true), showAdv = ref(false)
const structurePages = ref(''), assayPages = ref(''), mzPages = ref('')
const markushConfig = ref(''), markushPages = ref('')
const kinds = ref('example'), core = ref(''), activityCol = ref('')
const zoom = ref(3.0), textZoomMult = ref(2.0), mzTol = ref(1.5), vlm = ref('auto')
const job = ref(null), error = ref(''), busy = ref(false), selftest = ref(null)
const summaryMd = ref('')
let timer

onMounted(async () => {
  try { avail.value = await api.patentsarAvailability() } catch (_) {}
  try { example.value = await api.patentsarExample() } catch (_) {}
})
onUnmounted(() => clearInterval(timer))

const llmOn = computed(() => !!avail.value?.llm?.configured)
const pngArtifacts = computed(() => (job.value?.artifacts || []).filter(a => a.name.endsWith('.png')))

function loadExampleParams() {
  useExample.value = true
  if (!example.value) return
  const p = example.value.params || {}
  pathType.value = p.path_type || 'born_digital'
  structurePages.value = p.structure_pages || ''
  assayPages.value = p.assay_pages || ''
  mzPages.value = p.mz_pages || ''
  markushPages.value = p.markush_pages || ''
  if (!target.value || target.value === 'TARGET') target.value = p.target || 'TARGET'
  if (!name.value.trim()) name.value = example.value.name
}

async function submit() {
  error.value = ''
  if (!useExample.value && !pdfFile.value) { error.value = '请上传专利 PDF，或勾选使用内置示例'; return }
  busy.value = true
  try {
    const fd = new FormData()
    fd.append('name', name.value.trim())
    fd.append('example', useExample.value)
    fd.append('path_type', pathType.value)
    fd.append('target', target.value.trim() || 'TARGET')
    fd.append('with_sar', withSar.value)
    fd.append('with_summary', withSummary.value)
    fd.append('structure_pages', structurePages.value.trim())
    fd.append('assay_pages', assayPages.value.trim())
    fd.append('mz_pages', mzPages.value.trim())
    fd.append('markush_config', markushConfig.value)
    fd.append('markush_pages', markushPages.value.trim())
    fd.append('kinds', kinds.value.trim())
    fd.append('core', core.value.trim())
    fd.append('activity_col', activityCol.value.trim())
    fd.append('zoom', String(zoom.value))
    fd.append('text_zoom_mult', String(textZoomMult.value))
    fd.append('mz_tol', String(mzTol.value))
    fd.append('vlm', vlm.value)
    if (!useExample.value && pdfFile.value) fd.append('pdf', pdfFile.value)
    job.value = await api.patentsarMine(fd)
    poll()
  } catch (e) { error.value = e.message } finally { busy.value = false }
}

function poll() {
  clearInterval(timer)
  timer = setInterval(async () => {
    if (!job.value) return
    try {
      job.value = await api.patentsarJob(job.value.id)
      if (['completed', 'failed'].includes(job.value.status)) {
        clearInterval(timer)
        if (job.value.status === 'completed') {
          summaryMd.value = ''
          const hasSummary = (job.value.artifacts || []).some(a => a.name === 'sar/sar_summary.md')
          if (hasSummary) {
            try {
              const res = await fetch(api.artifactUrl(job.value.id, 'sar/sar_summary.md'))
              summaryMd.value = await res.text()
            } catch (_) {}
          }
        }
      }
    } catch (_) {}
  }, 2000)
}

async function runSelftest() {
  selftest.value = null
  try { selftest.value = await api.patentsarSelftest() } catch (e) { error.value = e.message }
}

// 迷你 markdown 渲染（标题/粗体/列表/表格/删除线/行内码，覆盖 sar_summary.md 子集）
function mdToHtml(md) {
  const esc = (s) => s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
  const inline = (s) => esc(s)
    .replace(/\*\*([^*]+)\*\*/g, '<b>$1</b>')
    .replace(/~~([^~]+)~~/g, '<s>$1</s>')
    .replace(/`([^`]+)`/g, '<code>$1</code>')
    .replace(/\*([^*]+)\*/g, '<i>$1</i>')
  let html = '', inTable = false, inList = false
  for (const raw of (md || '').split('\n')) {
    const line = raw.trimEnd()
    if (line.startsWith('|')) {
      const cells = line.split('|').slice(1, -1).map(c => c.trim())
      if (cells.length && cells.every(c => c === '' || /^:?-{2,}:?$/.test(c))) continue
      if (!inTable) { html += '<table class="tbl md-table"><tr>' + cells.map(c => `<th>${inline(c)}</th>`).join('') + '</tr>'; inTable = true }
      else html += '<tr>' + cells.map(c => `<td>${inline(c)}</td>`).join('') + '</tr>'
      continue
    }
    if (inTable) { html += '</table>'; inTable = false }
    const h = line.match(/^(#{1,4}) (.*)$/)
    if (h) { html += `<h${h[1].length + 1}>${inline(h[2])}</h${h[1].length + 1}>`; continue }
    if (/^[-*] /.test(line)) {
      if (!inList) { html += '<ul>'; inList = true }
      html += `<li>${inline(line.slice(2))}</li>`; continue
    }
    if (inList) { html += '</ul>'; inList = false }
    if (line === '---') { html += '<hr>'; continue }
    if (line.trim()) html += `<p>${inline(line)}</p>`
  }
  if (inTable) html += '</table>'
  if (inList) html += '</ul>'
  return html
}
</script>

<template>
  <h1 class="page-title">专利 SAR 挖掘</h1>
  <p class="page-sub">
    上传药物专利 PDF（有文本层或纯扫描件均可），自动提取<b>结构 + 活性数据集</b>，做
    R 基团分析 / 3D-QSAR，并生成带引用校验的叙事 SAR 总结。Markush R 表定义的实施例
    通过显式通式模板枚举。输出为<b>研究参考数据集</b>：结构经 OCR/枚举，引用前请对照
    check.html 人工复核；仅覆盖实施例层面，不含权利要求 Markush 空间。
  </p>
  <div v-if="avail && !avail.available" class="alert err">引擎不可用：{{ avail.hint }}</div>
  <div v-if="error" class="alert err">{{ error }}</div>

  <div class="card">
    <div class="grid2">
      <label class="field">任务名称（可选）<input v-model="name" placeholder="如 WO2025XXXXXX AMY3R 系列" /></label>
      <label class="field">专利 PDF
        <input type="file" accept=".pdf" :disabled="useExample"
          @change="(e) => pdfFile = e.target.files[0] || null" /></label>
    </div>
    <label class="field" style="margin-top:10px">靶点名称（用于活性列识别与报告标题）
      <input v-model="target" /></label>
    <div style="margin-top:10px; display:flex; gap:18px; flex-wrap:wrap; align-items:center">
      <label><input type="checkbox" v-model="withSar" /> 挖掘后自动跑 SAR 分析（RGA/3D-QSAR）</label>
      <label><input type="checkbox" v-model="withSummary" /> 生成叙事 SAR 总结</label>
    </div>
    <div style="margin-top:12px; display:flex; gap:10px; align-items:center; flex-wrap:wrap">
      <button class="btn" :disabled="busy" @click="submit">{{ busy ? '挖掘中…' : '开始挖掘' }}</button>
      <button class="btn ghost" @click="loadExampleParams">📋 使用内置示例（合成硬基准专利）</button>
      <span v-if="useExample" class="pill">示例模式：忽略上传文件</span>
      <button class="btn ghost sm" @click="showAdv = !showAdv">{{ showAdv ? '收起高级选项' : '高级选项（页码/Markush/VLM）' }}</button>
    </div>

    <div v-if="showAdv" style="margin-top:12px; border-top:1px dashed rgba(128,140,170,.4); padding-top:10px">
      <div class="grid2">
        <label class="field">路径
          <select v-model="pathType">
            <option value="auto">自动探测（有文本层 → born-digital）</option>
            <option value="born_digital">born-digital（有文本层）</option>
            <option value="scanned">扫描件（OCR + DECIMER）</option>
          </select></label>
        <label class="field">结构页范围（如 1:5，扫描件同）<input v-model="structurePages" placeholder="留空=全部页" /></label>
        <label class="field">活性表页（born: "6"；扫描: "6:amy3r+ctr"）<input v-model="assayPages" /></label>
        <label class="field">LCMS 文本页（如 7；born 与扫描通用）<input v-model="mzPages" /></label>
        <label class="field">Markush R 表页（如 5）<input v-model="markushPages" /></label>
        <label class="field">结构类型过滤（扫描件；逗号分隔 example,intermediate）<input v-model="kinds" /></label>
        <label class="field">OCR 倍率（扫描结构页，建议 ≥3）<input type="number" step="0.5" v-model.number="zoom" /></label>
        <label class="field">文本页倍率乘子（表格防粘列，默认 2）<input type="number" step="0.5" v-model.number="textZoomMult" /></label>
      </div>
      <label class="field" style="margin-top:8px">Markush 配置 JSON（core_template + positions；留空且选了示例时用内置配置）
        <textarea v-model="markushConfig" rows="5" placeholder='{"core_template": "{r1}NC(=O)c1ccc2{r2}ncc2c1", "positions": {"r1": {"Me": "C", "Et": "CC"}, "r2": {"H": "[nH]", "Me": "[n](C)"}}}'></textarea></label>
      <div class="grid2" style="margin-top:8px">
        <label class="field">母核 SMILES（可选，指定则不用 FMCS 自动找核）<input v-model="core" /></label>
        <label class="field">活性列名（默认自动探测第一个 ec50* 列）<input v-model="activityCol" /></label>
      </div>
      <div style="margin-top:8px; display:flex; gap:18px; flex-wrap:wrap; align-items:center">
        <label class="field">VLM 标签消歧
          <select v-model="vlm">
            <option value="auto">auto（配置了 API key 时启用）</option>
            <option value="off">off（仅几何启发式）</option>
          </select></label>
        <label class="field">m/z 容差 (Da)<input type="number" step="0.5" v-model.number="mzTol" style="width:90px" /></label>
      </div>
    </div>
  </div>

  <div v-if="job" style="margin-top:14px">
    <JobCard :job="job" @refresh="job = null" />
  </div>

  <!-- ============ 结果 ============ -->
  <template v-if="job && job.status === 'completed'">
    <div class="card" style="margin-top:14px">
      <div style="display:flex; gap:10px; flex-wrap:wrap; align-items:center">
        <h2 style="margin:0">结果 · {{ job.name }}</h2>
        <span class="pill">{{ job.stats?.path_type }}</span>
        <span v-if="job.stats?.rga_positions" class="pill">RGA: {{ (job.stats.rga_positions || []).join(' / ') }}</span>
        <span style="flex:1"></span>
        <a v-for="a in (job.artifacts || []).filter(x => ['check.html', 'sar_report.md', 'sar_summary.md', 'sar_summary.json'].includes(x.name) || x.name.endsWith('_check.html') || x.name.endsWith('_sar_report.md') || x.name.endsWith('_sar_summary.md'))"
           :key="a.name" class="btn ghost sm" :href="api.artifactUrl(job.id, a.name)" target="_blank">{{ a.name }}</a>
      </div>
      <div v-if="job.stats?.summary_headline" class="stat" style="margin-top:12px">
        <b>{{ job.stats.summary_headline }}</b><span>总结总览（每条结论均经引用校验）</span>
      </div>
      <!-- 叙事总结 -->
      <div v-if="summaryMd" style="margin-top:10px" v-html="mdToHtml(summaryMd)"></div>
      <div v-else class="muted" style="margin-top:10px; font-size:12.5px">本次未生成叙事总结。</div>
    </div>

    <div class="card" style="margin-top:14px">
      <h3>SAR 图表</h3>
      <div v-if="!pngArtifacts.length" class="muted">本次未生成图表（可在高级选项中开启 SAR 分析）。</div>
      <div v-else style="display:grid; grid-template-columns:repeat(auto-fill,minmax(340px,1fr)); gap:12px">
        <figure v-for="a in pngArtifacts" :key="a.name" style="margin:0; text-align:center">
          <img :src="api.artifactUrl(job.id, a.name)" style="max-width:100%; border:1px solid rgba(128,140,170,.3); border-radius:6px" />
          <figcaption class="muted" style="font-size:12px; margin-top:4px">{{ a.name }}</figcaption>
        </figure>
      </div>
    </div>

    <div class="card" style="margin-top:14px">
      <h3>产物下载</h3>
      <table class="tbl">
        <tr><th>文件</th><th></th></tr>
        <tr v-for="a in (job.artifacts || [])" :key="a.path">
          <td><code>{{ a.name }}</code></td>
          <td><a class="btn ghost sm" :href="api.artifactUrl(job.id, a.name)" target="_blank">打开/下载</a></td>
        </tr>
      </table>
      <div class="muted" style="margin-top:8px; font-size:12px">
        dataset.csv = 结构化数据集；_check.html = 人工复核表（异常行高亮）；
        sar/ 下为 RGA 图表 + 3D-QSAR + 叙事总结（markdown/JSON）。
      </div>
    </div>
  </template>

  <!-- ============ 状态与自检 ============ -->
  <div v-if="avail" class="card" style="margin-top:14px">
    <div style="display:flex; gap:10px; flex-wrap:wrap; align-items:center">
      <h3 style="margin:0">引擎状态</h3>
      <span class="pill">{{ avail.available ? '✅ 可用' : '❌ 不可用' }}</span>
      <span class="pill">born-digital {{ avail.paths?.born_digital ? '✅' : '❌' }}</span>
      <span class="pill">扫描件 {{ avail.paths?.scanned ? '✅' : '❌' }}</span>
      <span class="pill">{{ llmOn ? `LLM 已配置（${avail.llm.model}）` : 'LLM 未配置 → 总结走 digest 模式，VLM 消歧关闭' }}</span>
      <span style="flex:1"></span>
      <button class="btn ghost sm" @click="runSelftest">运行离线自检</button>
    </div>
    <div v-if="selftest" class="muted" style="font-size:12.5px; margin-top:8px">
      自检：{{ selftest.ok ? '✅ 全部通过' : '❌ 存在失败' }} ·
      VLM 解析用例失败 {{ selftest.vlm_parse_fails }} · 引用校验门用例失败 {{ selftest.citation_gate_fails }}
    </div>
    <div v-if="avail.missing_dependencies?.length" class="muted" style="font-size:12.5px; margin-top:6px">
      缺失依赖：{{ avail.missing_dependencies.join('、') }}
      （pip install pymupdf decimer opencv-python easyocr；DECIMER 权重缓存见 $PYSTOW_HOME）
    </div>
  </div>
</template>

<style scoped>
.md-table { margin: 8px 0; }
.md-table code { font-size: 12px; }
</style>

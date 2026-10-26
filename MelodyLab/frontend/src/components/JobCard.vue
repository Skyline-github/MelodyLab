<script setup>
// 可复用：任务进度条 + 控制按钮
import { computed } from 'vue'
import { api } from '../api'

const props = defineProps({ job: { type: Object, required: true } })
const emit = defineEmits(['refresh'])

const running = computed(() => ['running', 'queued', 'stopping'].includes(props.job.status))

async function act(fn) {
  try { await fn(); emit('refresh') } catch (e) { alert(e.message) }
}
</script>

<template>
  <div class="card">
    <div style="display:flex; align-items:center; gap:10px; flex-wrap:wrap">
      <b>{{ job.name }}</b>
      <span class="status" :class="job.status">{{ job.status }}</span>
      <span class="pill">{{ job.kind }}</span>
      <span class="muted" style="font-size:12.5px">{{ job.created_at }}</span>
      <span style="flex:1"></span>
      <button v-if="running" class="btn ghost sm" @click="act(() => api.pauseJob(job.id))">暂停</button>
      <button v-if="job.status === 'paused'" class="btn sm" @click="act(() => job.kind === 'triage' ? api.triageResume(job.id) : job.kind === 'degrader' ? api.degraderResume(job.id) : job.kind === 'md' ? api.mdResume(job.id) : api.resumeRun(job.id))">恢复</button>
      <button v-if="running || job.status === 'paused'" class="btn danger sm" @click="act(() => api.stopJob(job.id))">停止</button>
      <button class="btn danger sm" @click="act(() => api.deleteJob(job.id))">删除</button>
    </div>
    <div class="progress" style="margin:12px 0 8px"><div :style="{ width: job.progress + '%' }"></div></div>
    <div class="muted" style="font-size:13px">
      {{ job.phase || '—' }} · {{ job.message || '' }}
      <template v-if="job.error"> · <span style="color:var(--warn)">{{ job.error.slice(0, 120) }}</span></template>
    </div>
    <div v-if="job.artifacts?.length" style="margin-top:10px; display:flex; gap:8px; flex-wrap:wrap">
      <a v-for="a in job.artifacts" :key="a.path" class="btn ghost sm" :href="api.artifactUrl(job.id, a.name)">
        ⬇ {{ a.name }}
      </a>
    </div>
    <details v-if="job.kind === 'docking'" style="margin-top:8px">
      <summary class="muted" style="cursor:pointer; font-size:13px">运行详情</summary>
      <pre class="mono muted" style="font-size:12px; white-space:pre-wrap">{{ JSON.stringify(job.stats, null, 2) }}</pre>
    </details>
  </div>
</template>

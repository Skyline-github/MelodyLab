<script setup>
// 任务中心：全部任务列表 + 轮询刷新
import { onMounted, onUnmounted, ref } from 'vue'
import { api } from '../api'
import JobCard from '../components/JobCard.vue'

const jobs = ref([])
const kind = ref('')
let timer

async function refresh() {
  try { jobs.value = (await api.jobs(kind.value || undefined)).jobs } catch (_) {}
}
onMounted(() => { refresh(); timer = setInterval(refresh, 3000) })
onUnmounted(() => clearInterval(timer))
</script>

<template>
  <h1 class="page-title">任务中心</h1>
  <p class="page-sub">所有计算任务的统一视图 · job.json 断点恢复 · 结果下载</p>

  <div style="display:flex; gap:8px; margin-bottom:14px">
    <button v-for="k in ['', 'screening', 'docking', 'generation', 'linker']" :key="k"
      class="btn sm" :class="kind === k ? '' : 'ghost'" @click="kind = k; refresh()">
      {{ k === '' ? '全部' : ({ screening: '筛选', docking: '对接', generation: '生成', linker: 'Linker' })[k] }}
    </button>
    <span style="flex:1"></span>
    <button class="btn ghost sm" @click="refresh">刷新</button>
  </div>

  <div v-if="!jobs.length" class="card muted">暂无任务</div>
  <JobCard v-for="j in jobs" :key="j.id" :job="j" @refresh="refresh" />
</template>

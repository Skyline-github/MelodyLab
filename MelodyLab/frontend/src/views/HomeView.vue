<script setup>
// 门户首页：环境状态 + 工作流概览 + 最近任务
import { onMounted, onUnmounted, ref } from 'vue'
import { api } from '../api'

const env = ref(null)
const jobs = ref([])
let timer

onMounted(async () => {
  await refresh()
  timer = setInterval(refresh, 5000)
})
onUnmounted(() => clearInterval(timer))

async function refresh() {
  try {
    env.value = await api.environment()
    jobs.value = (await api.jobs()).jobs.slice(0, 5)
  } catch (_) {}
}

const steps = [
  { title: '靶点立项评估', desc: '立项五问证据筛查 · 9 大公开数据库 · 红/黄/绿风险信号', to: '/triage' },
  { title: '分子准入筛选', desc: '标准化 · Ro5/QED/SA · PAINS 警报 · 多样性选择', to: '/screening' },
  { title: '批量分子对接', desc: 'GPU 优先引擎 · 口袋检测 · 断点续跑', to: '/docking' },
  { title: '3D 分子生成', desc: 'FlashDiff 扩散模型口袋内从头生成', to: '/generation' },
  { title: 'Linker 生成', desc: 'DiffPROTACs 双片段连接子设计', to: '/linker' },
]
</script>

<template>
  <h1 class="page-title">MelodyLab</h1>
  <p class="page-sub">开源 AI 药物发现工作台 · 从化合物库到结构证据的一体化流水线</p>

  <div class="cards3" v-if="env">
    <div class="stat"><b>{{ env.status === 'ready' ? '✅ 就绪' : '⚠️ 未就绪' }}</b><span>运行环境</span></div>
    <div class="stat"><b>{{ env.gpu?.available ? env.gpu.name || 'NVIDIA GPU' : 'CPU' }}</b><span>计算设备</span></div>
    <div class="stat"><b>{{ (env.docking_modes || []).join(' / ') || '无可用引擎' }}</b><span>对接引擎</span></div>
  </div>
  <div v-if="env && env.status !== 'ready'" class="alert warn">
    尚未检测到可用对接引擎。请在 assets/tools/ 放置 vina 可执行文件（或运行 scripts/setup_assets.py），然后重启后端。
  </div>

  <div class="grid2" style="margin-bottom:18px">
    <router-link v-for="s in steps" :key="s.title" :to="s.to" class="card" style="text-decoration:none;color:inherit">
      <h3>{{ s.title }}</h3>
      <div class="muted" style="font-size:13.5px">{{ s.desc }}</div>
    </router-link>
  </div>

  <div class="card">
    <h3>最近任务</h3>
    <div v-if="!jobs.length" class="muted">还没有任务，从上面的功能卡片开始。</div>
    <table v-else class="tbl">
      <tr><th>任务</th><th>类型</th><th>状态</th><th>进度</th></tr>
      <tr v-for="j in jobs" :key="j.id">
        <td>{{ j.name }}</td><td>{{ j.kind }}</td>
        <td><span class="status" :class="j.status">{{ j.status }}</span></td>
        <td style="min-width:160px"><div class="progress"><div :style="{width: j.progress+'%'}"></div></div></td>
      </tr>
    </table>
    <router-link to="/jobs" class="muted" style="font-size:13px">查看全部 →</router-link>
  </div>
</template>

<script setup>
import { onMounted, ref } from 'vue'
import { api } from './api'

const dark = ref(localStorage.getItem('ml-theme') === 'dark')

function toggleTheme() {
  dark.value = !dark.value
  document.documentElement.classList.toggle('dark', dark.value)
  localStorage.setItem('ml-theme', dark.value ? 'dark' : 'light')
}
onMounted(() => document.documentElement.classList.toggle('dark', dark.value))

const env = ref(null)
onMounted(async () => { try { env.value = await api.environment() } catch (_) {} })
</script>

<template>
  <div class="layout">
    <aside class="sidebar">
      <div class="logo"><span class="mark"></span> MelodyLab</div>
      <router-link class="nav-item" to="/">🏠 门户首页</router-link>
      <router-link class="nav-item" to="/triage">🧭 靶点立项评估</router-link>
      <router-link class="nav-item" to="/screening">🧹 分子准入筛选</router-link>
      <router-link class="nav-item" to="/docking">🎯 批量分子对接</router-link>
      <router-link class="nav-item" to="/generation">✨ 3D 分子生成</router-link>
      <router-link class="nav-item" to="/linker">🔗 Linker 生成</router-link>
      <router-link class="nav-item" to="/degrader">🧬 降解剂评估</router-link>
      <router-link class="nav-item" to="/ternary">🔬 三元复合物预测</router-link>
      <div class="nav-sep" style="margin:10px 0; border-top:1px solid rgba(128,140,170,0.25)"></div>
      <router-link class="nav-item" to="/md">🌀 分子动力学</router-link>
      <router-link class="nav-item" to="/knowledge">🗄 数据库发现</router-link>
      <router-link class="nav-item" to="/ppi">🕸 蛋白互作网络</router-link>
      <router-link class="nav-item" to="/geo">📈 GEO 分析</router-link>
      <router-link class="nav-item" to="/retro">🧪 逆合成规划</router-link>
      <router-link class="nav-item" to="/patentsar">📑 专利 SAR 挖掘</router-link>
      <router-link class="nav-item" to="/settings">⚙ 设置</router-link>
      <router-link class="nav-item" to="/jobs">📋 任务中心</router-link>
      <div class="side-footer">
        <button class="btn ghost sm" @click="toggleTheme">{{ dark ? '☀️ 亮色' : '🌙 暗色' }}</button>
        <div style="margin-top:8px">开源 AI 药物发现工作台<br />v0.2.0</div>
      </div>
    </aside>
    <main class="main">
      <router-view />
    </main>
  </div>
</template>

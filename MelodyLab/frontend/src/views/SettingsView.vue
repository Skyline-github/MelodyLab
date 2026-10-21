<script setup>
import { onMounted, ref } from 'vue'
import { api } from '../api'

const settings = ref(null), newPath = ref(''), error = ref(''), ok = ref(''), busy = ref(false)

onMounted(load)
async function load() {
  try { settings.value = await api.getSettings() } catch (e) { error.value = e.message }
}
async function save() {
  error.value = ''; ok.value = ''
  if (!newPath.value.trim()) { error.value = '请输入数据保存路径'; return }
  busy.value = true
  try {
    const r = await api.setDataRoot({ path: newPath.value.trim() })
    ok.value = `已保存：${r.data_root}（重启后端后生效）`
    load()
  } catch (e) { error.value = e.message } finally { busy.value = false }
}
</script>

<template>
  <h1 class="page-title">设置</h1>
  <p class="page-sub">数据保存位置、引擎与运行信息</p>

  <div v-if="error" class="alert err">{{ error }}</div>
  <div v-if="ok" class="alert ok">{{ ok }}</div>

  <div class="card">
    <h3>数据保存路径</h3>
    <p style="color:var(--muted)">任务数据、分析报告与结果文件（jobs/outputs）的保存根目录。</p>
    <label class="field">当前路径
      <input :value="settings?.data_root" readonly />
    </label>
    <p style="color:var(--muted); font-size:12px">
      来源：{{ settings?.source === 'env' ? '环境变量 MELODYLAB_DATA_ROOT（优先级最高）'
              : settings?.source === 'config' ? '配置文件 data_root.txt'
              : '默认（用户主目录）'}}
      <template v-if="settings?.marker_file">｜配置文件：{{ settings?.marker_file }}</template>
    </p>
    <label class="field" style="margin-top:12px">新路径<input v-model="newPath" placeholder="例如 D:\MelodyLabData" /></label>
    <div style="margin-top:14px"><button class="btn" :disabled="busy" @click="save">{{ busy ? '保存中…' : '保存并应用' }}</button></div>
  </div>

  <div class="card">
    <h3>说明</h3>
    <ul style="line-height:1.9">
      <li>目录结构：<code>jobs/</code>（全部任务的数据与产物）、<code>logs/</code>（运行日志）</li>
      <li>修改路径后<b>需重启后端</b>才能生效；已有任务数据不会自动迁移，可手动整目录拷贝</li>
      <li>MD / GEO / 对接等产物（报告、轨迹、CSV/Excel、ZIP）均保存在该目录对应任务文件夹下</li>
    </ul>
  </div>
</template>

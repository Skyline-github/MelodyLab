<script setup>
import { ref } from 'vue'
import { api } from '../api'

const smiles = ref(''), engine = ref('hybrid'), maxRoutes = ref(5), maxDepth = ref(2)
const result = ref(null), error = ref(''), busy = ref(false)

async function submit() {
  error.value = ''; result.value = null
  if (!smiles.value.trim()) { error.value = '请输入目标分子 SMILES'; return }
  busy.value = true
  try {
    const fd = new FormData()
    fd.append('smiles', smiles.value.trim()); fd.append('engine', engine.value)
    fd.append('max_routes', maxRoutes.value); fd.append('max_depth', maxDepth.value)
    result.value = await api.retroPlan(fd)
  } catch (e) { error.value = e.message } finally { busy.value = false }
}
</script>

<template>
  <h1 class="page-title">逆合成路线规划</h1>
  <p class="page-sub">RDKit BRICS + RECAP 混合断键 → 递归分层路线树 + 启发式置信度（纯本地计算）</p>
  <div v-if="error" class="alert err">{{ error }}</div>

  <div class="card">
    <label class="field">目标分子 SMILES<input v-model="smiles" placeholder="CC(=O)NC1=CC=C(C=C1)O" /></label>
    <div class="grid3" style="margin-top:12px">
      <label class="field">引擎
        <select v-model="engine"><option value="hybrid">hybrid（BRICS+RECAP）</option><option value="brics">BRICS</option><option value="recap">RECAP</option></select></label>
      <label class="field">最大路线数<input v-model="maxRoutes" type="number" min="1" max="20" /></label>
      <label class="field">递归深度<input v-model="maxDepth" type="number" min="1" max="4" /></label>
    </div>
    <div style="margin-top:16px"><button class="btn" :disabled="busy" @click="submit">{{ busy ? '规划中…' : '规划路线' }}</button></div>
  </div>

  <div v-if="result" class="card">
    <h3>目标：{{ result.target }}（{{ result.routes.length }} 条路线）</h3>
    <div v-for="r in result.routes" :key="r.route_id" class="card" style="margin-top:12px; background:rgba(255,255,255,0.03)">
      <b>路线 {{ r.route_id }}</b>　置信度 <span class="pill">{{ (r.confidence * 100).toFixed(0) }}%</span>
      <div v-for="(s, i) in r.steps" :key="i" style="margin-top:10px">
        <div style="font-family:var(--mono, monospace); font-size:13px">
          {{ s.smiles }} <span style="color:var(--muted)">--键型 {{ s.bond_type }}--&gt;</span>
          <span v-for="(p, j) in s.precursors" :key="j">{{ p }}<template v-if="j < s.precursors.length - 1"> + </template></span>
        </div>
      </div>
    </div>
    <p style="color:var(--muted); margin-top:12px; font-size:12px">{{ result.note }}</p>
  </div>
</template>

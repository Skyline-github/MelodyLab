// 后端 API 客户端
const BASE = ''

async function request(path, options = {}) {
  const res = await fetch(BASE + path, options)
  if (!res.ok) {
    let detail = res.statusText
    try { detail = (await res.json()).detail || detail } catch (_) {}
    throw new Error(typeof detail === 'string' ? detail : JSON.stringify(detail))
  }
  return res.json()
}

export const api = {
  // 系统
  health: () => request('/api/health'),
  environment: () => request('/api/environment'),
  getTools: () => request('/api/tools'),
  saveTools: (payload) => request('/api/tools', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(payload) }),

  // 任务中心
  jobs: (kind) => request(`/api/jobs${kind ? `?kind=${kind}` : ''}`),
  job: (id) => request(`/api/jobs/${id}`),
  jobLog: (id) => request(`/api/jobs/${id}/log`),
  pauseJob: (id) => request(`/api/jobs/${id}/pause`, { method: 'POST' }),
  stopJob: (id) => request(`/api/jobs/${id}/stop`, { method: 'POST' }),
  deleteJob: (id) => request(`/api/jobs/${id}`, { method: 'DELETE' }),
  artifactUrl: (id, name) => `${BASE}/api/jobs/${id}/artifacts/${encodeURIComponent(name)}`,

  // 筛选
  screening: (form) => request('/api/screening', { method: 'POST', body: form }),

  // 对接
  docking: (form) => request('/api/docking', { method: 'POST', body: form }),
  dockingResults: (id) => request(`/api/docking/runs/${id}/results`),
  resumeRun: (id) => request(`/api/docking/runs/${id}/resume`, { method: 'POST' }),
  pockets: (form) => request('/api/pockets', { method: 'POST', body: form }),

  // 生成
  generationAvailability: () => request('/api/generation/availability'),
  generation: (form) => request('/api/generation', { method: 'POST', body: form }),

  // Linker
  linkerAvailability: () => request('/api/linker/availability'),
  linker: (form) => request('/api/linker', { method: 'POST', body: form }),

  // MD
  mdAvailability: () => request('/api/md/availability'),
  mdJobs: () => request('/api/md/jobs'),
  mdJob: (form) => request('/api/md/jobs', { method: 'POST', body: form }),
  mdResults: (id) => request(`/api/md/jobs/${id}/results`),
  mdTrajectoryPreview: (id) => `${BASE}/api/md/jobs/${id}/trajectory-preview`,

  // 知识库
  knowledgeActivity: (form) => request('/api/knowledge/activity/search', { method: 'POST', body: form }),
  knowledgeStructures: (form) => request('/api/knowledge/structures/search', { method: 'POST', body: form }),
  knowledgeZinc: (form) => request('/api/knowledge/zinc/search', { method: 'POST', body: form }),

  // PPI
  ppiCreate: (form) => request('/api/ppi/networks', { method: 'POST', body: form }),
  ppiResults: (id) => request(`/api/ppi/networks/${id}/results`),

  // GEO
  geoAnalyze: (form) => request('/api/geo/analyze', { method: 'POST', body: form }),

  // 逆合成
  retroPlan: (form) => request('/api/retro/plan', { method: 'POST', body: form }),

  // 靶点立项评估
  triageAvailability: () => request('/api/triage/availability'),
  triageAssess: (form) => request('/api/triage/assess', { method: 'POST', body: form }),
  triageJob: (id) => request(`/api/triage/jobs/${id}`),
  triageResume: (id) => request(`/api/triage/jobs/${id}/resume`, { method: 'POST' }),
  triageReports: () => request('/api/triage/reports'),
  triageReport: (id) => request(`/api/triage/reports/${id}`),

  // 设置
  getSettings: () => request('/api/settings'),
  setDataRoot: (payload) => request('/api/settings/data-root', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(payload) }),
}

export function fmtTime(s) {
  return s || '—'
}

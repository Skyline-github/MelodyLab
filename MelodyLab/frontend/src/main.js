import { createApp } from 'vue'
import { createRouter, createWebHashHistory } from 'vue-router'
import App from './App.vue'
import HomeView from './views/HomeView.vue'
import TriageView from './views/TriageView.vue'
import ScreeningView from './views/ScreeningView.vue'
import DockingView from './views/DockingView.vue'
import GenerationView from './views/GenerationView.vue'
import LinkerView from './views/LinkerView.vue'
import DegraderView from './views/DegraderView.vue'
import TernaryView from './views/TernaryView.vue'
import MdView from './views/MdView.vue'
import KnowledgeView from './views/KnowledgeView.vue'
import PpiView from './views/PpiView.vue'
import GeoView from './views/GeoView.vue'
import RetroView from './views/RetroView.vue'
import SettingsView from './views/SettingsView.vue'
import JobsView from './views/JobsView.vue'
import './styles.css'

const router = createRouter({
  history: createWebHashHistory(),
  routes: [
    { path: '/', name: 'home', component: HomeView },
    { path: '/triage', name: 'triage', component: TriageView },
    { path: '/screening', name: 'screening', component: ScreeningView },
    { path: '/docking', name: 'docking', component: DockingView },
    { path: '/generation', name: 'generation', component: GenerationView },
    { path: '/linker', name: 'linker', component: LinkerView },
    { path: '/degrader', name: 'degrader', component: DegraderView },
    { path: '/ternary', name: 'ternary', component: TernaryView },
    { path: '/md', name: 'md', component: MdView },
    { path: '/knowledge', name: 'knowledge', component: KnowledgeView },
    { path: '/ppi', name: 'ppi', component: PpiView },
    { path: '/geo', name: 'geo', component: GeoView },
    { path: '/retro', name: 'retro', component: RetroView },
    { path: '/settings', name: 'settings', component: SettingsView },
    { path: '/jobs', name: 'jobs', component: JobsView },
  ],
})

createApp(App).use(router).mount('#app')

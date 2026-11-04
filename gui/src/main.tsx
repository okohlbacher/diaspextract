import { createRoot } from 'react-dom/client'
import App from './App'
import { tauriBridge, type Bridge } from './bridge'
import './index.css'

// Inside the app: the Rust backend. In a plain browser (npm run dev, the UI tests): the mock.
const bridge: Bridge =
  '__TAURI_INTERNALS__' in window ? await tauriBridge() : (await import('./mockBridge')).mockBridge()

createRoot(document.getElementById('root')!).render(<App bridge={bridge} />)

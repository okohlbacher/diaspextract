import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

// Tauri expects the dev server on a fixed port.
export default defineConfig({
  plugins: [react()],
  clearScreen: false,
  build: { target: 'es2022' },
  server: { port: 1420, strictPort: true, watch: { ignored: ['**/src-tauri/**'] } }
})

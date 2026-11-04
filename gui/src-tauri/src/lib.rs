// The backend: find the bundled diaspextract, run it with arguments built from an allowlist, stream its output to the
// window as events, and kill it on Cancel or when the window closes.
//
// The CLI is spawned with an argv array, never a shell string, and build_args() emits only the options the CLI
// declares in params.generated.json (its own -write_ini), so a crafted invoke cannot pass anything else.

use std::collections::{BTreeMap, HashMap};
use std::io::{BufRead, BufReader, Read};
use std::path::PathBuf;
use std::process::{Child, Command, Stdio};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use serde::Serialize;
use tauri::{AppHandle, Emitter, Manager, State, WindowEvent};

const MANIFEST: &str = include_str!("../../src/params.generated.json");
const EXE: &str = if cfg!(windows) { "diaspextract.exe" } else { "diaspextract" };

#[derive(Default)]
struct Running(Mutex<Option<Arc<Mutex<Child>>>>);

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct BinaryInfo {
    bin: String,
    source: &'static str,
    version: Option<String>,
    detail: String,
}

#[derive(Serialize, Clone)]
#[serde(rename_all = "camelCase")]
struct Done {
    code: Option<i32>,
    cancelled: bool,
    elapsed_ms: u64,
}

/// name -> type, from the manifest.
fn declared() -> HashMap<String, String> {
    let v: serde_json::Value = serde_json::from_str(MANIFEST).expect("params.generated.json is valid JSON");
    v["params"]
        .as_array()
        .expect("params array")
        .iter()
        .map(|p| (p["name"].as_str().unwrap().to_string(), p["type"].as_str().unwrap().to_string()))
        .collect()
}

/// The argv after the program name. `in` and `out` are required; `threads` and every declared option are optional.
/// A TOPP flag (type bool) is passed bare when true and left out when false.
pub fn build_args(values: &BTreeMap<String, String>) -> Result<Vec<String>, String> {
    let decl = declared();
    let mut args = Vec::new();
    for key in ["in", "out"] {
        match values.get(key) {
            Some(v) if !v.trim().is_empty() => args.extend([format!("-{key}"), v.clone()]),
            _ => return Err(format!("-{key} is required")),
        }
    }
    if let Some(t) = values.get("threads") {
        t.parse::<u32>().map_err(|_| format!("-threads must be a whole number, not '{t}'"))?;
        args.extend(["-threads".into(), t.clone()]);
    }
    for (k, v) in values {
        if matches!(k.as_str(), "in" | "out" | "threads") {
            continue;
        }
        match decl.get(k).map(String::as_str) {
            None => return Err(format!("unknown option -{k}")),
            Some("bool") => {
                if v == "true" {
                    args.push(format!("-{k}"))
                }
            }
            Some(_) => args.extend([format!("-{k}"), v.clone()]),
        }
    }
    Ok(args)
}

/// DIASPEXTRACT_BIN, else the bundled CLI (the app's resources, or src-tauri/resources in a dev build), else PATH.
fn resolve(app: &AppHandle) -> (PathBuf, &'static str) {
    if let Some(b) = std::env::var_os("DIASPEXTRACT_BIN").map(PathBuf::from).filter(|b| b.exists()) {
        return (b, "DIASPEXTRACT_BIN");
    }
    let mut roots = Vec::new();
    if let Ok(rd) = app.path().resource_dir() {
        roots.push(rd.join("diaspextract"));
    }
    roots.push(PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/resources/diaspextract")));
    for r in roots {
        let b = r.join("bin").join(EXE);
        if b.exists() {
            return (b, "bundled");
        }
    }
    (PathBuf::from(EXE), "PATH")
}

fn command(bin: &PathBuf) -> Command {
    let mut c = Command::new(bin);
    // libOpenMS loads OpenBLAS, which starts a spinning thread per core; the tool never calls BLAS.
    if std::env::var_os("OPENBLAS_NUM_THREADS").is_none() {
        c.env("OPENBLAS_NUM_THREADS", "1");
    }
    // The release app has no console, so Windows would open one for every child.
    #[cfg(windows)]
    {
        use std::os::windows::process::CommandExt;
        c.creation_flags(0x0800_0000); // CREATE_NO_WINDOW
    }
    c
}

/// "Version: 1.9.0 (OpenMS ...)" from --help: proof the binary starts and finds its libraries here.
fn parse_version(text: &str) -> Option<String> {
    let rest = &text[text.find("Version: ")? + 9..];
    Some(rest.split_whitespace().next()?.to_string())
}

#[tauri::command]
async fn probe(app: AppHandle) -> BinaryInfo {
    // off the main thread: a cold first start (Gatekeeper, antivirus) takes seconds
    tauri::async_runtime::spawn_blocking(move || {
        let (bin, source) = resolve(&app);
        let out = command(&bin).arg("--help").stdin(Stdio::null()).output();
        let (version, detail) = match out {
            Ok(o) => {
                let text = format!("{}{}", String::from_utf8_lossy(&o.stdout), String::from_utf8_lossy(&o.stderr));
                match parse_version(&text) {
                    Some(v) => (Some(v), String::new()),
                    None => (None, text.lines().rev().take(5).collect::<Vec<_>>().join("\n")),
                }
            }
            Err(e) => (None, e.to_string()),
        };
        BinaryInfo { bin: bin.display().to_string(), source, version, detail }
    })
    .await
    .unwrap()
}

#[tauri::command]
fn command_line(values: BTreeMap<String, String>) -> Result<Vec<String>, String> {
    build_args(&values)
}

#[tauri::command]
fn run(app: AppHandle, running: State<Running>, values: BTreeMap<String, String>) -> Result<(), String> {
    let mut slot = running.0.lock().unwrap();
    if slot.is_some() {
        return Err("a run is already in progress".into());
    }
    let args = build_args(&values)?;
    let (bin, _) = resolve(&app);
    CANCELLED.store(false, std::sync::atomic::Ordering::SeqCst);
    let mut child = command(&bin)
        .args(&args)
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| format!("cannot start {}: {e}", bin.display()))?;
    let readers = [
        pump(app.clone(), child.stdout.take().unwrap()),
        pump(app.clone(), child.stderr.take().unwrap()),
    ];
    let child = Arc::new(Mutex::new(child));
    *slot = Some(child.clone());
    let started = Instant::now();
    std::thread::spawn(move || {
        let status = loop {
            if let Ok(Some(s)) = child.lock().unwrap().try_wait() {
                break s;
            }
            std::thread::sleep(Duration::from_millis(100));
        };
        for r in readers {
            let _ = r.join(); // every line is out before "done"
        }
        let running = app.state::<Running>();
        *running.0.lock().unwrap() = None;
        let cancelled = CANCELLED.swap(false, std::sync::atomic::Ordering::SeqCst);
        let elapsed_ms = started.elapsed().as_millis() as u64;
        let _ = app.emit("dx:done", Done { code: status.code(), cancelled, elapsed_ms });
    });
    Ok(())
}

static CANCELLED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

fn pump(app: AppHandle, stream: impl Read + Send + 'static) -> std::thread::JoinHandle<()> {
    std::thread::spawn(move || {
        for line in BufReader::new(stream).lines().map_while(Result::ok) {
            let _ = app.emit("dx:log", line);
        }
    })
}

#[tauri::command]
fn cancel(running: State<Running>) -> bool {
    match running.0.lock().unwrap().as_ref() {
        Some(c) => {
            CANCELLED.store(true, std::sync::atomic::Ordering::SeqCst);
            c.lock().unwrap().kill().is_ok()
        }
        None => false,
    }
}

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run_app() {
    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_opener::init())
        .manage(Running::default())
        .invoke_handler(tauri::generate_handler![probe, command_line, run, cancel])
        // closing the window must not leave the CLI running
        .on_window_event(|w, e| {
            if let WindowEvent::Destroyed = e {
                if let Some(c) = w.app_handle().state::<Running>().0.lock().unwrap().as_ref() {
                    let _ = c.lock().unwrap().kill();
                }
            }
        })
        .run(tauri::generate_context!())
        .expect("error while running the DIAspeXtract app");
}

#[cfg(test)]
mod tests {
    use super::*;

    fn m(kv: &[(&str, &str)]) -> BTreeMap<String, String> {
        kv.iter().map(|(k, v)| (k.to_string(), v.to_string())).collect()
    }

    #[test]
    fn required_and_order() {
        let a = build_args(&m(&[("out", "o.mzML"), ("in", "r.d"), ("threads", "0"), ("gate:delta_rt", "5")])).unwrap();
        assert_eq!(a, ["-in", "r.d", "-out", "o.mzML", "-threads", "0", "-gate:delta_rt", "5"]);
        assert!(build_args(&m(&[("in", "r.d")])).is_err());
        assert!(build_args(&m(&[("in", " "), ("out", "o")])).is_err());
    }

    #[test]
    fn allowlist_and_flags() {
        let base = [("in", "r.d"), ("out", "o.mzML")];
        let mut v = m(&base);
        v.insert("ini".into(), "evil.ini".into());
        assert!(build_args(&v).unwrap_err().contains("unknown option"));
        let mut v = m(&base);
        v.insert("threads".into(), "4; rm -rf".into());
        assert!(build_args(&v).is_err());
        let mut v = m(&base);
        v.insert("diag:spectrum_meta".into(), "true".into());
        assert!(build_args(&v).unwrap().ends_with(&["-diag:spectrum_meta".to_string()]));
        v.insert("diag:spectrum_meta".into(), "false".into());
        assert_eq!(build_args(&v).unwrap().len(), 4);
    }

    #[test]
    fn version() {
        assert_eq!(parse_version("x\nVersion: 1.9.0 (OpenMS 3.6)\n").as_deref(), Some("1.9.0"));
        assert_eq!(parse_version("no"), None);
    }
}

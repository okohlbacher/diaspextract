// No console window next to the app on Windows in a release build.
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

fn main() {
    diaspextract_gui_lib::run_app()
}

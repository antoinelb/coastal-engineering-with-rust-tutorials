use colored::*;
use std::io::Write;

pub fn load_print(msg: &str) {
    let symbol = "[*]".bold();
    print!("\r{} {:<pad$}", symbol, msg, pad = 80);
    std::io::stdout().flush().unwrap();
}

pub fn done_print(msg: &str) {
    let symbol = "[+]".bright_green().bold();
    println!("\r{} {:<pad$}", symbol, msg, pad = 80);
}

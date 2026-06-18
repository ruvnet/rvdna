//! `dna-storage` — DNA Storage Codec Simulator CLI.
//!
//! Turn any file into synthetic DNA strands, inject sequencing errors, and
//! reconstruct the original — demonstrating Reed–Solomon + fountain-code error
//! correction with no wet lab.
//!
//! ```text
//! dna-storage encode   <input> [--out archive.json] [--fasta strands.fasta]
//!                              [--block-size N] [--parity N] [--overhead F]
//! dna-storage simulate <input> [--out recovered.bin] [--sub P] [--ins P]
//!                              [--del P] [--drop P] [--coverage N] [...encode opts]
//! dna-storage decode   <archive.json> <reads.fasta> [--out recovered.bin]
//! dna-storage stats    <archive.json>
//! ```

use std::collections::HashMap;
use std::fs;
use std::process::exit;

use rvdna::storage::{channel::ErrorModel, DnaArchive, DnaStorageCodec, EncodeParams};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        usage();
        exit(2);
    }

    let cmd = args[1].as_str();
    let rest = &args[2..];
    let result = match cmd {
        "encode" => cmd_encode(rest),
        "simulate" => cmd_simulate(rest),
        "decode" => cmd_decode(rest),
        "stats" => cmd_stats(rest),
        "-h" | "--help" | "help" => {
            usage();
            Ok(())
        }
        other => {
            eprintln!("unknown command: {other}\n");
            usage();
            exit(2);
        }
    };

    if let Err(e) = result {
        eprintln!("error: {e}");
        exit(1);
    }
}

fn usage() {
    eprintln!(
        "dna-storage — DNA Storage Codec Simulator\n\n\
         USAGE:\n\
         \x20 dna-storage encode   <input> [--out archive.json] [--fasta strands.fasta]\n\
         \x20                              [--block-size N] [--parity N] [--overhead F] [--seed N]\n\
         \x20 dna-storage simulate <input> [--out recovered.bin] [--sub P] [--ins P] [--del P]\n\
         \x20                              [--drop P] [--coverage N] [encode opts]\n\
         \x20 dna-storage decode   <archive.json> <reads.fasta> [--out recovered.bin]\n\
         \x20 dna-storage stats    <archive.json>\n"
    );
}

// ---------------------------------------------------------------------------
// argument parsing helpers
// ---------------------------------------------------------------------------

struct Opts {
    positionals: Vec<String>,
    flags: HashMap<String, String>,
}

fn parse_opts(args: &[String]) -> Opts {
    let mut positionals = Vec::new();
    let mut flags = HashMap::new();
    let mut i = 0;
    while i < args.len() {
        let a = &args[i];
        if let Some(key) = a.strip_prefix("--") {
            let val = if i + 1 < args.len() && !args[i + 1].starts_with("--") {
                i += 1;
                args[i].clone()
            } else {
                "true".to_string()
            };
            flags.insert(key.to_string(), val);
        } else {
            positionals.push(a.clone());
        }
        i += 1;
    }
    Opts { positionals, flags }
}

impl Opts {
    fn get_usize(&self, key: &str, default: usize) -> usize {
        self.flags
            .get(key)
            .and_then(|v| v.parse().ok())
            .unwrap_or(default)
    }
    fn get_f64(&self, key: &str, default: f64) -> f64 {
        self.flags
            .get(key)
            .and_then(|v| v.parse().ok())
            .unwrap_or(default)
    }
    fn get_u64(&self, key: &str, default: u64) -> u64 {
        self.flags
            .get(key)
            .and_then(|v| v.parse().ok())
            .unwrap_or(default)
    }
    fn get_str(&self, key: &str) -> Option<&str> {
        self.flags.get(key).map(|s| s.as_str())
    }
}

fn encode_params(o: &Opts) -> EncodeParams {
    let d = EncodeParams::default();
    EncodeParams {
        block_size: o.get_usize("block-size", d.block_size),
        rs_parity: o.get_usize("parity", d.rs_parity),
        overhead: o.get_f64("overhead", d.overhead),
        max_homopolymer: d.max_homopolymer,
        seed: o.get_u64("seed", d.seed),
    }
}

// ---------------------------------------------------------------------------
// FASTA I/O
// ---------------------------------------------------------------------------

fn write_fasta(path: &str, archive: &DnaArchive) -> std::io::Result<()> {
    let mut out = String::new();
    for s in &archive.strands {
        out.push_str(&format!(">strand_{}\n{}\n", s.index, s.sequence));
    }
    fs::write(path, out)
}

fn read_fasta_sequences(path: &str) -> std::io::Result<Vec<String>> {
    let text = fs::read_to_string(path)?;
    let mut reads = Vec::new();
    let mut cur = String::new();
    for line in text.lines() {
        if line.starts_with('>') {
            if !cur.is_empty() {
                reads.push(std::mem::take(&mut cur));
            }
        } else {
            cur.push_str(line.trim());
        }
    }
    if !cur.is_empty() {
        reads.push(cur);
    }
    Ok(reads)
}

// ---------------------------------------------------------------------------
// commands
// ---------------------------------------------------------------------------

fn cmd_encode(args: &[String]) -> Result<(), String> {
    let o = parse_opts(args);
    let input = o.positionals.first().ok_or("encode: missing <input>")?;
    let data = fs::read(input).map_err(|e| format!("read {input}: {e}"))?;

    let codec = DnaStorageCodec::new(encode_params(&o));
    let archive = codec.encode(input, &data).map_err(|e| e.to_string())?;

    print_archive_stats(&archive);

    let out = o.get_str("out").unwrap_or("archive.json");
    let json = serde_json::to_string_pretty(&archive).map_err(|e| e.to_string())?;
    fs::write(out, json).map_err(|e| format!("write {out}: {e}"))?;
    println!("\n  archive written  -> {out}");

    if let Some(fasta) = o.get_str("fasta") {
        write_fasta(fasta, &archive).map_err(|e| format!("write {fasta}: {e}"))?;
        println!("  strands (FASTA)  -> {fasta}");
    }
    Ok(())
}

fn cmd_simulate(args: &[String]) -> Result<(), String> {
    let o = parse_opts(args);
    let input = o.positionals.first().ok_or("simulate: missing <input>")?;
    let data = fs::read(input).map_err(|e| format!("read {input}: {e}"))?;

    let codec = DnaStorageCodec::new(encode_params(&o));
    let model = ErrorModel {
        p_sub: o.get_f64("sub", 0.01),
        p_ins: o.get_f64("ins", 0.0),
        p_del: o.get_f64("del", 0.0),
        p_drop: o.get_f64("drop", 0.10),
        coverage: o.get_usize("coverage", 3),
    };
    let channel_seed = o.get_u64("channel-seed", 1);

    let (archive, report) = codec
        .simulate(input, &data, &model, channel_seed)
        .map_err(|e| e.to_string())?;

    print_archive_stats(&archive);
    println!("\n  ── noisy channel ─────────────────────────────");
    println!("  substitution rate  {:.3}", model.p_sub);
    println!("  insertion rate     {:.3}", model.p_ins);
    println!("  deletion rate      {:.3}", model.p_del);
    println!("  strand dropout     {:.3}", model.p_drop);
    println!("  coverage           {}x", model.coverage);

    println!("\n  ── recovery ──────────────────────────────────");
    println!("  reads sequenced    {}", report.reads_in);
    println!(
        "  strands recovered  {} (RS-corrected)",
        report.strands_recovered
    );
    println!(
        "  blocks recovered   {}/{} (fountain peel)",
        report.blocks_recovered, archive.num_blocks
    );
    let exact = report.bytes == data;
    println!(
        "  CRC32 verify       {}",
        if report.crc_ok {
            "PASS ✅"
        } else {
            "FAIL ❌"
        }
    );
    println!("  byte-exact match   {}", if exact { "yes" } else { "no" });

    if report.crc_ok {
        let out = o.get_str("out").unwrap_or("recovered.bin");
        fs::write(out, &report.bytes).map_err(|e| format!("write {out}: {e}"))?;
        println!("\n  recovered file     -> {out}");
        println!(
            "\n  🧬 stored {} bytes in DNA and got them back.",
            data.len()
        );
        Ok(())
    } else {
        Err("recovery failed — increase --overhead/--coverage or lower error rates".into())
    }
}

fn cmd_decode(args: &[String]) -> Result<(), String> {
    let o = parse_opts(args);
    let archive_path = o
        .positionals
        .first()
        .ok_or("decode: missing <archive.json>")?;
    let reads_path = o
        .positionals
        .get(1)
        .ok_or("decode: missing <reads.fasta>")?;

    let json = fs::read_to_string(archive_path).map_err(|e| format!("read {archive_path}: {e}"))?;
    let archive: DnaArchive = serde_json::from_str(&json).map_err(|e| e.to_string())?;
    let reads = read_fasta_sequences(reads_path).map_err(|e| format!("read {reads_path}: {e}"))?;

    let codec = DnaStorageCodec::new(archive.params.clone());
    let report = codec.decode(&archive, &reads).map_err(|e| e.to_string())?;

    println!("  reads in           {}", report.reads_in);
    println!("  strands recovered  {}", report.strands_recovered);
    println!(
        "  blocks recovered   {}/{}",
        report.blocks_recovered, archive.num_blocks
    );
    println!(
        "  CRC32 verify       {}",
        if report.crc_ok { "PASS" } else { "FAIL" }
    );

    if report.crc_ok {
        let out = o.get_str("out").unwrap_or("recovered.bin");
        fs::write(out, &report.bytes).map_err(|e| format!("write {out}: {e}"))?;
        println!("  recovered file     -> {out}");
        Ok(())
    } else {
        Err("recovery failed".into())
    }
}

fn cmd_stats(args: &[String]) -> Result<(), String> {
    let o = parse_opts(args);
    let archive_path = o
        .positionals
        .first()
        .ok_or("stats: missing <archive.json>")?;
    let json = fs::read_to_string(archive_path).map_err(|e| format!("read {archive_path}: {e}"))?;
    let archive: DnaArchive = serde_json::from_str(&json).map_err(|e| e.to_string())?;
    print_archive_stats(&archive);
    Ok(())
}

fn print_archive_stats(archive: &DnaArchive) {
    let strand_len = archive
        .strands
        .first()
        .map(|s| s.sequence.len())
        .unwrap_or(0);
    println!("  ── DNA archive ───────────────────────────────");
    println!("  file               {}", archive.filename);
    println!("  original size      {} bytes", archive.byte_len);
    println!("  source blocks      {}", archive.num_blocks);
    println!("  strands (oligos)   {}", archive.strands.len());
    println!("  strand length      {} bases", strand_len);
    println!("  total bases        {}", archive.total_bases());
    println!(
        "  density            {:.3} bits/base",
        archive.bits_per_base()
    );
    println!("  mean GC content    {:.1}%", archive.mean_gc() * 100.0);
    if let Some(s) = archive.strands.first() {
        let run = rvdna::storage::constraints::max_homopolymer_run(&s.sequence);
        println!("  max homopolymer    {} (run length)", run);
    }
}

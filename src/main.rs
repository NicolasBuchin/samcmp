use clap::{Arg, Command};
use flate2::read::GzDecoder;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

#[derive(Debug, Clone)]
struct AlignmentRecord {
    qname: String,
    flag: u16,
    rname: String,
    pos: u32,
    mapq: u8,
    cigar: String,
    alignment_score: Option<i32>,
}

#[derive(Debug, Default)]
struct ComparisonStats {
    unmapped_to_mapped: usize,
    mapped_to_unmapped: usize,
    both_mapped_better: usize,
    both_mapped_worse: usize,
    both_mapped_same: usize,
    both_mapped: usize,
    both_unmapped: usize,
    only_in_first: usize,
    only_in_second: usize,
    total_score_first: i64,
    total_score_second: i64,
    total_reads_first: usize,
    total_reads_second: usize,
    mapped_reads_first: usize,
    mapped_reads_second: usize,
    min_score_first: Option<i32>,
    max_score_first: Option<i32>,
    min_score_second: Option<i32>,
    max_score_second: Option<i32>,
    median_first: Option<f64>,
    median_second: Option<f64>,
    q1_first: Option<f64>,
    q1_second: Option<f64>,
    q3_first: Option<f64>,
    q3_second: Option<f64>,
    scored_reads_first: usize,
    scored_reads_second: usize,
    better_perfect_to_imperfect: usize,
    better_imperfect_to_perfect: usize,
    better_both_perfect: usize,
    better_both_imperfect: usize,
    worse_perfect_to_imperfect: usize,
    worse_imperfect_to_perfect: usize,
    worse_both_perfect: usize,
    worse_both_imperfect: usize,
}

fn main() {
    let matches = Command::new("sam-compare")
        .version("1.0")
        .about("Compare SAM files to analyze alignment improvements")
        .arg(
            Arg::new("file1")
                .help("First SAM file (baseline)")
                .required(true)
                .index(1),
        )
        .arg(
            Arg::new("file2")
                .help("Second SAM file (comparison)")
                .required(true)
                .index(2),
        )
        .arg(
            Arg::new("extract")
                .short('x')
                .value_name("FASTQ")
                .help("Extract reads with worse scores from this FASTQ file (supports .gz)"),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .value_name("OUTPUT")
                .help("Output file for extracted reads (required if -x is used)"),
        )
        .arg(
            Arg::new("limit")
                .short('n')
                .value_name("NUMBER")
                .value_parser(clap::value_parser!(usize))
                .help("Limit the number of reads to extract (only valid with -x)"),
        )
        .arg(
            Arg::new("sim0")
                .long("sim0")
                .action(clap::ArgAction::SetTrue)
                .help("Track perfect CIGAR strings (e.g., '150=') and filter extraction to perfect→imperfect reads"),
        )
        .get_matches();

    let file1_path = matches.get_one::<String>("file1").unwrap();
    let file2_path = matches.get_one::<String>("file2").unwrap();
    let extract_fastq = matches.get_one::<String>("extract");
    let output_file = matches.get_one::<String>("output");
    let extract_limit = matches.get_one::<usize>("limit");
    let sim0_mode = matches.get_flag("sim0");

    if extract_fastq.is_some() && output_file.is_none() {
        eprintln!("Error: Output file (-o) is required when using extraction (-x)");
        std::process::exit(1);
    }

    if extract_limit.is_some() && extract_fastq.is_none() {
        eprintln!("Error: Limit (-n) can only be used with extraction (-x)");
        std::process::exit(1);
    }

    println!("Comparing SAM files:");
    println!("  Baseline: {}", file1_path);
    println!("  Comparison: {}", file2_path);
    if sim0_mode {
        println!("  Mode: Perfect CIGAR tracking enabled");
    }
    if let Some(fastq_path) = extract_fastq {
        println!("  FASTQ for extraction: {}", fastq_path);
        if let Some(output_path) = output_file {
            println!("  Output file: {}", output_path);
        }
        if let Some(limit) = extract_limit {
            println!("  Extraction limit: {} reads", limit);
        }
        if sim0_mode {
            println!("  Extraction filter: Only reads going from perfect to imperfect CIGAR");
        }
    }
    println!();

    match compare_sam_files(file1_path, file2_path, sim0_mode) {
        Ok(stats) => {
            print_comparison_results(&stats, sim0_mode);

            if let (Some(fastq_path), Some(output_path)) = (extract_fastq, output_file) {
                println!("\n=== Extracting Reads with Worse Scores ===");
                match extract_worse_reads(
                    file1_path,
                    file2_path,
                    fastq_path,
                    output_path,
                    extract_limit.copied(),
                    sim0_mode,
                ) {
                    Ok(extracted_count) => {
                        println!(
                            "Successfully extracted {} reads to {}",
                            extracted_count, output_path
                        );
                    }
                    Err(e) => {
                        eprintln!("Error during extraction: {}", e);
                        std::process::exit(1);
                    }
                }
            }
        }
        Err(e) => {
            eprintln!("Error: {}", e);
            std::process::exit(1);
        }
    }
}

fn parse_sam_line(line: &str) -> Result<AlignmentRecord, Box<dyn std::error::Error>> {
    let fields: Vec<&str> = line.split('\t').collect();

    if fields.len() < 11 {
        return Err("Invalid SAM line: insufficient fields".into());
    }

    let qname = fields[0].to_string();
    let flag = fields[1].parse::<u16>()?;
    let rname = fields[2].to_string();
    let pos = fields[3].parse::<u32>()?;
    let mapq = fields[4].parse::<u8>()?;
    let cigar = fields[5].to_string();

    let mut alignment_score = None;
    for field in &fields[11..] {
        if field.starts_with("AS:i:") {
            if let Ok(score) = field[5..].parse::<i32>() {
                alignment_score = Some(score);
                break;
            }
        }
    }

    Ok(AlignmentRecord {
        qname,
        flag,
        rname,
        pos,
        mapq,
        cigar,
        alignment_score,
    })
}

fn is_mapped(record: &AlignmentRecord) -> bool {
    (record.flag & 0x4) == 0 && record.rname != "*"
}

fn is_perfect_cigar(cigar: &str) -> bool {
    if cigar == "*" {
        return false;
    }

    let trimmed = cigar.trim();
    if trimmed.is_empty() {
        return false;
    }

    let mut has_content = false;
    for c in trimmed.chars() {
        if c.is_ascii_digit() {
            has_content = true;
        } else if c == '=' {
            if !has_content {
                return false;
            }
        } else {
            return false;
        }
    }

    trimmed.ends_with('=') && has_content
}

fn read_sam_file(
    file_path: &str,
) -> Result<HashMap<String, AlignmentRecord>, Box<dyn std::error::Error>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    let mut alignments = HashMap::new();
    let mut line_count = 0;

    for line in reader.lines() {
        let line = line?;
        line_count += 1;

        if line.starts_with('@') {
            continue;
        }

        if line.trim().is_empty() {
            continue;
        }

        match parse_sam_line(&line) {
            Ok(record) => {
                alignments.insert(record.qname.clone(), record);
            }
            Err(e) => {
                eprintln!("Warning: Failed to parse line {}: {}", line_count, e);
                continue;
            }
        }
    }

    println!("Loaded {} alignments from {}", alignments.len(), file_path);
    Ok(alignments)
}

fn get_worse_read_names(
    file1_path: &str,
    file2_path: &str,
    sim0_mode: bool,
) -> Result<HashSet<String>, Box<dyn std::error::Error>> {
    let alignments1 = read_sam_file(file1_path)?;
    let alignments2 = read_sam_file(file2_path)?;

    let worse_reads: HashSet<String> = alignments1
        .par_iter()
        .filter_map(|(read_name, record1)| {
            if let Some(record2) = alignments2.get(read_name) {
                let mapped1 = is_mapped(record1);
                let mapped2 = is_mapped(record2);

                if mapped1 && !mapped2 {
                    if sim0_mode {
                        if is_perfect_cigar(&record1.cigar) {
                            return Some(read_name.clone());
                        }
                    } else {
                        return Some(read_name.clone());
                    }
                }

                if mapped1 && mapped2 {
                    if let (Some(score1), Some(score2)) =
                        (record1.alignment_score, record2.alignment_score)
                    {
                        if score2 < score1 {
                            if sim0_mode {
                                let perfect1 = is_perfect_cigar(&record1.cigar);
                                let perfect2 = is_perfect_cigar(&record2.cigar);
                                if perfect1 && !perfect2 {
                                    return Some(read_name.clone());
                                }
                            } else {
                                return Some(read_name.clone());
                            }
                        }
                    }
                }
            }
            None
        })
        .collect();

    Ok(worse_reads)
}

fn create_fastq_reader(file_path: &str) -> Result<Box<dyn BufRead>, Box<dyn std::error::Error>> {
    let file = File::open(file_path)?;

    if file_path.ends_with(".gz") {
        let decoder = GzDecoder::new(file);
        Ok(Box::new(BufReader::new(decoder)))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn extract_fastq_read_name(header: &str) -> String {
    let name = if header.starts_with('@') {
        &header[1..]
    } else {
        header
    };

    let base_name = name.split_whitespace().next().unwrap_or(name);

    if let Some(slash_pos) = base_name.rfind('/') {
        if let Some(suffix) = base_name.get(slash_pos + 1..) {
            if suffix == "1" || suffix == "2" {
                return base_name[..slash_pos].to_string();
            }
        }
    }

    base_name.to_string()
}

fn extract_worse_reads(
    file1_path: &str,
    file2_path: &str,
    fastq_path: &str,
    output_path: &str,
    limit: Option<usize>,
    sim0_mode: bool,
) -> Result<usize, Box<dyn std::error::Error>> {
    let worse_read_names = get_worse_read_names(file1_path, file2_path, sim0_mode)?;

    if worse_read_names.is_empty() {
        println!("No reads with worse scores found - creating empty output file");
        File::create(output_path)?;
        return Ok(0);
    }

    let reader = create_fastq_reader(fastq_path)?;
    let output_file = File::create(output_path)?;
    let mut writer = BufWriter::new(output_file);

    let mut lines = reader.lines();
    let mut extracted_count = 0;

    while let Some(header_line) = lines.next() {
        let header = header_line?;

        if !header.starts_with('@') {
            continue;
        }

        let sequence = lines.next().ok_or("Unexpected end of file")??;
        let plus = lines.next().ok_or("Unexpected end of file")??;
        let quality = lines.next().ok_or("Unexpected end of file")??;

        let read_name = extract_fastq_read_name(&header);

        if worse_read_names.contains(&read_name) {
            if let Some(limit_val) = limit {
                if extracted_count >= limit_val {
                    break;
                }
            }

            writeln!(writer, "{}", header)?;
            writeln!(writer, "{}", sequence)?;
            writeln!(writer, "{}", plus)?;
            writeln!(writer, "{}", quality)?;

            extracted_count += 1;
        }
    }

    writer.flush()?;
    Ok(extracted_count)
}

fn compare_sam_files(
    file1_path: &str,
    file2_path: &str,
    sim0_mode: bool,
) -> Result<ComparisonStats, Box<dyn std::error::Error>> {
    let alignments1 = read_sam_file(file1_path)?;
    let alignments2 = read_sam_file(file2_path)?;

    let mut stats = ComparisonStats::default();
    stats.total_reads_first = alignments1.len();
    stats.total_reads_second = alignments2.len();

    let (score_stats1, mapped_count1) = calculate_score_stats(&alignments1);
    let (score_stats2, mapped_count2) = calculate_score_stats(&alignments2);

    stats.mapped_reads_first = mapped_count1;
    stats.mapped_reads_second = mapped_count2;
    stats.total_score_first = score_stats1.total_score;
    stats.total_score_second = score_stats2.total_score;
    stats.min_score_first = score_stats1.min_score;
    stats.max_score_first = score_stats1.max_score;
    stats.median_first = score_stats1.median;
    stats.q1_first = score_stats1.q1;
    stats.q3_first = score_stats1.q3;
    stats.min_score_second = score_stats2.min_score;
    stats.max_score_second = score_stats2.max_score;
    stats.median_second = score_stats2.median;
    stats.q1_second = score_stats2.q1;
    stats.q3_second = score_stats2.q3;
    stats.scored_reads_first = score_stats1.scored_count;
    stats.scored_reads_second = score_stats2.scored_count;

    let mut all_reads: std::collections::HashSet<String> = alignments1.keys().cloned().collect();
    all_reads.extend(alignments2.keys().cloned());

    let comparison_results: Vec<_> = all_reads
        .par_iter()
        .map(|read_name| {
            let record1 = alignments1.get(read_name);
            let record2 = alignments2.get(read_name);

            match (record1, record2) {
                (Some(r1), Some(r2)) => {
                    let mapped1 = is_mapped(r1);
                    let mapped2 = is_mapped(r2);

                    match (mapped1, mapped2) {
                        (false, true) => ComparisonResult::UnmappedToMapped,
                        (true, false) => ComparisonResult::MappedToUnmapped,
                        (true, true) => {
                            let perfect1 = is_perfect_cigar(&r1.cigar);
                            let perfect2 = is_perfect_cigar(&r2.cigar);

                            match (r1.alignment_score, r2.alignment_score) {
                                (Some(score1), Some(score2)) => {
                                    if score2 > score1 {
                                        ComparisonResult::BothMappedBetter(perfect1, perfect2)
                                    } else if score2 < score1 {
                                        ComparisonResult::BothMappedWorse(perfect1, perfect2)
                                    } else {
                                        ComparisonResult::BothMappedSame
                                    }
                                }
                                _ => ComparisonResult::BothMappedSame,
                            }
                        }
                        (false, false) => ComparisonResult::BothUnmapped,
                    }
                }
                (Some(_), None) => ComparisonResult::OnlyInFirst,
                (None, Some(_)) => ComparisonResult::OnlyInSecond,
                (None, None) => unreachable!(),
            }
        })
        .collect();

    for result in comparison_results {
        match result {
            ComparisonResult::UnmappedToMapped => stats.unmapped_to_mapped += 1,
            ComparisonResult::MappedToUnmapped => stats.mapped_to_unmapped += 1,
            ComparisonResult::BothMappedBetter(perfect1, perfect2) => {
                stats.both_mapped_better += 1;
                if sim0_mode {
                    match (perfect1, perfect2) {
                        (true, false) => stats.better_perfect_to_imperfect += 1,
                        (false, true) => stats.better_imperfect_to_perfect += 1,
                        (true, true) => stats.better_both_perfect += 1,
                        (false, false) => stats.better_both_imperfect += 1,
                    }
                }
            }
            ComparisonResult::BothMappedWorse(perfect1, perfect2) => {
                stats.both_mapped_worse += 1;
                if sim0_mode {
                    match (perfect1, perfect2) {
                        (true, false) => stats.worse_perfect_to_imperfect += 1,
                        (false, true) => stats.worse_imperfect_to_perfect += 1,
                        (true, true) => stats.worse_both_perfect += 1,
                        (false, false) => stats.worse_both_imperfect += 1,
                    }
                }
            }
            ComparisonResult::BothMappedSame => stats.both_mapped_same += 1,
            ComparisonResult::OnlyInFirst => stats.only_in_first += 1,
            ComparisonResult::OnlyInSecond => stats.only_in_second += 1,
            ComparisonResult::BothUnmapped => stats.both_unmapped += 1,
        }
    }

    stats.both_mapped = stats.both_mapped_better + stats.both_mapped_worse + stats.both_mapped_same;

    Ok(stats)
}

#[derive(Debug)]
struct ScoreStats {
    total_score: i64,
    min_score: Option<i32>,
    max_score: Option<i32>,
    scored_count: usize,
    median: Option<f64>,
    q1: Option<f64>,
    q3: Option<f64>,
}

fn calculate_score_stats(alignments: &HashMap<String, AlignmentRecord>) -> (ScoreStats, usize) {
    let mut total_score = 0i64;
    let mut min_score: Option<i32> = None;
    let mut max_score: Option<i32> = None;
    let mut scored_count = 0usize;
    let mut mapped_count = 0usize;
    let mut scores: Vec<i32> = Vec::new();

    for record in alignments.values() {
        if is_mapped(record) {
            mapped_count += 1;

            if let Some(score) = record.alignment_score {
                total_score += score as i64;
                scored_count += 1;
                scores.push(score);

                min_score = Some(min_score.map_or(score, |min| min.min(score)));
                max_score = Some(max_score.map_or(score, |max| max.max(score)));
            }
        }
    }

    let (median, q1, q3) = if scores.is_empty() {
        (None, None, None)
    } else {
        scores.sort_unstable();
        let median = calculate_median(&scores);
        let q1 = calculate_quartile(&scores, 0.25);
        let q3 = calculate_quartile(&scores, 0.75);
        (Some(median), Some(q1), Some(q3))
    };

    (
        ScoreStats {
            total_score,
            min_score,
            max_score,
            scored_count,
            median,
            q1,
            q3,
        },
        mapped_count,
    )
}

fn calculate_median(sorted_scores: &[i32]) -> f64 {
    let len = sorted_scores.len();
    if len % 2 == 0 {
        (sorted_scores[len / 2 - 1] + sorted_scores[len / 2]) as f64 / 2.0
    } else {
        sorted_scores[len / 2] as f64
    }
}

fn calculate_quartile(sorted_scores: &[i32], quantile: f64) -> f64 {
    let len = sorted_scores.len();
    let index = quantile * (len - 1) as f64;
    let lower_index = index.floor() as usize;
    let upper_index = index.ceil() as usize;

    if lower_index == upper_index {
        sorted_scores[lower_index] as f64
    } else {
        let weight = index - lower_index as f64;
        sorted_scores[lower_index] as f64 * (1.0 - weight)
            + sorted_scores[upper_index] as f64 * weight
    }
}

#[derive(Debug)]
enum ComparisonResult {
    UnmappedToMapped,
    MappedToUnmapped,
    BothMappedBetter(bool, bool),
    BothMappedWorse(bool, bool),
    BothMappedSame,
    OnlyInFirst,
    OnlyInSecond,
    BothUnmapped,
}

fn print_comparison_results(stats: &ComparisonStats, sim0_mode: bool) {
    println!("=== SAM File Comparison Results ===");
    println!();

    println!("Read Count Summary:");
    println!(
        "  Baseline file:        {:>8} reads",
        stats.total_reads_first
    );
    println!(
        "  Comparison file:      {:>8} reads",
        stats.total_reads_second
    );
    println!(
        "  Only in baseline:     {:>8} reads ({:.1}%)",
        stats.only_in_first,
        stats.only_in_first as f64 / stats.total_reads_first as f64 * 100.0
    );
    println!(
        "  Only in comparison:   {:>8} reads ({:.1}%)",
        stats.only_in_second,
        stats.only_in_second as f64 / stats.total_reads_second as f64 * 100.0
    );
    println!();

    println!("Mapping Status:");
    println!(
        "  Mapped in baseline:   {:>8} reads ({:.1}%)",
        stats.mapped_reads_first,
        stats.mapped_reads_first as f64 / stats.total_reads_first as f64 * 100.0
    );
    println!(
        "  Mapped in comparison: {:>8} reads ({:.1}%)",
        stats.mapped_reads_second,
        stats.mapped_reads_second as f64 / stats.total_reads_second as f64 * 100.0
    );

    let common_reads = stats.both_mapped
        + stats.both_unmapped
        + stats.unmapped_to_mapped
        + stats.mapped_to_unmapped;

    println!(
        "  Mapped in both:       {:>8} reads ({:.1}%)",
        stats.both_mapped,
        if common_reads > 0 {
            stats.both_mapped as f64 / common_reads as f64 * 100.0
        } else {
            0.0
        }
    );
    println!(
        "  Unmapped in both:     {:>8} reads ({:.1}%)",
        stats.both_unmapped,
        if common_reads > 0 {
            stats.both_unmapped as f64 / common_reads as f64 * 100.0
        } else {
            0.0
        }
    );
    println!();

    println!("Mapping Status Changes:");
    println!(
        "  Unmapped → Mapped:    {:>8} reads ({:.1}%)",
        stats.unmapped_to_mapped,
        if common_reads > 0 {
            stats.unmapped_to_mapped as f64 / common_reads as f64 * 100.0
        } else {
            0.0
        }
    );
    println!(
        "  Mapped → Unmapped:    {:>8} reads ({:.1}%)",
        stats.mapped_to_unmapped,
        if common_reads > 0 {
            stats.mapped_to_unmapped as f64 / common_reads as f64 * 100.0
        } else {
            0.0
        }
    );
    println!();

    println!("Alignment Score Changes (for reads mapped in both files):");
    println!(
        "  Improved scores:      {:>8} reads ({:.1}%)",
        stats.both_mapped_better,
        if stats.both_mapped > 0 {
            stats.both_mapped_better as f64 / stats.both_mapped as f64 * 100.0
        } else {
            0.0
        }
    );
    println!(
        "  Worse scores:         {:>8} reads ({:.1}%)",
        stats.both_mapped_worse,
        if stats.both_mapped > 0 {
            stats.both_mapped_worse as f64 / stats.both_mapped as f64 * 100.0
        } else {
            0.0
        }
    );
    println!(
        "  Same scores:          {:>8} reads ({:.1}%)",
        stats.both_mapped_same,
        if stats.both_mapped > 0 {
            stats.both_mapped_same as f64 / stats.both_mapped as f64 * 100.0
        } else {
            0.0
        }
    );
    println!();

    if sim0_mode {
        println!("Perfect CIGAR Analysis:");
        println!("  Improved scores breakdown:");
        println!(
            "    Imperfect → Perfect:  {:>8} reads ({:.1}%)",
            stats.better_imperfect_to_perfect,
            if stats.both_mapped_better > 0 {
                stats.better_imperfect_to_perfect as f64 / stats.both_mapped_better as f64 * 100.0
            } else {
                0.0
            }
        );
        println!(
            "    Both perfect:         {:>8} reads ({:.1}%)",
            stats.better_both_perfect,
            if stats.both_mapped_better > 0 {
                stats.better_both_perfect as f64 / stats.both_mapped_better as f64 * 100.0
            } else {
                0.0
            }
        );
        println!(
            "    Both imperfect:       {:>8} reads ({:.1}%)",
            stats.better_both_imperfect,
            if stats.both_mapped_better > 0 {
                stats.better_both_imperfect as f64 / stats.both_mapped_better as f64 * 100.0
            } else {
                0.0
            }
        );
        println!();
        println!("  Worse scores breakdown:");
        println!(
            "    Perfect → Imperfect:  {:>8} reads ({:.1}%)",
            stats.worse_perfect_to_imperfect,
            if stats.both_mapped_worse > 0 {
                stats.worse_perfect_to_imperfect as f64 / stats.both_mapped_worse as f64 * 100.0
            } else {
                0.0
            }
        );
        println!(
            "    Both perfect:         {:>8} reads ({:.1}%)",
            stats.worse_both_perfect,
            if stats.both_mapped_worse > 0 {
                stats.worse_both_perfect as f64 / stats.both_mapped_worse as f64 * 100.0
            } else {
                0.0
            }
        );
        println!(
            "    Both imperfect:       {:>8} reads ({:.1}%)",
            stats.worse_both_imperfect,
            if stats.both_mapped_worse > 0 {
                stats.worse_both_imperfect as f64 / stats.both_mapped_worse as f64 * 100.0
            } else {
                0.0
            }
        );
        println!();
    }

    println!("Alignment Score Statistics:");
    println!("  Baseline file:");
    if stats.scored_reads_first > 0 {
        let avg_first = stats.total_score_first as f64 / stats.scored_reads_first as f64;
        println!("    Total score:        {:>12}", stats.total_score_first);
        println!("    Average score:      {:>12.2}", avg_first);
        if let Some(min) = stats.min_score_first {
            println!("    Minimum score:      {:>12}", min);
        }
        if let Some(q1) = stats.q1_first {
            println!("    1st quartile:       {:>12.2}", q1);
        }
        if let Some(median) = stats.median_first {
            println!("    Median score:       {:>12.2}", median);
        }
        if let Some(q3) = stats.q3_first {
            println!("    3rd quartile:       {:>12.2}", q3);
        }
        if let Some(max) = stats.max_score_first {
            println!("    Maximum score:      {:>12}", max);
        }
    } else {
        println!("    No alignment scores found");
    }

    println!("  Comparison file:");
    if stats.scored_reads_second > 0 {
        let avg_second = stats.total_score_second as f64 / stats.scored_reads_second as f64;
        println!("    Total score:        {:>12}", stats.total_score_second);
        println!("    Average score:      {:>12.2}", avg_second);
        if let Some(min) = stats.min_score_second {
            println!("    Minimum score:      {:>12}", min);
        }
        if let Some(q1) = stats.q1_second {
            println!("    1st quartile:       {:>12.2}", q1);
        }
        if let Some(median) = stats.median_second {
            println!("    Median score:       {:>12.2}", median);
        }
        if let Some(q3) = stats.q3_second {
            println!("    3rd quartile:       {:>12.2}", q3);
        }
        if let Some(max) = stats.max_score_second {
            println!("    Maximum score:      {:>12}", max);
        }
    } else {
        println!("    No alignment scores found");
    }

    if stats.scored_reads_first > 0 && stats.scored_reads_second > 0 {
        let score_diff = stats.total_score_second - stats.total_score_first;
        let avg_first = stats.total_score_first as f64 / stats.scored_reads_first as f64;
        let avg_second = stats.total_score_second as f64 / stats.scored_reads_second as f64;
        let avg_diff = avg_second - avg_first;

        println!("  Score differences:");
        println!("    Total difference:   {:>+12}", score_diff);
        println!("    Average difference: {:>+12.2}", avg_diff);

        if stats.total_score_first != 0 {
            let percent_change = (score_diff as f64 / stats.total_score_first as f64) * 100.0;
            println!("    Percent change:     {:>+11.2}%", percent_change);
        }
    }
}

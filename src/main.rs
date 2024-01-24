// use clap::{Arg, Command, ArgAction};
use clap::Parser;

use std::fs::File;
use std::path::Path;
use std::collections::HashMap;
use serde_json::Value;
use std::io::{self, BufRead, BufReader};
use bzip2::read::BzDecoder;

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

fn tajd_cstes(nb: u64) -> (f64, f64) {
	let mut a1: f64 = 0.0;
	let mut a2: f64 = 0.0;
	for n in 1..nb {
		a1 += 1.0 / n as f64;
		a2 += 1.0 / (n.pow(2) as f64);
	}

	let b1: f64 = (nb + 1) as f64 / ((3 * (nb - 1)) as f64);
	let b2: f64 = ((2 * (nb.pow(2) + nb + 3)) as f64 )/ ((9 * nb * (nb - 1)) as f64);
	let c1: f64 = b1 - (1.0 / a1);
	let c2: f64 = b2 - ((nb + 2) as f64 / (a1 * nb as f64)) + (a2 / a1.powi(2));
	let e1: f64 = c1 / a1;
	let e2: f64 = c2 / (a1.powi(2) + a2);
	return (e1, e2);
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// gtf file
    #[arg(short, long)]
    gtf: String,

    /// dbsnp json file
    #[arg(short, long)]
    dbsnp: String,
}

#[derive(Debug, Eq, PartialEq, Hash)]
struct PiTmp {
    total: u64,
    refcount: u64,
}
fn main() {
  let args = Args::parse();
  let mut exon_posset: HashMap<String, HashMap<u32, bool>> = HashMap::new();

  // Read in pi values
  if let Ok(lines) = read_lines(args.gtf) {
	for line in lines {
        if let Ok(ip) = line {
			let spl: Vec<&str> = ip.split("\t").collect();
            if spl[2] != "exon" {
                continue;
            }
			let startp = spl[3].parse::<u32>().unwrap();
			let stopp = spl[4].parse::<u32>().unwrap();
            let mut currexon: HashMap<u32, bool> = HashMap::new();
            for i in startp..(stopp+1) {
                currexon.insert(i, true);
            }
            let desc = spl[8];
            let desc_parts: Vec<&str> = desc.split("; ").collect();
            for part in desc_parts {
                let part_parts: Vec<&str> = part.split(" ").collect();
                if part_parts[0] == "gene_name" {
                    let gene_id = part_parts[1].replace("\"", "");
                    if !exon_posset.contains_key(&gene_id) {
                        exon_posset.insert(gene_id, currexon.clone());
                    }
                    else {
                        let tmp = exon_posset.get(&gene_id).unwrap();
                        currexon.extend(tmp.clone().into_iter());
                        exon_posset.insert(gene_id, currexon.clone());
                    }
                }
            }
		}        
	}
  }
  
    let mut gene_tmppi: HashMap<String, HashMap<PiTmp, u16>> = HashMap::new();
    let file = File::open(args.dbsnp).unwrap();
    let reader = BufReader::new(BzDecoder::new(file));
    for line in reader.lines() {
        let l = line.unwrap();
        let v: Value = serde_json::from_str(&l).unwrap();

        let mut a_count: u64 = 0;
        let mut c_count: u64 = 0;
        let mut freq_idx: usize = 127;
        let mut tmp_pos: u32 = 0;
        let mut on_exon: bool = false;

        'marker: for idx in 0..v["primary_snapshot_data"]["placements_with_allele"].as_array().unwrap().len() {
            if v["primary_snapshot_data"]["placements_with_allele"][idx]["is_ptlp"].as_bool().unwrap() {
                for idx1 in 0..v["primary_snapshot_data"]["placements_with_allele"][idx]["alleles"].as_array().unwrap().len() {
                    if v["primary_snapshot_data"]["placements_with_allele"][idx]["alleles"][idx1]["allele"]["spdi"]["deleted_sequence"].as_str().unwrap() == v["primary_snapshot_data"]["placements_with_allele"][idx]["alleles"][idx1]["allele"]["spdi"]["inserted_sequence"].as_str().unwrap() {
                        tmp_pos = v["primary_snapshot_data"]["placements_with_allele"][idx]["alleles"][idx1]["allele"]["spdi"]["position"].as_u64().unwrap() as u32;
                        freq_idx = idx1;
                        break 'marker;
                    }
                }
            }
        }

        let gene_symbol_option = v["primary_snapshot_data"]["allele_annotations"][freq_idx]["assembly_annotation"][0]["genes"][0]["locus"].as_str();
        let gene_symbol = v["primary_snapshot_data"]["allele_annotations"][freq_idx]["assembly_annotation"][0]["genes"][0]["locus"].as_str().unwrap();

        if let Some(gene_symbol) = gene_symbol_option {
            // gene_symbol is not None, you can use it here
            if exon_posset.contains_key(&gene_symbol.to_string()) {
                let tmp = exon_posset.get(&gene_symbol.to_string()).unwrap();
                if tmp.contains_key(&tmp_pos) {
                    on_exon = true;
                }
                else {
                    continue;
                }
            } 
        } else {
            println!("failed to get gene symbol: {}", v["refsnp_id"]);
            continue;
        }        

        if freq_idx == 127 {
            println!("failed to get freq_idx: {}", v["refsnp_id"]);
            continue;
        }
        // let vt1 = v["primary_snapshot_data"]["allele_annotations"][freq_idx]["assembly_annotation"][0]["genes"][0]["rnas"][0]["sequence_ontology"][0]["name"].to_string();
        for idx in 0..v["primary_snapshot_data"]["allele_annotations"][freq_idx]["frequency"].as_array().unwrap().len() {
            a_count += v["primary_snapshot_data"]["allele_annotations"][freq_idx]["frequency"][idx]["allele_count"].as_u64().unwrap() as u64;
            c_count += v["primary_snapshot_data"]["allele_annotations"][freq_idx]["frequency"][idx]["total_count"].as_u64().unwrap() as u64;
        }
        // let vt2: String = v["primary_snapshot_data"]["variant_type"].to_string();
        if a_count == 0 || c_count == 0 || a_count == c_count || !on_exon {
            continue;
        }
        let tmp_pos: PiTmp = PiTmp {total: c_count, refcount: a_count};

        if let Some(tmp) = gene_tmppi.get_mut(&gene_symbol.to_string()) {
            if tmp.contains_key(&tmp_pos) {
                let tmp2 = tmp.get(&tmp_pos).unwrap() + 1;
                tmp.insert(tmp_pos, tmp2);
                println!("test1: {}\t{}\t{}", gene_symbol, a_count, c_count);
            }
            else {
                tmp.insert(tmp_pos, 1);
                println!("test2: {}\t{}\t{}", gene_symbol, a_count, c_count);
            }
        }
        else {
            gene_tmppi.insert(gene_symbol.to_string(), HashMap::from([(tmp_pos, 1),]));
            println!("test3: {}\t{}\t{}", gene_symbol, a_count, c_count);
        }

    };
    

    for gene in gene_tmppi.keys() {
        let tmp = gene_tmppi.get(gene).unwrap();
  
        let mut all_nb: u64 = 0;
        let det: f64 = exon_posset.get(gene).unwrap().len() as f64;
        let mut pi_tmp: u64 = 0;
        let mut pi_min_tmp: u64 = 0;
        let mut h_tmp: u64 = 0;
        let mut seg: u64 = 0;
        let mut tajd: f64 = 0.0;
        let mut tajd_min: f64 = 0.0;
        for pos in tmp.keys() {
            if pos.total > all_nb {
                all_nb = pos.total;
            }
            pi_tmp += pos.refcount * (pos.total - pos.refcount) * (tmp.get(pos).unwrap().clone() as u64);
            pi_min_tmp += pos.total * (pos.total - 1) * (tmp.get(pos).unwrap().clone() as u64);
            h_tmp +=  pos.refcount.pow(2) * (tmp.get(pos).unwrap().clone() as u64);
            seg += tmp.get(pos).unwrap().clone() as u64;
            println!("{}\t{}\t{}\t{}", gene, pos.refcount, pos.total, tmp.get(pos).unwrap());
        }

        // calculate a_1
        let mut a_1: f64 = 0.0;
        for m in 1..(all_nb+1) {
            a_1 += 1.0 / m as f64;
        }

         let tipi: u64 = all_nb;
         let pi = pi_tmp as f64 / ((tipi as f64 * (tipi as f64 - 1.0)) / 2.0) / det ;
          let pi_min = pi_min_tmp as f64 / ((tipi as f64 * (tipi as f64 - 1.0)) / 2.0) / det ;
          let theta = seg as f64 / a_1 / det;

          let mut d: f64 = 0.0;
         let mut dmin: f64 = 0.0;
          let (e1, e2) = tajd_cstes(all_nb);
          let p1p2: f64 = e1 * seg as f64 + e2 * seg as f64 * (seg as f64 - 1.0);
          if pi > 0.0 && theta > 0.0 {
          	d = pi - theta;
          	dmin = pi_min - theta;
            dmin = dmin.abs();
          }

          if p1p2 > 0.0 {
            tajd = d * det / p1p2.sqrt();
            tajd_min = dmin * det / p1p2.sqrt();
          }

          let hache = h_tmp as f64/ ((tipi as f64 * (tipi as f64 - 1.0)) / 2.0);
          let final_h = pi - hache;

          println!("pi:\tthetaW:\tTajima's D:\tTajima's D normalized:\tH:");
          println!("{}\t{}\t{}\t{}\t{}", pi, theta, tajd, tajd/tajd_min, final_h);
    }
}
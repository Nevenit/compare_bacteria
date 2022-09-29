mod bad_profiler;
use crate::bad_profiler::bad_profiler::Profiler;

//use rustc_hash::FxHashMap;
use std::fs;
use std::str;
use std::env;
use std::time::{Instant, SystemTime};

const CODE: [i32; 26] = [0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3];
const LEN: usize = 6;
const  AA_NUMBER: usize = 20;
const EPSILON: f64 = 1e-010;
const M: usize = usize::pow(AA_NUMBER, LEN as u32);
const M1: usize = usize::pow(AA_NUMBER, (LEN - 1) as u32);
const M2: usize = usize::pow(AA_NUMBER, (LEN - 2) as u32);

pub struct Vars {
    bacteria_count: i32,
    bacteria_path: Vec<String>
}

impl Vars {
    pub fn add_bacteria_path(&mut self, name: String) {
        self.bacteria_path.push(name)
    }
    pub fn set_bacteria_count(&mut self, num: i32) {
        self.bacteria_count = num;
    }
}

// Stores the public variables of Bacteria
pub struct Bacteria {
    pub count: i64,
    pub tv: Vec<f64>,
    pub ti: Vec<i64>
}

// Stores the private variables of Bacteria
pub struct BacteriaCounters {
    vector: Box<[i64]>,
    second: Box<[i64]>,
    one_l: Box<[i64]>,
    indexes: i64,
    total: i64,
    total_l: i64,
    complement: i64
}

// Implements the Bacteria class
impl Bacteria {
    // This function initialises the bacteria struct, reads the file and calculates the kmers
    pub fn init(&mut self, file: &str, profiler: &mut Profiler) {
        profiler.start("read_file");

        // Read the file
        let bacteria_file = fs::read_to_string(&file).unwrap();

        // convert it to bytes
        let file_bytes = bacteria_file.as_bytes();

        profiler.end("read_file");

        profiler.start("init_counters");

        // Initialise the BacteriaCounters which is used as the private variables for the struct
        let mut bc = BacteriaCounters {
            vector: vec![0; M].into_boxed_slice(),
            second: vec![0; M1].into_boxed_slice(),
            one_l: vec![0; AA_NUMBER].into_boxed_slice(),
            indexes: 0,
            total: 0,
            total_l: 0,
            complement: 0
        };
        profiler.end("init_counters");

        profiler.start("extract_kmers");

        // Loop through all bytes of the file
        let mut i = 0;
        while i < file_bytes.len(){
            if file_bytes[i] as char == '>' {
                profiler.start("compute_buffer");

                // Skip all bytes until new lines because this is a comment line
                while file_bytes[i] as char != '\n' {
                    i += 1;
                }

                // Initialise the buffer array
                let mut buffer: [u8; LEN - 1] = [0 as u8; LEN - 1];

                // Loop for the length of our kmer and add the bytes to the buffer
                let mut j = 0;
                while j < LEN - 1 {
                    i += 1;
                    buffer[j] = file_bytes[i];
                    j += 1;
                }

                profiler.end("compute_buffer");
                profiler.start("init_buffer");

                // Add buffer to vectors
                self.init_buffer(&mut bc, buffer);
                profiler.end("init_buffer");
            }// If char isn't new line or end of file
            else if file_bytes[i] as char != '\n' && file_bytes[i] != 13 {
                profiler.start("cont_buffer");

                // Add char to kmers
                self.cont_buffer(&mut bc, file_bytes[i] as char);
                profiler.end("cont_buffer");
            }
            i += 1;
        }
        profiler.end("extract_kmers");
        profiler.start("fill_arrays");

        // Initialise variables
        let total_plus_complement: f64 = bc.total as f64 + bc.complement as f64 ;
        let total_div_2: f64 = bc.total as f64 * 0.5;
        let mut i_mod_aa_number: usize = 0;
        let mut i_div_aa_number: usize = 0;
        let mut i_mod_m1: usize = 0;
        let mut i_div_m1: usize = 0;

        // Fill the one_l_div_total array
        let mut one_l_div_total = [0.0; AA_NUMBER];
        for i in 0..AA_NUMBER {
            one_l_div_total[i] = bc.one_l[i] as f64 / bc.total_l as f64;
        }

        // Fill the second_div_total array
        let mut second_div_total = vec![0.0; M1].into_boxed_slice();
        for i in 0..M1 {
            second_div_total[i] = bc.second[i] as f64 / total_plus_complement;
        }

        // Initialise variables
        //let mut t: Box<[f64]> = vec![0.0; M].into_boxed_slice();

        profiler.end("fill_arrays");
        profiler.start("calculate_t");

        self.tv = vec![];
        self.ti = vec![];

        // Loop Through all possible kmers
        for i in 0..M {
            let p1: f64 = second_div_total[i_div_aa_number];
            let p2: f64 = one_l_div_total[i_mod_aa_number];
            let p3: f64 = second_div_total[i_mod_m1];
            let p4: f64 = one_l_div_total[i_div_m1];
            let stochastic: f64 = (p1 * p2 + p3 * p4) * total_div_2;

            if i_mod_aa_number == AA_NUMBER - 1 {
                i_mod_aa_number = 0;
                i_div_aa_number += 1;
            }
            else {
                i_mod_aa_number += 1;
            }

            if i_mod_m1 == M1 - 1 {
                i_mod_m1 = 0;
                i_div_m1 += 1;
            }
            else {
                i_mod_m1 += 1;
            }

            if stochastic > EPSILON {
                //t[i] = (bc.vector[i] as f64 - stochastic) / stochastic;
                self.tv.push((bc.vector[i] as f64 - stochastic) / stochastic);
                self.ti.push(i as i64);
                self.count += 1;
            }
        }
        profiler.end("calculate_t");
    }

    fn init_buffer(&mut self, bc: &mut BacteriaCounters, buffer: [u8; LEN - 1]) {
        bc.complement += 1;
        bc.indexes = 0;
        for i in 0..buffer.len() {
            let enc = encode(buffer[i] as char);
            bc.one_l[enc] += 1;
            bc.total_l += 1;
            bc.indexes = bc.indexes * AA_NUMBER as i64 + enc as i64;
        }
        bc.second[bc.indexes as usize] += 1;
    }

    fn cont_buffer(&mut self, bc: &mut BacteriaCounters, ch: char) {
        let enc = encode(ch);
        bc.one_l[enc] += 1;
        bc.total_l += 1;
        let index = bc.indexes * AA_NUMBER as i64 + enc as i64;
        bc.vector[index as usize] += 1;
        bc.total += 1;
        bc.indexes = (bc.indexes % M2 as i64) * AA_NUMBER as i64 + enc as i64;
        bc.second[bc.indexes as usize] += 1;
    }
}

// Encode the character as in int 'A' as 0, 'B' as 1 and so on
fn encode(ch: char) -> usize {
    CODE[(ch as u8 - 'A' as u8) as usize] as usize
}

// Reads the file that contains the list of bacteria name and writes them in the "global" struct
fn read_input_file(file: &str, program_vars: &mut Vars) {
    // Read file to string
    let file_contents = fs::read_to_string(&file).unwrap();

    // Split the string into an iterator line array
    let mut file_lines = file_contents.lines();

    // Parse the first line into a number
    let num: i32 = file_lines.next().unwrap().parse::<i32>().unwrap();

    // Update the bacteria count in the "global" struct
    program_vars.set_bacteria_count(num);

    // Iterate through the remaining lines adding them to the bacteria name array as path strings
    for string in file_lines {
        program_vars.add_bacteria_path(format!("data/{}.faa", string));
    }
}

fn compare_bacteria(b1: &Bacteria, b2: &Bacteria, profiler: &mut Profiler) -> f64 {
    let mut correlation: f64 = 0.0;
    let mut vector_len1: f64 = 0.0;
    let mut vector_len2: f64 = 0.0;
    let mut p1: usize = 0;
    let mut p2: usize = 0;

    profiler.start("p1_and_p2");
    while p1 < b1.count as usize && p2 < b2.count as usize {
        let _n1 = b1.ti[p1];
        let _n2 = b2.ti[p2];

        if _n1 < _n2 {
            let t1 = b1.tv[p1];
            vector_len1 += t1 * t1;
            p1 += 1;
        } else if _n2 < _n1 {
            let t2 = b2.tv[p2];
            p2 += 1;
            vector_len2 += t2 *  t2;
        } else {
            let t1 = b1.tv[p1];
            let t2 = b2.tv[p2];
            p1 += 1;
            p2 += 1;
            vector_len1 += t1 * t1;
            vector_len2 += t2 * t2;
            correlation += t1 * t2;
        }
    }
    profiler.end("p1_and_p2");
    profiler.start("p1");
    while p1 < b1.count as usize {
        let n1 = b1.ti[p1];
        let t1 = b1.tv[p1];
        p1 += 1;
        vector_len1 += t1 * t1;
    }
    profiler.end("p1");
    profiler.start("p2");
    while p2 < b2.count as usize {
        let n2 = b2.ti[p2];
        let t2 = b2.tv[p2];
        p2 += 1;
        vector_len2 += t2 * t2;
    }
    profiler.end("p2");
    return correlation / (vector_len1.sqrt() * vector_len2.sqrt());
}

/*
This function is where most of the work happens.
It starts by initialising a new class for each bacteria, reading the files, loading all kmers
and calculating the statistics. Then for every bacteria class it compares it with all of the
other bacteria and prints the result
 */
fn compare_all_bacteria(program_vars: &mut Vars, profiler: &mut Profiler) {
    profiler.start("init_bacteria");

    let mut bacteria_array = vec![];

    // Loop through every bacteria
    for i in 0..program_vars.bacteria_count {
        println!("load {} of {}", i + 1, program_vars.bacteria_count);

        // Initialise the struct
        bacteria_array.push(Bacteria {
            count: 0,
            tv: vec![],
            ti: vec![]
        });

        profiler.start(format!("init_bacteria {}", i).as_str());

        // Call the bacteria initialization function
        bacteria_array[i as usize].init(&program_vars.bacteria_path[i as usize], profiler);

        profiler.end(format!("init_bacteria {}", i).as_str());
    }

    profiler.end("init_bacteria");
    profiler.start("compare_bacteria");

    // Loop through every bacteria
    for i in 0..program_vars.bacteria_count {
        // This skips the comparisons that have already been done
        for j in i+1..program_vars.bacteria_count {
            print!("{} {} -> ", i, j);

            // Calculate and print the correlation between the two bacteria
            let correlation = compare_bacteria(&bacteria_array[i as usize], &bacteria_array[j as usize], profiler);
            println!("{:.20}", correlation);
        }
    }
    profiler.end("compare_bacteria");
}

fn main() {
    // Read get the arguments passed to the program
    let args: Vec<String> = env::args().collect();

    // Create a new instance of the profiler
    let mut profiler = Profiler::new();

    profiler.start("main");

    // Create a struct to be used as a global variables
    let mut program_vars = Vars {
        bacteria_count: 0,
        bacteria_path: vec![]
    };

    // Start measuring time
    let start_time = Instant::now();

    // Read the list file to get the count and names of bacteria to compare
    read_input_file(&args[1], &mut program_vars);

    // Analise each bacteria file and compare them
    compare_all_bacteria(&mut program_vars, &mut profiler);

    profiler.end("main");

    let elapsed_time = start_time.elapsed();

    println!("time elapsed: {0} millis", elapsed_time.as_millis());

    profiler.print_results();
}

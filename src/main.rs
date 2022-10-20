mod bad_profiler;
use crate::bad_profiler::bad_profiler::Profiler;

//use std::time::Instant;
//use rustc_hash::FxHashMap;
use std::fs;
use std::str;
use std::env;
use std::time::Instant;
use std::sync::{Arc, Mutex};
use std::thread;
use threadpool::ThreadPool;
use std::sync::mpsc::channel;

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

pub struct Bacteria {
    //vector: Box<[i64]>,
    //second: Box<[i64]>,
    //one_l: Box<[i64]>,
    //indexes: i64,
    //total: i64,
    //total_l: i64,
    //complement: i64,
    pub count: usize,
    pub tv: Vec<f64>,
    pub ti: Vec<i64>
}

pub struct BacteriaCounters {
    vector: Box<[i32]>,
    second: Box<[i32]>,
    one_l: Box<[i32]>,
    indexes: i64,
    total: i64,
    total_l: i64,
    complement: i64
}

impl Bacteria {
    // This function initialises the bacteria struct, reads the file and calculates the kmers
    pub fn init(&mut self, file: &str) {
        //profiler.start("read_file");

        // Read the file
        let bacteria_file = fs::read_to_string(&file).unwrap();

        // convert it to bytes
        let file_bytes = bacteria_file.as_bytes();

        //profiler.end("read_file");

        //profiler.start("init_counters");

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
        //profiler.end("init_counters");

        //profiler.start("extract_kmers");

        // Loop through all bytes of the file
        let mut i = 0;
        while i < file_bytes.len(){
            if file_bytes[i] as char == '>' {
                //profiler.start("compute_buffer");

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

                //profiler.end("compute_buffer");
                //profiler.start("init_buffer");

                // Add buffer to vectors
                self.init_buffer(&mut bc, buffer);
                //profiler.end("init_buffer");
            }// If char isn't new line or end of file
            else if file_bytes[i] as char != '\n' && file_bytes[i] != 13 {
                //profiler.start("cont_buffer");

                // Add char to kmers
                self.cont_buffer(&mut bc, file_bytes[i] as char);
                //profiler.end("cont_buffer");
            }
            i += 1;
        }
        //profiler.end("extract_kmers");
        //profiler.start("fill_arrays");

        // Initialise variables
        let total_plus_complement: f64 = bc.total as f64 + bc.complement as f64 ;
        let total_div_2: f64 = bc.total as f64 * 0.5;

        // Fill the one_l_div_total array
        let mut one_l_div_total: Vec<(usize, f64)> = vec![];
        for i in 0..AA_NUMBER {
            if bc.one_l[i] == 0 {continue}
            one_l_div_total.push((i, bc.one_l[i] as f64 / bc.total_l as f64));
        }

        // Fill the second_div_total array
        let mut second_div_total: Vec<(usize, f64)> = vec![];
        for i in 0..M1 {
            if bc.second[i] == 0 {continue}
            second_div_total.push((i, bc.second[i] as f64 / total_plus_complement));
        }

        //profiler.end("fill_arrays");
        //profiler.start("calculate p1 p2 p3 p4");

        // Calculate and fill the p1p2 array
        let mut p1p2: Vec<(usize, f64)> = vec![];
        for div_aa in &second_div_total {
            for mod_aa in &one_l_div_total {
                p1p2.push((div_aa.0 * AA_NUMBER + mod_aa.0, div_aa.1 * mod_aa.1))
            }
        }

        // Calculate and fill the 0304 array
        let mut p3p4: Vec<(usize, f64)> = vec![];
        for div_m1 in &one_l_div_total {
            for mod_m1 in &second_div_total {
                p3p4.push((div_m1.0 * M1 + mod_m1.0, div_m1.1 * mod_m1.1))
            }
        }

        //profiler.end("calculate p1 p2 p3 p4");
        //profiler.start("calculate_t");

        self.tv = vec![];
        self.ti = vec![];

        let mut p1p2_i = 0;
        let mut p3p4_i = 0;

        while p1p2_i < p1p2.len() && p3p4_i < p3p4.len() {
            let index: usize;
            let stochastic: f64;

            if p1p2[p1p2_i].0 < p3p4[p3p4_i].0 {
                stochastic = p1p2[p1p2_i].1 * total_div_2;
                index = p1p2[p1p2_i].0;
                p1p2_i += 1;
            } else if p1p2[p1p2_i].0 > p3p4[p3p4_i].0 {
                stochastic = p3p4[p3p4_i].1 * total_div_2;
                index = p3p4[p3p4_i].0;
                p3p4_i += 1;
            } else {
                stochastic = (p1p2[p1p2_i].1 + p3p4[p3p4_i].1) * total_div_2;
                index = p1p2[p1p2_i].0;
                p1p2_i += 1;
                p3p4_i += 1;
            }

            if stochastic > EPSILON {
                self.tv.push((bc.vector[index] as f64 - stochastic) / stochastic);
                self.ti.push(index as i64);
                self.count += 1;
            }
        }

        while p1p2_i < p1p2.len() {
            let stochastic = p1p2[p1p2_i].1 * total_div_2;
            if stochastic > EPSILON {
                self.tv.push((bc.vector[p1p2[p1p2_i].0] as f64 - stochastic) / stochastic);
                self.ti.push(p1p2[p1p2_i].0 as i64);
                self.count += 1;
            }
            p1p2_i += 1;
        }

        while p3p4_i < p3p4.len() {
            let stochastic = p3p4[p3p4_i].1 * total_div_2;
            if stochastic > EPSILON {
                self.tv.push((bc.vector[p3p4[p3p4_i].0] as f64 - stochastic) / stochastic);
                self.ti.push(p3p4[p3p4_i].0 as i64);
                self.count += 1;
            }
            p3p4_i += 1;
        }

        //profiler.end("calculate_t");
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

fn encode(ch: char) -> usize {
    CODE[(ch as u8 - 'A' as u8) as usize] as usize
}

fn read_input_file(file: &str, program_vars: &mut Vars) {
    let file_contents = fs::read_to_string(&file).unwrap();
    let mut file_lines = file_contents.lines();

    let num: i32 = file_lines.next().unwrap().parse::<i32>().unwrap();

    program_vars.set_bacteria_count(num);

    for string in file_lines {
        program_vars.add_bacteria_path(format!("data/{}.faa", string));
    }
}

fn compare_bacteria(b1: &Bacteria, b2: &Bacteria, /*profiler: &mut Profiler*/) -> f64 {
    let mut correlation: f64 = 0.0;
    let mut vector_len1: f64 = 0.0;
    let mut vector_len2: f64 = 0.0;
    let mut p1: usize = 0;
    let mut p2: usize = 0;

    //profiler.start("p1_and_p2");
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
    //profiler.end("p1_and_p2");
    //profiler.start("p1");
    while p1 < b1.count as usize {
        let t1 = b1.tv[p1];
        p1 += 1;
        vector_len1 += t1 * t1;
    }
    //profiler.end("p1");
    //profiler.start("p2");
    while p2 < b2.count as usize {
        let t2 = b2.tv[p2];
        p2 += 1;
        vector_len2 += t2 * t2;
    }
    //profiler.end("p2");
    return correlation / (vector_len1.sqrt() * vector_len2.sqrt());
}

fn triangular_number(n: usize) -> usize{
    let mut nn = n;
    let mut i = n;
    loop {
        if i == 0 {return 0}
        if i == 1 {return nn}
        i -= 1;
        nn += i;
    }
}


fn compare_all_bacteria(program_vars: &mut Vars, profiler: &mut Profiler, thread_count: usize) -> String {
    profiler.start("init_bacteria");
    let mut bacteria_array = vec![];

    let pool = ThreadPool::new(thread_count);
    let (tx, rx) = channel();
    for i in 0..program_vars.bacteria_count {
        println!("load {} of {}", i + 1, program_vars.bacteria_count);
        bacteria_array.push(Arc::new(Mutex::new(Bacteria {
            count: 0,
            tv: vec![],
            ti: vec![]
        })));
        //profiler.start(format!("init_bacteria {}", i).as_str());

        let tx = tx.clone();
        let path = program_vars.bacteria_path[i as usize].clone();
        let thread_bacteria_array = Arc::clone(&bacteria_array[i as usize]);
        pool.execute(move || {
            let mut thread_bacteria_array_local = thread_bacteria_array.lock().unwrap();
            thread_bacteria_array_local.init(&path);
            tx.send(1);
        });

        //profiler.end(format!("init_bacteria {}", i).as_str());
    }
    profiler.end("init_bacteria");
    profiler.start("wait for pool to finish");
    for _ in 0..program_vars.bacteria_count {
        rx.recv().unwrap();
    }
    profiler.end("wait for pool to finish");

    profiler.start("clone array");
    let mut bacteria_array_arc = vec![];
    for bacteria in bacteria_array {
        let test = bacteria.lock().unwrap();
        bacteria_array_arc.push(Arc::new(Bacteria {
            count: test.count.clone(),
            tv: test.tv.clone(),
            ti: test.ti.clone()
        }));
    }
    profiler.end("clone array");
    profiler.start("compare_bacteria");
    let mut output = String::new();
    let result_array_arc = Arc::new(Mutex::new(vec![String::new(); ((program_vars.bacteria_count.pow(2) - program_vars.bacteria_count) / 2) as usize]));

    let pool = ThreadPool::new(thread_count);

    let (tx, rx) = channel();
    for i in 0..program_vars.bacteria_count as usize {
        for j in i+1..program_vars.bacteria_count as usize {
            //print!("{} {} -> ", i, j);
            let b1 = Arc::clone(&bacteria_array_arc[i as usize]);
            let b2 = Arc::clone(&bacteria_array_arc[j as usize]);
            let result_array = Arc::clone(&result_array_arc);
            let _i = i.clone();
            let _j = j.clone();
            let bacteria_count = program_vars.bacteria_count.clone();
            let tx = tx.clone();
            let compared_indexes = format!("{} {} ->", i, j);
            pool.execute(move || {
                let correlation = compare_bacteria(&b1, &b2/*, profiler*/);
                let index = ((_i * (bacteria_count as usize - 1)) - triangular_number(_i as usize)) + _j - 1;
                let mut result_array_lock = result_array.lock().unwrap();
                //let result_array_ptr = (result_array.as_ptr() as *mut String).offset(index as isize);
                //*result_array_ptr = format!("{} {:.20}\n", compared_indexes, correlation);
                result_array_lock[index] = format!("{} {:.20}\n", compared_indexes, correlation);
                drop(result_array_lock);
                tx.send(1);
           });
        }
    }


    for _ in 0..((program_vars.bacteria_count.pow(2) - program_vars.bacteria_count) / 2) {
        rx.recv().unwrap();
    }

    let results = result_array_arc.clone().lock().unwrap().clone();
    //let results = Arc::clone(&result_array_arc);
    for i in results {
        output.push_str(&*i);
    }
    profiler.end("compare_bacteria");
    print!("{}", output);

    return output;
}

fn verify_output(mut output: String) {
    let mut file_contents = fs::read_to_string("validation.txt").unwrap();

    // Remove new lines because different operating systems use different line endings
    output = output.replace("\r\n", "");
    file_contents = file_contents.replace("\r\n", "");
    output = output.replace('\n', "");
    file_contents = file_contents.replace('\n', "");

    if output == file_contents {
        println!("Validation successful");
    } else {
        println!("Validation failed");
    }
    //fs::write("output_data.txt", output);
    //fs::write("true_data.txt", file_contents);
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let mut profiler = Profiler::new();
    /*
        profiler.start("profiler profiler");
        for i in 0..640000000 {
            SystemTime::now();
        }
        profiler.end("profiler profiler");
        profiler.print_results();
     */
    //return;
    profiler.start("main");

    let mut program_vars = Vars {
        bacteria_count: 0,
        bacteria_path: vec![]
    };


    let start_time = Instant::now();
    read_input_file(&args[1], &mut program_vars);

    // Analise each bacteria file and compare them
    let output = compare_all_bacteria(&mut program_vars, &mut profiler, 12);

    // Verify output is correct, used to make sure functionality isnt broken
    verify_output(output);

    profiler.end("main");

    let elapsed_time = start_time.elapsed();
    println!("time elapsed: {0} millis", elapsed_time.as_millis());
    profiler.print_results();

}

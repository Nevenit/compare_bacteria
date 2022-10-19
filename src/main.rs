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
    pub tv: Box<[f64]>,
    pub ti: Box<[i64]>
}

pub struct BacteriaCounters {
    vector: Box<[i64]>,
    second: Box<[i64]>,
    one_l: Box<[i64]>,
    indexes: i64,
    total: i64,
    total_l: i64,
    complement: i64
}

impl Bacteria {
    pub fn init(&mut self, file: &str, profiler: &mut Profiler, thread_count: usize) {
        profiler.start("read_file");
        let bacteria_file = fs::read_to_string(&file).unwrap();
        let file_bytes = bacteria_file.as_bytes();
        profiler.end("read_file");

        profiler.start("init_counters");
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
        let mut i = 0;
        while i < file_bytes.len(){
            if file_bytes[i] as char == '>' {
                profiler.start("compute_buffer");
                while file_bytes[i] as char != '\n' {
                    i += 1;
                }

                let mut buffer: [u8; LEN - 1] = [0 as u8; LEN - 1];

                let mut j = 0;
                while j < LEN - 1 {
                    i += 1;
                    buffer[j] = file_bytes[i];
                    j += 1;
                }
                profiler.end("compute_buffer");
                profiler.start("init_buffer");
                self.init_buffer(&mut bc, buffer);
                profiler.end("init_buffer");
            }
            else if file_bytes[i] as char != '\n' && file_bytes[i] != 13 {
                profiler.start("cont_buffer");
                self.cont_buffer(&mut bc, file_bytes[i] as char);
                profiler.end("cont_buffer");
            }
            i += 1;
        }
        profiler.end("extract_kmers");

        profiler.start("fill_arrays");
        let total_plus_complement: f64 = bc.total as f64 + bc.complement as f64 ;
        let total_div_2: f64 = bc.total as f64 * 0.5;

        let mut one_l_div_total =[0.0; AA_NUMBER];
        for i in 0..AA_NUMBER {
            one_l_div_total[i] = bc.one_l[i] as f64 / bc.total_l as f64;
        }

        let mut second_div_total = vec![0.0; M1].into_boxed_slice();
        for i in 0..M1 {
            second_div_total[i] = bc.second[i] as f64 / total_plus_complement;
        }

        profiler.end("fill_arrays");

        profiler.start("calculate_t");

        // Calculate the starting and ending index for each thread
        // we have to do that before starting the threads because we need this information
        // to declare the array length used for the thread output
        let mut index_bounds_arc= vec![];
        for i in 0..thread_count {
            let starting_index: usize = M / thread_count * i;

            let ending_index: usize;
            if i == thread_count - 1 {
                ending_index = M;
            } else {
                ending_index = M / thread_count * (i + 1);
            }
            index_bounds_arc.push(Arc::new((starting_index, ending_index)));
        }

        // These next few lines convert the variables needed in each thread to a thread safe pointer
        //let thread_id_arc_pointer = Arc::new(Mutex::new(0));
        let count_arc = Arc::new(Mutex::new(0));
        let one_l_div_total_arc = Arc::new(one_l_div_total);
        let second_div_total_arc = Arc::new(second_div_total);
        let bc_vector_arc = Arc::new(bc.vector);
        //let arc_index_bounds_arc = Arc::new(index_bounds.clone());
        let mut t_arc = vec![];
        for i in 0..thread_count {
            t_arc.push(Arc::new(Mutex::new(vec![0.0; index_bounds_arc[i].1 - index_bounds_arc[i].0].into_boxed_slice())));
        }

        let mut handles = vec![];

        for thread_index in 0..thread_count {
            // Create a copy of the thread safe pointer for each thread
            //let thread_id = Arc::clone(&thread_id_arc_pointer);
            let _count = Arc::clone(&count_arc);
            let _one_l_div_total = Arc::clone(&one_l_div_total_arc);
            let _second_div_total = Arc::clone(&second_div_total_arc);
            let _bc_vector = Arc::clone(&bc_vector_arc);
            let _index_bounds = Arc::clone(&index_bounds_arc[thread_index]);
            let _t = Arc::clone(&t_arc[thread_index]);

            let handle = thread::spawn(move || {
                //let mut current_thread_id = thread_id.lock().unwrap();
                //let my_id = current_thread_id.clone();
                //*current_thread_id += 1;
                //drop(current_thread_id);
                //let my_id = thread_index;

                let mut counter = 0;

                // Because only this thread is using this lock there is no performance impact
                let mut t_lock = _t.lock().unwrap();


                let starting_index: usize = _index_bounds.0;
                let mut i_mod_aa_number: usize = starting_index % AA_NUMBER;
                let mut i_div_aa_number: usize = starting_index / AA_NUMBER;
                let mut i_mod_m1: usize = starting_index % M1;
                let mut i_div_m1: usize = starting_index / M1;
                let ending_index: usize = _index_bounds.1;
                //println!("id: {}  {} {} {} {}",my_id, i_mod_aa_number, i_div_aa_number, i_mod_m1, i_div_m1);

                //println!("id: {}  start: {}  end: {}",my_id, starting_index, ending_index);
                for i in starting_index..ending_index {
                    let p1: f64 = _second_div_total[i_div_aa_number];
                    let p2: f64 = _one_l_div_total[i_mod_aa_number];
                    let p3: f64 = _second_div_total[i_mod_m1];
                    let p4: f64 = _one_l_div_total[i_div_m1];
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
                        t_lock[i - starting_index] = (_bc_vector[i] as f64 - stochastic) / stochastic;
                        counter += 1;
                    } else {
                        t_lock[i - starting_index] = 0.0;
                    }
                }

                let mut count_lock = _count.lock().unwrap();
                *count_lock += counter;
                drop(count_lock);

                //*thread_id += 1;
            });
            handles.push(handle);
        }

        // Wait for all threads to finish
        for handle in handles {
            handle.join().unwrap();
        }

        self.count = *count_arc.lock().unwrap();

        profiler.end("calculate_t");

        profiler.start("calculate_tv_and_ti");
        //let thread_id_arc = Arc::new(Mutex::new(0));

        // Because we dont know the distribution of values between threads we have to assume the worst
        // and give each thread an array with the length of count

        let mut thread_tv_arc = vec![];
        let mut thread_ti_arc = vec![];
        for i in 0..thread_count {
            thread_tv_arc.push(Arc::new(Mutex::new(vec![0.0; self.count].into_boxed_slice())));
            thread_ti_arc.push(Arc::new(Mutex::new(vec![0 as i64; self.count].into_boxed_slice())));
        }

        let mut handles = vec![];

        for thread_index in 0..thread_count {
            let _thread_t = Arc::clone(&t_arc[thread_index]);
            let _thread_tv = Arc::clone(&thread_tv_arc[thread_index]);
            let _thread_ti = Arc::clone(&thread_ti_arc[thread_index]);
            let _thread_index_bounds = Arc::clone(&index_bounds_arc[thread_index]);
            let _thread_array_length = Arc::new(self.count);

            let handle = thread::spawn(move || {
                //let my_id = thread_index;
                let t_lock = _thread_t.lock().unwrap();
                let mut tv_lock = _thread_tv.lock().unwrap();
                let mut ti_lock = _thread_ti.lock().unwrap();

                let mut pos = 0;
                for i in 0.._thread_index_bounds.1 - _thread_index_bounds.0 {
                    if t_lock[i] != 0.0 {
                        tv_lock[pos] = t_lock[i];
                        ti_lock[pos] = (_thread_index_bounds.0 + i) as i64;
                        pos += 1;
                    }
                }
            });

            handles.push(handle);
        }

        for handle in handles {
            handle.join().unwrap();
        }

        self.tv = vec![0.0; self.count].into_boxed_slice();
        self.ti = vec![0; self.count].into_boxed_slice();

        let mut pos = 0;
        for thread_index in 0..thread_count {
            let _thread_tv = Arc::clone(&thread_tv_arc[thread_index]);
            let _thread_ti = Arc::clone(&thread_ti_arc[thread_index]);
            let thread_tv_array = _thread_tv.lock().unwrap();
            let thread_ti_array = _thread_ti.lock().unwrap();
            for i in 0..self.count {
                if thread_tv_array[i] != 0.0 {
                    self.tv[pos] = thread_tv_array[i];
                    self.ti[pos] = thread_ti_array[i];
                    pos += 1;
                }
            }
        }



        /*
        self.tv = vec![0.0; self.count].into_boxed_slice();
        self.ti = vec![0; self.count].into_boxed_slice();
        let mut pos = 0;
        let mut index_counter: usize = 0;
        for i in 0..thread_count {
            let _t = t_arc[i].lock().unwrap();
            for j in 0..index_bounds_arc[i].1 - index_bounds_arc[i].0 {
                if _t[j] != 0.0 {
                    self.tv[pos] = _t[j];
                    self.ti[pos] = (index_counter + j) as i64;
                    pos += 1;
                }
            }
            index_counter += index_bounds_arc[i].1 - index_bounds_arc[i].0;
        }

         */

        /*
        for i in 0..M {
            if t[i] != 0.0 {
                self.tv[pos] = t[i];
                self.ti[pos] = i as i64;
                pos += 1;
            }
        }

         */
        profiler.end("calculate_tv_and_ti");

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

fn compare_all_bacteria(program_vars: &mut Vars, profiler: &mut Profiler) -> String {
    profiler.start("init_bacteria");
    let mut bacteria_array = vec![];

    for i in 0..program_vars.bacteria_count {
        println!("load {} of {}", i + 1, program_vars.bacteria_count);
        bacteria_array.push(Bacteria {
            count: 0,
            tv: vec![].into_boxed_slice(),
            ti: vec![].into_boxed_slice()
        });
        profiler.start(format!("init_bacteria {}", i).as_str());
        bacteria_array[i as usize].init(&program_vars.bacteria_path[i as usize], profiler, 12);
        profiler.end(format!("init_bacteria {}", i).as_str());
    }
    profiler.end("init_bacteria");

    profiler.start("compare_bacteria");
    //let mut output = String::new();
    for i in 0..program_vars.bacteria_count {
        for j in i+1..program_vars.bacteria_count {
            print!("{} {} -> ", i, j);
            //output += &format!("{} {} -> ", i, j).to_string();
            let correlation = compare_bacteria(&bacteria_array[i as usize], &bacteria_array[j as usize], profiler);
            //output += &format!("{:.20}\n", correlation).to_string();
            println!("{:.20}", correlation);
        }
    }
    profiler.end("compare_bacteria");
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
    let output = compare_all_bacteria(&mut program_vars, &mut profiler);

    // Verify output is correct, used to make sure functionality isnt broken
    verify_output(output);

    profiler.end("main");

    let elapsed_time = start_time.elapsed();
    println!("time elapsed: {0} millis", elapsed_time.as_millis());
    profiler.print_results();

}

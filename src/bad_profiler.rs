pub mod bad_profiler {
    //use std::alloc::System;
    use std::time::SystemTime;
    use colored::*;
    pub struct Profiler {
        enabled: bool,
        indent: i32,
        name: Vec<String>,
        start: Vec<SystemTime>,
        end: Vec<SystemTime>,
        elapsed: Vec<i64>,
        status: Vec<bool>,
        iterations: Vec<i64>,
        indents: Vec<i32>,
    }

    impl Profiler {
        pub fn new() -> Profiler {
            Profiler {
                enabled: true,
                indent: 0,
                name: vec![],
                start: vec![],
                elapsed: vec![],
                end: vec![],
                status: vec![],
                iterations: vec![],
                indents: vec![]
            }
        }

        pub fn start(&mut self, name: &str) {
            if !self.enabled {return};
            let index = self.get_index(name);

            match index {
                Some(i) => {
                    if !self.status[i] {
                        self.status[i] = true;
                        self.start[i] = SystemTime::now();
                    } else {
                        panic!("Profiler exception: Can't start a new iteration before the previous one has finished: {}", name)
                    }
                }
                None => {
                    self.name.push(name.to_owned());
                    self.end.push(SystemTime::UNIX_EPOCH);
                    self.elapsed.push(0);
                    self.status.push(true);
                    self.iterations.push(0);
                    self.indents.push(self.indent);
                    self.start.push(SystemTime::now());
                }
            }

            self.indent += 1;
        }

        pub fn end(&mut self, name: &str) {
            if !self.enabled {return};
            let index = self.get_index(name);
            match index {
                Some(i) => {
                    if self.status[i] {
                        self.end[i] = SystemTime::now();
                        self.elapsed[i] += self.end[i].duration_since(self.start[i]).expect("Clock may have gone backwards").as_nanos() as i64;
                        self.status[i] = false;
                        self.iterations[i] += 1;
                    } else {
                        println!("Profiler error: Can't stop a profiler that hasn't started: {}", name)
                    }
                },
                None => println!("Profiler error: Can't stop a profiler that doesn't exist: {}", name)
            }
            self.indent -= 1;
        }

        pub fn print_results(&mut self) {
            if !self.enabled {return};
            for i in 0..self.name.len(){
                let mut indent: String = String::new();
                /*if i > 0 {
                    for j in 0..i {
                        let duration =  self.end[i].duration_since(self.end[j]);
                        match duration {
                            Ok(duration) => {
                                if duration.as_nanos() != 0 {
                                    continue;
                                }
                            },
                            Err(..) => ()
                        }
                        if indent.len() == 0 {
                            indent = format!("{}{}", "└", indent);
                        } else {
                            indent = format!("{}{}", " ", indent);
                        }
                    }
                }
                */
                //indent = format!("{}{}", "└", indent);
                for _j in 0..self.indents[i] {
                    indent = indent + " ";
                }
                indent = indent + "└";

                if !self.status[i] {
                    let lindent = indent.truecolor(203, 119, 47);
                    let name = self.name[i].truecolor(255, 198, 109);
                    let elapsed = (self.elapsed[i] / 1000000).to_string().truecolor(100, 255, 100);
                    let iterations = self.iterations[i].to_string().truecolor(100, 255, 100);
                    if self.iterations[i] > 1 {
                        println!("{3}{4}{0} {5}{1} {6} {2}", ":".truecolor(203, 119, 47), "ms over".truecolor(120, 166, 100), "iterations".truecolor(120, 166, 100), lindent, name, elapsed, iterations)
                    } else {
                        println!("{2}{3}{0} {4}{1}", ":".truecolor(203, 119, 47), "ms".truecolor(120, 166, 100), lindent, name, elapsed);
                    }
                } else {
                    println!("{}{}: hasn't finished", indent, self.name[i])
                }
            }
        }

        fn get_index(&mut self, name: &str) -> Option<usize> {
            self.name.iter().position(|r| r == name)
        }
    }

}

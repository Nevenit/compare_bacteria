use criterion::{black_box, criterion_group, criterion_main, Criterion};
use std::{thread, time};

fn do_for(n: usize) {
    let mut total:i64 = 0;

    for i in 1..n {
        total += (i * 2) as i64;
    }
    thread::sleep(time::Duration::from_millis(1));
}

fn do_while(n: usize) {
    let mut total:i64 = 0;

    let mut i = 0;
    while i < n {
        total += i as i64;
        i += 1;
    }
    thread::sleep(time::Duration::from_millis(1));
}

fn test_loops(c: &mut Criterion) {
    let n = 1000000;
    let mut group = c.benchmark_group("for_vs_while");
    group.bench_function("for", |b| {
       b.iter(|| black_box(do_for(n)))
    });
    group.bench_function("while", |b| {
        b.iter(|| black_box(do_while(n)))
    });
    group.finish();
}

criterion_group!(benches, test_loops);
criterion_main!(benches);
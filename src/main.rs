use std::collections::HashMap;

use rand::{Rng, prelude::ThreadRng};

const MUTATE_RATE: f64 = 0.01;
const CROSS_RATE: f64 = 0.7;
const N_GENES: usize = 12;
const N_CHROMOS: usize = 128;

#[derive(Clone, Debug)]
struct Chromosome {
    genes: Vec<u8>,
    value: f64,
    fit: f64,
}

/// Takes in vector of chromosomes, and returns vector of n selected ones
fn roulette(chromosomes: &Vec<Chromosome>, r: &mut ThreadRng) -> Chromosome {
    let n_chromo = chromosomes.len();
    let end = r.gen_range(1..3) as f64;
    let mut ticker: f64 = 0.;
    let mut idx = 0;

    while ticker < end {
        ticker += chromosomes[idx % n_chromo].fit;
        idx += 1;
    }

    // &chromosomes[idx % n_chromo]
    chromosomes[idx % n_chromo].clone()
}

fn fitness(target: f64, res: f64) -> f64 {
    if target == res {
        return 1.
    }
    1. / (target - res)
}

// TODO: clean this garbage
fn eval_gene(genes: &Vec<u8>, gene_map: &mut HashMap<u8, char>) -> f64 {
    let e: Vec<&char> = genes.into_iter().map(|gene| gene_map.get(gene).unwrap_or(&' ')).collect();
    let mut eqn = String::from("");
    for op in e {
        if *op != ' ' {
            eqn.push(*op);
        }
    }
    // println!("'{}', {}", eqn, eqn.len());
    meval::eval_str(eqn).unwrap_or(0.)
}

fn mutate(chromo: &mut Chromosome, r: &mut ThreadRng) {
    if r.gen_bool(MUTATE_RATE) {
        let m = 1u8 << r.gen_range(0..8);
        let idx = r.gen_range(0..N_GENES);
        chromo.genes[idx] = chromo.genes[idx] ^ m;
    }
}

fn crossover(c1: &mut Chromosome, c2: &mut Chromosome, r: &mut ThreadRng) {
    if r.gen_bool(CROSS_RATE) {
        let cross = r.gen_range(0..N_GENES-1);
        // TODO: len 0?
        let g1 = c1.genes.clone();
        let g2 = c2.genes.clone();
        for i in cross..N_GENES {
            c1.genes[i] = g2[i];
            c2.genes[i] = g1[i];
        }
    }
}

fn gen_chromo(r: &mut ThreadRng) -> Chromosome {
    // try different length genes
    let mut genes = Vec::<u8>::with_capacity(N_GENES);
    for _ in 0..N_GENES {
        genes.push(r.gen_range(0..=14) as u8);
    }
    Chromosome {
        genes,
        value: 0.,
        fit: 0.
    }
    // match genes.try_fill(r) {
    //     Ok(_) => 
    //         Chromosome {
    //             genes,
    //             value: 0.,
    //             fit: 0.
    //         },
    //     Err(err) => panic!("Error generating chromosome: {}", err),
    // }
}

fn main() {
    let mut gene_map: HashMap<u8, char> = HashMap::from([
        (0u8, '0'),
        (1u8, '1'),
        (2u8, '2'),
        (3u8, '3'),
        (4u8, '4'),
        (5u8, '5'),
        (6u8, '6'),
        (7u8, '7'),
        (8u8, '8'),
        (9u8, '9'),
        (10u8, '+'),
        (11u8, '-'),
        (12u8, '*'),
        (13u8, '/'),
        (14u8, '%'),
    ]);

    let mut r = rand::thread_rng();
    let mut chromosomes = Vec::<Chromosome>::with_capacity(N_CHROMOS);
    for _ in 0..128 {
        chromosomes.push(gen_chromo(&mut r));
    }

    let mut next_chromos = Vec::<Chromosome>::with_capacity(N_CHROMOS);
    let mut change_list = Vec::<(i32, Chromosome)>::new();

    let target = 2964.37;

    // assign fitness
    for mut chromo in &mut chromosomes {
        chromo.value = eval_gene(&chromo.genes, &mut gene_map);
        chromo.fit = fitness(target, chromo.value);
    }

    let mut iters = 0;
    let mut best_so_far = chromosomes[0].clone();

    while best_so_far.fit < 0.95 {
        // one run through chromosomes
        for i in 0..(N_CHROMOS/2) {
            let mut c1 = roulette(&chromosomes, &mut r);
            let mut c2 = roulette(&chromosomes, &mut r);

            crossover(&mut c1, &mut c2, &mut r);
            mutate(&mut c1, &mut r);
            mutate(&mut c2, &mut r);

            // TODO: make closure
            let c1_val = eval_gene(&c1.genes, &mut gene_map);
            if c1_val == target {
                println!("{:?}", c1.genes);
                break;
            }

            next_chromos.push(Chromosome {
                genes: c1.genes,
                value: c1_val,
                fit: fitness(target, c1_val)
            });

            if next_chromos[i*2].fit > best_so_far.fit {
                best_so_far = next_chromos[i*2].clone();
                change_list.push((iters, best_so_far.clone()));
                println!("New best {:?}", best_so_far);
            }

            let c2_val = eval_gene(&c2.genes, &mut gene_map);
            if c2_val == target {
                println!("{:?}", c2.genes);
                break;
            }

            next_chromos.push(Chromosome {
                genes: c2.genes,
                value: c2_val,
                fit: fitness(target, c2_val)
            });

            if next_chromos[i*2+1].fit > best_so_far.fit {
                best_so_far = next_chromos[i*2+1].clone();
                change_list.push((iters, best_so_far.clone()));
                println!("New best {:?}", best_so_far);
            }
        }
        
        iters += 1;
        if iters % 100 == 0 {
            println!("Gen {}", iters);
            println!("Best so far: {:?}", best_so_far);
        }

        chromosomes = next_chromos;
        next_chromos = Vec::<Chromosome>::with_capacity(N_CHROMOS);
    }

    println!("In {} generations!", iters);
    println!("Changes:\n {:#?}", change_list);
}

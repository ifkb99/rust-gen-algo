use std::collections::HashMap;

use rand::{prelude::ThreadRng, Rng};

const MUTATE_RATE: f64 = 0.0025;
const CROSS_RATE: f64 = 0.7;
const N_GENES: usize = 12;
const N_CHROMOS: usize = 512;
const NOT_IT: f64 = -9999999999.99999;

#[derive(Clone, Debug)]
struct Chromosome {
    genes: Vec<u8>,
    value: f64,
    fit: f64,
}

/// Takes in vector of chromosomes, and returns vector of n selected ones
fn roulette(chromosomes: &Vec<Chromosome>, avg: f64, r: &mut ThreadRng) -> Chromosome {
    let end = avg * r.gen_range(1.0..2.5) * N_CHROMOS as f64;
    let mut ticker: f64 = 0.;
    let mut idx = 0;

    loop {
        ticker += chromosomes[idx % N_CHROMOS].fit;
        if ticker > end {
            break;
        }
        idx += 1;
    }

    // &chromosomes[idx % n_chromo]
    chromosomes[idx % N_CHROMOS].clone()
}

fn fitness(target: f64, res: f64) -> f64 {
    1. / (target - res).abs()
}

fn eval_gene(genes: &Vec<u8>, gene_map: &mut HashMap<u8, char>) -> f64 {
    // let eqn: String = genes
    //     .into_iter()
    //     .map(|gene| gene_map.get(gene).unwrap_or(&'#'))
    //     .filter(|gene| **gene != '#')
    //     .collect::<String>();
    let eqn: String = genes
        .into_iter()
        .map(|gene| gene_map.get(gene).unwrap_or(&'#'))
        .collect::<String>();

    // on error, returns effectively dead gene
    meval::eval_str(eqn).unwrap_or(NOT_IT)
}

fn mutate(chromo: &mut Chromosome, r: &mut ThreadRng) {
    if r.gen_bool(MUTATE_RATE) {
        let m = 1u8 << r.gen_range(0..5);
        // let m = 1u8 << r.gen_range(0..8);
        let idx = r.gen_range(0..N_GENES);
        chromo.genes[idx] = chromo.genes[idx] ^ m;
    }
}

fn crossover(c1: &mut Chromosome, c2: &mut Chromosome, r: &mut ThreadRng) {
    if r.gen_bool(CROSS_RATE) {
        let cross = r.gen_range(0..N_GENES);
        let g1 = c1.genes.clone();
        let g2 = c2.genes.clone();
        for i in cross..N_GENES {
            c1.genes[i] = g2[i];
            c2.genes[i] = g1[i];
        }
    }
}

fn gen_chromo(target: f64, r: &mut ThreadRng, gene_map: &mut HashMap<u8, char>) -> Chromosome {
    // try different length genes
    let mut genes = Vec::<u8>::with_capacity(N_GENES);
    for _ in 0..N_GENES {
        genes.push(r.gen_range(0..=14) as u8);
    }
    let value = eval_gene(&genes, gene_map);
    Chromosome {
        genes,
        value,
        fit: fitness(target, value),
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

    let target = 2976.37;

    let mut r = rand::thread_rng();
    let mut chromosomes = Vec::<Chromosome>::with_capacity(N_CHROMOS);
    for _ in 0..N_CHROMOS {
        chromosomes.push(gen_chromo(target, &mut r, &mut gene_map));
    }

    let mut next_chromos = Vec::<Chromosome>::with_capacity(N_CHROMOS);
    let mut change_list = Vec::<(i32, Chromosome)>::new();

    let mut iters = 0;
    let mut best_so_far = chromosomes
        .iter()
        .max_by_key(|x| (x.fit * 10e11) as u64)
        .unwrap()
        .clone();

    println!("Starting...");
    let mut avg = 0.001; //chromosomes.iter().fold(0., |acc, cur| acc + cur.fit) / N_CHROMOS as f64;
    let mut sum_fit = 0.;
    println!("avg: {}", avg);
    println!("bsf: {:?}", best_so_far);
    // TODO: why NaN?

    while best_so_far.fit < 0.95 {
        let mut best_now = Chromosome {
            genes: Vec::new(),
            value: NOT_IT,
            fit: NOT_IT,
        };
        // one run through chromosomes
        for i in 0..(N_CHROMOS / 2) {
            let mut c1 = roulette(&chromosomes, avg, &mut r);
            let mut c2 = roulette(&chromosomes, avg, &mut r);

            crossover(&mut c1, &mut c2, &mut r);
            mutate(&mut c1, &mut r);
            mutate(&mut c2, &mut r);

            // TODO: make closure
            let c1_val = eval_gene(&c1.genes, &mut gene_map);
            // if c1_val == target {
            //     println!("{:?}", c1.genes);
            //     break;
            // }
            let c1_fit = fitness(target, c1_val);
            next_chromos.push(Chromosome {
                genes: c1.genes,
                value: c1_val,
                fit: c1_fit,
            });
            sum_fit += c1_fit;

            if next_chromos[i * 2].fit > best_now.fit {
                best_now = next_chromos[i * 2].clone();
                // println!("Gen {} Best: {:?}", iters, best_now);
                if best_now.fit > best_so_far.fit {
                    best_so_far = best_now.clone();
                    change_list.push((iters, best_so_far.clone()));
                    println!("New best!");
                }
            }

            let c2_val = eval_gene(&c2.genes, &mut gene_map);
            // if c2_val == target {
            //     println!("{:?}", c2.genes);
            //     break;
            // }

            let c2_fit = fitness(target, c2_val);
            next_chromos.push(Chromosome {
                genes: c2.genes,
                value: c2_val,
                fit: c2_fit,
            });
            sum_fit += c2_fit;

            // if next_chromos[i * 2 + 1].fit > best_so_far.fit {
            //     best_so_far = next_chromos[i * 2 + 1].clone();
            //     change_list.push((iters, best_so_far.clone()));
            //     println!("Gen {}: New best {:?}", iters, best_so_far);
            // }
            if next_chromos[i * 2 + 1].fit > best_now.fit {
                best_now = next_chromos[i * 2 + 1].clone();
                // println!("Gen {} Best: {:?}", iters, best_now);
                if best_now.fit > best_so_far.fit {
                    best_so_far = best_now.clone();
                    change_list.push((iters, best_so_far.clone()));
                    println!("New best! {}", best_now.value);
                }
            }
        }

        if iters % 1 == 0 {
            println!("Gen {}, Avg fitness: {}", iters, avg);
            println!("Best now: {:?}", best_now);
            println!("Bst totl: {:?}", best_so_far);
        }
        iters += 1;

        avg = sum_fit / N_CHROMOS as f64;
        sum_fit = 0.;
        chromosomes = next_chromos;
        next_chromos = Vec::<Chromosome>::with_capacity(N_CHROMOS);
    }

    println!("Changes:\n{:?}", change_list);
    println!("In {} generations!", iters);
}

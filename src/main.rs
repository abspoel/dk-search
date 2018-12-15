extern crate argparse;
extern crate num;
extern crate num_bigint;
extern crate num_traits;
extern crate primal;
extern crate time;

use argparse::{ArgumentParser, Store, StoreFalse, StoreTrue};
use num::Integer;
use num_bigint::{BigInt, ToBigInt};
use std::collections::{HashMap, HashSet};
use std::fmt;
use std::ops::ShrAssign;
use time::precise_time_s;

#[derive(Clone)]
struct Sieve {
    modulus: usize,
    qr: Vec<u8>,
}

trait NumAssignOps<Rhs = Self>: num::traits::NumAssignOps<Rhs> + ShrAssign<Rhs> {}

impl<T, Rhs> NumAssignOps<Rhs> for T where T: num::traits::NumAssignOps<Rhs> + ShrAssign<Rhs> {}

trait Num<Base = Self>:
    num::traits::Num
    + NumAssignOps
    + From<u8>
    + Clone
    + Copy
    + PartialEq<Base>
    + PartialOrd<Base>
    + ToBigInt
{
}

impl<T, Base> Num<Base> for T where
    T: num::traits::Num
        + NumAssignOps
        + From<u8>
        + Clone
        + Copy
        + PartialEq<Base>
        + PartialOrd<Base>
        + ToBigInt
{
}

/// Calculate the Legendre symbol (a/p), which is returned as an element
/// in {0,1} = F_2. It is 0 iff `a` is a quadratic residue modulo `p`.
fn legendre<T>(a: T, p: T) -> u8
where
    T: Num,
{
    let e = ((p - T::from(1)) / T::from(2)).to_bigint().unwrap();
    (a.to_bigint().unwrap().modpow(&e, &(p.to_bigint().unwrap())) == BigInt::from(1u8)) as u8
}

/// Simple Eratosthenes sieve to get all primes up to bound
struct Primes {
    cur: usize,
    bound: usize,
    sieve: Vec<u8>,
}

impl Primes {
    fn new(bound: usize) -> Self {
        Primes {
            cur: 1,
            bound: bound,
            sieve: vec![0; bound],
        }
    }

    fn next_prime(&mut self, after: usize) -> Option<usize> {
        let mut search = after;
        while search < self.cur {
            search += 1;
            if self.sieve[search] == 0 {
                return Some(search);
            }
        }

        while self.cur < self.bound - 1 {
            self.cur += 1;
            if self.sieve[self.cur] == 0 {
                let mut k = 2 * self.cur;
                while k < self.bound {
                    self.sieve[k] = 1;
                    k += self.cur;
                }
                if self.cur > after {
                    return Some(self.cur);
                }
            }
        }
        None
    }
}

/// PrimesIterator allows iterating over Primes, instead of having to call Primes::next_prime()
struct PrimesIterator<'a> {
    cur: usize,
    primes: &'a mut Primes,
}

impl<'a> IntoIterator for &'a mut Primes {
    type Item = usize;
    type IntoIter = PrimesIterator<'a>;

    fn into_iter(self) -> Self::IntoIter {
        PrimesIterator {
            cur: 1,
            primes: self,
        }
    }
}

impl<'a> Iterator for PrimesIterator<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        let r = self.primes.next_prime(self.cur);
        if let Some(x) = r {
            self.cur = x;
        }
        r
    }
}

fn lcm<T>(a: T, b: T) -> Result<T, ()>
where
    T: num::CheckedMul + num::Integer,
{
    match a.checked_mul(&b) {
        Some(x) => Ok(x / a.gcd(&b)),
        None => Err(()),
    }
}

impl Sieve {
    fn new(modulus: usize) -> Sieve {
        let mut qr: Vec<u8> = Vec::with_capacity(modulus);
        let minus_one = legendre(modulus - 1, modulus);

        qr.push(0);

        for a in 1..modulus {
            qr.push((1 + minus_one + legendre(a, modulus)) % 2);
        }

        Sieve {
            qr: qr,
            modulus: modulus,
        }
    }

    fn combine(sieves: &[Sieve]) -> Result<Sieve, ()> {
        let modulus = try!(sieves
            .iter()
            .map(|s| s.modulus)
            .fold(Ok(1), |res, x| res.and_then(|acc| lcm(acc, x))));
        let mut qr: Vec<u8> = Vec::with_capacity(modulus);
        for a in 0..modulus {
            qr.push(sieves.iter().map(|s| s.qr[a % s.modulus]).all(|x| x == 1) as u8);
        }
        Ok(Sieve {
            modulus: modulus,
            qr: qr,
        })
    }

    fn take_first(
        sorted: &mut Vec<Sieve>,
        group: &mut Vec<Sieve>,
        max_size: &mut i64,
    ) -> Result<(), ()> {
        if sorted.is_empty() {
            Err(())
        } else {
            let modulus = sorted[0].modulus as i64;
            if modulus <= *max_size {
                group.push(sorted.remove(0));
                *max_size /= modulus;
                Ok(())
            } else {
                Err(())
            }
        }
    }

    fn combine_all(sieves: Vec<Sieve>, bound: usize) -> Vec<Sieve> {
        // This assumes all sieves are relatively prime
        let mut sorted = sieves;
        //sorted.sort_unstable_by_key(|s| (s.modulus as i64));
        let mut combined: Vec<Sieve> = Vec::new();
        let mut group: Vec<Sieve> = Vec::new();

        let mut max_size: i64;

        loop {
            if sorted.len() == 0 {
                break;
            }

            let sieve = sorted.remove(0);
            max_size = (bound / sieve.modulus) as i64;
            group.push(sieve);

            loop {
                if let Err(_) = Sieve::take_first(&mut sorted, &mut group, &mut max_size) {
                    break;
                }
            }

            combined.push(Sieve::combine(&group).unwrap());
            group.clear();
        }
        combined
    }

    fn from_pair(pair: RelatedPair) -> Self {
        let (p, q) = pair;
        let modulus = p * q;
        let s_p = Sieve::new(p);
        let s_q = Sieve::new(q);
        let mut qr = Vec::with_capacity(modulus);
        let mut a_p = 0;
        let mut a_q = 0;
        for _ in 0..modulus {
            qr.push(s_p[a_p] | s_q[a_q]);
            a_p += 1;
            if (a_p) >= p {
                a_p -= p;
            }
            a_q += 1;
            if (a_q) >= q {
                a_q -= q;
            }
        }
        Sieve {
            qr: qr,
            modulus: modulus,
        }
    }

    #[allow(dead_code)]
    fn from_triple(triple: RelatedTriple) -> Self {
        let (p, q, r) = triple;
        let modulus = p * q * r;
        let s_p = Sieve::new(p);
        let s_q = Sieve::new(q);
        let s_r = Sieve::new(r);
        let mut qr = Vec::with_capacity(modulus);
        let mut a_p = 0;
        let mut a_q = 0;
        let mut a_r = 0;
        for _ in 0..modulus {
            qr.push(s_p[a_p] | s_q[a_q] | s_r[a_r]);
            a_p += 1;
            if (a_p) >= p {
                a_p -= p;
            }
            a_q += 1;
            if (a_q) >= q {
                a_q -= q;
            }
            a_r += 1;
            if (a_r) >= r {
                a_r -= r;
            }
        }
        Sieve {
            qr: qr,
            modulus: modulus,
        }
    }

    fn from_pairs(pairs: Vec<RelatedPair>) -> Result<Self, ()> {
        let sieves = pairs
            .iter()
            .map(|pair| Sieve::from_pair(*pair))
            .collect::<Vec<_>>();
        Sieve::combine(&sieves)
    }
}

// TODO: make this into HashSet
type RelatedPair = (usize, usize);

fn get_related_pairs(d0: usize, d1: usize, dk_star: bool) -> HashSet<RelatedPair> {
    assert!(
        d0 * d0 >= d1,
        "For this function to work, d0^2 should be greater than d1"
    );
    assert!(d0 >= 4, "d0 should be greater than max(difference)");

    let mut primes = Primes::new(d1);
    let mut pairs = HashSet::new();

    let mut p = d0;

    let differences = if dk_star {
        [-2i64, 2, -4, 4]
    } else {
        [-1i64, 1, -2, 2]
    };

    loop {
        let mut marked = vec![0u8; d1 + 5]; // d1+1 + max(difference)
        if let Some(prime) = primes.next_prime(p) {
            p = prime;
        } else {
            break;
        }

        // Mark multiples of p : d0 < p <= d1
        {
            let mut k = 1;
            while k * p <= d1 {
                marked[k * p] = 1;
                k += if dk_star { 2 } else { 1 };
            }
        }

        // Test whether we have k*q+s is marked, for d0 < q < p < d1
        let mut q = d0;
        loop {
            if let Some(prime) = primes.next_prime(q) {
                q = prime;
            } else {
                break;
            }
            if q >= p {
                break;
            }

            let mut k = 1;
            while k * q <= d1 {
                for s in &differences {
                    if marked[((k * q) as i64 + s) as usize] != 0 {
                        // println!("Add edge: ({}, {}) since ?*{} == {}*{}+{}", q, p, p, k, q, s);
                        pairs.insert((q, p));
                    }
                }
                k += if dk_star { 2 } else { 1 };
            }
        }
    }
    pairs
}

#[allow(dead_code)]
type RelatedTriple = (usize, usize, usize);

#[allow(dead_code)]
fn get_related_triples(d0: usize, d2: usize, dk_star: bool) -> HashSet<RelatedTriple> {
    assert!(
        d0 * d0 >= d2,
        "For this function to work, d0^2 should be greater than d2"
    );
    assert!(d0 >= 6, "d0 should be greater than max(difference)");

    let mut primes = Primes::new(d2);
    let mut triples = HashSet::new();

    let mut p = d0;

    let differences = if dk_star {
        [-2i64, 2, -4, 4, -6, 6, -8, 8]
    } else {
        [-1i64, 1, -2, 2, -3, 3, -4, 4]
    };

    let star_factor = if dk_star { 2 } else { 1 };

    loop {
        let mut marked = vec![0u8; d2 + 9]; // d2+1 + max(s)
        if let Some(prime) = primes.next_prime(p) {
            p = prime;
        } else {
            break;
        }

        // Mark multiples of p : d0 < p <= d2
        {
            let mut k = 1;
            while k * p <= d2 {
                marked[k * p] = 1;
                k += if dk_star { 2 } else { 1 };
            }
        }

        // Test whether we have k*q+s is marked, for d0 < r < q < p < d1
        let mut q = d0;
        loop {
            if let Some(prime) = primes.next_prime(q) {
                q = prime;
            } else {
                break;
            }
            if q >= p {
                break;
            }

            let mut r = d0;
            loop {
                if let Some(prime) = primes.next_prime(r) {
                    r = prime;
                } else {
                    break;
                }
                if r >= q {
                    break;
                }

                let mut k = 1;
                while k * r <= d2 {
                    for s1 in &differences {
                        if marked[((k * r) as i64 + s1) as usize] != 0 {
                            for s2 in &differences {
                                if (s2 - s1).abs() <= star_factor * 4
                                    && marked[((k * r) as i64 + s2) as usize] != 0
                                {
                                    triples.insert((r, q, p));
                                }
                            }
                        }
                    }
                    k += if dk_star { 2 } else { 1 };
                }
            }
        }
    }
    triples
}

/// Get a list of combined sieves (at large as possible, bounded in size by `sieve_bound`) for
/// the related prime pairs `pairs`.
fn related_pair_sieves(pairs: HashSet<RelatedPair>, sieve_bound: usize) -> Vec<Sieve> {
    // To construct a list of sieves as compact as possible, we see the primes of the related pairs
    // as the vertices of a graph, with the edges between all related primes. We want to partition
    // the edges of this graph such that the set of primes that are incident to each edge in the
    // part is minimal.
    let mut combined: Vec<Sieve> = Vec::new();
    let mut edges = pairs;

    let mut degree: HashMap<usize, usize> = HashMap::new();

    for edge in &edges {
        for vertex in &[edge.0, edge.1] {
            *degree.entry(*vertex).or_insert(0) += 1;
        }
    }

    while !edges.is_empty() {
        let mut group_edges: Vec<RelatedPair> = Vec::new();
        let mut group_vertices: HashSet<usize> = HashSet::new();
        while !edges.is_empty() {
            let group_size;
            match group_vertices
                .iter()
                .fold(Ok(1usize), |res, x| res.and_then(|acc| lcm(acc, *x)))
            {
                Ok(x) => {
                    group_size = x;
                }
                Err(_) => {
                    break;
                }
            }

            let max_size = sieve_bound / group_size;
            // The next edge is chosen such that:
            // 1. the resulting group_size does not exceed sieve_bound
            // 2. we prioritize edges that have incident vertices already in our group
            // 3. there are edges remaining that have common vertices incident to the chosen edge
            // 4. the growth of group_size is small
            let edge = edges
                .iter()
                .max_by_key(|&&edge| {
                    let (v1, v2) = edge;
                    let has_v1 = group_vertices.contains(&v1) as usize;
                    let has_v2 = group_vertices.contains(&v2) as usize;
                    let mut size = 1;
                    if has_v1 == 0 {
                        size *= v1;
                    }
                    if has_v2 == 0 {
                        size *= v2;
                    }
                    let size_score = if size > max_size { 0 } else { -(size as i64) };
                    let eligible = 1 - ((size > max_size) as i8);
                    let vertex_score = (has_v1 + has_v2) as i64;
                    let v1_degree = degree.get(&v1).unwrap();
                    let v2_degree = degree.get(&v2).unwrap();
                    let max_degree = v1_degree.max(v2_degree);
                    let min_degree = v1_degree.min(v2_degree);
                    let score = (
                        eligible as i8,
                        vertex_score,
                        *max_degree,
                        *min_degree,
                        size_score,
                    );
                    score
                })
                .unwrap()
                .clone();
            let (v1, v2) = edge;
            let mut size = 1;
            if !group_vertices.contains(&v1) {
                size *= v1;
            }
            if !group_vertices.contains(&v2) {
                size *= v2;
            }
            if size > max_size {
                break;
            }
            edges.take(&edge);

            // Update `degree` map
            *degree.entry(edge.0).or_insert(0) -= 1;
            *degree.entry(edge.1).or_insert(0) -= 1;

            group_vertices.insert(edge.0);
            group_vertices.insert(edge.1);
            group_edges.push(edge);
        }

        // let modulus = group_vertices.iter().fold(1, |acc, x| lcm(acc, *x).unwrap());
        let sieve = Sieve::from_pairs(group_edges).unwrap();
        combined.push(sieve);
    }

    combined
}

#[allow(dead_code)]
fn related_triple_sieves(triples: HashSet<RelatedTriple>, sieve_bound: usize) -> Vec<Sieve> {
    let mut combined: Vec<Sieve> = Vec::new();

    for triple in &triples {
        let modulus = triple.0 * triple.1 * triple.2;
        assert!(
            modulus <= sieve_bound,
            "Modulus {} is greater than sieve_bound {}",
            modulus,
            sieve_bound
        );
        let sieve = Sieve::from_triple(*triple);
        combined.push(sieve);
    }

    combined
}

impl std::ops::Index<usize> for Sieve {
    type Output = u8;

    fn index(&self, index: usize) -> &u8 {
        &self.qr[index]
    }
}

impl fmt::Display for Sieve {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        try!(write!(f, "Sieve(modulus={}, ", self.modulus));
        for a in &self.qr {
            try!(write!(f, "{}", a));
        }
        write!(f, ")")
    }
}

struct JumpList {
    modulus: usize,
    first: usize,
    jl: Vec<usize>,
}

impl From<Sieve> for JumpList {
    fn from(sieve: Sieve) -> Self {
        let mut a = 0usize;
        let modulus = sieve.modulus;
        let mut jl: Vec<usize> = Vec::new();
        while a < modulus {
            if sieve[a] == 1 {
                break;
            }
            a += 1;
        }
        let first = a;
        let mut prev = a;
        a += 1;
        while a < modulus {
            if sieve[a] == 1 {
                jl.push(a - prev);
                prev = a;
            }
            a += 1;
        }

        jl.push(first + modulus - prev);
        JumpList {
            modulus: modulus,
            first: first,
            jl: jl,
        }
    }
}

impl fmt::Display for JumpList {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        try!(write!(
            f,
            "Sieve(modulus={}, first={}, [",
            self.modulus, self.first
        ));
        for a in &self.jl {
            try!(write!(f, "{}, ", a));
        }
        write!(f, "])")
    }
}

impl JumpList {
    /// Return the first integer matching the jump list after `start`,
    /// together with the next jump index.
    fn first_after(&self, start: u64) -> (u64, usize) {
        let mut a = start - (start % (self.modulus as u64)) + (self.first as u64);
        for index in 0..self.jl.len() {
            if a >= start {
                return (a, index);
            }
            a += self.jl[index] as u64;
        }
        if a < start {
            if self.jl.len() == 0 {
                panic!("JumpList is empty");
            } else {
                panic!("Unexpected behaviour");
            }
        }
        (a, 0)
    }
}

struct Wheel {
    modulus: usize,
    qr: Vec<u8>,
    jump: Vec<u64>,
    jump_to_ix: Vec<usize>,
    consecutive: Vec<u64>,
}

fn to_wheel(modulus: usize, qr: &Vec<u8>, step: u64) -> Result<Wheel, ()> {
    if step.gcd(&(modulus as u64)) != 1 {
        return Err(());
    }
    let mut jumps: Vec<u64> = Vec::with_capacity(modulus);
    let mut index = (step % (modulus as u64)) as usize;
    let mut consecutive: Vec<u64> = Vec::new();
    let mut jump_to_ix: Vec<usize> = Vec::new();
    let mut cs_map = vec![modulus; modulus];
    let step_add = index;

    let mut start_index = 0;
    while qr[start_index] == 0 {
        start_index += 1;
    }

    {
        let mut index = start_index;
        let mut jump;
        let mut first = true;
        while (index != start_index) || first {
            first = false;
            jump = step;
            index += step_add;
            if index >= modulus {
                index -= modulus;
            }
            while qr[index] == 0 {
                jump += step;
                index += step_add;
                if index >= modulus {
                    index -= modulus;
                }
            }
            consecutive.push(jump);
            cs_map[index] = consecutive.len() - 1;
        }
    }
    assert_eq!(
        consecutive.len(),
        qr.iter().fold(0usize, |sum, val| sum + (*val as usize))
    );

    for i in 0..(modulus as u64) {
        let mut jump = step;
        // This should terminate since gcd(step, modulus) == 1, and we're
        // assuming at least one qr[i] != 0
        index = ((i + step) % (modulus as u64)) as usize;
        while qr[index] == 0 {
            jump += step;
            index += step_add;
            if index >= modulus {
                index -= modulus;
            }
        }
        jumps.push(jump);
        assert_ne!(cs_map[index], modulus);
        jump_to_ix.push(cs_map[index]);
    }

    Ok(Wheel {
        modulus: modulus,
        qr: qr.to_vec(),
        jump: jumps,
        jump_to_ix: jump_to_ix,
        consecutive: consecutive,
    })
}

/// Convert sieves into a list of wheels, and possibly remaining sieves.
fn to_wheels(
    mut sieves: Vec<Sieve>,
    zeroth_step: u64,
    step_bound: u64,
) -> (Vec<Wheel>, Vec<Sieve>) {
    let mut step = zeroth_step;
    let mut wheels: Vec<Wheel> = Vec::with_capacity(sieves.len());

    let mut ix = 0;
    while ix < sieves.len() {
        let modulus = sieves[ix].modulus;
        if let Ok(wheel) = to_wheel(modulus, &sieves[ix].qr, step) {
            wheels.push(wheel);
            sieves.remove(ix);

            match step.checked_mul(modulus as u64) {
                Some(x) => {
                    if x < step_bound {
                        step = x
                    } else {
                        break;
                    }
                }
                None => break,
            }
        } else {
            ix += 1;
        }
    }
    (wheels, sieves)
}

fn get_dk(p: u64, k: usize, lower_bound: usize, dk_star: bool) -> u64 {
    let mut window: Vec<i32> = Vec::new();
    for i in ((lower_bound as i64) + 1) - (k as i64)..((lower_bound as i64) + 2) + (k as i64) {
        if dk_star {
            window.push(if legendre(2 * (i as i128) + 1, p as i128) == 1 {
                1
            } else {
                -1
            });
        } else {
            window.push(if legendre(i as i128, p as i128) == 1 {
                1
            } else {
                -1
            });
        }
    }
    let mut d = (lower_bound as u64) + 1;
    while window.iter().fold(0i32, |sum, val| sum + *val) > 0 {
        d += 1;
        window.remove(0);
        if dk_star {
            let pushval = if legendre(2 * (d + (k as u64)) + 1, p) == 1 {
                1
            } else {
                -1
            };
            window.push(pushval);
        } else {
            window.push(if legendre(d + (k as u64), p) == 1 {
                1
            } else {
                -1
            });
        }
    }
    d - 1
}

struct HighDkPrime {
    p: u64,
    df: u64,
    dk: u64,
}

fn output_num(
    x: u64,
    start_time: f64,
    dk_star: bool,
    bests: &mut Vec<HighDkPrime>,
    progress: f64,
    k: usize,
    df_min: usize,
    dk_min: usize,
) {
    if primal::is_prime(x) {
        let elapsed = precise_time_s() - start_time;
        let star = if dk_star { "*" } else { "" };
        let bits = (x as f64).log2().ceil();
        let f = k - 1;
        let df = get_dk(x, f, df_min, dk_star);
        let dk = get_dk(x, k, dk_min, dk_star);
        let mut insert;
        let mut ix = 0usize;

        // bests should for each (z0, z1) contain the minimal p such that d0(p) >= z0, dk(p) >= z1
        if bests.is_empty() {
            insert = true;
        } else if x < bests[0].p {
            insert = true;
        } else if let Err(found_ix) = bests.binary_search_by_key(&x, |hdp| hdp.p) {
            insert = true;
            for best in bests[..found_ix].iter() {
                if best.df >= df && best.dk >= dk {
                    insert = false;
                }
            }
            ix = found_ix;
        } else {
            // Element already in there
            return;
        }

        if insert {
            bests.retain(|hdp| hdp.p < x || hdp.df > df || hdp.dk > dk);
            bests.insert(
                ix,
                HighDkPrime {
                    p: x,
                    df: df,
                    dk: dk,
                },
            );
            println!(
                "[{:.3} s elapsed, {:.3}%] {} ({:.0} bits), d{}{}={}, d{}{}={}",
                elapsed,
                progress * 100.0,
                x,
                bits,
                f,
                star,
                df,
                k,
                star,
                dk
            );
        }
    }
}

fn search(
    jumplist: &JumpList,
    sieves: &[Sieve],
    start: u64,
    stop: u64,
    dk_star: bool,
    k: usize,
    df_min: usize,
    dk_min: usize,
) -> Vec<HighDkPrime> {
    let (mut p, mut i) = jumplist.first_after(start);
    let len = jumplist.jl.len();
    let start_time = precise_time_s();
    let mut bests: Vec<HighDkPrime> = Vec::new();
    while p < stop {
        let mut b = 1;
        for sieve in sieves {
            if sieve.qr[(p % sieve.modulus as u64) as usize] == 0 {
                b = 0;
                break;
            }
        }

        if b == 1 {
            let progress = (p - start) as f64 / (stop - start) as f64;
            output_num(
                p, start_time, dk_star, &mut bests, progress, k, df_min, dk_min,
            );
        }

        p += jumplist.jl[i] as u64;
        i += 1;
        if i >= len {
            i -= len;
        }
    }
    let elapsed = precise_time_s() - start_time;
    println!("=== Done in {:.3} s ===", elapsed);
    bests
}

fn search_with_wheels(
    wheels: &[Wheel],
    sieves: &[Sieve],
    start: u64,
    stop: u64,
    dk_star: bool,
    k: usize,
    df_min: usize,
    dk_min: usize,
) -> Vec<HighDkPrime> {
    let len = wheels.len();
    let mut x = vec![0u64; len];
    let mut bests: Vec<HighDkPrime> = Vec::new();

    // This is the number of times a jump has occurred in the i-th wheel. We start counting at 1:
    // the zeroth jump is the initial one (which is skipped if qr[index] == 1). We terminate and
    // reduce i (except for the last wheel, so i==len-1) when this counter reaches
    // sum(wheels[i].qr[j] for j in 0..wheels[i].modulus), which is equal to
    // wheels[i].consecutive.len().
    let mut cnt = vec![0usize; len];

    let mut cs_index = vec![0usize; len];
    let mut is_indexed = vec![false; len];
    let mut i = 0;
    x[0] = start;
    let mut index = (x[i] % (wheels[i].modulus as u64)) as usize;

    let mut do_jump = wheels[i].qr[index] == 0;
    cnt[0] = wheels[i].qr[index] as usize;

    let start_time = precise_time_s();
    let mut overflow = false;
    loop {
        if do_jump {
            if !is_indexed[i] {
                index = (x[i] % (wheels[i].modulus as u64)) as usize;

                // Rust should fix tuple assignment
                let (xi, of) = x[i].overflowing_add(wheels[i].jump[index]);
                x[i] = xi;
                overflow = of;

                cs_index[i] = wheels[i].jump_to_ix[index];
                is_indexed[i] = true;
            } else {
                cs_index[i] += 1;
                if cs_index[i] >= wheels[i].consecutive.len() {
                    cs_index[i] -= wheels[i].consecutive.len();
                }

                // Rust should fix tuple assignment
                let (xi, of) = x[i].overflowing_add(wheels[i].consecutive[cs_index[i]]);
                x[i] = xi;
                overflow = of;
            }

            cnt[i] += 1;
        }

        if do_jump
            && (((x[i] >= stop) || (i < len - 1 && cnt[i] > wheels[i].consecutive.len()))
                || overflow)
        {
            overflow = false;
            if i > 0 {
                i -= 1;
                do_jump = true;
            } else {
                break;
            }
        } else {
            if i < len - 1 {
                i += 1;
                x[i] = x[i - 1];

                index = (x[i] % (wheels[i].modulus as u64)) as usize;
                cnt[i] = wheels[i].qr[index] as usize;
                do_jump = wheels[i].qr[index] == 0;
                is_indexed[i] = false;
            } else {
                // check all sieves
                let mut b = 1;
                for sieve in sieves {
                    if sieve.qr[(x[i] % sieve.modulus as u64) as usize] == 0 {
                        b = 0;
                        break;
                    }
                }

                if b == 1 {
                    let mut progress = 0.0;
                    let mut fraction = 1.0;
                    for j in 0..len - 1 {
                        progress +=
                            fraction * ((cnt[j] - 1) as f64 / wheels[j].consecutive.len() as f64);
                        fraction /= wheels[j].consecutive.len() as f64;
                    }
                    output_num(
                        x[i], start_time, dk_star, &mut bests, progress, k, df_min, dk_min,
                    );
                }
                do_jump = true;
            }
        }
    }
    let elapsed = precise_time_s() - start_time;
    println!("=== Done in {:.3} s ===", elapsed);
    bests
}

/// Set jumplist_bound to 0 to disable JumpList
fn preprocess(
    d0: usize,
    d1: usize,
    jumplist_bound: usize,
    sieve_bound: usize,
    dk_star: bool,
) -> (Option<JumpList>, Vec<Sieve>) {
    let primes_bound = if dk_star { 2 * d0 + 1 } else { d0 };
    let fixed_sieve = if dk_star {
        Sieve {
            modulus: 4,
            qr: vec![0, 0, 0, 1],
        }
    } else {
        Sieve {
            modulus: 8,
            qr: vec![0, 0, 0, 0, 0, 0, 0, 1],
        }
    };
    // println!("Fixed sieve={}", fixed_sieve);

    let mut primes = Primes::new(primes_bound);
    let iter = PrimesIterator {
        cur: 2,
        primes: &mut primes,
    };
    let prime_sieves: Vec<Sieve> = iter.map(|p| Sieve::new(p)).collect();

    let mut sorted = prime_sieves;
    // This should already be sorted by ascending modulus. We want the sparsest sieves first,
    // since then each jump (what we do most) is as far as possible.

    let mut max_size = if jumplist_bound > 0 {
        (jumplist_bound / fixed_sieve.modulus) as i64
    } else {
        (sieve_bound / fixed_sieve.modulus) as i64
    };

    let mut group = vec![fixed_sieve];
    loop {
        if let Err(_) = Sieve::take_first(&mut sorted, &mut group, &mut max_size) {
            break;
        }
    }

    let first_sieve = Sieve::combine(&group).unwrap();
    let mut combined: Vec<Sieve> = vec![];
    let jumplist;
    if jumplist_bound > 0 {
        jumplist = Some(JumpList::from(first_sieve));
    } else {
        jumplist = None;
        combined.push(first_sieve);
    };

    combined.append(&mut Sieve::combine_all(sorted, sieve_bound));

    let d1_bound = if dk_star { 2 * d1 + 1 } else { d1 };
    let pairs = get_related_pairs(primes_bound, d1_bound, dk_star);
    let mut pair_sieves = related_pair_sieves(pairs, sieve_bound);
    combined.append(&mut pair_sieves);

    (jumplist, combined)
}

fn main() {
    let mut d0 = 0usize;
    let mut d1 = 0usize;
    let mut k = 1usize;
    let mut jumplist_bound = 5000000usize;
    let mut sieve_bound = 4000000usize;
    let mut start = 1u64;
    let mut stop = 0u64;
    let mut with_wheel = true;
    let mut dk_star = false;
    {
        let mut parser = ArgumentParser::new();
        parser.refer(&mut d0).required().add_argument(
            "d0",
            Store,
            "Minimum value of d_0 to sieve for",
        );
        parser.refer(&mut d1).required().add_argument(
            "d1",
            Store,
            "Minimum value of d_1 to sieve for",
        );
        parser.refer(&mut stop).required().add_argument(
            "stop",
            Store,
            "Maximum value of p to search for",
        );
        parser.refer(&mut k).add_option(
            &["-k"],
            Store,
            "The parameter k, either 1 or 2 (default: k=1)",
        );
        parser.refer(&mut jumplist_bound).add_option(
            &["--jumplist-bound"],
            Store,
            "Maximum size (modulus) of the jumplist. Set to 0 to disable.",
        );
        parser.refer(&mut sieve_bound).add_option(
            &["--sieve-bound"],
            Store,
            "Maximum size (modulus) of each sieve",
        );
        parser.refer(&mut with_wheel).add_option(
            &["--no-wheel"],
            StoreFalse,
            "Use traditional sieves and a jumplist instead of wheels",
        );
        parser.refer(&mut start).add_option(
            &["--start"],
            Store,
            "Value of p to start searching at",
        );
        parser.refer(&mut dk_star).add_option(
            &["-s", "--dk-star"],
            StoreTrue,
            "Find d_k^* instead of d_k, i.e. only look at odd indices for the Legendre symbol",
        );
        parser.parse_args_or_exit();
    }

    if k != 1 && k != 2 {
        println!("k should be either 1 or 2");
        std::process::exit(1);
    }

    if with_wheel {
        jumplist_bound = 0;
    };
    println!("Preprocessing...");
    let (jumplist, sieves): (Option<JumpList>, Vec<Sieve>) =
        preprocess(d0, d1, jumplist_bound, sieve_bound, dk_star);

    for sieve in &sieves {
        println!("Sieve has modulus {}", sieve.modulus);
    }

    let df_min;
    let dk_min;
    // f = k-1
    if k == 1 {
        df_min = d0;
        dk_min = if d1 >= d0 { d1 } else { d0 };
    } else {
        df_min = d1;
        dk_min = d1;
    }

    let bests;
    if with_wheel {
        let wheel_start = precise_time_s();
        println!("Instantiating wheels...");
        let (wheels, sieves) = to_wheels(sieves, 1u64, stop - start);
        println!(
            "Preprocessing done in {:.1} seconds",
            precise_time_s() - wheel_start
        );
        println!("{} wheels, {} sieves", wheels.len(), sieves.len());
        for wheel in &wheels {
            println!("Wheel modulus: {}", wheel.modulus);
        }

        println!("Searching with wheels...");
        bests = search_with_wheels(&wheels, &sieves, start, stop, dk_star, k, df_min, dk_min);
    } else {
        println!("Searching with {} sieves and a jumplist...", sieves.len());
        bests = search(
            &jumplist.unwrap(),
            &sieves,
            start,
            stop,
            dk_star,
            k,
            df_min,
            dk_min,
        );
    }

    for hdp in &bests {
        let bits = (hdp.p as f64).log2().ceil();
        let star = if dk_star { "*" } else { "" };
        let f = k - 1;
        println!(
            "{} ({:.0} bits), d{}{}={}, d{}{}={}",
            hdp.p, bits, f, star, hdp.df, k, star, hdp.dk
        );
    }
}

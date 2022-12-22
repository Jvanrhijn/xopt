use rand::Rng;


fn generate_positions<R: Rng>(n: usize, rng: &mut R) -> Vec<Vec<f64>> {
    let mut out = vec![];
    for _ in 0..n {
        out.push(vec![rng.gen(), rng.gen()]);
    }
    out
}


fn subtr(vec_a: &Vec<f64>, vec_b: &Vec<f64>) -> Vec<f64> {
    vec_a.iter().zip(vec_b).map(|(a, b)| a - b).collect()
}

fn norm(x: &Vec<f64>) -> f64 {
    f64::sqrt(x.iter().map(|a| a*a).sum())
}



fn ccw(a: &Vec<f64>, b: &Vec<f64>, c: &Vec<f64>) -> bool {
    (c[1] - a[1]) * (b[0] - a[0]) > (b[1] - a[1]) * (c[0] - a[0])
}

fn intersect(a: &Vec<f64>, b: &Vec<f64>, c: &Vec<f64>, d: &Vec<f64>) -> bool {
    ccw(a, c, d) != ccw(b, c, d) && ccw(a, b, c) != ccw(a, b, d)
}

fn length(tour: &Vec<usize>, xs: &Vec<Vec<f64>>) -> f64 {
    let mut l = 0.0;
    for (i, j) in tour.iter().zip(&tour[1..]) {
        l += norm(&subtr(&xs[*i], &xs[*j]))
    }
    l + norm(&subtr(&xs[0], &xs[xs.len()-1]))
}

fn find_tangle(tour: &Vec<usize>, x: &Vec<Vec<f64>>) -> Option<(usize, usize)> {
    let n = tour.len();
    for (v, i) in tour.iter().enumerate() {
        for (u, k) in tour[i+1..].iter().enumerate() {
            let j = k + i;
            let e = [tour[i % n], tour[(i+1) % n]];
            let f = [tour[j % n], tour[(j+1) % n]];
            if e == f {
                continue
            }
            if intersect(&x[e[0]], &x[e[1]], &x[f[0]], &x[f[1]]) {
                let idxe = i % n;
                let idxf = j % n;
                if f.contains(&e[0]) || f.contains(&e[1]) {
                    continue
                }
                return Some((idxe, idxf))
            }
        }
    }
    None
}


fn untangle_iteration(i: usize, j: usize, tour: &Vec<usize>) -> Vec<usize> {
    let n = tour.len();
    let (i1, i2) = (usize::min(i, j), usize::max(i, j));
    let mut p1: Vec<_> = tour[..i1].into();
    let mut p2: Vec<_> = tour[i1+2..i2].into();
    p2.reverse();
    let mut out: Vec<usize> = vec![];
    if i2 == n-1 {
        out.append(&mut p1[1..].into());
        out.push(tour[i1]);
        out.push(tour[i2]);
        out.append(&mut p2);
        out.push(tour[i1+1]);
        out.push(tour[(i2+1) % n]);
    } else if i2 == n-2 {
        out.append(&mut p1);
        out.push(tour[i1]);
        out.push(tour[i2]);
        out.append(&mut p2);
        out.push(tour[i1+1]);
        out.push(tour[(i2+1) % n]);
    } else {
        let mut p3: Vec<usize> = tour[i2+2..].into();
        out.append(&mut p1);
        out.push(tour[i1]); out.push(tour[i2]);
        out.append(&mut p2);
        out.push(tour[i1+1]); out.push(tour[i2+1]);
        out.append(&mut p3);
    }

    return out;
}


fn untangle_tour(mut tour: Vec<usize>, x: &Vec<Vec<f64>>, iter_max: usize) -> (Vec<usize>, f64) {
    for _ in 0..iter_max {
        if let Some((i, j)) = find_tangle(&tour, x) {
            let new_tour = untangle_iteration(i, j, &tour);
            tour = new_tour;
        } else {
            break;
        }
    }
    let l = length(&tour, &x);
    (tour, l)
}


fn main() {
    let mut ns: Vec<_> = (100..=1000).step_by(100).collect();
    let mut ns2: Vec<_> = (10000..=10000).step_by(1000).collect();
    ns.append(&mut ns2);
    let nsamples: usize = 1000;

    let mut rng = rand::thread_rng();


    for n in ns {
        let mut lav = 0.0;
        for _ in 0..nsamples {
            let x = generate_positions(n, &mut rng);
            let tour: Vec<usize> = (0..n).collect();
            let (_, length) = untangle_tour(tour, &x, 1000000);
            lav += length / nsamples as f64;
        }
        println!("{} {:.10}", n, lav);
    }
}
use primal::is_prime;
fn main() {
    let n =4;
    let a = vec![5,6,7,8]; // z.B. G(x) = 1 + 2x + 3x^2 + 4x^3
    let q = 7681;
    let res = find_primitive_nth_roots(n,q);
    if is_ntt_compatible(n, q) {
        let transformed = ntt_based_on_omega(&a, res[0], q);
        println!("NTT(a) = {:?}", transformed);
    }
    println!("omega = {:?}", res);
    if let Some(inv) = inveres_erw_euklid(23, 17) {
        println!("Das Inverse ist: {}", inv);
    } else {
        println!("Kein Inverses vorhanden");
    }

}
//q muss primzahl sein und n muss ein teiler von q-1 sein weil es sonst keine einheitswurzel gibt (Satz von euler)
fn is_ntt_compatible(n: u64, q: u64) -> bool {
    is_prime(q) && (q - 1) % n == 0
}



//seite 5 ntt for beginners
fn ntt_based_on_omega(a: &[u64], omega: u64, q: u64) -> Vec<u64> { // O(n²)
    let n = a.len();
    let mut result = vec![0u64; n];

    for j in 0..n {
        let mut sum = 0u64;
        for i in 0..n {
            let power = mod_pow(omega, (i * j) as u64, q);
            sum = (sum + a[i] * power) % q;
        }
        result[j] = sum;
    }

    result
}
//Square & multiply algo
fn mod_pow(mut base: u64, mut exp: u64, modulus: u64) -> u64 {
    let mut result = 1;
    base = base % modulus;//mann darf vorher schon die basis kürzen damit man später leichter rechnen kann
    while exp > 0 {
        if exp % 2 == 1 {
            result = result * base % modulus;
        }
        base = base * base % modulus;
        exp /= 2;
    }
    result
}
//seite 4 von ntt for beginners
fn is_primitive_root(root: u64, n: u64, q: u64) -> bool {
    if mod_pow(root, n, q) != 1 {
        return false;
    }
    for k in 1..n {
        if mod_pow(root, k, q) == 1 {
            return false;
        }
    }
    true
}

fn find_nth_roots_of_unity(n: u64, q: u64) -> Vec<u64> {
    let mut roots = Vec::new();
    for omega in 2..q {//omega 0,1 sind Trivialfälle ,das wegen starten wir bei der 2 
        if mod_pow(omega, n, q) == 1 {
            roots.push(omega);
        }
    }
    roots
}

fn find_primitive_nth_roots(n: u64, q: u64) -> Vec<u64> {
    let all_roots = find_nth_roots_of_unity(n, q);
    all_roots
        .into_iter()
        .filter(|&omega| is_primitive_root(omega, n, q))
        .collect()
}
fn inveres_erw_euklid(a: u64, b: u64) -> Option<u64> {
    let mut t = 0;
    let mut new_t =1;
    let mut r = b;
    let mut new_r = a;
    while new_r!=0{
        let q = r/ new_r;
        let temp_t = t - q * new_t;
        t = new_t;
        new_t = temp_t;
        let temp_r = r -q * new_r;
        r = new_r;
        new_r = temp_r;
    }
    if r >1 {
        return None;
    }

    Some(t)
}
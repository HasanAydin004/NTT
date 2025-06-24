use primal::is_prime;

fn main() {
    let inputs: Vec<i128> = vec![7, 3, 5, 1, 4, 9, 2, 6]; // Beispiel-Vektor (Länge muss Potenz von 2 sein)

    // Schritt 1: Bestimme ein geeignetes Modulus
    let max_val = *inputs.iter().max().unwrap();
    let len = inputs.len();
    let min_mod = max_val * max_val * len as i128 + 1;
    let modulus = find_modulus(len, min_mod);

    // Schritt 2: Finde eine primitive n-te Einheitswurzel modulo modulus
    let root = find_primitive_root(len as i128, modulus - 1, modulus);

    println!("Input vector:      {:?}", inputs);
    println!("Modulus:           {}", modulus);
    println!("Primitive root:    {}", root);

    // Schritt 3: NTT berechnen
    let ntt = ntt_recursive(&inputs, root, modulus);
    println!("NTT result:        {:?}", ntt);

    // Schritt 4: INTT berechnen
    let intt = intt_recursive(&ntt, root, modulus);
    println!("Reconstructed:     {:?}", intt);

    // Schritt 5: Vergleich
    if intt == inputs {
        println!("Test bestanden: INTT entspricht dem Original!");
    } else {
        println!("Test fehlgeschlagen!");
    }
    // Kyber beispiel
    let inputs: Vec<i128> = (0..256).collect(); // Testdaten: [0, 1, 2, ..., 255]

    // Parameter aus Kyber
    let modulus: i128 = 3329;
    let root: i128 = 17;

    println!("Input vector (Kyber-like):   {:?}", &inputs[..]);
    println!("Modulus (Kyber):             {}", modulus);
    println!("Primitive root (Kyber):      {}", root);

    // NTT anwenden
    let ntt = ntt_recursive(&inputs, root, modulus);
    println!("NTT result (first ):        {:?}", &ntt[..]);
    
    // INTT anwenden
    let intt = intt_recursive(&ntt, root, modulus);
    println!("Reconstructed (first ):     {:?}", &intt[..]);

    if intt == inputs {
        println!("Test bestanden: INTT entspricht dem Original!");
    } else {
        println!("Test fehlgeschlagen!");
    }
}


// Computes the forward Number-Theoretic Transform (NTT) of the input vector `invec`,
// using the provided primitive n-th root of unity `root`, modulo `modulus`.
fn ntt_recursive(invec: &Vec<i128>, root: i128, modulus: i128) -> Vec<i128> {
    let n = invec.len();
    if n == 1 {
        return vec![invec[0] % modulus];
    }

    let half = n / 2;

    // Split into even and odd indices
    let even: Vec<i128> = invec.iter().step_by(2).cloned().collect();
    let odd: Vec<i128> = invec.iter().skip(1).step_by(2).cloned().collect();

    // Square the root to pass down recursively
    let root_squared = pow_mod(root, 2, modulus);
    let even_ntt = ntt_recursive(&even, root_squared, modulus);
    let odd_ntt = ntt_recursive(&odd, root_squared, modulus);

    // Combine step
    let mut out = vec![0; n];
    let mut omega = 1;

    for i in 0..half {
        let t = (omega * odd_ntt[i]) % modulus;
        out[i] = (even_ntt[i] + t + modulus) % modulus;
        out[i + half] = (even_ntt[i] - t + modulus) % modulus;
        omega = (omega * root) % modulus;
    }

    out
}

// Computes the inverse Number-Theoretic Transform (INTT) of the input vector `invec`,
// using the provided primitive n-th root of unity `root`, modulo `modulus`.
//
// The inverse is computed by applying a forward NTT with the reciprocal root
// and then multiplying each coefficient by the modular inverse of the vector length.
fn intt_recursive(invec: &Vec<i128>, root: i128, modulus: i128) -> Vec<i128> {
    let n = invec.len() as i128;
    let inv_root = reciprocal_mod(root, modulus);
    let mut outvec = ntt_recursive(invec, inv_root, modulus);

    let scaler = reciprocal_mod(n, modulus);
    for i in 0..outvec.len() {
        outvec[i] = mul(outvec[i], scaler, modulus);
    }

    outvec
}



// Gibt die kleinste Primzahl `mod` zurück, sodass:
//   mod = i * len + 1
// für ein i ≥ 1 gilt und mod > len ist.
//
// Die Funktion erhöht i schrittweise, bis sie eine solche Primzahl findet.
// Auch wenn die Schleife theoretisch lange laufen und große Zahlen erzeugen kann,
// garantiert Dirichlets Satz über arithmetische Folgen,
// dass es unendlich viele Primzahlen dieser Form gibt.
// Daher wird die Funktion immer ein Ergebnis liefern.
fn find_modulus(len: usize, min: i128) -> i128 {
    let len_i128 = len as i128;

    // Finde das kleinste k, sodass k * len + 1 >= min
    let mut k = (min - 1 + len_i128 - 1) / len_i128;
    if k < 1 {
        k = 1;
    }

    // Berechne den Startwert für n
    let mut candidate = k * len_i128 + 1;

    // Suche die nächste Primzahl, die die Form k * len + 1 hat
    loop {
        if is_prime(candidate as u64) {
            return candidate;
        }
        candidate += len_i128;
    }
}

// Returns a primitive `degree`-th root of unity modulo `modulus`.
// Requires that `totient` (Euler's totient of modulus) is divisible by `degree`.
// If `modulus` is prime, such a root is guaranteed to exist.
fn find_primitive_root(degree: i128, totient: i128, modulus: i128) -> i128 {
    let generator = find_generator(totient, modulus);
    pow_mod(generator, totient / degree, modulus)
}

/// Finds a generator (primitive root) of the multiplicative group modulo `modulus`.
/// Assumes `totient` = φ(modulus). For prime `modulus`, such a generator always exists.
fn find_generator(totient: i128, modulus: i128) -> i128 {
    for candidate in 1..modulus {
        if is_primitive_root(candidate, totient, modulus) {
            return candidate;
        }
    }
    panic!("No primitive root found. Are you sure the modulus is correct?");
}

/// Checks whether `val` is a primitive root modulo `modulus`.
/// That means: val^totient ≡ 1 (mod modulus), but not for any smaller exponent dividing totient.
fn is_primitive_root(val: i128, totient: i128, modulus: i128) -> bool {
    // Must satisfy Fermat's little theorem
    if pow_mod(val, totient, modulus) != 1 {
        return false;
    }

    // Must NOT satisfy for any smaller divisor totient/p (p prime)
    for p in unique_prime_factors(totient) {
        if pow_mod(val, totient / p, modulus) == 1 {
            return false;
        }
    }

    true
}

// Returns the unique prime factors of `n` in ascending order.
// Example: unique_prime_factors(60) = [2, 3, 5]
fn unique_prime_factors(mut n: i128) -> Vec<i128> {
    let mut result = Vec::new();
    let mut i = 2;

    while i * i <= n {
        if n % i == 0 {
            result.push(i);
            while n % i == 0 {
                n /= i;
            }
        }
        i += 1;
    }

    if n > 1 {
        result.push(n);
    }

    result
}
fn pow_mod(base: i128, exp: i128, modulus: i128) -> i128 {
    let mut result = 1;
    let mut base = base % modulus;//mann darf vorher schon die basis kürzen damit man später leichter rechnen kann
    let mut exp = exp;
    while exp > 0 {
        if exp % 2 == 1 {
            result = result * base % modulus;
        }
        base = base * base % modulus;
        exp /= 2;
    }
    result
}
fn add(lhs: i128, rhs: i128, modulo: i128) -> i128 {
    ((lhs % modulo + rhs % modulo) + modulo) % modulo
}

fn sub(lhs: i128, rhs: i128, modulo: i128) -> i128 {
    ((lhs % modulo - rhs % modulo) + modulo) % modulo
}

fn mul(lhs: i128, rhs: i128, modulo: i128) -> i128 {
    ((lhs % modulo) * (rhs % modulo)) % modulo
}

// Computes the modular inverse of `a` modulo `modulus`, i.e., x such that (a * x) % modulus == 1.
//Panics if no inverse exists.
pub fn reciprocal_mod(a: i128, modulus: i128) -> i128 {
    let (g, x, _) = extended_gcd(a, modulus);
    if g != 1 {
        panic!("No modular inverse exists for {} mod {}", a, modulus);
    }
    (x % modulus + modulus) % modulus
}

//Extended Euclidean Algorithm: returns gcd(a, b) and Bézout coefficients (x, y)
fn extended_gcd(a: i128, b: i128) -> (i128, i128, i128) {
    if b == 0 {
        (a, 1, 0)
    } else {
        let (g, x1, y1) = extended_gcd(b, a % b);
        (g, y1, x1 - (a / b) * y1)
    }
}

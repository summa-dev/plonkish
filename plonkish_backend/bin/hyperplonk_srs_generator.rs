use halo2_curves::{pairing::Engine, serde::SerdeObject};
use itertools::Itertools;
use plonkish_backend::{
    halo2_curves::{
        bn256::Bn256,
        group::ff::Field,
    },
    util::{
        arithmetic::{batch_projective_to_affine, fixed_base_msm, window_size, window_table},
        izip,
        parallel::parallelize,
    },
};
use std::{
    env,
    iter,
    fs::File,
    io::Write,
};
use rand::rngs::OsRng;

// Some of code and logic are referenced from `https://github.com/han0110/halo2-kzg-srs`
fn main() {
    let dst_prefix = env::args()
        .nth(1)
        .expect("Please specify destination file path prefix (will be appended with suffix k)");
    let desired_k = env::args().nth(2).and_then(|s| s.parse::<u32>().ok()).expect("Please specify the number of K");

    // Generate destination file
    //
    // The logic is referenced from the `src/pcs/multilinear/kzg.rs` file
    let num_vars = desired_k as usize;
    let ss: Vec<<Bn256 as Engine>::Scalar> = iter::repeat_with(|| <Bn256 as Engine>::Scalar::random(OsRng))
            .take(num_vars)
            .collect_vec();

    let g1 = <Bn256 as Engine>::G1Affine::generator();
    let eqs = {
        let mut eqs = Vec::with_capacity(1 << (num_vars + 1));
        eqs.push(vec![<Bn256 as Engine>::Scalar::ONE]);

        for s_i in ss.iter() {
            let last_evals = eqs.last().unwrap();
            let mut evals = vec![<Bn256 as Engine>::Scalar::ZERO; 2 * last_evals.len()];

            let (evals_lo, evals_hi) = evals.split_at_mut(last_evals.len());

            parallelize(evals_hi, |(evals_hi, start)| {
                izip!(evals_hi, &last_evals[start..])
                    .for_each(|(eval_hi, last_eval)| *eval_hi = *s_i * last_eval);
            });
            parallelize(evals_lo, |(evals_lo, start)| {
                izip!(evals_lo, &evals_hi[start..], &last_evals[start..])
                    .for_each(|(eval_lo, eval_hi, last_eval)| *eval_lo = *last_eval - eval_hi);
            });

            eqs.push(evals)
        }

        let window_size = window_size((2 << num_vars) - 2);
        let window_table = window_table(window_size, g1);

        let mut eqs: Vec<<Bn256 as Engine>::G1Affine> =
            batch_projective_to_affine(&fixed_base_msm(
                window_size,
                &window_table,
                eqs.iter().flat_map(|evals| evals.iter()),
            ));

        let eqs = &mut eqs.drain(..);
        (0..num_vars + 1)
            .map(move |idx| eqs.take(1 << idx).collect_vec())
            .collect_vec()
    };

    let g2 = <Bn256 as Engine>::G2Affine::generator();
    let ss: Vec<<Bn256 as Engine>::G2Affine> = {
        let window_size = window_size(num_vars as usize);
        let window_table = window_table(window_size, <Bn256 as Engine>::G2Affine::generator());
        batch_projective_to_affine(&fixed_base_msm(window_size, &window_table, &ss))
    };

    // Exports the SRS to a file
    let mut writer = File::create(format!("{}{}", dst_prefix, num_vars)).unwrap();
    writer.write_all(&desired_k.to_le_bytes()).unwrap();
    g1.write_raw(&mut writer).unwrap();
    // Flatten eqs and writes them to the file
    for e in eqs.iter().flat_map(|e| e.iter()) {
        e.write_raw(&mut writer).unwrap();
    }
    g2.write_raw(&mut writer).unwrap();
    for s in ss.iter() {
        s.write_raw(&mut writer).unwrap();
    }

    println!("SRS generated successfully");
}

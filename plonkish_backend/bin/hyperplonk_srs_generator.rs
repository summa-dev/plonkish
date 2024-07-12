use byteorder::{LittleEndian, ReadBytesExt};
use halo2_curves::{pairing::Engine, serde::SerdeObject};
use itertools::Itertools;
use plonkish_backend::{
    halo2_curves::{
        bn256::{Bn256},
        group::ff::{Field, PrimeField},
    },
    util::{
        arithmetic::{batch_projective_to_affine, fixed_base_msm, window_size, window_table},
        izip,
        parallel::parallelize,
    },
};
use std::{
    env,
    fs::File,
    io::{self, Read, Seek, SeekFrom, Write},
};

// Some of code and logic are referenced from `https://github.com/han0110/halo2-kzg-srs`
pub const HEADER_SIZE_OFFSET: u64 = 16;
pub const HEADER_OFFSET: u64 = HEADER_SIZE_OFFSET + 8;

fn main() {
    // Read environment variables
    let src = env::args()
        .nth(1)
        .expect("Please specify source file path to convert");
    let dst_prefix = env::args()
        .nth(2)
        .expect("Please specify destination file path prefix (will be appended with suffix k)");
    let desired_k = env::args().nth(3).and_then(|s| s.parse::<u32>().ok());
    
    // Read randomness from source file
    //
    // We are using `g1` points, which are elliptic curve points in group 1, as a random seed in the HyperPlonk SRS.
    // If need more details about the source file, please refer to https://github.com/iden3/snarkjs#7-prepare-phase-2
    let reader = &mut File::open(src).expect("Couldn't open file");
    reader.seek(SeekFrom::Start(HEADER_SIZE_OFFSET)).unwrap();
    let header_size = reader.read_u64::<LittleEndian>().unwrap();

    let k_offset = HEADER_OFFSET + header_size - 8;
    reader.seek(io::SeekFrom::Start(k_offset)).unwrap();
    let k = reader.read_u32::<LittleEndian>().unwrap();

    // The size of the source, the srs file, should be sufficient even if it's smaller than the `desired_k` value.
    // becuase the srs file for HyperPlonk backend can be generated with `k` random points.
    let n = 1 << k;
    let num_vars = desired_k.unwrap();
    assert!(
        n >= num_vars,
        "{}",
        format!("Source file size not enough for desired k: {:?}", desired_k)
    );

    // G1 offset
    let _ = reader.seek(SeekFrom::Start(HEADER_OFFSET + header_size + 12));

    let mut reprs = vec![[0u8; 32]; 2 * n as usize];
    for repr in reprs.iter_mut() {
        reader.read_exact(repr.as_mut()).unwrap();
    }

    let mut random_source = Vec::with_capacity(num_vars as usize);

    for i in 0..(2 * num_vars) {
        // Update `random_source` at even indices, which represent, which is x coordinate of G1
        let repr = reprs[i as usize];
        if i % 2 == 0 {
            random_source.push(<Bn256 as Engine>::Scalar::from_repr(repr).unwrap());
        }
    }

    // Generate destination file
    //
   // The logic is referenced from the `src/pcs/multilinear/kzg.rs` 
    let g1 = <Bn256 as Engine>::G1Affine::generator();
    let eqs = {
        let mut eqs = Vec::with_capacity(1 << (num_vars + 1));
        eqs.push(vec![<Bn256 as Engine>::Scalar::ONE]);

        for s_i in random_source.iter() {
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
        batch_projective_to_affine(&fixed_base_msm(window_size, &window_table, &random_source))
    };

    // Exports the SRS to a file
    let mut writer = File::create(format!("{}{}", dst_prefix, num_vars)).unwrap();
    writer.write_all(&num_vars.to_le_bytes()).unwrap();
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

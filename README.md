This code can be used to search for primes p that have large values of d_k(p) and d_k^\*(p), as defined in our paper:

> Mark Abspoel and Niek J. Bouman and Berry Schoenmakers and Niels de Vreede. Fast Secure Comparison for Medium-Sized Integers and Its Application in Binarized Neural Networks. Cryptographers Track RSA Conference, 2019. To appear.

To use, install Rust using [rustup](https://rustup.rs/). Clone the repository and use `cargo run --release -- [args]` to execute.

Two helper scripts `run.py` and `table.py` are provided to start multiple processes to speed up the computation, and aggregate the resulting log files, respectively.

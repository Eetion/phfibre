
# How to install

* Install Rust
* Clone this repo and `SOLAR`
* This repo contains a file called `Cargo.toml`.  In your cloned copy, update the file paths for `SOLAR` so that Rust will know where to look in your computer for this code.

# How to run

* To run a calculation, you will typically write a program file, save it to the `bin` folder, and then run it via the following command (it is **very important** to include the `--release` at the end; if you exclude it the program will run in debug mode, which is much slower):

```
cd path/to/this/repo
cargo run --bin file_name --release
```

* **Example** Try the following for a simplicial complex:

```
cd path/to/this/repo
cargo run --bin B_m_n --release
```

and the following for a CW complex:

```
cd path/to/this/repo
cargo run --bin T_2 --release
```

# Notes for the developers

## Unit tests

* tests have been written in tandem with the code
* to check correctness of the main (top-level) algorithm, we
    * checked results with computations from the original paper on PH fibre (this is the `B_m_n.rs` file)
    * added code to most of the binary files to automatically check that each output facet (i) engenders the correct barcode, and (ii) is distinct from all other output facets.  therefore the list of output polyhedra for each input should at least form a subset of the actual whole
    * applied two different "permutation tests" to a number of examples; these are designed to check that the algorithm returns equivalent results no matter how we permute the vertices of the underlying simplicial complex (which is to be filtered)


## To-d0

- add test to check whether a given cycle must have infinite or finite lifespan
- try exporting to Poly-make
- convert all references to BiMapSequential to BiMapSequential



## Examples (try soon)

- delta and CW complex representations of 
    - the torus
    - rp2
    - klein bottle
    - circle with n edges

## Examples (try later)

- homotopy collapse
- performance for antitranspose calculation
- reproduce examples from theory paper
- tree + sphere
- variants for a single space
    * how do fibres for different barcodes compare? 
    * effect of gluing two complexes?
    * is fibre null-homologous?
- restrict to lower star filtrations
- (a seemingly reliable option) report some statistics for parametrized or random family of space

## Statistics to consider

- betti numbers
- number of cells in each dimension
- computation time / memory
- is it a manifold?

## Paper

- consider moving away from notation \sigma_b
- don't pre-assign a class to be critical versus non-critical
- merge barcode search
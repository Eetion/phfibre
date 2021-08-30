
# How to install

* Install Rust
* Clone this repo and `SOLAR`
* This repo contains a file called `Cargo.toml`.  In your cloned copy, update the file paths for `SOLAR` so that Rust will know where to look in your computer for this code.

# How to run

* **IMPORTANT** Functions that take combinatorial simplices (in the form of vectors of inteters) assume that the entries of these vectorsa appear in **sorted order**.

* To run a calculation, you will typically write a program file, save it to the `bin` folder, and then run it via the following command (it is **very important** to include the `--release` at the end; if you exclude it the program will run in debug mode, which is much slower):

```
cd path/to/this/repo
cargo run --bin file_name --release
```

* **Example** Try the following examples:

The 1-skeleton of the 2-simplex (from paper by Leygonie and Tillmann)

```
cd path/to/this/repo
cargo run --bin B_m_n --release
```

A CW complex:

```
cd path/to/this/repo
cargo run --bin T_2 --release
```

Lower-star and lower-edge filtrations over the the unit interval, subdivided into 5 1-simplces and 6 0-simplices

```
cd path/to/this/repo
cargo run --bin I_m --release
```

# Features of this library

* aribtrary boundary matrices with order constraints
* lower-star and lower-edge filtrations
* short-circuit optimizations using previously stored information

#  Important lessons from implementation

* search order matters in this depth-first paradigm; say that you sit at a given node of the search tree, and must add some positive and negative cells; since we look for facets of the polyhedral fibre complex, we generally want our level sets to be as small as possible; exploring all "death additions" first fills our stockpile of results with things that we can use to short-circuit explorations of later branches; by contrast, adding births first essentially guarantees that we build all the largest level sets (thus ths minimal polytopes that are less useful for short-circuiting future branches)
    * i've found that placing death-additions ahead of birth-additions is **synergistic** with methods that use prior results to short-circuit branch exploration; in particular, for the simplex_d3 examples, using the `certifies_prior_exploration` short-circuiting method alone achieved almost no speed-up (and in fact seemed to slow things a little, from 12sec to 13sec); just placing death-additions ahead of birth-additions dropped computation time to 6.75 seconds; but *both* prioritizing death-additions *and* using using `certifies_prior_exploration` dropped run time down to 3.75sec


# Notes for the developers

## Unit tests

* tests have been written in tandem with the code
* to check correctness of the main (top-level) algorithm, we
    * checked results with computations from the original paper on PH fibre (this is the `B_m_n.rs` file)
    * added code to most of the binary files to automatically check that each output facet (i) engenders the correct barcode, and (ii) is distinct from all other output facets.  therefore the list of output polyhedra for each input should at least form a subset of the actual whole
    * applied two different "permutation tests" to a number of examples; these are designed to check that the algorithm returns equivalent results no matter how we permute the vertices of the underlying simplicial complex (which is to be filtered)


## To-do

- turn B_m_n.rs into a unit test
- typically the number of vertices is very small, compared to the number of cells in the complex; one could either try to enumerate vertices first, then either use this information to dramatically cut down the search space for higher dimensinoal polytopes, or attempt to evaluate the truth/falsehood of the statement (this set of vertices is contained within a polytope) directly, in order to construct the dowker dual to the nerve complex
- parallelization
- add test to check whether a given cycle must have infinite or finite lifespan
- try exporting to Poly-make
- convert all references to BiMapSequential to BiMapSequential
- (NOW, UPON REFLECTION, I THINK WE HAVE TO BE CAREFUL ABOUT WHETHER THIS IS TRULY WELL-FOUNDED, MATHEMATICALLY) ADD A SCREENER FOR WHETHER A NEW POSTIVE CELL SHOULD HAVE FINITE OR INFINITE LIFE



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
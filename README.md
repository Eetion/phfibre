Welcome!  This is the git repository for the paper  [Algorithmic Reconstruction of the Fiber of Persistent Homology on Cell Complexes](https://arxiv.org/abs/2110.14676).  Here you can find source code for the experiments reported in that work.


# How to install


* Install [Rust](https://www.rust-lang.org)
* Obtain a copy of [this repository](https://github.com/Eetion/phfibre.git)
* Obtain a copy of the version of the `SOLAR` repository *specifically associated with this repository*.  It is important to use the correct version.  There are two ways to do this.
  * obtain a copy from `https://github.com/ExHACT/SOLAR.git`.  Then run the following in a command shell to check out the correct version of the repo

  ```
  cd path/to/the/solar/repo
  git checkout 241343ee54bfe63624eb5a07b95da97f55832603
  ```
  * alternatively, you can find a copy of the correct version of  `SOLAR` in the present repository, under `dependencies/SOLAR` (it is a git repo with the main branch checked out; that is the correct branch and version)).
* The present repo contains a file called `Cargo.toml`.  In your cloned copy, update the file path for `SOLAR` so that Rust will know where to find this code.


# How to run


* **IMPORTANT** If you write your own code, keep in mind that functions which take combinatorial simplices (in the form of vectors of inteters) as input typically depend on the assumption that simplies are represented as vectors whose entries are **sorted in ascending order**.

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

# Experiments reported in the paper


Each experiment reported in the paper has a corresponding source file in the current repository, under `srce/bin/`.  For example, 

- the interval subdivided into several 0- and 1-simplies (`I_m.rs`)
- the torus (`SimplicialTorus.rs`)
- rp2 (`SimplicialProjectivePlane.rs`)
- klein bottle (`SimplicialKlein`)
- B_m_n, examples from Appendix A of Leygonie, Tillmann: *The Fiber of Persistent Homology for simplicial complexes* (`B_m_n.rs`)
- circle with n edges


# Features of this library


* aribtrary boundary matrices 
  * one can specify any boundary matrix via the `boundary` field of the `Node` object defined in `phfibre.rs`
  * however, there are utility functions for generating the boundary matrices of simplicial complexes; see the examples files in `src/bin/` for illustration

* special order constraints on the addition of cells
  * one can impose arbitrary contraints of form "before adding cell A, each cell in the set S must be added" via the `cell_id_to_prereqs` field of the `Node` object defined in `phfibre.rs`

* lower-star and lower-edge filtrations
  * One can restrict the calculation to points in the fiber which are lower-star or lower-edge filtrations by passing a certain struct to the main solver.  See `phfibre.rs` for further details.

* short-circuit optimizations using previously stored information
  * We use several heuristics to short-circuit branches of the search tree which are gauranteed not to add new cells to the fibre.


#  Important lessons learned from this implementation


* search order matters in this depth-first paradigm; say that you sit at a given node of the search tree, and must add some positive and negative cells; since we look for facets of the polyhedral fibre complex, we generally want our level sets to be as small as possible; exploring all "death additions" first fills our stockpile of results with things that we can use to short-circuit explorations of later branches; by contrast, adding births first essentially guarantees that we build all the largest level sets (thus ths minimal polytopes that are less useful for short-circuiting future branches)
    * i've found that placing death-additions ahead of birth-additions is **synergistic** with methods that use prior results to short-circuit branch exploration; in particular, for the simplex_d3 examples, using the `certifies_prior_exploration` short-circuiting method alone achieved almost no speed-up (and in fact seemed to slow things a little, from 12sec to 13sec); just placing death-additions ahead of birth-additions dropped computation time to 6.75 seconds; but *both* prioritizing death-additions *and* using using `certifies_prior_exploration` dropped run time down to 3.75sec


# Unit tests


* tests have been written in tandem with the code
* to check correctness of the main (top-level) algorithm, we
    * checked results with computations from the original paper on PH fibre (this is the `B_m_n.rs` file)
    * added code to most of the binary files to automatically check that each output facet (i) engenders the correct barcode, and (ii) is distinct from all other output facets.  therefore the list of output polyhedra for each input should at least form a subset of the actual whole
    * applied two different "permutation tests" to a number of examples; these are designed to check that the algorithm returns equivalent results no matter how we permute the vertices of the underlying simplicial complex (which is to be filtered)


# Notes for the developers

See the file `README.md` in the branch `dev_greg` for notes which may be relevant to researchers in this area and code developers.
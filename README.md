





# To-d0

* try exporting to Poly-make



# Debug / verification

* write function to compute barcode, given equiv classes
* write a simple program to test every permutation
    + generate all permutations
    + for each permutation:
        + reduce the boundary matrix
        + cycle over (all partitions of n) x (all selections of k critical classes)
            + if resulting barcode matches, keep this leaf


# Examples

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
- # cells in each dimension
- computation time / memory
- is it a manifold?

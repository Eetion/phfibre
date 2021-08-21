//  DESCRIPTION

//  A square decomposed into two 2-simplices.


//  DEPENDENCIES

use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use phfibre::pipelines::simplex_pipeline;
use solar::utilities::sequences_and_ordinals::ordinate_unique_vals;
use solar::cell_complexes::simplices_unweighted::facets::ordered_subsimplices_up_thru_dim_concatenated_vec; 

fn main() {

    //  ----------------------------------------------------------------------------------------------    
    //  NO FINITE BARS
    //  ----------------------------------------------------------------------------------------------

    //  Define the base space, barcode, and ring

    let complex_facets      =   vec![  vec![0,1,2], vec![1,2,3] ];

    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 2); 

    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0} ],
                                    fin: vec![ ],
                                    ordinal: ordinate_unique_vals( & vec![ 0 ] ),
                                };

    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

    simplex_pipeline(
        &   simplex_sequence,
        &   barcode,
        &   ring
    );

    //  RESULTS
    //  -------

    // BASE SPACE
    // simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [1, 2], [1, 3], [2, 3], [0, 1, 2], [1, 2, 3]]
    // POLYHEDRAL COMPLEX FACETS
    // number of facets (total): 160
    // number of facets (binned by dimension): [0, 0, 0, 0, 0, 160]
    // number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 160]
    // POLYHEDRAL COMPLEX CELLS
    // number of polytopes (total): 2731
    // number of polytopes (binned by dimension): [32, 241, 702, 964, 632, 160]
    // betti numbers (of polyhedral complex): [1, 0, 0, 0, 0, 0]
    // INTERSECTIONS
    // number of pairs of intersecting facets, binned by the dimension of the intersection polytope: [5472, 3872, 2120, 928, 328, 0]
    // DOWKER NERVE COMPLEX FACETS
    // number of dowker nerve complex facets (total): 160
    // number of dowker nerve complex facets (binned by dimension): [0, 0, 0, 0, 0, 160]
    // DOWKER NERVE COMPLEX CELLS
    // number of nerve dowker complex cells (total): 2731
    // number of nerve dowker complex cells (binned by dimension): [32, 241, 702, 964, 632, 160]
    // betti numbers (of dowker nerve): [1, 0, 0, 0, 0, 0]
}
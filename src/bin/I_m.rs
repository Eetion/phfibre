//  DESCRIPTION

//  The closed interavl, decomposed into 1-simplices.


//  DEPENDENCIES

use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use phfibre::pipelines::simplex_pipeline;
use solar::utilities::sequences_and_ordinals::ordinate_unique_vals;
use solar::cell_complexes::simplices_unweighted::facets::ordered_subsimplices_up_thru_dim_concatenated_vec; 
use std::iter::FromIterator;

//  CODE

fn main() {

    //  ----------------------------------------------------------------------------------------------    
    //  5-VERTEX INTERVAL WITH NO FINITE BARS
    //  ----------------------------------------------------------------------------------------------


    //  Define the base space, barcode, and ring

    let complex_facets      =   Vec::from_iter(     // this command is equivalent to writing vec![ vec![0,1], vec![1,2], vec![2,3], vec![3,4] ]
                                    ( 0 .. 5 ).map(|x| vec![x, x+1] )
                                ); 
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 1);    // build all faces up to dimension 1

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

    // POLYHEDRAL COMPLEX FACETS
    // number of facets (total): 32
    // number of facets (binned by dimension): [0, 0, 0, 0, 0, 32]
    // POLYHEDRAL COMPLEX CELLS
    // number of polytopes (total): 791
    // number of polytopes (binned by dimension): [21, 105, 231, 258, 144, 32]
    // betti numbers (of polyhedral complex): [1, 0, 0, 0, 0, 0]
    // INTERSECTIONS
    // number of pairs of intersecting facets, binned by the dimension of the intersection polytope: [126, 126, 112, 84, 48, 0]
    // DOWKER NERVE COMPLEX FACETS
    // number of dowker nerve complex facets (total): 791
    // number of dowker nerve complex facets (binned by dimension): [21, 105, 231, 258, 144, 32]
    // DOWKER NERVE COMPLEX CELLS
    // number of nerve dowker complex cells (total): 791
    // number of nerve dowker complex cells (binned by dimension): [21, 105, 231, 258, 144, 32]
    // betti numbers (of dowker nerve): [1, 0, 0, 0, 0, 0]  

    

    //  ----------------------------------------------------------------------------------------------    
    //  5-VERTEX INTERVAL WITH ONE FINITE BAR
    //  ----------------------------------------------------------------------------------------------


    //  Define the base space, barcode, and ring

    let complex_facets      =   Vec::from_iter(     // this command is equivalent to writing vec![ vec![0,1], vec![1,2], vec![2,3], vec![3,4] ]
                                    ( 0 .. 5 ).map(|x| vec![x, x+1] )
                                ); 
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 1);    // build all faces up to dimension 1

    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0} ],
                                    fin: vec![ BarFinite{dim:0,birth:1,death:2} ],
                                    ordinal: ordinate_unique_vals( & vec![ 0, 1, 2 ] ),
                                };

    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

    simplex_pipeline(
        &   simplex_sequence,
        &   barcode,
        &   ring
    );

    //  RESULTS
    //  -------

    // POLYHEDRAL COMPLEX FACETS
    // number of facets (total): 1824
    // number of facets (binned by dimension): [0, 0, 0, 0, 1824]
    // POLYHEDRAL COMPLEX CELLS
    // number of polytopes (total): 22054
    // number of polytopes (binned by dimension): [990, 4674, 8214, 6352, 1824]
    // betti numbers (of polyhedral complex): [2, 0, 0, 0, 0]
    // INTERSECTIONS
    // number of pairs of intersecting facets, binned by the dimension of the intersection polytope: [38870, 27344, 12548, 4272, 0]
    // DOWKER NERVE COMPLEX FACETS
    // number of dowker nerve complex facets (total): 1180834
    // number of dowker nerve complex facets (binned by dimension): [990, 13654, 64756, 160186, 248672, 269712, 216320, 129712, 56608, 16896, 3072, 256]
    // DOWKER NERVE COMPLEX CELLS
    // number of nerve dowker complex cells (total): 1180834
    // number of nerve dowker complex cells (binned by dimension): [990, 13654, 64756, 160186, 248672, 269712, 216320, 129712, 56608, 16896, 3072, 256]
    // betti numbers (of dowker nerve): [2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  


}
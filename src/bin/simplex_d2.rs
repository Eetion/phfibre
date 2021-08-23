//  DESCRIPTION

//  The 2-simplex.


//  DEPENDENCIES

use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use phfibre::pipelines::simplex_pipeline;
use phfibre::phfibre::{NoExtraCondition, LowerStarCondition, LowerEdgeCondition};
use solar::utilities::sequences_and_ordinals::ordinate_unique_vals;
use solar::cell_complexes::simplices_unweighted::facets::ordered_subsimplices_up_thru_dim_concatenated_vec; 
use std::iter::FromIterator;

//  CODE

fn main() {

    //  ----------------------------------------------------------------------------------------------    
    //  NO FINITE BARS
    //  ----------------------------------------------------------------------------------------------


    //  Define the base space, barcode, and ring

    let complex_facets      =   vec![  vec![0,1,2] ];
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 1);    // build all faces up to dimension 1

    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:1} ],
                                    fin: vec![ ],
                                    ordinal: ordinate_unique_vals( & vec![ 0, 1 ] ),
                                };

    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

    let precondition_to_make_new_lev_set_lower_none     =   NoExtraCondition{};      // this struct won't impose any extra conditions on the filtrations we build    

    simplex_pipeline(
        &   simplex_sequence,
        &   barcode,
        &   ring,
        &   precondition_to_make_new_lev_set_lower_none,
    );

    //  RESULTS
    //  -------

    // BASE SPACE
    // simplices of the base space: [[0], [1], [2], [0, 1], [0, 2], [1, 2]]
    // POLYHEDRAL COMPLEX FACETS
    // number of facets (total): 12
    // number of facets (binned by dimension): [0, 0, 12]
    // number of facets (binned by number of vertices): [0, 0, 12]
    // POLYHEDRAL COMPLEX CELLS
    // number of polytopes (total): 42
    // number of polytopes (binned by dimension): [9, 21, 12]
    // betti numbers (of polyhedral complex): [1, 1, 0]
    // INTERSECTIONS
    // number of pairs of intersecting facets, binned by the dimension of the intersection polytope: [24, 15, 0]
    // DOWKER NERVE COMPLEX FACETS
    // number of dowker nerve complex facets (total): 12
    // number of dowker nerve complex facets (binned by dimension): [0, 0, 12]
    // DOWKER NERVE COMPLEX CELLS
    // number of nerve dowker complex cells (total): 42
    // number of nerve dowker complex cells (binned by dimension): [9, 21, 12]
    // betti numbers (of dowker nerve): [1, 1, 0]


}
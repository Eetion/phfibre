//  DESCRIPTION

//  A "musical symbol" consisting of a triangle and a dangling edge


//  DEPENDENCIES

use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use phfibre::pipelines::simplex_pipeline;
use phfibre::phfibre::{NoExtraCondition, LowerStarCondition, LowerEdgeCondition};
use solar::utilities::sequences_and_ordinals::ordinate_unique_vals;
use solar::cell_complexes::simplices_unweighted::facets::ordered_subsimplices_up_thru_dim_concatenated_vec; 


//  CODE

fn main() {
        

    //  ----------------------------------------------------------------------------------------------    
    //  EMPTY TRIANGLE WITH EDGE ATTACHED
    //  ----------------------------------------------------------------------------------------------


    //  Define the base space, barcode, and ring
    let complex_facets      =   vec![  vec![0,1,2], vec![2, 3] ]; // triangle with edge attached
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 1);    // build all faces up to dimension 1

    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:2} ],
                                    fin: vec![ BarFinite{dim:0,birth:0,death:1} ],
                                    ordinal: ordinate_unique_vals( & vec![0, 1, 2] ),
                                };

    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new(); // this just defines the rational numbers

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
    // simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [1, 2], [2, 3]]
    // POLYHEDRAL COMPLEX FACETS
    // number of facets (total): 104
    // number of facets (binned by dimension): [0, 0, 104]
    // number of facets (binned by number of vertices): [0, 0, 104]
    // POLYHEDRAL COMPLEX CELLS
    // number of polytopes (total): 393
    // number of polytopes (binned by dimension): [92, 197, 104]
    // betti numbers (of polyhedral complex): [1, 2, 0]
    // INTERSECTIONS
    // number of pairs of intersecting facets, binned by the dimension of the intersection polytope: [271, 161, 0]
    // DOWKER NERVE COMPLEX FACETS
    // number of dowker nerve complex facets (total): 104
    // number of dowker nerve complex facets (binned by dimension): [0, 0, 58, 46]
    // DOWKER NERVE COMPLEX CELLS
    // number of nerve dowker complex cells (total): 669
    // number of nerve dowker complex cells (binned by dimension): [92, 289, 242, 46]
    // betti numbers (of dowker nerve): [1, 2, 0, 0]


}
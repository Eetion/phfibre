//  DESCRIPTION

//  The 3-simplex.


//  DEPENDENCIES

use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use phfibre::pipelines::simplex_pipeline;
use phfibre::phfibre::{ConditionNone, ConditionLowerStar, ConditionLowerEdge};
use solar::utilities::sequences_and_ordinals::ordinate_unique_vals;
use solar::cell_complexes::simplices_unweighted::facets::ordered_subsimplices_up_thru_dim_concatenated_vec; 
use std::iter::FromIterator;

//  CODE

fn main() {

    //  ----------------------------------------------------------------------------------------------    
    //  2-SKELETON, NO FINITE BARS
    //  ----------------------------------------------------------------------------------------------

    //  Define the base space, barcode, and ring

    let complex_facets      =   vec![  vec![0, 1, 2, 3] ];
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 2); 

    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:2,birth:1} ],
                                    fin: vec![ ],
                                    ordinal: ordinate_unique_vals( & vec![ 0, 1 ] ),
                                };

    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};      // this struct won't impose any extra conditions on the filtrations we build    

    let analyze_dowker_dual =   true;

    simplex_pipeline(
        &   simplex_sequence,
        &   barcode,
        &   ring,
        &   precondition_to_make_new_lev_set_lower_none,
            analyze_dowker_dual,
    );


    //  ----------------------------------------------------------------------------------------------    
    //  3-SKELETON, NO FINITE BARS
    //  ----------------------------------------------------------------------------------------------

    //  Define the base space, barcode, and ring

    let complex_facets      =   vec![  vec![0, 1, 2, 3] ];
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 3); 

    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0} ],
                                    fin: vec![ ],
                                    ordinal: ordinate_unique_vals( & vec![ 0 ] ),
                                };

    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};      // this struct won't impose any extra conditions on the filtrations we build    

    simplex_pipeline(
        &   simplex_sequence,
        &   barcode,
        &   ring,
        &   precondition_to_make_new_lev_set_lower_none,
            analyze_dowker_dual,
    );

}

// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 2, birth: 1 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0, 1], val_to_ord: {1: 1, 0: 0} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 11.073037962s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 1920
// number of facets (binned by dimension): [0, 0, 0, 0, 0, 0, 1920]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 0, 1920]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 522.300995ms"
// number of polytopes (total): 39566
// number of polytopes (binned by dimension): [64, 858, 4480, 10860, 13320, 8064, 1920]
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 0, 1, 0, 0, 0, 0]
// DOWKER NERVE COMPLEX FACETS
// number of dowker nerve complex facets (total): 1920
// number of dowker nerve complex facets (binned by dimension): [0, 0, 0, 0, 0, 0, 1920]
// DOWKER NERVE COMPLEX CELLS
// number of nerve dowker complex cells (total): 39566
// number of nerve dowker complex cells (binned by dimension): [64, 858, 4480, 10860, 13320, 8064, 1920]
// betti numbers (of dowker nerve, user-specified ring coefficients): [1, 0, 1, 0, 0, 0, 0]
// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3], [0, 1, 2, 3]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0], val_to_ord: {0: 0} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 11.822844753s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 1920
// number of facets (binned by dimension): [0, 0, 0, 0, 0, 0, 0, 1920]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 0, 0, 1920]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 981.149253ms"
// number of polytopes (total): 79133
// number of polytopes (binned by dimension): [65, 922, 5338, 15340, 24180, 21384, 9984, 1920]
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 0, 0, 0, 0, 0, 0, 0]
// DOWKER NERVE COMPLEX FACETS
// number of dowker nerve complex facets (total): 1920
// number of dowker nerve complex facets (binned by dimension): [0, 0, 0, 0, 0, 0, 0, 1920]
// DOWKER NERVE COMPLEX CELLS
// number of nerve dowker complex cells (total): 79133
// number of nerve dowker complex cells (binned by dimension): [65, 922, 5338, 15340, 24180, 21384, 9984, 1920]
// betti numbers (of dowker nerve, user-specified ring coefficients): [1, 0, 0, 0, 0, 0, 0, 0]
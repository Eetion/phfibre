//  DESCRIPTION

//  The 3-simplex.


//  DEPENDENCIES

use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use phfibre::pipelines::{simplex_pipeline, analyze_fibre};
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

    let complex_facets      =   vec![  vec![0, 1, 2, 3], vec![3, 4] ];
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 2); 

    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:2,birth:1} ],
                                    fin: vec![ ],
                                    ordinal: ordinate_unique_vals( & vec![ 0, 1 ] ),
                                };

    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};      // this struct won't impose any extra conditions on the filtrations we build    
    let analyze_dowker_dual =   true;

    // simplex_pipeline(
    //     &   simplex_sequence,
    //     &   barcode,
    //     &   ring,
    //     &   precondition_to_make_new_lev_set_lower_none,
    //         analyze_dowker_dual,
    // );
    let poly_complex_facets =       simplex_pipeline(
        &   simplex_sequence,
        &   barcode,
        &   ring,
        &   precondition_to_make_new_lev_set_lower_none,
            false, // do not analyze the dowker dual to the nerve complex
    );  

    //  ANALYZE    
    let analyze_dowker_dual =   false;
    analyze_fibre( 
        &   poly_complex_facets,
            ring.clone(),
            analyze_dowker_dual,
    );      
    

}

// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [4], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [3, 4], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 2, birth: 1 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0, 1], val_to_ord: {0: 0, 1: 1} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 10069.209142083s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 12816
// number of facets (binned by dimension): [0, 0, 0, 0, 0, 0, 0, 12816]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 0, 0, 12816]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 7.760420012s"
// number of polytopes (total): 423452
// number of polytopes (binned by dimension): [183, 3262, 22283, 72710, 127173, 122937, 62088, 12816]
// Time elapsed to compute binary-coeff polyhedral betti numbers 942.068145ms
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 0, 1, 0, 0, 0, 0, 0]
// DOWKER NERVE COMPLEX FACETS
// number of dowker nerve complex facets (total): 12816
// number of dowker nerve complex facets (binned by dimension): [0, 0, 0, 0, 0, 0, 0, 10896, 0, 0, 0, 0, 0, 1920]
// DOWKER NERVE COMPLEX CELLS
// Time elapsed to compute nerve dowker dual betti numbers (user specified coeff) 109.320907342s
// number of nerve dowker complex cells (total): 14627744
// number of nerve dowker complex cells (binned by dimension): [183, 4978, 51737, 274888, 890313, 1940305, 2996712, 3350076, 2713680, 1572360, 634368, 169344, 26880, 1920]
// betti numbers (of dowker nerve, user-specified ring coefficients): [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]



// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [4], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [3, 4], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 2, birth: 1 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0, 1], val_to_ord: {0: 0, 1: 1} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 3890.206579721s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.



// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [4], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [3, 4], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 2, birth: 1 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0, 1], val_to_ord: {1: 1, 0: 0} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 3709.737526693s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 12816
// number of facets (binned by dimension): [0, 0, 0, 0, 0, 0, 0, 12816]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 0, 0, 12816]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 7.998989202s"
// number of polytopes (total): 423452
// number of polytopes (binned by dimension): [183, 3262, 22283, 72710, 127173, 122937, 62088, 12816]
// Time elapsed to compute binary-coeff polyhedral betti numbers 968.328079ms
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 0, 1, 0, 0, 0, 0, 0]
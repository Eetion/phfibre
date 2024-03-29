//  DESCRIPTION

//  The 3-simplex.


//  DEPENDENCIES

use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use phfibre::pipelines::{simplex_pipeline, analyze_fibre};
use phfibre::phfibre::{ConditionNone, ConditionLowerStar, ConditionLowerEdge};
use solar::utilities::sequences_and_ordinals::ordinate_unique_vals;
use solar::utilities::cell_complexes::simplices_unweighted::facets::ordered_subsimplices_up_thru_dim_concatenated_vec;
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

    println!("\n\n2-SKELETON, NO FINITE BARS");
    let save_dir_opt        =   None;
    let poly_complex_facets =   simplex_pipeline(
                                    &   simplex_sequence,
                                    &   barcode,
                                    &   ring,
                                    &   precondition_to_make_new_lev_set_lower_none,
                                        analyze_dowker_dual,
                                        save_dir_opt,
                                ); 

    let analyze_dowker_dual =   true;
    analyze_fibre( 
        &   poly_complex_facets,
            ring.clone(),
            analyze_dowker_dual,
            save_dir_opt,
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

    let analyze_dowker_dual =   true;

    println!("\n\n3-SKELETON, NO FINITE BARS");
    let save_dir_opt        =   None;
    let poly_complex_facets =   simplex_pipeline(
                                    &   simplex_sequence,
                                    &   barcode,
                                    &   ring,
                                    &   precondition_to_make_new_lev_set_lower_none,
                                        analyze_dowker_dual,
                                        save_dir_opt,
                                ); 
                                    
    let analyze_dowker_dual =   true;
    analyze_fibre( 
        &   poly_complex_facets,
            ring.clone(),
            analyze_dowker_dual,
            save_dir_opt,
    );    

}



// WITH NEW, CORRECTED, (SLOWER) CODE


// 2-SKELETON, NO FINITE BARS
// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 2, birth: 1 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0, 1], val_to_ord: {0: 0, 1: 1} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 20.244175473s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 1920
// number of facets (binned by dimension): [0, 0, 0, 0, 0, 0, 1920]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 0, 1920]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 571.232966ms"
// number of polytopes (total): 39566
// number of polytopes (binned by dimension): [64, 858, 4480, 10860, 13320, 8064, 1920]
// Time elapsed to compute binary-coeff polyhedral betti numbers 60.936718ms
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 0, 1, 0, 0, 0, 0]
// DOWKER NERVE COMPLEX FACETS
// number of dowker nerve complex facets (total): 1920
// number of dowker nerve complex facets (binned by dimension): [0, 0, 0, 0, 0, 0, 1920]
// DOWKER NERVE COMPLEX CELLS
// Time elapsed to compute nerve dowker dual betti numbers (user specified coeff) 109.563696ms
// number of nerve dowker complex cells (total): 39566
// number of nerve dowker complex cells (binned by dimension): [64, 858, 4480, 10860, 13320, 8064, 1920]
// betti numbers (of dowker nerve, user-specified ring coefficients): [1, 0, 1, 0, 0, 0, 0]


// 3-SKELETON, NO FINITE BARS
// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3], [0, 1, 2, 3]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0], val_to_ord: {0: 0} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 31.434274409s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 1920
// number of facets (binned by dimension): [0, 0, 0, 0, 0, 0, 0, 1920]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 0, 0, 1920]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 1.073394944s"
// number of polytopes (total): 79133
// number of polytopes (binned by dimension): [65, 922, 5338, 15340, 24180, 21384, 9984, 1920]
// Time elapsed to compute binary-coeff polyhedral betti numbers 82.959582ms
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 0, 0, 0, 0, 0, 0, 0]
// DOWKER NERVE COMPLEX FACETS
// number of dowker nerve complex facets (total): 1920
// number of dowker nerve complex facets (binned by dimension): [0, 0, 0, 0, 0, 0, 0, 1920]
// DOWKER NERVE COMPLEX CELLS
// Time elapsed to compute nerve dowker dual betti numbers (user specified coeff) 142.42424ms
// number of nerve dowker complex cells (total): 79133
// number of nerve dowker complex cells (binned by dimension): [65, 922, 5338, 15340, 24180, 21384, 9984, 1920]
// betti numbers (of dowker nerve, user-specified ring coefficients): [1, 0, 0, 0, 0, 0, 0, 0]







// WITH OLD, BUGGY CODE

// 2-SKELETON, NO FINITE BARS
// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 2, birth: 1 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0, 1], val_to_ord: {1: 1, 0: 0} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 3.437952036s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 1920
// number of facets (binned by dimension): [0, 0, 0, 0, 0, 0, 1920]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 0, 1920]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 601.688788ms"
// number of polytopes (total): 39566
// number of polytopes (binned by dimension): [64, 858, 4480, 10860, 13320, 8064, 1920]
// Time elapsed to compute binary-coeff polyhedral betti numbers 65.030835ms
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 0, 1, 0, 0, 0, 0]
// DOWKER NERVE COMPLEX FACETS
// number of dowker nerve complex facets (total): 1920
// number of dowker nerve complex facets (binned by dimension): [0, 0, 0, 0, 0, 0, 1920]
// DOWKER NERVE COMPLEX CELLS
// Time elapsed to compute nerve dowker dual betti numbers (user specified coeff) 113.248706ms
// number of nerve dowker complex cells (total): 39566
// number of nerve dowker complex cells (binned by dimension): [64, 858, 4480, 10860, 13320, 8064, 1920]
// betti numbers (of dowker nerve, user-specified ring coefficients): [1, 0, 1, 0, 0, 0, 0]


// 3-SKELETON, NO FINITE BARS
// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3], [0, 1, 2, 3]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0], val_to_ord: {0: 0} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 3.731470418s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 1920
// number of facets (binned by dimension): [0, 0, 0, 0, 0, 0, 0, 1920]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 0, 0, 1920]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 1.12540889s"
// number of polytopes (total): 79133
// number of polytopes (binned by dimension): [65, 922, 5338, 15340, 24180, 21384, 9984, 1920]
// Time elapsed to compute binary-coeff polyhedral betti numbers 87.270717ms
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 0, 0, 0, 0, 0, 0, 0]
// DOWKER NERVE COMPLEX FACETS
// number of dowker nerve complex facets (total): 1920
// number of dowker nerve complex facets (binned by dimension): [0, 0, 0, 0, 0, 0, 0, 1920]
// DOWKER NERVE COMPLEX CELLS
// Time elapsed to compute nerve dowker dual betti numbers (user specified coeff) 150.460046ms
// number of nerve dowker complex cells (total): 79133
// number of nerve dowker complex cells (binned by dimension): [65, 922, 5338, 15340, 24180, 21384, 9984, 1920]
// betti numbers (of dowker nerve, user-specified ring coefficients): [1, 0, 0, 0, 0, 0, 0, 0]
//  DESCRIPTION

//  The 2-sphere


//  DEPENDENCIES

use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use phfibre::pipelines::{simplex_pipeline, analyze_fibre};
use phfibre::phfibre::{ConditionNone, ConditionLowerStar, ConditionLowerEdge};
use solar::utilities::sequences_and_ordinals::{ordinate_unique_vals, BiMapSequential};
use solar::cell_complexes::simplices_unweighted::facets::ordered_subsimplices_up_thru_dim_concatenated_vec; 
use std::iter::FromIterator;

//  CODE // !!!!! SEE BELOW FOR RESULTS

fn main() {

    //  ----------------------------------------------------------------------------------------------    
    //  NO FINITE BARS, S^1 WITH 3 EDGES
    //  ----------------------------------------------------------------------------------------------

    println!("\nBARCODE WITH 0 FINITE BARS // S^1 WITH 3 EDGES");
    println!("===============================================================================================\n");  


    //  Define the base space, barcode, and ring

    let complex_facets      =   vec![  vec![0,1,2] ];
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 1);    // build all faces up to dimension 1

    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:1} ],
                                    fin: vec![ ],
                                    ordinal: ordinate_unique_vals( & vec![ 0, 1 ] ),
                                };

    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

    let simplex_bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence.clone() );

    //  Define any extra conditions on starting a new level set

    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};        
    let precondition_to_make_new_lev_set_lower_star     =   ConditionLowerStar{ simplex_bimap_sequential: & simplex_bimap_sequential };
    let precondition_to_make_new_lev_set_lower_edge     =   ConditionLowerEdge{ simplex_bimap_sequential: & simplex_bimap_sequential };    
    
    let analyze_dowker_dual                             =   true;

    //  UNCONSTRAINED FILTRATION
    //  -------------------------------------------------------------------------------------------

    println!("\nUNCONSTRAINED FILTRATIONS");
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

    //  LOWER EDGE
    //  -------------------------------------------------------------------------------------------

    println!("\nLOWER EDGE");
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

    //  LOWER STAR
    //  -------------------------------------------------------------------------------------------

    println!("\nLOWER STAR");
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


    println!("\n\nBARCODE WITH 0 FINITE BARS // S^1 WITH 5 EDGES");
    println!("===============================================================================================\n");  


    //  Define the base space, barcode, and ring

    let complex_facets      =   vec![ vec![0,4], vec![0,1], vec![1,2], vec![2,3], vec![3,4] ];
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 1);    // build all faces up to dimension 1

    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:1} ],
                                    fin: vec![ ],
                                    ordinal: ordinate_unique_vals( & vec![ 0, 1 ] ),
                                };

    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

    let simplex_bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence.clone() );

    //  Define any extra conditions on starting a new level set

    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};        
    let precondition_to_make_new_lev_set_lower_star     =   ConditionLowerStar{ simplex_bimap_sequential: & simplex_bimap_sequential };
    let precondition_to_make_new_lev_set_lower_edge     =   ConditionLowerEdge{ simplex_bimap_sequential: & simplex_bimap_sequential };    


    //  UNCONSTRAINED FILTRATION
    //  -------------------------------------------------------------------------------------------

    println!("\nUNCONSTRAINED FILTRATIONS");
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

    //  LOWER EDGE
    //  -------------------------------------------------------------------------------------------

    println!("\nLOWER EDGE");
    simplex_pipeline(
        &   simplex_sequence,
        &   barcode,
        &   ring,
        &   precondition_to_make_new_lev_set_lower_edge,
            analyze_dowker_dual,
            save_dir_opt,
    );        

    //  LOWER STAR
    //  -------------------------------------------------------------------------------------------

    println!("\nLOWER STAR");
    simplex_pipeline(
        &   simplex_sequence,
        &   barcode,
        &   ring,
        &   precondition_to_make_new_lev_set_lower_star,
            analyze_dowker_dual,
            save_dir_opt,
    );             

}

// BARCODE WITH 0 FINITE BARS // S^1 WITH 3 EDGES
// ===============================================================================================


// UNCONSTRAINED FILTRATIONS
// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [0, 1], [0, 2], [1, 2]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 1 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0, 1], val_to_ord: {0: 0, 1: 1} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 1.031775ms

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 12
// number of facets (binned by dimension): [0, 0, 12]
// number of facets (binned by number of vertices): [0, 0, 12]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 586.029µs"
// number of polytopes (total): 42
// number of polytopes (binned by dimension): [9, 21, 12]
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 1, 0]
// DOWKER NERVE COMPLEX FACETS
// number of dowker nerve complex facets (total): 12
// number of dowker nerve complex facets (binned by dimension): [0, 0, 12]
// DOWKER NERVE COMPLEX CELLS
// number of nerve dowker complex cells (total): 42
// number of nerve dowker complex cells (binned by dimension): [9, 21, 12]
// betti numbers (of dowker nerve, user-specified ring coefficients): [1, 1, 0]

// LOWER EDGE
// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [0, 1], [0, 2], [1, 2]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 1 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0, 1], val_to_ord: {0: 0, 1: 1} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 1.672419ms

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 12
// number of facets (binned by dimension): [0, 0, 12]
// number of facets (binned by number of vertices): [0, 0, 12]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 640.725µs"
// number of polytopes (total): 42
// number of polytopes (binned by dimension): [9, 21, 12]
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 1, 0]
// DOWKER NERVE COMPLEX FACETS
// number of dowker nerve complex facets (total): 12
// number of dowker nerve complex facets (binned by dimension): [0, 0, 12]
// DOWKER NERVE COMPLEX CELLS
// number of nerve dowker complex cells (total): 42
// number of nerve dowker complex cells (binned by dimension): [9, 21, 12]
// betti numbers (of dowker nerve, user-specified ring coefficients): [1, 1, 0]

// LOWER STAR
// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [0, 1], [0, 2], [1, 2]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 1 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0, 1], val_to_ord: {0: 0, 1: 1} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 1.142797ms

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 6
// number of facets (binned by dimension): [0, 6]
// number of facets (binned by number of vertices): [0, 6]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 93.701µs"
// number of polytopes (total): 12
// number of polytopes (binned by dimension): [6, 6]
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 1]
// DOWKER NERVE COMPLEX FACETS
// number of dowker nerve complex facets (total): 6
// number of dowker nerve complex facets (binned by dimension): [0, 6]
// DOWKER NERVE COMPLEX CELLS
// number of nerve dowker complex cells (total): 12
// number of nerve dowker complex cells (binned by dimension): [6, 6]
// betti numbers (of dowker nerve, user-specified ring coefficients): [1, 1]


// BARCODE WITH 0 FINITE BARS // S^1 WITH 5 EDGES
// ===============================================================================================


// UNCONSTRAINED FILTRATIONS
// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [4], [0, 1], [0, 4], [1, 2], [2, 3], [3, 4]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 1 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0, 1], val_to_ord: {0: 0, 1: 1} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 531.699656ms

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 80
// number of facets (binned by dimension): [0, 0, 0, 0, 80]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 80]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 9.509912ms"
// number of polytopes (total): 820
// number of polytopes (binned by dimension): [25, 150, 305, 260, 80]
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 1, 0, 0, 0]
// DOWKER NERVE COMPLEX FACETS
// number of dowker nerve complex facets (total): 80
// number of dowker nerve complex facets (binned by dimension): [0, 0, 0, 0, 80]
// DOWKER NERVE COMPLEX CELLS
// number of nerve dowker complex cells (total): 820
// number of nerve dowker complex cells (binned by dimension): [25, 150, 305, 260, 80]
// betti numbers (of dowker nerve, user-specified ring coefficients): [1, 1, 0, 0, 0]

// LOWER EDGE
// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [4], [0, 1], [0, 4], [1, 2], [2, 3], [3, 4]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 1 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0, 1], val_to_ord: {0: 0, 1: 1} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 795.001837ms

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 80
// number of facets (binned by dimension): [0, 0, 0, 0, 80]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 80]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 8.614727ms"
// number of polytopes (total): 820
// number of polytopes (binned by dimension): [25, 150, 305, 260, 80]
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 1, 0, 0, 0]
// DOWKER NERVE COMPLEX FACETS
// number of dowker nerve complex facets (total): 80
// number of dowker nerve complex facets (binned by dimension): [0, 0, 0, 0, 80]
// DOWKER NERVE COMPLEX CELLS
// number of nerve dowker complex cells (total): 820
// number of nerve dowker complex cells (binned by dimension): [25, 150, 305, 260, 80]
// betti numbers (of dowker nerve, user-specified ring coefficients): [1, 1, 0, 0, 0]

// LOWER STAR
// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [4], [0, 1], [0, 4], [1, 2], [2, 3], [3, 4]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 1 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0, 1], val_to_ord: {0: 0, 1: 1} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 655.330807ms

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 40
// number of facets (binned by dimension): [0, 0, 0, 40]
// number of facets (binned by number of vertices): [0, 0, 0, 40]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 2.134443ms"
// number of polytopes (total): 240
// number of polytopes (binned by dimension): [20, 80, 100, 40]
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 1, 0, 0]
// DOWKER NERVE COMPLEX FACETS
// number of dowker nerve complex facets (total): 40
// number of dowker nerve complex facets (binned by dimension): [0, 0, 0, 40]
// DOWKER NERVE COMPLEX CELLS
// number of nerve dowker complex cells (total): 240
// number of nerve dowker complex cells (binned by dimension): [20, 80, 100, 40]
// betti numbers (of dowker nerve, user-specified ring coefficients): [1, 1, 0, 0]
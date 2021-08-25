//  DESCRIPTION

//  The closed interavl, decomposed into 1-simplices.


//  DEPENDENCIES

use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use phfibre::pipelines::{simplex_pipeline, analyze_fibre};
use phfibre::phfibre::{ConditionNone, ConditionLowerStar, ConditionLowerEdge};
use solar::utilities::sequences_and_ordinals::{ordinate_unique_vals, BiMapSequential};
use solar::cell_complexes::simplices_unweighted::facets::ordered_subsimplices_up_thru_dim_concatenated_vec; 
use std::iter::FromIterator;

//  CODE  ( !!!!  SEE BELOW FOR RESULTS  )

fn main() {

    //  -------------------------------------------------------------------------------------------    
    //  5-VERTEX INTERVAL WITH NO FINITE BARS
    //  -------------------------------------------------------------------------------------------

    println!("\nBARCODE WITH 0 FINITE BARS");
    println!("===============================================================================================\n");    


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

    let simplex_bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence.clone() );

    //  Define any extra conditions on starting a new level set

    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};        
    let precondition_to_make_new_lev_set_lower_star     =   ConditionLowerStar{ simplex_bimap_sequential: & simplex_bimap_sequential };
    let precondition_to_make_new_lev_set_lower_edge     =   ConditionLowerEdge{ simplex_bimap_sequential: & simplex_bimap_sequential };    

    //  This example is small enough to allow analysis of the nerve
    let analyze_dowker_dual                             =   true;

    //  UNCONSTRAINED FILTRATION
    //  -------------------------------------------------------------------------------------------

    println!("\nUNCONSTRAINED FILTRATIONS");
    simplex_pipeline(
        &   simplex_sequence,
        &   barcode,
        &   ring,
        &   precondition_to_make_new_lev_set_lower_none,
            analyze_dowker_dual,
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
    );    

    

    //  ----------------------------------------------------------------------------------------------    
    //  5-VERTEX INTERVAL WITH ONE FINITE BAR
    //  ----------------------------------------------------------------------------------------------

    println!("\n\nBARCODE WITH 1 FINITE BAR");
    println!("===============================================================================================\n");    

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

    //  UNCONSTRAINED FILTRATION
    //  -------------------------------------------------------------------------------------------

    println!("\nUNCONSTRAINED FILTRATIONS");
    simplex_pipeline(
        &   simplex_sequence,
        &   barcode,
        &   ring,
        &   precondition_to_make_new_lev_set_lower_none,
            analyze_dowker_dual,
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
    );      


}



// BARCODE WITH 0 FINITE BARS
// ===============================================================================================


// UNCONSTRAINED FILTRATIONS
// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [4], [5], [0, 1], [1, 2], [2, 3], [3, 4], [4, 5]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0], val_to_ord: {0: 0} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 7.120471372s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 32
// number of facets (binned by dimension): [0, 0, 0, 0, 0, 32]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 32]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 5.768158ms"
// number of polytopes (total): 791
// number of polytopes (binned by dimension): [21, 105, 231, 258, 144, 32]
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 0, 0, 0, 0, 0]
// DOWKER NERVE COMPLEX FACETS
// number of dowker nerve complex facets (total): 32
// number of dowker nerve complex facets (binned by dimension): [0, 0, 0, 0, 0, 32]
// DOWKER NERVE COMPLEX CELLS
// number of nerve dowker complex cells (total): 791
// number of nerve dowker complex cells (binned by dimension): [21, 105, 231, 258, 144, 32]
// betti numbers (of dowker nerve, user-specified ring coefficients): [1, 0, 0, 0, 0, 0]

// LOWER EDGE
// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [4], [5], [0, 1], [1, 2], [2, 3], [3, 4], [4, 5]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0], val_to_ord: {0: 0} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 10.700447188s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 32
// number of facets (binned by dimension): [0, 0, 0, 0, 0, 32]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 32]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 6.259558ms"
// number of polytopes (total): 791
// number of polytopes (binned by dimension): [21, 105, 231, 258, 144, 32]
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 0, 0, 0, 0, 0]
// DOWKER NERVE COMPLEX FACETS
// number of dowker nerve complex facets (total): 32
// number of dowker nerve complex facets (binned by dimension): [0, 0, 0, 0, 0, 32]
// DOWKER NERVE COMPLEX CELLS
// number of nerve dowker complex cells (total): 791
// number of nerve dowker complex cells (binned by dimension): [21, 105, 231, 258, 144, 32]
// betti numbers (of dowker nerve, user-specified ring coefficients): [1, 0, 0, 0, 0, 0]

// LOWER STAR
// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [4], [5], [0, 1], [1, 2], [2, 3], [3, 4], [4, 5]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0], val_to_ord: {0: 0} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 7.29180338s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 32
// number of facets (binned by dimension): [0, 0, 0, 0, 0, 32]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 32]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 5.758572ms"
// number of polytopes (total): 791
// number of polytopes (binned by dimension): [21, 105, 231, 258, 144, 32]
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 0, 0, 0, 0, 0]
// DOWKER NERVE COMPLEX FACETS
// number of dowker nerve complex facets (total): 32
// number of dowker nerve complex facets (binned by dimension): [0, 0, 0, 0, 0, 32]
// DOWKER NERVE COMPLEX CELLS
// number of nerve dowker complex cells (total): 791
// number of nerve dowker complex cells (binned by dimension): [21, 105, 231, 258, 144, 32]
// betti numbers (of dowker nerve, user-specified ring coefficients): [1, 0, 0, 0, 0, 0]


// BARCODE WITH 1 FINITE BAR
// ===============================================================================================


// UNCONSTRAINED FILTRATIONS
// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [4], [5], [0, 1], [1, 2], [2, 3], [3, 4], [4, 5]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2], val_to_ord: {0: 0, 1: 1, 2: 2} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 13.493090743s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 1824
// number of facets (binned by dimension): [0, 0, 0, 0, 1824]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 1824]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 286.361809ms"
// number of polytopes (total): 22054
// number of polytopes (binned by dimension): [990, 4674, 8214, 6352, 1824]
// betti numbers (of polyhedral complex, Z2 coefficients): [2, 0, 0, 0, 0]
// DOWKER NERVE COMPLEX FACETS
// number of dowker nerve complex facets (total): 1824
// number of dowker nerve complex facets (binned by dimension): [0, 0, 0, 0, 576, 0, 0, 704, 288, 0, 0, 256]
// DOWKER NERVE COMPLEX CELLS
// number of nerve dowker complex cells (total): 1180834
// number of nerve dowker complex cells (binned by dimension): [990, 13654, 64756, 160186, 248672, 269712, 216320, 129712, 56608, 16896, 3072, 256]
// betti numbers (of dowker nerve, user-specified ring coefficients): [2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

// LOWER EDGE
// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [4], [5], [0, 1], [1, 2], [2, 3], [3, 4], [4, 5]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2], val_to_ord: {0: 0, 1: 1, 2: 2} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 16.112812525s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 1824
// number of facets (binned by dimension): [0, 0, 0, 0, 1824]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 1824]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 286.005101ms"
// number of polytopes (total): 22054
// number of polytopes (binned by dimension): [990, 4674, 8214, 6352, 1824]
// betti numbers (of polyhedral complex, Z2 coefficients): [2, 0, 0, 0, 0]
// DOWKER NERVE COMPLEX FACETS
// number of dowker nerve complex facets (total): 1824
// number of dowker nerve complex facets (binned by dimension): [0, 0, 0, 0, 576, 0, 0, 704, 288, 0, 0, 256]
// DOWKER NERVE COMPLEX CELLS
// number of nerve dowker complex cells (total): 1180834
// number of nerve dowker complex cells (binned by dimension): [990, 13654, 64756, 160186, 248672, 269712, 216320, 129712, 56608, 16896, 3072, 256]
// betti numbers (of dowker nerve, user-specified ring coefficients): [2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

// LOWER STAR
// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [4], [5], [0, 1], [1, 2], [2, 3], [3, 4], [4, 5]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2], val_to_ord: {0: 0, 1: 1, 2: 2} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 6.261546729s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 416
// number of facets (binned by dimension): [0, 0, 0, 416]
// number of facets (binned by number of vertices): [0, 0, 0, 416]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 33.228035ms"
// number of polytopes (total): 2962
// number of polytopes (binned by dimension): [330, 1064, 1152, 416]
// betti numbers (of polyhedral complex, Z2 coefficients): [2, 0, 0, 0]
// DOWKER NERVE COMPLEX FACETS
// number of dowker nerve complex facets (total): 416
// number of dowker nerve complex facets (binned by dimension): [0, 0, 0, 160, 0, 224, 0, 32]
// DOWKER NERVE COMPLEX CELLS
// number of nerve dowker complex cells (total): 17554
// number of nerve dowker complex cells (binned by dimension): [330, 2216, 5056, 5408, 3136, 1120, 256, 32]
// betti numbers (of dowker nerve, user-specified ring coefficients): [2, 0, 0, 0, 0, 0, 0, 0]
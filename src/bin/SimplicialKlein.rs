//  DESCRIPTION

//  Simplicial Klein, Only Lower-Stars
// BECAUSE WE TAKE RATIONAL COEFFCIENTS, I THINK HOMOLOGY OF KLEIN IS H_0=Q, H_1=Q, H_n=0 for n>1. SO IN THE ASSOCIATED BARCODE THERE IS JUST 2 INFINITE BARS. If not true we need to modify infinite bars in barcodes below


     ///////               0---------1---------2--------0
     ///                   | \       |       /    \     |
     ///                   |    \    |    /         \   |
     ///                   |       \ | /              \ |
     ///                   3---------4---------–--------5
     ///                   | \       |  \             / |
     ///                   |    \    |    \         /   |
     ///                   |       \ |       \   /      |
     ///                   5---------6---------7--------3
     ///                   |       / |       /   _ /  / |
     ///                   |    /    |    /   _ /   /   |
     ///                   | /       | /___ /    /      |
     ///                   0---------1---------2--------0


//  DEPENDENCIES

use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use phfibre::pipelines::{simplex_pipeline, analyze_fibre};
use phfibre::phfibre::{ConditionNone, ConditionLowerStar, ConditionLowerEdge};
use solar::utilities::sequences_and_ordinals::{ordinate_unique_vals, BiMapSequential};
use solar::utilities::cell_complexes::simplices_unweighted::facets::ordered_subsimplices_up_thru_dim_concatenated_vec; 
use std::iter::FromIterator;

//  CODE

fn main() {

    let save_dir_opt        =   None; // we will not save any files

    println!("\n BECAUSE WE TAKE RATIONAL COEFFCIENTS, I THINK HOMOLOGY OF KLEIN IS Q for H0 AND H1, and null for higher degrees. SO IN THE ASSOCIATED BARCODE THERE IS JUST TWO INFINITE BARS. If not true we need to modify infinite bars in barcodes below");
    //  ----------------------------------------------------------------------------------------------    
    //  DEFINE COEFFICIENT RING
    //  ----------------------------------------------------------------------------------------------    

    //  This is the ring of rationals.
    type CeofficientRing    =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >; // the TYPE of ring object
    let ring                =   CeofficientRing::new(); // the object itself
    
    
    
    // -------------------- SIMPLICIAL KLEIN, BARCODE WITH NO FINITE BARS------------

    //  Define the base space, barcode, and ring
    println!("\n SIMPLICIAL KLEIN, BARCODE WITH NO FINITE BARS");
    println!("===============================================================================================\n");  
    let complex_facets      =   vec![    vec![0,3,4], vec![0,1,4], vec![1,2,4],vec![2,4,5],vec![0,2,5], vec![3,5,6],vec![3,4,6],vec![4,6,7],vec![4,5,7],vec![3,5,7],vec![0,5,6],vec![0,1,6], vec![1,6,7],vec![1,3,7],vec![1,2,3],vec![0,2,3]];
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 2);    // build all faces up to dimension 1

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:1}],
        fin: vec![ ],
        ordinal: ordinate_unique_vals( & vec![0, 1] ),
    };
    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

    let simplex_bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence.clone() );

    //  Define any extra conditions on starting a new level set

    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};        
    let precondition_to_make_new_lev_set_lower_star     =   ConditionLowerStar{ simplex_bimap_sequential: & simplex_bimap_sequential };
    let precondition_to_make_new_lev_set_lower_edge     =   ConditionLowerEdge{ simplex_bimap_sequential: & simplex_bimap_sequential };    
    
    let analyze_dowker_dual                             =   false;

    //  LOWER STAR
    //  -------------------------------------------------------------------------------------------

    println!("\nLOWER STAR");
     let poly_complex_facets =       simplex_pipeline(
                                         &   simplex_sequence,
                                         &   barcode,
                                        &   ring,
                                        &   precondition_to_make_new_lev_set_lower_star,
                                        analyze_dowker_dual,
                                        save_dir_opt,
                                    );                        

    let analyze_dowker_dual =   false;
    analyze_fibre( 
        &   poly_complex_facets,
            ring.clone(),
            analyze_dowker_dual,
            save_dir_opt,
    );  

    // -------------------- SIMPLICIAL KLEIN, BARCODE WITH ONE FINITE BAR ------------

    //  Define the base space, barcode, and ring
    println!("\n SIMPLICIAL KLEIN, BARCODE WITH ONE FINITE BAR");
    println!("===============================================================================================\n");  
    let complex_facets      =   vec![    vec![0,3,4], vec![0,1,4], vec![1,2,4],vec![2,4,5],vec![0,2,5], vec![3,5,6],vec![3,4,6],vec![4,6,7],vec![4,5,7],vec![3,5,7],vec![0,5,6],vec![0,1,6], vec![1,6,7],vec![1,3,7],vec![1,2,3],vec![0,2,3]];

    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 2);    // build all faces up to dimension 1

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:3}],
        fin: vec![ BarFinite{dim:0,birth:1,death:2} ],
        ordinal: ordinate_unique_vals( & vec![0, 1, 2, 3] ),
    };
    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

    let simplex_bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence.clone() );

    //  Define any extra conditions on starting a new level set

    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};        
    let precondition_to_make_new_lev_set_lower_star     =   ConditionLowerStar{ simplex_bimap_sequential: & simplex_bimap_sequential };
    let precondition_to_make_new_lev_set_lower_edge     =   ConditionLowerEdge{ simplex_bimap_sequential: & simplex_bimap_sequential };    
    
    let analyze_dowker_dual                             =   false;

    //  LOWER STAR
    //  -------------------------------------------------------------------------------------------

    println!("\nLOWER STAR");
     let poly_complex_facets =       simplex_pipeline(
                                         &   simplex_sequence,
                                         &   barcode,
                                        &   ring,
                                        &   precondition_to_make_new_lev_set_lower_star,
                                        analyze_dowker_dual,
                                        save_dir_opt,
                                    );

    let analyze_dowker_dual =   false;
    analyze_fibre( 
        &   poly_complex_facets,
            ring.clone(),
            analyze_dowker_dual,
            save_dir_opt,
    );  


    // -------------------- SIMPLICIAL KLEIN, BARCODE WITH TWO DISJOINT FINITE BARS ------------

    //  Define the base space, barcode, and ring
    println!("\n SIMPLICIAL KLEIN, BARCODE WITH ONE FINITE BAR");
    println!("===============================================================================================\n");  
    let complex_facets      =   vec![    vec![0,3,4], vec![0,1,4], vec![1,2,4],vec![2,4,5],vec![0,2,5], vec![3,5,6],vec![3,4,6],vec![4,6,7],vec![4,5,7],vec![3,5,7],vec![0,5,6],vec![0,1,6], vec![1,6,7],vec![1,3,7],vec![1,2,3],vec![0,2,3]];

    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 2);    // build all faces up to dimension 1

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:5}],
        fin: vec![ BarFinite{dim:0,birth:1,death:2}, BarFinite{dim:0,birth:3,death:4}  ],
        ordinal: ordinate_unique_vals( & vec![0, 1, 2, 3,4,5 ] ),
    };
    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

    let simplex_bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence.clone() );

    //  Define any extra conditions on starting a new level set

    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};        
    let precondition_to_make_new_lev_set_lower_star     =   ConditionLowerStar{ simplex_bimap_sequential: & simplex_bimap_sequential };
    let precondition_to_make_new_lev_set_lower_edge     =   ConditionLowerEdge{ simplex_bimap_sequential: & simplex_bimap_sequential };    
    
    let analyze_dowker_dual                             =   false;

    //  LOWER STAR
    //  -------------------------------------------------------------------------------------------

    println!("\nLOWER STAR");
     let poly_complex_facets =       simplex_pipeline(
                                         &   simplex_sequence,
                                         &   barcode,
                                        &   ring,
                                        &   precondition_to_make_new_lev_set_lower_star,
                                        analyze_dowker_dual,
                                        save_dir_opt,
                                    );

    let analyze_dowker_dual =   false;
    analyze_fibre( 
        &   poly_complex_facets,
            ring.clone(),
            analyze_dowker_dual,
            save_dir_opt,
    );  

        // -------------------- SIMPLICIAL KLEIN, BARCODE WITH TWO OVERLAPPING FINITE BARS ------------

    //  Define the base space, barcode, and ring
    println!("\n SIMPLICIAL KLEIN, BARCODE WITH TWO OVERLAPPING FINITE BARS");
    println!("===============================================================================================\n");  
    let complex_facets      =   vec![    vec![0,3,4], vec![0,1,4], vec![1,2,4],vec![2,4,5],vec![0,2,5], vec![3,5,6],vec![3,4,6],vec![4,6,7],vec![4,5,7],vec![3,5,7],vec![0,5,6],vec![0,1,6], vec![1,6,7],vec![1,3,7],vec![1,2,3],vec![0,2,3]];

    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 2);    // build all faces up to dimension 1

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:5}],
        fin: vec![ BarFinite{dim:0,birth:1,death:3}, BarFinite{dim:0,birth:2,death:4}  ],
        ordinal: ordinate_unique_vals( & vec![0, 1, 2, 3,4,5 ] ),
    };
    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

    let simplex_bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence.clone() );

    //  Define any extra conditions on starting a new level set

    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};        
    let precondition_to_make_new_lev_set_lower_star     =   ConditionLowerStar{ simplex_bimap_sequential: & simplex_bimap_sequential };
    let precondition_to_make_new_lev_set_lower_edge     =   ConditionLowerEdge{ simplex_bimap_sequential: & simplex_bimap_sequential };    
    
    let analyze_dowker_dual                             =   false;

    //  LOWER STAR
    //  -------------------------------------------------------------------------------------------

    println!("\nLOWER STAR");
     let poly_complex_facets =       simplex_pipeline(
                                         &   simplex_sequence,
                                         &   barcode,
                                        &   ring,
                                        &   precondition_to_make_new_lev_set_lower_star,
                                        analyze_dowker_dual,
                                        save_dir_opt,
                                    );

    let analyze_dowker_dual =   false;
    analyze_fibre( 
        &   poly_complex_facets,
            ring.clone(),
            analyze_dowker_dual,
            save_dir_opt,
    );  

        // -------------------- SIMPLICIAL KLEIN, BARCODE WITH 1 FINITE BAR IN ANOTHER ------------

    //  Define the base space, barcode, and ring
    println!("\n SIMPLICIAL KLEIN,BARCODE WITH 1 FINITE BAR IN ANOTHER");
    println!("===============================================================================================\n");  
    let complex_facets      =   vec![    vec![0,3,4], vec![0,1,4], vec![1,2,4],vec![2,4,5],vec![0,2,5], vec![3,5,6],vec![3,4,6],vec![4,6,7],vec![4,5,7],vec![3,5,7],vec![0,5,6],vec![0,1,6], vec![1,6,7],vec![1,3,7],vec![1,2,3],vec![0,2,3]];

    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 2);    // build all faces up to dimension 1

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:5}],
        fin: vec![ BarFinite{dim:0,birth:2,death:3}, BarFinite{dim:0,birth:1,death:4}  ],
        ordinal: ordinate_unique_vals( & vec![0, 1, 2, 3,4,5 ] ),
    };
    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

    let simplex_bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence.clone() );

    //  Define any extra conditions on starting a new level set

    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};        
    let precondition_to_make_new_lev_set_lower_star     =   ConditionLowerStar{ simplex_bimap_sequential: & simplex_bimap_sequential };
    let precondition_to_make_new_lev_set_lower_edge     =   ConditionLowerEdge{ simplex_bimap_sequential: & simplex_bimap_sequential };    
    
    let analyze_dowker_dual                             =   false;

    //  LOWER STAR
    //  -------------------------------------------------------------------------------------------

    println!("\nLOWER STAR");
     let poly_complex_facets =       simplex_pipeline(
                                         &   simplex_sequence,
                                         &   barcode,
                                        &   ring,
                                        &   precondition_to_make_new_lev_set_lower_star,
                                        analyze_dowker_dual,
                                        save_dir_opt,
                                    );

    let analyze_dowker_dual =   false;
    analyze_fibre( 
        &   poly_complex_facets,
            ring.clone(),
            analyze_dowker_dual,
            save_dir_opt,
    );  
}
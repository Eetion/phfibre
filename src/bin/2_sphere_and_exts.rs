//  DESCRIPTION

//  The 2-simplex.

//  This file was shared by Jacob with Greg in August 2021; greg modified by switching the 1 to a 2 in
//  each copy of line
//  ```
// let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 2);    // build all faces up to dimension 1
// ```
// The results of the computation are coppied at the end of this file.


//  DEPENDENCIES

use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use phfibre::pipelines::{simplex_pipeline, analyze_fibre};
use phfibre::phfibre::{ConditionNone, ConditionLowerStar, ConditionLowerEdge};
use solar::utilities::sequences_and_ordinals::{ordinate_unique_vals, BiMapSequential};
use solar::cell_complexes::simplices_unweighted::facets::ordered_subsimplices_up_thru_dim_concatenated_vec; 
use std::iter::FromIterator;



//  ============================================================================================
//  ============================================================================================

//  SEE THE BOTTOM OF THIS FILE FOR OUTPUT FROM GREG'S EXECUTION

//  ============================================================================================
//  ============================================================================================




//  CODE

//   2-SPHERE (BOUNDARY OF 3-SIMPLEX) AND TRIVIAL ELEMENTARY EXTENSIONS (E.G. ADDING EDGES ON THE SIDE)
fn main() {

let save_dir_opt        =   None; // we will not save any files

 
// ------------------------------------------- CASE 1: 2-sphere and just 2 infinite bars ----------------------------

println!("\n 2-sphere and just 2 infinite bars");
println!("===============================================================================================\n");  


//  Define the base space, barcode, and ring

let complex_facets      =  vec![  vec![0, 1, 2, 3] ];
let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 2);    // build all faces up to dimension 1

let barcode             =   Barcode{
                                inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:2,birth:1} ],
                                fin: vec![],
                                ordinal: ordinate_unique_vals( & vec![ 0,1 ] ),
                            };

let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

let simplex_bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence.clone() );

//  Define any extra conditions on starting a new level set

let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};        
let precondition_to_make_new_lev_set_lower_star     =   ConditionLowerStar{ simplex_bimap_sequential: & simplex_bimap_sequential };
let precondition_to_make_new_lev_set_lower_edge     =   ConditionLowerEdge{ simplex_bimap_sequential: & simplex_bimap_sequential };    

let analyze_dowker_dual                             =   false;
                 
let poly_complex_facets =   simplex_pipeline(
                                &   simplex_sequence,
                                &   barcode,
                                &   ring,
                                &   precondition_to_make_new_lev_set_lower_none,
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

    
    //  RESULTS
    //  -------

    // BASE SPACE
    // simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
    // POLYHEDRAL COMPLEX FACETS
    // number of facets (total): 1920
    // number of facets (binned by dimension): [0, 0, 0, 0, 0, 0, 1920]
    // number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 0, 1920]
    // POLYHEDRAL COMPLEX CELLS
    // number of polytopes (total): 39566
    // number of polytopes (binned by dimension): [64, 858, 4480, 10860, 13320, 8064, 1920]
    // betti numbers (of polyhedral complex): [1, 0, 1, 0, 0, 0, 0]
    // INTERSECTIONS
    // number of pairs of intersecting facets, binned by the dimension of the intersection polytope: [555120, 310848, 139392, 55248, 19152, 5376, 0]
    // DOWKER NERVE COMPLEX FACETS
    // number of dowker nerve complex facets (total): 1920
    // number of dowker nerve complex facets (binned by dimension): [0, 0, 0, 0, 0, 0, 1920]
    // DOWKER NERVE COMPLEX CELLS
    // number of nerve dowker complex cells (total): 39566
    // number of nerve dowker complex cells (binned by dimension): [64, 858, 4480, 10860, 13320, 8064, 1920]
    // betti numbers (of dowker nerve): [1, 0, 1, 0, 0, 0, 0]

// ------------------------------------------- CASE 2: 2-sphere and 1 deg-0 bounded bar ----------------------------

    //  Define the base space, barcode, and ring

    println!("\n 2-sphere and 1 deg-0 bounded bar ");
println!("===============================================================================================\n");  


//  Define the base space, barcode, and ring

let complex_facets      =  vec![  vec![0, 1, 2, 3] ];
let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 2);    // build all faces up to dimension 1

let barcode             =   Barcode{
                                inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:2,birth:3} ],
                                fin: vec![ BarFinite{dim:0,birth:1,death:2}],
                                ordinal: ordinate_unique_vals( & vec![ 0,1,2,3] ),
                            };

let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

let simplex_bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence.clone() );

//  Define any extra conditions on starting a new level set

let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};        
let precondition_to_make_new_lev_set_lower_star     =   ConditionLowerStar{ simplex_bimap_sequential: & simplex_bimap_sequential };
let precondition_to_make_new_lev_set_lower_edge     =   ConditionLowerEdge{ simplex_bimap_sequential: & simplex_bimap_sequential };    

let analyze_dowker_dual                             =   false;

let save_dir_opt        =   None;                      
let poly_complex_facets =   simplex_pipeline(
                                &   simplex_sequence,
                                &   barcode,
                                &   ring,
                                &   precondition_to_make_new_lev_set_lower_none,
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
    // BASE SPACE
    // simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
    // POLYHEDRAL COMPLEX FACETS
    // number of facets (total): 18240
    // number of facets (binned by dimension): [0, 0, 0, 0, 0, 18240]
    // number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 18240]
    // POLYHEDRAL COMPLEX CELLS
    // number of polytopes (total): 322824
    // number of polytopes (binned by dimension): [3936, 29928, 84516, 113244, 72960, 18240]
    // betti numbers (of polyhedral complex): [1, 1, 1, 1, 0, 0]
    // INTERSECTIONS
    // number of pairs of intersecting facets, binned by the dimension of the intersection polytope: [1702584, 1221960, 565056, 205176, 56064, 0]
    // RUN INTERRUPTED (probably runtime error) "Killed: 9"

    
// ------------------------------------------- CASE 3: 2-sphere and 1 deg-1 bounded bar ----------------------------

      //  Define the base space, barcode, and ring

      println!("\n 2-sphere and 1 deg-1 bounded bar ");
      println!("===============================================================================================\n");  
      
      
      //  Define the base space, barcode, and ring
      
      let complex_facets      =  vec![  vec![0, 1, 2, 3] ];
      let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 2);    // build all faces up to dimension 1
      
      let barcode             =   Barcode{
                                      inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:2,birth:3} ],
                                      fin: vec![ BarFinite{dim:1,birth:1,death:2}],
                                      ordinal: ordinate_unique_vals( & vec![ 0,1,2,3 ] ),
                                  };
      
      let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();
      
      let simplex_bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence.clone() );
      
      //  Define any extra conditions on starting a new level set
      
      let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};        
      let precondition_to_make_new_lev_set_lower_star     =   ConditionLowerStar{ simplex_bimap_sequential: & simplex_bimap_sequential };
      let precondition_to_make_new_lev_set_lower_edge     =   ConditionLowerEdge{ simplex_bimap_sequential: & simplex_bimap_sequential };    
      
      let analyze_dowker_dual                             =   false;

      let save_dir_opt        =   None;     
      let poly_complex_facets =       simplex_pipeline(
                                           &   simplex_sequence,
                                          &   barcode,
                                        &   ring,
                                     &   precondition_to_make_new_lev_set_lower_none,
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
    
    // BASE SPACE
    // simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
    // POLYHEDRAL COMPLEX FACETS
    // number of facets (total): 18240
    // number of facets (binned by dimension): [0, 0, 0, 0, 0, 18240]
    // number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 18240]
    // POLYHEDRAL COMPLEX CELLS
    // number of polytopes (total): 322824
    // number of polytopes (binned by dimension): [3936, 29928, 84516, 113244, 72960, 18240]
    // betti numbers (of polyhedral complex): [1, 1, 1, 1, 0, 0]
    // INTERSECTIONS
    // number of pairs of intersecting facets, binned by the dimension of the intersection polytope: [1702584, 1221960, 565056, 205176, 56064, 0]
    // Killed: 9

// ------------------------------------------- CASE 4: 2-sphere, 1 deg-0 and 1 deg-1 bounded bars ----------------------------

   //  Define the base space, barcode, and ring

   println!("\n 2-sphere, 1 deg-0 and 1 deg-1 bounded bars ");
   println!("===============================================================================================\n");  
   
   
   //  Define the base space, barcode, and ring
   
   let complex_facets      =  vec![  vec![0, 1, 2, 3] ];
   let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 2);    // build all faces up to dimension 1
   
   let barcode             =   Barcode{
                                   inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:2,birth:5} ],
                                   fin: vec![ BarFinite{dim:0,birth:1,death:2}, BarFinite{dim:1,birth:3,death:4}],
                                   ordinal: ordinate_unique_vals( & vec![ 0,1,2,3,4,5 ] ),
                               };
   
   let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();
   
   let simplex_bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence.clone() );
   
   //  Define any extra conditions on starting a new level set
   
   let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};        
   let precondition_to_make_new_lev_set_lower_star     =   ConditionLowerStar{ simplex_bimap_sequential: & simplex_bimap_sequential };
   let precondition_to_make_new_lev_set_lower_edge     =   ConditionLowerEdge{ simplex_bimap_sequential: & simplex_bimap_sequential };    
   
   let analyze_dowker_dual                             =   false;
   
   let save_dir_opt        =   None;     
   let poly_complex_facets =       simplex_pipeline(
                                        &   simplex_sequence,
                                       &   barcode,
                                     &   ring,
                                  &   precondition_to_make_new_lev_set_lower_none,
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

// NEW COMPUTATION
//--------------------------------------------------------------------------------------------- 
//BASE SPACE
//simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
//BARCODE
//barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 2, birth: 5 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }, BarFinite { dim: 1, birth: 3, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4, 5], val_to_ord: {0: 0, 3: 3, 4: 4, 5: 5, 2: 2, 1: 1} } }
//TIME TO COMPUTE FIBRE FACETS
//Time elapsed to compute facets of PH fibre: 4965.946629583s

//ANALYSIS
//Each polytope facet has been checked for compatiblity with the given barcode.
//POLYHEDRAL COMPLEX FACETS
//number of facets (total): 162048
//number of facets (binned by dimension): [0, 0, 0, 0, 162048]
//number of facets (binned by number of vertices): [0, 0, 0, 0, 162048]
//POLYHEDRAL COMPLEX CELLS
//"Time elapsed to compute facets of PH fibre: 37.709009654s"
//number of polytopes (total): 1987032
//number of polytopes (binned by dimension): [87588, 422112, 743892, 571392, 162048]
//Killed: 9``
// --------------------------------------------------------------------------------------------- 
//------------ GREG's Computation
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 162048
// number of facets (binned by dimension): [0, 0, 0, 0, 162048]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 162048]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 39.787473234s"
// number of polytopes (total): 1987032
// number of polytopes (binned by dimension): [87588, 422112, 743892, 571392, 162048]
// Time elapsed to compute binary-coeff polyhedral betti numbers 50.653351456s
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 2, 27, 2, 0]


// ------------------------------------------- CASE 5: 2-sphere, 1 deg-1 and 1 deg-0 bounded bars ----------------------------

    //  Define the base space, barcode, and ring

    println!("\n 2-sphere, 1 deg-1 and 1 deg-0 bounded bars ");
    println!("===============================================================================================\n");  
    
    
    //  Define the base space, barcode, and ring
    
    let complex_facets      =  vec![  vec![0, 1, 2, 3] ];
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 2);    // build all faces up to dimension 1
    
    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:2,birth:5} ],
                                    fin: vec![ BarFinite{dim:1,birth:1,death:2}, BarFinite{dim:0,birth:3,death:4}],
                                    ordinal: ordinate_unique_vals( & vec![ 0,1,2,3,4,5 ] ),
                                };
    
    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();
    
    let simplex_bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence.clone() );
    
    //  Define any extra conditions on starting a new level set
    
    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};        
    let precondition_to_make_new_lev_set_lower_star     =   ConditionLowerStar{ simplex_bimap_sequential: & simplex_bimap_sequential };
    let precondition_to_make_new_lev_set_lower_edge     =   ConditionLowerEdge{ simplex_bimap_sequential: & simplex_bimap_sequential };    
    
    let analyze_dowker_dual                             =   false;
    
    let save_dir_opt        =   None;     
    let poly_complex_facets =       simplex_pipeline(
                                         &   simplex_sequence,
                                        &   barcode,
                                      &   ring,
                                   &   precondition_to_make_new_lev_set_lower_none,
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
 

// ------------------------------------------- CASE 6: 2-sphere, 2 consecutive deg-1 bounded bars ----------------------------

   //  Define the base space, barcode, and ring

   println!("\n 2-sphere, 2 consecutive deg-1 bounded bars ");
   println!("===============================================================================================\n");  
   
   
   //  Define the base space, barcode, and ring
   
   let complex_facets      =  vec![  vec![0, 1, 2, 3] ];
   let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 2);    // build all faces up to dimension 1
   
   let barcode             =   Barcode{
                                   inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:2,birth:5} ],
                                   fin: vec![ BarFinite{dim:1,birth:1,death:2}, BarFinite{dim:1,birth:3,death:4}],
                                   ordinal: ordinate_unique_vals( & vec![ 0,1,2,3,4,5 ] ),
                               };
   
   let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();
   
   let simplex_bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence.clone() );
   
   //  Define any extra conditions on starting a new level set
   
   let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};        
   let precondition_to_make_new_lev_set_lower_star     =   ConditionLowerStar{ simplex_bimap_sequential: & simplex_bimap_sequential };
   let precondition_to_make_new_lev_set_lower_edge     =   ConditionLowerEdge{ simplex_bimap_sequential: & simplex_bimap_sequential };    
   
   let analyze_dowker_dual                             =   false;
   
   let save_dir_opt        =   None;     
   let poly_complex_facets =       simplex_pipeline(
                                        &   simplex_sequence,
                                       &   barcode,
                                     &   ring,
                                  &   precondition_to_make_new_lev_set_lower_none,
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

// ------------------------------------------- CASE 7: 2-sphere, 2 overlapping deg-1 bounded bars ----------------------------

     //  Define the base space, barcode, and ring

     println!("\n2-sphere, 2 overlapping deg-1 bounded bars ");
     println!("===============================================================================================\n");  
     
     
     //  Define the base space, barcode, and ring
     
     let complex_facets      =  vec![  vec![0, 1, 2, 3] ];
     let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 2);    // build all faces up to dimension 1
     
     let barcode             =   Barcode{
                                     inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:2,birth:5} ],
                                     fin: vec![ BarFinite{dim:1,birth:1,death:3}, BarFinite{dim:1,birth:2,death:4}],
                                     ordinal: ordinate_unique_vals( & vec![ 0,1,2,3,4,5 ] ),
                                 };
     
     let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();
     
     let simplex_bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence.clone() );
     
     //  Define any extra conditions on starting a new level set
     
     let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};        
     let precondition_to_make_new_lev_set_lower_star     =   ConditionLowerStar{ simplex_bimap_sequential: & simplex_bimap_sequential };
     let precondition_to_make_new_lev_set_lower_edge     =   ConditionLowerEdge{ simplex_bimap_sequential: & simplex_bimap_sequential };    
     
     let analyze_dowker_dual                             =   false;

     let save_dir_opt        =   None;     
     let poly_complex_facets =       simplex_pipeline(
                                          &   simplex_sequence,
                                         &   barcode,
                                       &   ring,
                                    &   precondition_to_make_new_lev_set_lower_none,
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
// --------------------------------------- CASE 8: 2-sphere, 2 deg-1 bounded bars, one incldued into the other --------------

    //  Define the base space, barcode, and ring

    println!("\n 2-sphere, 2 deg-1 bounded bars, one incldued into the other ");
    println!("===============================================================================================\n");  
    
    
    //  Define the base space, barcode, and ring
    
    let complex_facets      =  vec![  vec![0, 1, 2, 3] ];
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 2);    // build all faces up to dimension 1
    
    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:2,birth:5} ],
                                    fin: vec![ BarFinite{dim:1,birth:2,death:3}, BarFinite{dim:1,birth:1,death:4}],
                                    ordinal: ordinate_unique_vals( & vec![ 0,1,2,3,4,5 ] ),
                                };
    
    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();
    
    let simplex_bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence.clone() );
    
    //  Define any extra conditions on starting a new level set
    
    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};        
    let precondition_to_make_new_lev_set_lower_star     =   ConditionLowerStar{ simplex_bimap_sequential: & simplex_bimap_sequential };
    let precondition_to_make_new_lev_set_lower_edge     =   ConditionLowerEdge{ simplex_bimap_sequential: & simplex_bimap_sequential };    
    
    let analyze_dowker_dual                             =   false;

    let save_dir_opt        =   None;     
    let poly_complex_facets =       simplex_pipeline(
                                         &   simplex_sequence,
                                        &   barcode,
                                      &   ring,
                                   &   precondition_to_make_new_lev_set_lower_none,
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

// ---------------------------------------CASE 9: 2-sphere and danggling edge: just 2 infinite bars ----------------------------

  //  Define the base space, barcode, and ring

    println!("\n 2-sphere and danggling edge: just 2 infinite bars ");
    println!("===============================================================================================\n");  
    
    
    //  Define the base space, barcode, and ring
    
    let complex_facets      =  vec![  vec![0, 1, 2, 3],vec![0, 4]   ];
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 2);    // build all faces up to dimension 1
    
    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:2,birth:1} ],
                                    fin: vec![ ],
                                    ordinal: ordinate_unique_vals( & vec![ 0,1 ] ),
                                };
    
    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();
    
    let simplex_bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence.clone() );
    
    //  Define any extra conditions on starting a new level set
    
    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};        
    let precondition_to_make_new_lev_set_lower_star     =   ConditionLowerStar{ simplex_bimap_sequential: & simplex_bimap_sequential };
    let precondition_to_make_new_lev_set_lower_edge     =   ConditionLowerEdge{ simplex_bimap_sequential: & simplex_bimap_sequential };    
    
    let analyze_dowker_dual                             =   false;
    
    let save_dir_opt        =   None;     
    let poly_complex_facets =       simplex_pipeline(
                                         &   simplex_sequence,
                                        &   barcode,
                                      &   ring,
                                   &   precondition_to_make_new_lev_set_lower_none,
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
    
    // NEW COMPUTATION
    //--------------------------------------------------------------------------------------------- 
    //BASE SPACE
    //simplices of the base space: [[0], [1], [2], [3], [4], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [3, 4], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
    //BARCODE
    //barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 2, birth: 1 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0, 1], val_to_ord: {1: 1, 0: 0} } }
    //TIME TO COMPUTE FIBRE FACETS
    //Time elapsed to compute facets of PH fibre: 10992.052942728s
    
    //ANALYSIS
    //Each polytope facet has been checked for compatiblity with the given barcode.
    //POLYHEDRAL COMPLEX FACETS
    //number of facets (total): 12816
    //number of facets (binned by dimension): [0, 0, 0, 0, 0, 0, 0, 12816]
    //number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 0, 0, 12816]
    //POLYHEDRAL COMPLEX CELLS
    //"Time elapsed to compute facets of PH fibre: 8.303524417s"
    //number of polytopes (total): 423452
    //number of polytopes (binned by dimension): [183, 3262, 22283, 72710, 127173, 122937, 62088, 12816]
    //betti numbers (of polyhedral complex, Z2 coefficients): [1, 0, 1, 0, 0, 0, 0, 0]
    //DOWKER NERVE COMPLEX FACETS
    //number of dowker nerve complex facets (total): 12816
    //number of dowker nerve complex facets (binned by dimension): [0, 0, 0, 0, 0, 0, 0, 10896, 0, 0, 0, 0, 0, 1920]
    //DOWKER NERVE COMPLEX CELLS
    //number of nerve dowker complex cells (total): 14627744
    //number of nerve dowker complex cells (binned by dimension): [183, 4978, 51737, 274888, 890313, 1940305, 2996712, 3350076, 2713680, 1572360, 634368, 169344, 26880, 1920]
    //betti numbers (of dowker nerve, user-specified ring coefficients): [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

/* 
    ------------------------------------- CASE 10: 2-sphere and danggling edge-edge: just 2 infinite bars ----------------------------

    //  Define the base space, barcode, and ring

    let complex_facets      =   vec![  vec![0, 1, 2, 3], vec![3,4],  vec![4,5]  ];
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 2);

    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:2,birth:1} ],
                                    fin: vec![ ],
                                    ordinal: ordinate_unique_vals( & vec![ 0, 1 ] ),
                                };

    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

    simplex_pipeline(
        &   simplex_sequence,
        &   barcode,
        &   ring
    );
    
// ------------------------------------- CASE 11: 2-sphere and 4 danggling edges: just 2 infinite bars ----------------------------

    //  Define the base space, barcode, and ring

    let complex_facets      =   vec![  vec![0, 1, 2, 3], vec![0,4],  vec![1,5], vec![2,6], vec![3,7]  ];
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 2);

    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:2,birth:1} ],
                                    fin: vec![ ],
                                    ordinal: ordinate_unique_vals( & vec![ 0, 1 ] ),
                                };

    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

    simplex_pipeline(
        &   simplex_sequence,
        &   barcode,
        &   ring
    );

// -------------------------- CASE 12: 2-sphere and a 3-simplex joined at a vertex: just 2 infinite bars --------------------------

    //  Define the base space, barcode, and ring

    let complex_facets      =   vec![  vec![0, 1, 2],vec![0, 1, 3], vec![0, 2, 3] ,  vec![1, 2, 3], vec![3,4,5,6] ];
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 3);

    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:2,birth:1} ],
                                    fin: vec![ ],
                                    ordinal: ordinate_unique_vals( & vec![ 0, 1 ] ),
                                };

    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

    simplex_pipeline(
        &   simplex_sequence,
        &   barcode,
        &   ring
    );
*/

}


//  ============================================================================================
//  ============================================================================================

//  OUTPUT FROM GREG'S EXECUTION OF THIS FILE

//  ============================================================================================
//  ============================================================================================



// 2-sphere and just 2 infinite bars
// ===============================================================================================

// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 2, birth: 1 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0, 1], val_to_ord: {0: 0, 1: 1} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 3.568017168s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 1920
// number of facets (binned by dimension): [0, 0, 0, 0, 0, 0, 1920]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 0, 1920]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 637.488565ms"
// number of polytopes (total): 39566
// number of polytopes (binned by dimension): [64, 858, 4480, 10860, 13320, 8064, 1920]
// Time elapsed to compute binary-coeff polyhedral betti numbers 68.463363ms
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 0, 1, 0, 0, 0, 0]

//  2-sphere and 1 deg-0 bounded bar 
// ===============================================================================================

// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 2, birth: 3 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3], val_to_ord: {2: 2, 1: 1, 3: 3, 0: 0} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 139.665238703s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 18240
// number of facets (binned by dimension): [0, 0, 0, 0, 0, 18240]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 18240]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 6.690194678s"
// number of polytopes (total): 322824
// number of polytopes (binned by dimension): [3936, 29928, 84516, 113244, 72960, 18240]
// Time elapsed to compute binary-coeff polyhedral betti numbers 3.138876354s
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 1, 1, 1, 0, 0]

//  2-sphere and 1 deg-1 bounded bar 
// ===============================================================================================

// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 2, birth: 3 }], fin: [BarFinite { dim: 1, birth: 1, death: 2 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3], val_to_ord: {3: 3, 0: 0, 1: 1, 2: 2} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 95.281256833s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 18240
// number of facets (binned by dimension): [0, 0, 0, 0, 0, 18240]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 18240]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 7.869258681s"
// number of polytopes (total): 322824
// number of polytopes (binned by dimension): [3936, 29928, 84516, 113244, 72960, 18240]
// Time elapsed to compute binary-coeff polyhedral betti numbers 1.174518212s
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 1, 1, 1, 0, 0]

//  2-sphere, 1 deg-0 and 1 deg-1 bounded bars 
// ===============================================================================================

// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 2, birth: 5 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }, BarFinite { dim: 1, birth: 3, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4, 5], val_to_ord: {0: 0, 3: 3, 4: 4, 5: 5, 2: 2, 1: 1} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 6560.816409761s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 162048
// number of facets (binned by dimension): [0, 0, 0, 0, 162048]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 162048]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 36.743753836s"
// number of polytopes (total): 1987032
// number of polytopes (binned by dimension): [87588, 422112, 743892, 571392, 162048]
// Time elapsed to compute binary-coeff polyhedral betti numbers 47.172512373s
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 2, 27, 2, 0]

//  2-sphere, 1 deg-1 and 1 deg-0 bounded bars 
// ===============================================================================================

// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 2, birth: 5 }], fin: [BarFinite { dim: 1, birth: 1, death: 2 }, BarFinite { dim: 0, birth: 3, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4, 5], val_to_ord: {0: 0, 1: 1, 2: 2, 4: 4, 3: 3, 5: 5} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 81.041571ms

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 576
// number of facets (binned by dimension): [0, 0, 0, 0, 576]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 576]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 103.193421ms"
// number of polytopes (total): 7056
// number of polytopes (binned by dimension): [324, 1512, 2628, 2016, 576]
// Time elapsed to compute binary-coeff polyhedral betti numbers 7.235522ms
// betti numbers (of polyhedral complex, Z2 coefficients): [4, 8, 4, 0, 0]

//  2-sphere, 2 consecutive deg-1 bounded bars 
// ===============================================================================================

// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 2, birth: 5 }], fin: [BarFinite { dim: 1, birth: 1, death: 2 }, BarFinite { dim: 1, birth: 3, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4, 5], val_to_ord: {2: 2, 4: 4, 0: 0, 3: 3, 1: 1, 5: 5} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 17.615275192s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 14592
// number of facets (binned by dimension): [0, 0, 0, 0, 14592]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 14592]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 2.579674965s"
// number of polytopes (total): 167616
// number of polytopes (binned by dimension): [6696, 33984, 62520, 49824, 14592]
// Time elapsed to compute binary-coeff polyhedral betti numbers 270.114223ms
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 14, 13, 0, 0]

// 2-sphere, 2 overlapping deg-1 bounded bars 
// ===============================================================================================

// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 2, birth: 5 }], fin: [BarFinite { dim: 1, birth: 1, death: 3 }, BarFinite { dim: 1, birth: 2, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4, 5], val_to_ord: {0: 0, 3: 3, 1: 1, 4: 4, 5: 5, 2: 2} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 21.400452311s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 14976
// number of facets (binned by dimension): [0, 0, 0, 0, 14976]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 14976]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 2.490736655s"
// number of polytopes (total): 169728
// number of polytopes (binned by dimension): [6552, 33984, 63336, 50880, 14976]
// Time elapsed to compute binary-coeff polyhedral betti numbers 254.198976ms
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 14, 13, 0, 0]

//  2-sphere, 2 deg-1 bounded bars, one incldued into the other 
// ===============================================================================================

// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 2, birth: 5 }], fin: [BarFinite { dim: 1, birth: 2, death: 3 }, BarFinite { dim: 1, birth: 1, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4, 5], val_to_ord: {5: 5, 1: 1, 4: 4, 0: 0, 2: 2, 3: 3} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 70.581630257s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 26496
// number of facets (binned by dimension): [0, 0, 0, 0, 26496]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 26496]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 4.571141822s"
// number of polytopes (total): 297696
// number of polytopes (binned by dimension): [11304, 59184, 111048, 89664, 26496]
// Time elapsed to compute binary-coeff polyhedral betti numbers 449.00208ms
// betti numbers (of polyhedral complex, Z2 coefficients): [2, 28, 26, 0, 0]

//  2-sphere and danggling edge: just 2 infinite bars 
// ===============================================================================================

// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [4], [0, 1], [0, 2], [0, 3], [0, 4], [1, 2], [1, 3], [2, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 2, birth: 1 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0, 1], val_to_ord: {0: 0, 1: 1} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 4612.49403266s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 12816
// number of facets (binned by dimension): [0, 0, 0, 0, 0, 0, 0, 12816]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 0, 0, 12816]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 8.417333535s"
// number of polytopes (total): 423452
// number of polytopes (binned by dimension): [183, 3262, 22283, 72710, 127173, 122937, 62088, 12816]
// Time elapsed to compute binary-coeff polyhedral betti numbers 1.006062869s
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 0, 1, 0, 0, 0, 0, 0]
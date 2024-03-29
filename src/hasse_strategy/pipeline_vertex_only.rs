
use crate::hasse_strategy::compute_vertex_only::{Node, verify_that_barcode_is_compatible, ExtraConditionToStartNewLevSet, explore_vertex_only};
use crate::intervals_and_ordinals::{Barcode, BarcodeInverse };
use crate::polytope::object_def::Polytope;
use crate::polytope::faces::{poly_complex_facets_to_whole_complex_bimapsequential};
use crate::polytope::intersection::{polytope_intersection};
use crate::polytope::differential::polyhedral_boundary_matrix_binary_coeff;
use crate::polytope::nerve::dowker_nerve_complex_facets;
use crate::rank_calculations::chain_cx_rank_nullity;
use solar::utilities::sequences_and_ordinals::{BiMapSequential};
use solar::utilities::statistics::histogram;
use solar::utilities::cell_complexes::simplices_unweighted::facets::{ ordered_subsimplices_up_thru_dim_concatenated_vec}; 
use solar::utilities::cell_complexes::simplices_unweighted::boundary_matrices::{boundary_matrix_from_complex_facets};   
use solar::rings::field_prime::GF2;
use solar::rings::ring::{Semiring, Ring, DivisionRing};
use std::collections::HashMap;
use std::iter::FromIterator;
use std::fmt::Debug;
use std::hash::Hash;

use itertools::Itertools;

use serde::{Deserialize, Serialize};
use serde_json::{json, Value};

use ordered_float::OrderedFloat;

type F  =   num::rational::Ratio<i64>;


// use solar::utilities::sequences_and_ordinals::BiMapSequential;
// use crate::phfibre::{Node, explore};
// use crate::intervals_and_ordinals::{Barcode, BarcodeInverse};
// use crate::polytope::object_def::Polytope;
// use solar::utilities::cell_complexes::simplices_unweighted::facets::{    
//     ordered_subsimplices_up_thru_dim_concatenated_vec 
// }; 
// use solar::utilities::cell_complexes::simplices_unweighted::boundary_matrices::{    
//     boundary_matrix_from_complex_facets 
// };   
// use ordered_float::OrderedFloat;


// type RingEltRational = OrderedFloat<f64>;
// type RingOpRational = solar::rings::ring_native::NativeDivisionRing<RingEltRational>;



//  --------------------------------------------------------------------------------------------
//  BOUNDARY MATRIX PIPELINE
//  --------------------------------------------------------------------------------------------


/// Compute the PH fibre of a given barcode for a given boundary matrix over a given ring.
/// 
/// The boundary matrix is stored in column-major vec-of-vec format.  Entries in each column should
/// appear in sorted order, according to row index.  The homological degree of
/// each basis vector (i.e. chain) is recorded in `cell_dims`.
pub fn  boundary_matrix_pipeline_vertex_only< FilRaw, RingOp, RingElt, PreconditionToMakeNewLevSet >(
            boundary:               &   Vec< Vec< ( usize, RingElt ) > >,
            cell_id_to_prereqs:         Option< & Vec< Vec< usize > > >,
            cell_dims:              &   Vec< usize >,            
            barcode:                &   Barcode< FilRaw >,
            ring:                   &   RingOp,           
            precondition_to_make_new_lev_set:   &PreconditionToMakeNewLevSet, 
            analyze_dowker_dual:        bool,
        )
    where   RingOp:     Ring<RingElt> + Semiring<RingElt> + DivisionRing<RingElt> + Clone,
            RingElt:    Clone + Debug + Ord,
            FilRaw:     Clone + Debug + Ord + Hash, 
            PreconditionToMakeNewLevSet:    ExtraConditionToStartNewLevSet + Clone + Debug, 
{
    let barcode_inverse     =   BarcodeInverse::from_barcode( & barcode );
    let root =  Node::make_root(
                        boundary.clone(),
                    &   barcode,
                    &   barcode_inverse,
                        cell_id_to_prereqs,
                    &   cell_dims,  
                        ring.clone(),
                        precondition_to_make_new_lev_set,
                );

    let mut poly_complex_facets     =   Vec::new();                                

    //  GATHER RESULTS
    //  --------------

    println!("\nNEW COMPUTATION\n--------------------------------------------------------------------------------------------- ");

    println!("BOUNDARY MATRIX");
    println!("boundary matrix sparse columns: {:?}", & boundary);
    println!("BARCODE");
    println!("barcode: {:?}", &barcode);

    let start = std::time::Instant::now();
    explore_vertex_only( & root, &mut poly_complex_facets );
    let duration = start.elapsed();
    println!("TIME TO COMPUTE FIBRE FACETS");
    println!("Time elapsed to compute facets of PH fibre: {:?}", duration);
    
    println!("\nANALYSIS");

    //  VERIFY COMPATIBILITY
    //  --------------------

    for poly in poly_complex_facets.iter() {
        verify_that_barcode_is_compatible(
            &   root,
            &   poly
        );
    }
    println!("Each polytope facet has been checked for compatiblity with the given barcode.");

    
    analyze_fibre_vertex_only( 
        &   poly_complex_facets,
            ring.clone(),
            analyze_dowker_dual,
            // save_dir_opt,
    ); 
        
}


//  --------------------------------------------------------------------------------------------
//  SIMPLEX PIPELINE
//  --------------------------------------------------------------------------------------------

//  -------------------------------------------------------------------
//  BEGIN:  THIS COMMENTED SECTION CAN BE DELETED; THE FUNCTION DEFS HAVE BEEN MOVED TO STRUCTS IN phfibre.rs 

// /// Determine whether all vertices in the 0-skeleton of a simplex have been assigned to level sets.
// /// 
// /// # Arguments
// /// 
// /// * poly: a polytope
// /// * simplex_bimap_sequential: assigns a unique id number to each simplex
// /// * simplex: the simplex in question
// pub fn  skel_0_added(
//             poly:                       & Polytope,
//             simplex_bimap_sequential:   & BiMapSequential< Vec< usize > >,
//             simplex:                    & Vec< usize >,
//         )
//         ->  
//         bool 
// {
//     for ( simplex_index, simplex ) in simplex_bimap_sequential.ord_to_val.iter().enumerate() {
//         if  ! poly.cell_id_to_has_lev_set( simplex_index )  // simplex has not been assinged to a level set
//             &&
//             simplex                                             // but all its vertices have
//                 .iter()
//                 .cloned()                                       
//                 .combinations( 1 )
//                 .all(   |x| 
//                         poly.cell_id_to_has_lev_set( 
//                             simplex_bimap_sequential.ord( &x ).unwrap()                                
//                         ) 
//                 ) 
//         {
//             return false
//         }
//     }
//     true 
// }

// /// Determine whether all vertices AND edges in the 1-skeleton of a simplex have been assigned to level sets.
// /// 
// /// # Arguments
// /// 
// /// * poly: a polytope
// /// * simplex_bimap_sequential: assigns a unique id number to each simplex
// /// * simplex: the simplex in question
// pub fn  skel_1_added(
//             poly:                       & Polytope,
//             simplex_bimap_sequential:   & BiMapSequential< Vec< usize > >,
//             simplex:                    & Vec< usize >,
//         )
//         ->  
//         bool 
// {
//     for ( simplex_index, simplex ) in simplex_bimap_sequential.ord_to_val.iter().enumerate() {
//         if  ! poly.cell_id_to_has_lev_set( simplex_index )  // simplex has not been assinged to a level set
//             &&
//             simplex                                         // but all its vertices have
//                 .iter() 
//                 .cloned()                                         
//                 .combinations( 1 )
//                 .all(   |x| 
//                         poly.cell_id_to_has_lev_set( 
//                             simplex_bimap_sequential.ord( &x ).unwrap()                                
//                         ) 
//                 ) 
//             &&
//             simplex                                         // but all its edges have, also
//                 .iter()
//                 .cloned()
//                 .combinations( 2 )
//                 .all(   |x| 
//                         poly.cell_id_to_has_lev_set( 
//                             simplex_bimap_sequential.ord( &x ).unwrap()                                
//                         ) 
//                 ) 
//         {
//             return false
//         }
//     }
//     true  // return true if no violations have been found
// }
//  END:  THIS COMMENTED SECTION CAN BE DELETED; THE FUNCTION DEFS HAVE BEEN MOVED TO STRUCTS IN phfibre.rs 
//  -------------------------------------------------------------------




/// Compute the PH fibre of a given barcode for a given simplicial complex over a given ring.
/// 
/// The `simplex_sequence` parameter should include all simplices in the complex.
pub fn  simplex_pipeline_vertex_only< FilRaw, RingOp, RingElt, PreconditionToMakeNewLevSet >(
            simplex_sequence:       &   Vec< Vec< usize > >,
            barcode:                &   Barcode< FilRaw >,
            ring:                   &   RingOp,     
            precondition_to_make_new_lev_set:   &   PreconditionToMakeNewLevSet,
            analyze_dowker_dual:        bool,
        )
        ->
        Vec< Polytope >
    where   RingOp:     Ring<RingElt> + Semiring<RingElt> + DivisionRing<RingElt> + Clone,
            RingElt:    Clone + Debug + Ord,
            FilRaw:     Clone + Debug + Ord + Hash, 
            PreconditionToMakeNewLevSet: ExtraConditionToStartNewLevSet + Clone + Debug, 
{

    //  DEFINE THE CELL DIMENSIONS + BOUNDARY MATRIX
    //  --------------------------------------------

    let cell_dims: Vec<_>   =   simplex_sequence.iter().map(|x| x.len()-1 ).collect();

    let bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence.clone() );
    let boundary            =   boundary_matrix_from_complex_facets( &bimap_sequential, ring.clone()); 


    //  DEFINE THE BARCODE + INVERSE BARCODE
    //  ------------------------------------
    
    let barcode_inverse     =   BarcodeInverse::from_barcode( & barcode );
            

    //  NO CELL HAS TO SATISFY A SPECIAL REQUIREMENT TO BE ADDED (ONLY "VANILLA" REQUIREMENTS)
    //  -----------------------------------------------------------------

    let cell_id_to_prereqs  =   None;

    //  FILTRATION CONSTRAINTS  <---- THIS SECTION SHOULD BE DEPRECATED NOW; PROBABLY OK TO DELETE
    //  ----------------------

    // // Panic if an invalid choice is given
    // if ! vec!["lower_star", "lower_edge", "none"].contains( &&filtration_constraints ) {
    //     panic!("Argument `filtration_constraints` must be \"lower_star\", \"lower_edge\", or \"none\"");
    // }

    // //  The vacuous default constraint that will be use if we don't use lower_star or lower_edge 
    // let precondition_to_make_new_lev_set    =               | poly: Polytope | -> bool { true } ;

    // if filtration_constraints == "lower_star" { 
    //     let precondition_to_make_new_lev_set    = | poly: Polytope | -> bool 
    //         { 
    //             return simplex_sequence.iter().all( |x| skel_0_added(&poly, &bimap_sequential, x) )
    //         };  
    // }

    // if filtration_constraints == "lower_edge" { 
    //     let precondition_to_make_new_lev_set    = | poly: Polytope | -> bool 
    //         { 
    //             return simplex_sequence.iter().all( |x| skel_1_added(&poly, &bimap_sequential, x) )
    //         };  
    // }    


    //  DEFINE THE ROOT NODE + RESULTS VECTOR
    //  -------------------------------------

    let root =  Node::make_root(
                        boundary.clone(),
                    &   barcode,
                    &   barcode_inverse,
                        cell_id_to_prereqs,
                    &   cell_dims,  
                        ring.clone(),
                        precondition_to_make_new_lev_set,
                );

    let mut poly_complex_facets     =   Vec::new();                                

    //  GATHER RESULTS
    //  --------------

    println!("--------------------------------------------------------------------------------------------- ");

    println!("BASE SPACE");
    println!("simplices of the base space: {:?}", & simplex_sequence);
    println!("BARCODE");
    println!("barcode: {:?}", &barcode);

    let start = std::time::Instant::now();
    explore_vertex_only( & root, &mut poly_complex_facets );
    let duration = start.elapsed();
    println!("TIME TO COMPUTE FIBRE FACETS");
    println!("Time elapsed to compute facets of PH fibre: {:?}", duration);
    
    println!("\nANALYSIS");

    //  VERIFY COMPATIBILITY
    //  --------------------

    for poly in poly_complex_facets.iter() {
        verify_that_barcode_is_compatible(
            &   root,
            &   poly
        );
    }
    println!("Each polytope facet has been checked for compatiblity with the given barcode.");

    poly_complex_facets

    
}                    




//  --------------------------------------------------------------------------------------------
//  OUTPUT ANALYSIS / DISPLAY
//  --------------------------------------------------------------------------------------------


pub fn analyze_fibre_vertex_only< RingOp, RingElt > (
            poly_complex_facets:    &   Vec< Polytope >, 
            ring:                       RingOp, 
            analyze_dowker_dual:        bool,           
        )
    where   RingOp:     Ring<RingElt> + Semiring<RingElt> + DivisionRing<RingElt> + Clone,
            RingElt:    Clone + Debug + Ord,
{
    //  REPORT
    //  ------

    println!("POLYHEDRAL COMPLEX FACETS");
    println!("number of facets (total): {:?}", poly_complex_facets.len() ); 
    println!("number of facets (binned by dimension): {:?}", histogram( poly_complex_facets.iter().map(|x| x.dim_cellagnostic().unwrap() ) ) );  
    println!("number of facets (binned by number of vertices): {:?}", histogram( poly_complex_facets.iter().map(|x| x.dim_cellagnostic().unwrap() ) ) );      

   
    println!("POLYHEDRAL COMPLEX CELLS"); 

    let start = std::time::Instant::now();    
    let poly_complex_bimapsequential    =   poly_complex_facets_to_whole_complex_bimapsequential( & poly_complex_facets );
    let duration = start.elapsed();
    let time_message_poly_faces = format!("Time elapsed to compute the faces of PH fibre (given the facets): {:?}", duration);   
    println!("{:?}", time_message_poly_faces); 

    // SQUARE EXAMPLE
    // (old face enumerator) Time elapsed to enumerate polytope faces: 109.550845ms
    // (new face enumerator) Time elapsed to enumerate polytope faces: 24.799084ms
    // 2-SKEL OF 3-SIMPLEX EXAMPLE
    // (old face enumerator) Time elapsed to enumerate polytope faces: 16.452426366s
    // (new face enumerator) Time elapsed to enumerate polytope faces: 605.019797ms

    let poly_complex_dims               =   Vec::from_iter( poly_complex_bimapsequential.ord_to_val.iter().map(|x| x.dim_cellagnostic().unwrap() ) );    
      
    println!("number of polytopes (total): {:?}", poly_complex_bimapsequential.ord_to_val.len() ); 
    println!("number of polytopes (binned by dimension): {:?}", histogram( poly_complex_dims.iter().cloned() ) );      

    let poly_complex_differential           =   polyhedral_boundary_matrix_binary_coeff( & poly_complex_bimapsequential );
   
    let start = std::time::Instant::now();      
    let poly_complex_rank_nullity           =   chain_cx_rank_nullity(
                                                    & poly_complex_differential,
                                                    & poly_complex_dims,
                                                    GF2{}
                                                );  
    let duration = start.elapsed();   
    println!("Time elapsed to compute binary-coeff polyhedral betti numbers {:?}", duration);                                                


    let poly_complex_betti_vec              =   poly_complex_rank_nullity.rank_homology_vec();         
    println!("betti numbers (of polyhedral complex, Z2 coefficients): {:?}", &poly_complex_betti_vec);


    // THIS WORKS FINE BUT IT TAKES TIME
    // let mut intersection_dim_bins  =   Vec::from_iter( std::iter::repeat(0).take( poly_complex_dim_top + 1) );    
    // for count_a in 0 .. num_facets{
    //     for count_b in count_a + 1 .. num_facets {
    //         if let Some( intersection_poly ) = polytope_intersection( &poly_complex_facets[count_a], &poly_complex_facets[count_b] ){
    //             dismat[ count_a ][ count_b ] = OrderedFloat( 0.1 );
    //             dismat[ count_b ][ count_a ] = OrderedFloat( 0.1 );  
    //             let dim                 =   intersection_poly.dim_cellagnostic().unwrap();                
    //             intersection_dim_bins[ dim ] += 1              
    //         }
    //     }
    // }    
    // println!("INTERSECTIONS");        
    // println!("number of pairs of intersecting facets, binned by the dimension of the intersection polytope: {:?}", &intersection_dim_bins );    

    
    // dowker complex
    if analyze_dowker_dual {
        println!("DOWKER NERVE COMPLEX FACETS");       

        let dowker_complex_facets               =       dowker_nerve_complex_facets(
                                                            & poly_complex_facets,
                                                            & poly_complex_bimapsequential.val_to_ord
                                                        );                                
                                                        
        println!("number of dowker nerve complex facets (total): {:?}", dowker_complex_facets.len() );
        println!("number of dowker nerve complex facets (binned by dimension): {:?}", histogram( dowker_complex_facets.iter().map(|x| x.len()-1 ) ) );    

        println!("DOWKER NERVE COMPLEX CELLS");   
        
        let dowker_complex_dim_top              =   dowker_complex_facets.iter().map(|x| x.len()-1 ).max().unwrap();

        let dowker_complex_simplex_sequence     =   ordered_subsimplices_up_thru_dim_concatenated_vec( 
                                                        & dowker_complex_facets, 
                                                        dowker_complex_dim_top
                                                    );    
        let dowker_complex_cell_dims: Vec<_>    =   dowker_complex_simplex_sequence.iter().map(|x| x.len()-1 ).collect();

        let dowker_complex_bimap_sequential     =   BiMapSequential::from_vec( dowker_complex_simplex_sequence.clone() );
        let dowker_complex_boundary             =   boundary_matrix_from_complex_facets( &dowker_complex_bimap_sequential, ring.clone()); 

        let start = std::time::Instant::now();       
        let dowker_complex_rank_nullity         =   chain_cx_rank_nullity(
                                                        & dowker_complex_boundary,
                                                        & dowker_complex_cell_dims,
                                                        ring.clone()
                                                    );
        let duration = start.elapsed();   
        println!("Time elapsed to compute nerve dowker dual betti numbers (user specified coeff) {:?}", duration);       
        
        println!("number of nerve dowker complex cells (total): {:?}", dowker_complex_cell_dims.len() );
        println!("number of nerve dowker complex cells (binned by dimension): {:?}", histogram( dowker_complex_cell_dims.iter().cloned() ));
        println!("betti numbers (of dowker nerve, user-specified ring coefficients): {:?}", &dowker_complex_rank_nullity.rank_homology_vec() );    
    }
    
}                    

// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 2, birth: 1 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0, 1], val_to_ord: {0: 0, 1: 1} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 11.626907392s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 1920
// number of facets (binned by dimension): [0, 0, 0, 0, 0, 0, 1920]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 0, 1920]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 604.425412ms"
// number of polytopes (total): 39566
// number of polytopes (binned by dimension): [64, 858, 4480, 10860, 13320, 8064, 1920]
// Time elapsed to compute binary-coeff polyhedral betti numbers 67.473867ms
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 0, 1, 0, 0, 0, 0]
// DOWKER NERVE COMPLEX FACETS
// number of dowker nerve complex facets (total): 1920
// number of dowker nerve complex facets (binned by dimension): [0, 0, 0, 0, 0, 0, 1920]
// DOWKER NERVE COMPLEX CELLS
// Time elapsed to compute nerve dowker dual betti numbers (user specified coeff) 102.759723ms
// number of nerve dowker complex cells (total): 39566
// number of nerve dowker complex cells (binned by dimension): [64, 858, 4480, 10860, 13320, 8064, 1920]
// betti numbers (of dowker nerve, user-specified ring coefficients): [1, 0, 1, 0, 0, 0, 0]
// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3], [0, 1, 2, 3]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0], val_to_ord: {0: 0} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 14.152683553s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 1920
// number of facets (binned by dimension): [0, 0, 0, 0, 0, 0, 0, 1920]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 0, 0, 1920]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 1.086333144s"
// number of polytopes (total): 79133
// number of polytopes (binned by dimension): [65, 922, 5338, 15340, 24180, 21384, 9984, 1920]
// Time elapsed to compute binary-coeff polyhedral betti numbers 85.649335ms
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 0, 0, 0, 0, 0, 0, 0]
// DOWKER NERVE COMPLEX FACETS
// number of dowker nerve complex facets (total): 1920
// number of dowker nerve complex facets (binned by dimension): [0, 0, 0, 0, 0, 0, 0, 1920]
// DOWKER NERVE COMPLEX CELLS
// Time elapsed to compute nerve dowker dual betti numbers (user specified coeff) 158.716406ms
// number of nerve dowker complex cells (total): 79133
// number of nerve dowker complex cells (binned by dimension): [65, 922, 5338, 15340, 24180, 21384, 9984, 1920]
// betti numbers (of dowker nerve, user-specified ring coefficients): [1, 0, 0, 0, 0, 0, 0, 0]
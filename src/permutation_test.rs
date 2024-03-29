
use crate::phfibre::{Node, explore, verify_that_barcode_is_compatible, ConditionNone};
use crate::intervals_and_ordinals::{Barcode, BarcodeInverse, BarInfinite, BarFinite, to_ordered_float};
use crate::polytope::faces::{polys_faces};
use crate::polytope::object_def::Polytope;
use crate::utilities::simplices_with_relabeled_vertices;
use solar::utilities::indexing_and_bijection::{compose_f_after_g, inverse_perm};
use solar::utilities::sequences_and_ordinals::BiMapSequential;
use solar::utilities::cell_complexes::simplices_unweighted::facets::{    
    ordered_subsimplices_up_thru_dim_concatenated_vec, 
}; 
use solar::utilities::cell_complexes::simplices_unweighted::boundary_matrices::{    
    boundary_matrix_from_complex_facets     }; 
use solar::utilities::cell_complexes::simplices_unweighted::simplex::{    
    simplex_perm_o2n_from_vertex_perm_o2n   }; 
use solar::rings::ring::{Semiring, Ring, DivisionRing};    
use num::rational::Ratio;
use ordered_float::OrderedFloat;
use std::collections::HashSet;
use std::iter::FromIterator;
use std::fmt::Debug;
use std::hash::Hash;
use itertools::Itertools;

type RingEltRational = OrderedFloat<f64>;
type RingOpRational = solar::rings::ring_native::NativeDivisionRing<RingEltRational>;



pub fn  superficial_permutation_test< FilRaw, RingOp, RingElt > (
            complex_facets:         & Vec< Vec< usize > >,
            barcode:                & Barcode< FilRaw >,
            ring:                   & RingOp,
            vertex_permutations:    & Vec< Vec< usize >>
        )
        where   RingOp:     Ring<RingElt> + Semiring<RingElt> + DivisionRing<RingElt> + Clone,
                RingElt:    Clone + Debug + Ord,
                FilRaw:     Clone + Debug + Ord + Hash,        
{ 



    //  COMPUTE FACETS
    //  --------------------------------------------                                    

    let mut fibre_facets    =   fibre_facets_from_complex_facets(
                                    & complex_facets,
                                    & barcode,
                                    ring,
                                );                                   

    //  GET THE ORIGINAL SEQUENCE OF SIMPLICES
    let max_dim             =   complex_facets.iter().map(|x| x.len()-1 ).max().unwrap();
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, max_dim); 
    println!("{:?}", &simplex_sequence);

    //  DEFINE A PERMUTAITON ON SIMPLICES INDUCED BY A PERMUTATION ON VERTICES
    for perm_v_o2n in vertex_permutations.iter() {
        let perm_s_o2n      =   simplex_perm_o2n_from_vertex_perm_o2n( &simplex_sequence, &perm_v_o2n );
        let perm_s_n2o      =   inverse_perm( & perm_s_o2n );    
    
        let mut fibre_facets_under_superficial_perm
               =    Vec::from_iter(
                        fibre_facets
                            .iter()
                            .map(   |x| 
                                    Polytope{
                                        data_l_to_fmin: x.data_l_to_fmin.clone(),
                                        data_c_to_l:    compose_f_after_g(
                                                            & x.data_c_to_l, 
                                                            & perm_s_n2o
                                                        )
                                    }
                            )
                    );
        fibre_facets_under_superficial_perm.sort();
        fibre_facets.sort();
        assert_eq!( &fibre_facets, &fibre_facets_under_superficial_perm);
    }
}  



pub fn  deep_permutation_test< FilRaw, RingOp, RingElt > (
            complex_facets:         & Vec< Vec< usize > >,
            barcode:                & Barcode< FilRaw >,
            ring:                   & RingOp,
            vertex_permutations:    & Vec< Vec< usize >>
        )
        where   RingOp:     Ring<RingElt> + Semiring<RingElt> + DivisionRing<RingElt> + Clone,
                RingElt:    Clone + Debug + Ord,
                FilRaw:     Clone + Debug + Ord + Hash,        
{ 

    //  COMPUTE FACETS
    //  --------------------------------------------                                    

    let mut fibre_facets    =   fibre_facets_from_complex_facets(
                                    & complex_facets,
                                    & barcode,
                                    ring,
                                );                                    

    //  GET THE ORIGINAL SEQUENCE OF SIMPLICES
    let max_dim             =   complex_facets.iter().map(|x| x.len()-1 ).max().unwrap();
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, max_dim);     

    //  DEFINE A PERMUTAITON ON SIMPLICES INDUCED BY A PERMUTATION ON VERTICES
    for perm_v_o2n in vertex_permutations.iter() {

        println!("{:?}", &perm_v_o2n );
        let complex_facets_permuted
            =   simplices_with_relabeled_vertices(
                    & complex_facets,
                    & perm_v_o2n
                );   
                
        println!("complex facets permuted: {:?}", complex_facets_permuted);
        println!("hard comp: starting");
        let fibre_facets_under_deep_perm
               =    fibre_facets_from_complex_facets(
                        & complex_facets_permuted,
                        & barcode,
                        ring,
                    );  
        println!("hard comp: finished");                    

        let perm_s_o2n      =   simplex_perm_o2n_from_vertex_perm_o2n( &simplex_sequence, &perm_v_o2n );
        let perm_s_n2o      =   inverse_perm( & perm_s_o2n );       
        
        let mut fibre_facets_under_deep_perm_internally_reindexed
               =    Vec::from_iter(
                        fibre_facets_under_deep_perm
                            .iter()
                            .map(   |x| 
                                    Polytope{
                                        data_l_to_fmin: x.data_l_to_fmin.clone(),
                                        data_c_to_l:    compose_f_after_g(
                                                            & x.data_c_to_l, 
                                                            & perm_s_n2o
                                                        )
                                    }
                            )
                    );        

        fibre_facets_under_deep_perm_internally_reindexed.sort();
        fibre_facets.sort();
        assert_eq!( &fibre_facets, &fibre_facets_under_deep_perm_internally_reindexed );
    }
}  

pub fn  deep_superficial_composite_permutation_test (
            complex_facets:         & Vec< Vec< usize > >,
            barcode:                & Barcode< usize >,
            vertex_permutations:    & Vec< Vec< usize >>
        )
{
        // define coefficient field
        let ring                =   RingOpRational::new();

        // run test
        superficial_permutation_test(
            & complex_facets,
            & barcode,
            & ring,
            & vertex_permutations,
        );

        deep_permutation_test(
            & complex_facets,
            & barcode,
            & ring,
            & vertex_permutations,
        );      
}  

pub fn  deep_permutation_test_comprehensive (
    complex_facets:         & Vec< Vec< usize > >,
    barcode:                & Barcode< usize >,
    )
{
    // define coefficient field
    let ring                    =   RingOpRational::new();

    let max_vertex              =   complex_facets.iter().flatten().max().unwrap();
    let vertex_permutations     =   Vec::from_iter(
                                        (0 .. max_vertex + 1).permutations( max_vertex + 1)
                                    );    

    deep_permutation_test(
        & complex_facets,
        & barcode,
        & ring,
        & vertex_permutations,
    );      
} 


pub fn  fibre_facets_from_complex_facets< FilRaw, RingOp, RingElt > (
            complex_facets:     & Vec< Vec< usize > >,
            barcode:            & Barcode< FilRaw >,
            ring:               & RingOp
        )
        ->
        Vec< Polytope >
        where   RingOp:     Ring<RingElt> + Semiring<RingElt> + DivisionRing<RingElt> + Clone,
                RingElt:    Clone + Debug + Ord,
                FilRaw:     Clone + Debug + Ord + Hash,
{

    //  DEFINE THE CELL DIMENSIONS + BOUNDARY MATRIX
    //  --------------------------------------------

    let max_dim             =   complex_facets.iter().map(|x| x.len()-1 ).max().unwrap();
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, max_dim);        
    let cell_dims: Vec<_>   =   simplex_sequence.iter().map(|x| x.len()-1 ).collect();

    let bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence.clone() );
    let boundary            =   boundary_matrix_from_complex_facets(&bimap_sequential, ring.clone());     

    
    //  DEFINE THE INVERSE BARCODE
    //  ------------------------------------ 
    
    let barcode_inverse     =   BarcodeInverse::from_barcode( & barcode );  

    //  WE DON'T HAVE TO SPECIFY ANY PREREQUISITES TO ADDING ANY CELL
    let cell_id_to_prereqs  =   None;

    //  ADDING A CERTAIN SET OF CELLS DOESN'T "FORCE" THE ADDITION OF ANY OTHER CELLS
    let precondition_to_make_new_lev_set = ConditionNone{};

    //  DEFINE THE ROOT NODE + RESULTS VECTOR
    //  -------------------------------------

    let root =  Node::make_root(
                        boundary.clone(),
                    &   barcode,
                    &   barcode_inverse,
                        cell_id_to_prereqs,                    
                    &   cell_dims,  
                        ring.clone(),
                    &   precondition_to_make_new_lev_set,
                        // last_must_be_crit,
                );

    let mut results         =   Vec::new();                                

    //  GATHER RESULTS
    //  --------------

    explore( & root, &mut results );    

    //  RETURN RESULTS
    //  -------------------

    results

}

// TEMPORARILY COMMENTED TO SAVE TIME ON COMPREHENSIVE TESTS

// #[cfg(test)]
// mod tests {
//     // Note this useful idiom: importing names from outer (for mod tests) scope.
//     use super::*;


//     #[test]
//     fn test_sim2_skel1_bc1inf() {

//         let complex_facets      =   vec![  vec![0,1], vec![1,2], vec![0,2] ];
//         let barcode             =   Barcode{
//                                         inf: vec![ BarInfinite{dim:0, birth:0}, BarInfinite{dim:1, birth:1} ],
//                                         fin: vec![],
//                                         ordinal: BiMapSequential::from_vec( vec![0, 1] ),
//                                     };
//         let perms               =   Vec::from_iter(
//                                         (0..3).permutations(3)
//                                     );
//         deep_superficial_composite_permutation_test( &complex_facets, &barcode, &perms);       
//     }  
    
    
//     #[test]
//     fn test_sqr_tri2_bc1inf() {

//         let complex_facets      =   vec![  vec![0,1,2], vec![1, 2, 3] ];
//         let barcode             =   Barcode{
//                                         inf: vec![ BarInfinite{dim:0, birth:0} ],
//                                         fin: vec![],
//                                         ordinal: BiMapSequential::from_vec( vec![0] ),
//                                     };
//         let perms               =   vec![
//                                         vec![ 0, 2, 1, 3 ],
//                                         vec![ 3, 1, 2, 0 ],
//                                         vec![ 3, 2, 1, 0 ]
//                                     ];
//         deep_superficial_composite_permutation_test( &complex_facets, &barcode, &perms);              
//     }    

//     #[test]
//     fn test_itv_5v_bc1inf() {

//         let complex_facets      =   vec![  vec![0,1], vec![1, 2], vec![2, 3], vec![3, 4] ];
//         let barcode             =   Barcode{
//                                         inf: vec![ BarInfinite{dim:0, birth:0} ],
//                                         fin: vec![],
//                                         ordinal: BiMapSequential::from_vec( vec![0] ),
//                                     };
//         let perms               =   vec![
//                                         vec![ 4, 3, 2, 1, 0 ],
//                                     ];
//         deep_superficial_composite_permutation_test( &complex_facets, &barcode, &perms);              
//     }  
    
//     #[test]
//     fn test_itv_5v_bc1inf_1fin() {

//         let complex_facets      =   vec![  vec![0,1], vec![1, 2], vec![2, 3], vec![3, 4] ];
//         let barcode             =   Barcode{
//                                         inf: vec![ BarInfinite{dim:0, birth:0} ],
//                                         fin: vec![ BarFinite{dim:0, birth:1, death: 2} ],
//                                         ordinal: BiMapSequential::from_vec( vec![0, 1, 2] ),
//                                     };
//         let perms               =   vec![
//                                         vec![ 4, 3, 2, 1, 0 ],
//                                     ];
//         deep_superficial_composite_permutation_test( &complex_facets, &barcode, &perms);              
//     } 
    
//     #[test]
//     fn test_itv_4v_bc1inf_1fin() {

//         let complex_facets      =   vec![  vec![0,1], vec![1, 2], vec![2, 3] ];
//         let barcode             =   Barcode{
//                                         inf: vec![ BarInfinite{dim:0, birth:0} ],
//                                         fin: vec![ BarFinite{dim:0, birth:1, death: 2} ],
//                                         ordinal: BiMapSequential::from_vec( vec![0, 1, 2] ),
//                                     };
//         let perms               =   vec![
//                                         vec![ 3, 2, 1, 0 ],
//                                     ];
//         deep_superficial_composite_permutation_test( &complex_facets, &barcode, &perms);      
//         deep_permutation_test_comprehensive( &complex_facets, &barcode );                  
//     }  
    
//     #[test]
//     fn test_symbol_bc2inf() {

//         let complex_facets      =   vec![  vec![0,1], vec![1, 2], vec![0, 2], vec![2, 3] ];
//         let barcode             =   Barcode{
//                                         inf: vec![ BarInfinite{dim:0, birth:0}, BarInfinite{dim:1, birth:1}, ],
//                                         fin: vec![ ],
//                                         ordinal: BiMapSequential::from_vec( vec![0, 1] ),
//                                     };

//         deep_permutation_test_comprehensive( &complex_facets, &barcode );              
//     }   

//     // TEST SEEMS TO PASS BUT IT HAS TO RUN CALCULATIONS FOR 60 DIFFERENT PERMUTATIONS SO I SHORTCIRCUITED IT
//     #[test]
//     fn test_alien_symbol_bc2inf() {

//         let complex_facets      =   vec![  vec![0,1], vec![1, 2], vec![0, 2], vec![2, 3], vec![2, 4] ];
//         let barcode             =   Barcode{
//                                         inf: vec![ BarInfinite{dim:0, birth:0}, BarInfinite{dim:1, birth:1}, ],
//                                         fin: vec![ ],
//                                         ordinal: BiMapSequential::from_vec( vec![0, 1] ),
//                                     };
//         deep_permutation_test(  &complex_facets, 
//                                 &barcode, 
//                                 &RingOpRational::new(), 
//                                 &Vec::from_iter(  (0..5).permutations(5).take(3)   )
//                             )
//     }       
    
    
// }
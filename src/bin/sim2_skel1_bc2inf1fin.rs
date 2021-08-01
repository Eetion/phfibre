
use solar::utilities::index::{BiMapSequential, compose_f_after_g, inverse_perm};
use phfibre::phfibre::{Node, explore, verify_that_barcode_is_compatible};
use phfibre::intervals_and_ordinals::{Barcode, BarcodeInverse, to_ordered_float};
use solar::cell_complexes::simplices_unweighted::maximal_cliques::{    
    ordered_subsimplices_up_thru_dim_concatenated_vec, 
}; 
use solar::cell_complexes::simplices_unweighted::boundary_matrices::{    
    boundary_matrix_from_complex_facets     }; 
use solar::cell_complexes::simplices_unweighted::simplex::{    
    simplex_perm_o2n_from_vertex_perm_o2n   }; 
use num::rational::Ratio;
use ordered_float::OrderedFloat;
use std::collections::HashSet;
use std::iter::FromIterator;

type RingEltRational = OrderedFloat<f64>;
type RingOpRational = solar::rings::ring_native::NativeDivisionRing<RingEltRational>;



fn main() { 

    //  DEFINE COEFFICIENT FIELD
    //  --------------------------------------------
    let ring                =   RingOpRational::new();

    //  DEFINE THE CELL DIMENSIONS + BOUNDARY MATRIX
    //  --------------------------------------------

    let complex_facets      =   vec![  vec![0,1,2] ];

    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 1);    
    let cell_dims: Vec<_>   =   simplex_sequence.iter().map(|x| x.len()-1 ).collect();

    let bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence.clone() );
    let boundary            =   boundary_matrix_from_complex_facets(bimap_sequential, ring.clone()); 


    //  DEFINE THE BARCODE + INVERSE BARCODE
    //  ------------------------------------

    let barcode_inf_dim     =   vec![0, 1];    
    let barcode_inf_brn     =   to_ordered_float( & vec![0., 3.] );

    let barcode_fin_dim     =   vec![ 0 ];
    let barcode_fin_brn     =   to_ordered_float( &vec![ 1. ]   );
    let barcode_fin_die     =   to_ordered_float( &vec![ 2. ]   );


    let barcode             =   Barcode::new(
                                    barcode_inf_dim,
                                    barcode_inf_brn,
                                    barcode_fin_dim,
                                    barcode_fin_brn,
                                    barcode_fin_die
                                );   
    
    let barcode_inverse     =   BarcodeInverse::from_barcode( & barcode );
            

    //  ALLOW FOR DEGENERATE CELLS AFTER THE LAST FINITE BARCODE ENDPOINT
    //  -----------------------------------------------------------------

    let last_must_be_crit   =   false;


    //  DEFINE THE ROOT NODE + RESULTS VECTOR
    //  -------------------------------------

    let root =  Node::make_root(
                        boundary.clone(),
                    &   barcode,
                    &   barcode_inverse,
                    &   cell_dims,  
                        RingOpRational::new(),
                        last_must_be_crit,
                );

    let mut results         =   Vec::new();                                

    //  GATHER RESULTS
    //  --------------


    explore( & root, &mut results );

    println!("\nRESULTS\n");

    // let mut result_vector_set  = HashSet::new();
    // let mut result_vector_vec  = Vec::new();    
    
    // for (result_count, result) in results.iter().cloned().enumerate() {
    //     let mut a   =   result.data_c_to_l.clone();
    //     let     b   =   result.data_l_to_fmin.clone();
    //     let mut c   =   solar::utilities::index::compose_f_after_g(&b, &a);
    //     a.append( &mut c );

    //     verify_that_barcode_is_compatible( 
    //         & root,
    //         & result
    //     );        
        
    //     println!("result number {:?}: {:?}", &result_count, &a );
    //     result_vector_set.insert( a.clone() );
    //     result_vector_vec.push( a.clone() );        
    // }     
    
    // let mut num_by_dim = vec![0, 0, 0, 0];
    // for (result_count, result) in results.iter().cloned().enumerate() {
    //     num_by_dim[ result.dim().unwrap() ] +=1;
    // }
    // println!("number of top dimensional cells: {:?}", num_by_dim );



    //-------------------------------------------------------------

    //  COLLECT AND SORT ALL VERTICES OF THE FIBRE
    let mut vertices = Vec::new();
    for result in results.iter().cloned() {
        if result.dim() == Some(0) {
            vertices.push( compose_f_after_g( & result.data_l_to_fmin, & result.data_c_to_l )   )
        }
    }    
    vertices.sort();

    for vertex in vertices.iter() { println!("vertex: {:?}", &vertex ) }

    
    //  DEFINE A PERMUTAITON ON SIMPLICES INDUCED BY A PERMUTATION ON VERTICES
    let perm_v_o2n      =   vec![2, 0, 1];
    let perm_s_o2n      =   simplex_perm_o2n_from_vertex_perm_o2n( &simplex_sequence, &perm_v_o2n );
    let perm_s_n2o      =   inverse_perm( & perm_s_o2n );

    //  OBTAIN A SEQUENCE OF FIBRE VERTICES CORRESPONDING TO THE PERMUTATION
    //  THIS WORKS B/C THE PERMUTATION DETERMINES AN ISOMORPHISM ON THE UNDERLYING SIMPLICIAL COMPLEX

    let mut vertices_under_group_action     =   vertices.clone();

    for i in 0 .. vertices_under_group_action.len() {
        vertices_under_group_action[ i ] = compose_f_after_g( &vertices_under_group_action[ i ], & perm_s_n2o);
    }

    //  COMPARE VERTEX SETS

    let hset_old    =   HashSet::<std::vec::Vec<usize>>::from_iter( vertices.iter().cloned() );
    let hset_new    =   HashSet::<std::vec::Vec<usize>>::from_iter( vertices_under_group_action );



    println!("OLD - NEW: {:?}", & hset_old.difference( & hset_new ) );
    println!("NEW - OLD: {:?}", & hset_new.difference( & hset_old ) );    

    // -----------------------------------------------------------------------------------------------------------
    //
    // NEW - OLD: [[2, 0, 1, 2, 2, 3], [2, 1, 0, 2, 2, 3]]

}  
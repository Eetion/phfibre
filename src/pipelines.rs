
use crate::phfibre::{Node, explore, verify_that_barcode_is_compatible};
use crate::intervals_and_ordinals::{Barcode, BarcodeInverse, BarFinite, BarInfinite};
use crate::polytope::object_def::Polytope;
use crate::polytope::faces::{poly_complex_facets_to_whole_complex_bimapsequential};
use crate::polytope::intersection::{polytope_intersection};
use crate::polytope::differential::polyhedral_boundary_matrix_binary_coeff;
use crate::polytope::nerve::dowker_nerve;
use crate::rank_calculations::chain_cx_rank_nullity;
use solar::utilities::sequences_and_ordinals::{BiMapSequential, ordinate_unique_vals};
use solar::utilities::statistics::histogram;
use solar::cell_complexes::simplices_unweighted::facets::{    
    ordered_subsimplices_up_thru_dim_concatenated_vec 
}; 
use solar::cell_complexes::simplices_unweighted::boundary_matrices::{    
    boundary_matrix_from_complex_facets 
};   
use solar::rings::field_prime::GF2;
use std::iter::FromIterator;
use ordered_float::OrderedFloat;

// use solar::utilities::sequences_and_ordinals::BiMapSequential;
// use crate::phfibre::{Node, explore};
// use crate::intervals_and_ordinals::{Barcode, BarcodeInverse};
// use crate::polytope::object_def::Polytope;
// use solar::cell_complexes::simplices_unweighted::facets::{    
//     ordered_subsimplices_up_thru_dim_concatenated_vec 
// }; 
// use solar::cell_complexes::simplices_unweighted::boundary_matrices::{    
//     boundary_matrix_from_complex_facets 
// };   
use solar::rings::ring::{Semiring, Ring, DivisionRing};
// use ordered_float::OrderedFloat;
use std::fmt::Debug;
use std::hash::Hash;

// type RingEltRational = OrderedFloat<f64>;
// type RingOpRational = solar::rings::ring_native::NativeDivisionRing<RingEltRational>;






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
    let boundary            =   boundary_matrix_from_complex_facets(bimap_sequential, ring.clone());     

    
    //  DEFINE THE INVERSE BARCODE
    //  ------------------------------------ 
    
    let barcode_inverse     =   BarcodeInverse::from_barcode( & barcode );  

    //  DEFINE THE ROOT NODE + RESULTS VECTOR
    //  -------------------------------------

    let root =  Node::make_root(
                        boundary.clone(),
                    &   barcode,
                    &   barcode_inverse,
                    &   cell_dims,  
                        ring.clone(),
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


pub fn  simplex_pipeline< FilRaw, RingOp, RingElt >(
            simplex_sequence:       &   Vec< Vec< usize > >,
            barcode:                &   Barcode< FilRaw >,
            ring:                   &   RingOp,            
        )
    where   RingOp:     Ring<RingElt> + Semiring<RingElt> + DivisionRing<RingElt> + Clone,
            RingElt:    Clone + Debug + Ord,
            FilRaw:     Clone + Debug + Ord + Hash,  
{
    //  DEFINE THE CELL DIMENSIONS + BOUNDARY MATRIX
    //  --------------------------------------------

    let cell_dims: Vec<_>   =   simplex_sequence.iter().map(|x| x.len()-1 ).collect();

    let bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence.clone() );
    let boundary            =   boundary_matrix_from_complex_facets(bimap_sequential, ring.clone()); 


    //  DEFINE THE BARCODE + INVERSE BARCODE
    //  ------------------------------------
    
    let barcode_inverse     =   BarcodeInverse::from_barcode( & barcode );
            

    //  ALLOW FOR DEGENERATE CELLS AFTER THE LAST FINITE BARCODE ENDPOINT
    //  -----------------------------------------------------------------

    // let last_must_be_crit   =   false;


    //  DEFINE THE ROOT NODE + RESULTS VECTOR
    //  -------------------------------------

    let root =  Node::make_root(
                        boundary.clone(),
                    &   barcode,
                    &   barcode_inverse,
                    &   cell_dims,  
                        ring.clone(),
                );

    let mut poly_complex_facets     =   Vec::new();                                

    //  GATHER RESULTS
    //  --------------

    // println!("{:?}", &root );

    explore( & root, &mut poly_complex_facets );

    //  VERIFY COMPATIBILITY
    //  --------------------

    for poly in poly_complex_facets.iter() {
        verify_that_barcode_is_compatible(
            &   root,
            &   poly
        );
    }
    

    //  REPORT
    //  ------

    println!("POLYHEDRAL COMPLEX FACETS");
    println!("number of facets (total): {:?}", poly_complex_facets.len() ); 
    println!("number of facets (binned by dimension): {:?}", histogram( poly_complex_facets.iter().map(|x| x.dim_cellagnostic().unwrap() ) ) );  

    let poly_complex_bimapsequential    =   poly_complex_facets_to_whole_complex_bimapsequential( & poly_complex_facets );
    let poly_complex_dims               =   Vec::from_iter( poly_complex_bimapsequential.ord_to_val.iter().map(|x| x.dim_cellagnostic().unwrap() ) );    
    let poly_complex_dim_top            =   poly_complex_dims.iter().max().unwrap();    
    println!("POLYHEDRAL COMPLEX CELLS");    
    println!("number of polytopes (total): {:?}", poly_complex_bimapsequential.ord_to_val.len() ); 
    println!("number of polytopes (binned by dimension): {:?}", histogram( poly_complex_dims.iter().cloned() ) );      

    let num_facets              =   poly_complex_facets.len();
    let mut dismat              =   Vec::from_iter(
                                        std::iter::repeat(
                                            Vec::from_iter(
                                                std::iter::repeat( OrderedFloat(1.) )
                                                .take( num_facets)                                                
                                            )
                                        )
                                        .take(num_facets)
                                    );
    for facet_id in 0 .. num_facets{ dismat[facet_id][facet_id] = OrderedFloat(0.) }

    let mut intersection_dim_bins  =   Vec::from_iter( std::iter::repeat(0).take( poly_complex_dim_top + 1) );

    for count_a in 0 .. num_facets{
        for count_b in count_a + 1 .. num_facets {
            if let Some( intersection_poly ) = polytope_intersection( &poly_complex_facets[count_a], &poly_complex_facets[count_b] ){
                dismat[ count_a ][ count_b ] = OrderedFloat( 0.1 );
                dismat[ count_b ][ count_a ] = OrderedFloat( 0.1 );  
                let dim                 =   intersection_poly.dim_cellagnostic().unwrap();                
                intersection_dim_bins[ dim ] += 1              
            }
        }
    }

    let poly_complex_differential           =   polyhedral_boundary_matrix_binary_coeff( & poly_complex_bimapsequential );
   
    let poly_complex_rank_nullity           =   chain_cx_rank_nullity(
                                                    & poly_complex_differential,
                                                    & poly_complex_dims,
                                                    GF2{}
                                                );                                                  
    let poly_complex_betti_vec              =   poly_complex_rank_nullity.rank_homology_vec();         
    println!("betti numbers (of polyhedral complex): {:?}", &poly_complex_betti_vec);


    println!("INTERSECTIONS");        
    println!("number of pairs of intersecting facets, binned by the dimension of the intersection polytope: {:?}", &intersection_dim_bins );    

    
    // dowker complex

    let dowker_complex_facet_bimapseq       =       dowker_nerve(
                                                        & poly_complex_facets,
                                                        & poly_complex_bimapsequential.val_to_ord
                                                    );
    
    let dowker_complex_dim_top              =   dowker_complex_facet_bimapseq.ord_to_val.iter().map(|x| x.len()-1 ).max().unwrap();

    let dowker_complex_simplex_sequence     =   ordered_subsimplices_up_thru_dim_concatenated_vec( 
                                                    & dowker_complex_facet_bimapseq.ord_to_val, 
                                                      dowker_complex_dim_top
                                                );    
    let dowker_complex_cell_dims: Vec<_>    =   dowker_complex_simplex_sequence.iter().map(|x| x.len()-1 ).collect();

    let dowker_complex_bimap_sequential     =   BiMapSequential::from_vec( dowker_complex_simplex_sequence.clone() );
    let dowker_complex_boundary             =   boundary_matrix_from_complex_facets(dowker_complex_bimap_sequential, ring.clone()); 

    let dowker_complex_rank_nullity         =   chain_cx_rank_nullity(
                                                    & dowker_complex_boundary,
                                                    & dowker_complex_cell_dims,
                                                      ring.clone()
                                                );

    println!("DOWKER NERVE COMPLEX FACETS");                                                        
    println!("number of dowker nerve complex facets (total): {:?}", dowker_complex_facet_bimapseq.ord_to_val.len() );
    println!("number of dowker nerve complex facets (binned by dimension): {:?}", histogram( dowker_complex_facet_bimapseq.ord_to_val.iter().map(|x| x.len()-1 ) ) );    

    println!("DOWKER NERVE COMPLEX CELLS");        
    println!("number of nerve dowker complex cells (total): {:?}", dowker_complex_cell_dims.len() );
    println!("number of nerve dowker complex cells (binned by dimension): {:?}", histogram( dowker_complex_cell_dims.iter().cloned() ));
    println!("betti numbers (of dowker nerve): {:?}", &dowker_complex_rank_nullity.rank_homology_vec() );    
}                    
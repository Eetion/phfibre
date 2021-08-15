
use solar::utilities::sequences_and_ordinals::BiMapSequential;
use crate::phfibre::{Node, explore};
use crate::intervals_and_ordinals::{Barcode, BarcodeInverse};
use crate::polytope::object_def::Polytope;
use solar::cell_complexes::simplices_unweighted::facets::{    
    ordered_subsimplices_up_thru_dim_concatenated_vec 
}; 
use solar::cell_complexes::simplices_unweighted::boundary_matrices::{    
    boundary_matrix_from_complex_facets 
};   
use solar::rings::ring::{Semiring, Ring, DivisionRing};
use ordered_float::OrderedFloat;
use std::fmt::Debug;
use std::hash::Hash;

type RingEltRational = OrderedFloat<f64>;
type RingOpRational = solar::rings::ring_native::NativeDivisionRing<RingEltRational>;






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

use solar::utilities::index::{BiMapSequential, histogram};
use phfibre::phfibre::{Node, explore, verify_that_barcode_is_compatible};
use phfibre::intervals_and_ordinals::{Barcode, BarcodeInverse, to_ordered_float};
use solar::cell_complexes::simplices_unweighted::maximal_cliques::{    
    ordered_subsimplices_up_thru_dim_concatenated_vec 
}; 
use solar::cell_complexes::simplices_unweighted::boundary_matrices::{    
    boundary_matrix_from_complex_facets 
};   
use num::rational::Ratio;
use ordered_float::OrderedFloat;


type RingEltRational = OrderedFloat<f64>;
type RingOpRational = solar::rings::ring_native::NativeDivisionRing<RingEltRational>;






pub fn  fibre_facets_from_complex_facets< RingOp, RingElt> (
            complex_facets:     & Vec< Vec< usize > >,
            barcode:            & Barcode,
            ring:               & RingOp
        )
{

    //  DEFINE THE CELL DIMENSIONS + BOUNDARY MATRIX
    //  --------------------------------------------

    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 1);    
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
                        RingOpRational::new(),
                        last_must_be_crit,
                );

    let mut results         =   Vec::new();                                

    //  GATHER RESULTS
    //  --------------

    explore( & root, &mut results );    

    //  RETURN FIBRE FACETS
    //  -------------------

    let fibre_facets    =   
        if let Some( d )    =   results.iter().map();

    Vec::from_iter(
        results
            .iter()
            .filter(|x| x.dim() == 0)
    )

}
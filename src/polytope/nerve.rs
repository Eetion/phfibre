use crate::polytope::object_def::Polytope;
use crate::polytope::faces::poly_faces;
use solar::cell_complexes::simplices_unweighted::facets::ordered_subsimplices_up_thru_dim_concatenated_vec;
use solar::utilities::sequences_and_ordinals::{BiMapSequential};


use std::iter::FromIterator;
use std::collections::HashMap;



/// 
pub fn  dowker_nerve_complex_facets(
            poly_complex_facets:            & Vec< Polytope >,
            poly_complex_vertex_order_map:  & HashMap< Polytope, usize >
        )
        -> 
        Vec< Vec< usize > >
{
    let mut dowker_nerve_facets         =   Vec::new();
    for poly_complex_facet in poly_complex_facets {
        let mut verts                   =   Vec::from_iter(
                                                poly_faces(
                                                    & poly_complex_facet,
                                                    0
                                                )
                                                .iter()
                                                .map(
                                                    |x|
                                                    poly_complex_vertex_order_map.get( x ).unwrap().clone()
                                                )
                                            );
        verts.sort();
        dowker_nerve_facets.push( verts );
        
    }

    // let cell_dim_max        =   dowker_nerve_facets.iter().map(|x| x.len() - 1 ).max().unwrap(); // this may DIFFER from the dimension of the corresponding polytope
    // let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &dowker_nerve_facets, cell_dim_max );    

    // BiMapSequential::from_vec( simplex_sequence )
    dowker_nerve_facets
    
}
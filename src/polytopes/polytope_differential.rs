
use crate::polytopes::polytope::Polytope;
use crate::polytopes::polytope_faces::{poly_faces_by_codim, polys_faces, poly_complex_facets_to_whole_complex_ordinal_data};
use crate::intervals_and_ordinals::{OrdinalData};
use solar::rings::ring::{Semiring, Ring, DivisionRing};
use std::cmp::{PartialOrd};
use std::fmt::Debug;
use std::iter::FromIterator;


type CellID     =   usize;





// !!!!!!!!!!!!!!!!!  
// TO DO
// !!!!!!!!!!!!!!!!!
// - code binary field (operats on bools)
// - calculate betti numbers homology 2 different ways
// - alternate main algorithm
// - subroutine to shortcircuit construction if you've already reached a similar subdivided state
// - instructions for jacob



//-------------------------------------------------------------------------
//  SIGNATURE FOR BOUNDARY MATRIX WITH ARBITRARY COEFFICIENTS;
// 
//  FILL THIS IN IF / WHEN SIGNS ARE WORKED OUT


// pub fn  polytope_facets_to_boundary_matrix< RingOp, RingElt > (
//             complex_facets:         Vec< Polytope >,
//             ring:                   RingOp,            
//         )
//         -> 
//         Vec< Vec< ( CellID, RingElt ) > >
        
// where   RingOp:     Clone + Semiring<RingElt> + Ring<RingElt> + DivisionRing<RingElt>,
//         RingElt:    Clone + Debug + PartialOrd  

//-------------------------------------------------------------------------



/// Given an OrdinalData object representing a bijection between {all polyhedra
/// in a complex} and {0, .., n}, generate the corresponding boundary matrix with
/// indices in {0, .., n} and coefficients in the 2-element field GF2.
/// 
/// The matrix has integer indices. 
pub fn  polyhedral_boundary_matrix_binary_coeff(
            polyhedra_ordinal_data: & OrdinalData< Polytope >,       
        )
        -> 
        Vec< Vec< ( CellID, bool ) > >
           
{

    // build boundary
    
    let mut boundary_matrix        =   Vec::new();
    for face in polyhedra_ordinal_data.ord_to_val.iter() {

        let mut col             =   Vec::from_iter(
                                        poly_faces_by_codim(
                                            face,
                                            1,
                                        )
                                        .iter()
                                        .map( 
                                            |x|
                                            ( polyhedra_ordinal_data.ord(x).unwrap(), true)
                                        )
                                    );
        col.sort();
        boundary_matrix.push( col );
    }
    
    boundary_matrix
} 

use solar::utilities::index::{BiMapSequential, compose_f_after_g, sort_perm};
use solar::utilities::ring::{MinusOneToPower};
use solar::rings::ring::{Semiring, Ring};
use std::cmp::Ordering;
use std::collections::{HashSet, HashMap};
use std::hash::Hash;
use itertools::Itertools;
use std::iter::FromIterator;
use std::fmt::Debug;

// type Vertex = u16;












//  ---------------------------------------------------------------------------
//  CONSTRUCT BOUNDARY MATRIX FROM FACE SEQUENCE
//  ---------------------------------------------------------------------------







#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    
    #[test]
    fn test_sequentialize_facets () {

        let ring                    =   solar::rings::ring_native::NativeDivisionRing::< f64 >::new();
        let complex_facets          =   vec![ vec![0,1,2] ];
        let bimap_sequential        =   BiMapSequential::from_vec(
                                            ordered_sequence_of_faces( complex_facets )
                                        );  
        let boundary                =   boundary_matrix_from_complex_facets( bimap_sequential, ring );

        for vec in boundary { println!("{:?}", vec) }
    }

    #[test]
    fn test_subsequences_up_to_card() {

        let facets                  =   vec![ vec![0, 1, 2] ];
        let subsequences            =   subsequences_up_to_card( facets, 3);
        println!("{:?}", & subsequences);

        assert_eq!(     &   subsequences,
                        &   vec![
                                vec![                                           ],
                                vec![  vec![0],     vec![1],    vec![2]         ],                                
                                vec![  vec![0,1],   vec![0,2],  vec![1,2]       ],
                                vec![  vec![0,1,2]                              ]
                            ]
        )

    }
    

}    
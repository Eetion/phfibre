
use solar::utilities::index::BiMapSequential;
use solar::utilities::ring::{MinusOneToPower};
use solar::rings::ring::{Semiring, Ring};
use std::cmp::Ordering;
use std::collections::{HashSet, HashMap};
use std::hash::Hash;
use itertools::Itertools;
use std::iter::FromIterator;
use std::fmt::Debug;



//  ---------------------------------------------------------------------------
//  FACES FROM FACETS-OF-THE-COMPLEX (NEW)
//  ---------------------------------------------------------------------------

pub fn  ordered_subsimplices_up_to_dim< Vertex >( 
    complex_facets: & Vec< Vec< Vertex >>, 
    max_dim: usize 
) 
-> 
Vec< Vec< Vec< Vertex >>> 
where Vertex: Ord + Clone
{
    let mut seq             =   Vec::with_capacity( max_dim );
    for card in 1 .. max_dim + 2 {
        let vec: Vec<_>     =   complex_facets
                                .iter()
                                .map( |x| x.iter().cloned().combinations(card)  )
                                .kmerge()
                                .dedup()
                                .collect();
        seq.push( vec );
    }
    seq
}


pub fn  ordered_subsimplices_up_to_dim_concatenated< Vertex >( 
    complex_facets: & Vec< Vec< Vertex >>, 
    max_dim: usize 
) 
-> 
Vec< Vec< Vertex >>
    where Vertex: Ord + Clone
{
    let mut a = ordered_subsimplices_up_to_dim( complex_facets, max_dim );
    let mut b = Vec::new();
    for i in 0 .. a.len() {
        b.append( &mut a[ i ] );
    }
    b
}



//  ---------------------------------------------------------------------------
//  FACET BIMAP TO BOUNDARY
//  ---------------------------------------------------------------------------


pub fn  boundary_matrix_from_complex_facets< Vertex, RingOp, RingElt >( 
            simplex_bimap:  BiMapSequential< Vec < Vertex > >,
            ring:           RingOp
        ) 
        ->
        Vec< Vec < (usize, RingElt) >>

        where   Vertex:    Ord + Hash + Clone + Debug,      
                RingOp:     Semiring< RingElt > + Ring< RingElt >,
{
    if simplex_bimap.ord_to_elt.is_empty() { return vec![] }

    let mut boundary            =   Vec::with_capacity( simplex_bimap.ord_to_elt.len() );  
    
    let mut simplex_dim         =   0;
    let mut simplex_num_verts   =   0;

    for simplex in simplex_bimap.ord_to_elt.iter().cloned() {

        simplex_num_verts       =   simplex.len();
        simplex_dim             =   simplex_num_verts - 1;

        // no need to calculate boundaries of dim-0 cells
        if simplex_dim == 0 {
            boundary.push( Vec::with_capacity(0) );
            continue;
        }  

        let mut vec             =   Vec::with_capacity( simplex_num_verts );    // num_vertices = NUMBER OF FACETS

        for (facet_count, facet)  in simplex.iter().cloned().combinations( simplex_dim ).enumerate() {
            vec.push( 
                (
                    simplex_bimap.ord( &facet ),
                    ring.minus_one_to_power( simplex_dim - facet_count )
                ) 
            )            
        }
        boundary.push( vec );
    }

    boundary

}










#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;




    #[test]
    fn test_ordered_subsimplices_up_to_dim() {

        let complex_facets          =   vec![ vec![0, 1, 2] ];

        assert_eq!(         ordered_subsimplices_up_to_dim( & complex_facets, 2),
                            vec![
                                vec![   vec![0],     vec![1],    vec![2]         ],                                
                                vec![   vec![0,1],   vec![0,2],  vec![1,2]       ],
                                vec![   vec![0,1,2]                              ]
                            ]
        );

        assert_eq!(         ordered_subsimplices_up_to_dim_concatenated( & complex_facets, 2),
                            vec![
                                        vec![0],     vec![1],    vec![2],                                         
                                        vec![0,1],   vec![0,2],  vec![1,2],       
                                        vec![0,1,2]                              
                            ]
        ) ;       


    }

    #[test]
    fn test_bimap_to_boundary () {

        let ring                    =   solar::rings::ring_native::NativeDivisionRing::< f64 >::new();
        let complex_facets          =   vec![ vec![0,1,2] ];


        let bimap_sequential        =   BiMapSequential::from_vec(
                                            ordered_subsimplices_up_to_dim_concatenated( & complex_facets, 2 )
                                        );  

        let boundary                =   boundary_matrix_from_complex_facets( bimap_sequential, ring );

        assert_eq!(     &   boundary,
                        &   vec![
                                    vec![],
                                    vec![],
                                    vec![],
                                    vec![(0, -1.0), (1, 1.0)],
                                    vec![(0, -1.0), (2, 1.0)],
                                    vec![(1, -1.0), (2, 1.0)],
                                    vec![(3, 1.0), (4, -1.0), (5, 1.0)]
                            ]
        )
    }    


}    
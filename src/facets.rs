
use solar::utilities::index::{BiMapSequential, compose_f_after_g, sort_perm, inverse_perm};
use solar::utilities::ring::{MinusOneToPower};
use solar::rings::ring::{Semiring, Ring};
use solar::cell_complexes::simplices::unweighted::Simplex;
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
//  PERMUTATION ON SIMPLICES INDUCED BY PERMUTATION ON VERTICES
//  ---------------------------------------------------------------------------


/// Given an vector f = [f0, .., fn] representing a function of form 
/// f: old_vertex_number -> new_vertex_number, obtain the vector 
/// g: old_simplex_number -> new_simplex_number.
pub fn  simplex_perm_o2n_from_vertex_perm_o2n( 
    simplex_sequence:           &   Vec< Vec< usize >>,
    vertex_perm_old_to_new:     &   Vec< usize >
    ) 
    ->
    Vec< usize >
{
    // Create vector of new simplices
    let mut new_simplex_sequence =  Vec::from_iter(
                                        simplex_sequence
                                            .iter()
                                            .cloned()
                                            .map(
                                                |x|
                                                Simplex{ 
                                                    vertices:  compose_f_after_g( &vertex_perm_old_to_new, &x )
                                                }
                                            )
                                    );

    // We must remember to sort the new vertices                                    
    for simplex in new_simplex_sequence.iter_mut() { simplex.vertices.sort()}

    // Obtain the sort permutation
    sort_perm( &new_simplex_sequence )
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


    #[test]
    fn test_simplex_perm_o2n_from_vertex_perm_o2n() {

        // sequence_old:          [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2]]
        // sequence_old_permuted: [[0], [1], [3], [2], [0, 1], [0, 3], [0, 2], [1, 2]]
        // new_sequence:          [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 3]]
        // permutation: simplex old -> new 
        //                        [ 0,   1,   3,   2,   4,      6,      5,      7]


        let complex_facets          =   vec![  vec![0,1,2], vec![0, 3] ];
        let simplex_sequence_old    =   ordered_subsimplices_up_to_dim_concatenated( &complex_facets, 1);   
        let perm_v_o2n              =   vec![0, 1, 3, 2];

        let mut simplex_sequence_new =  Vec::from_iter(
            simplex_sequence_old
                .iter()
                .cloned()
                .map(
                    |x|
                    Simplex{ 
                        vertices:  compose_f_after_g( &perm_v_o2n, &x )
                    }
                )
        );        

        // We must remember to sort the new vertices                                    
        for simplex in simplex_sequence_new.iter_mut() { simplex.vertices.sort() }        

        // perm: simplex OLD -> NEW
        let perm_s_o2n              =   simplex_perm_o2n_from_vertex_perm_o2n( &simplex_sequence_old, &perm_v_o2n );
        // perm: simplex NEW -> OLD
        let perm_s_n2o               =   inverse_perm( &perm_s_o2n );

        let mut simplex_sequence_permuted    =   compose_f_after_g( &simplex_sequence_old, &perm_s_n2o );

        let mut simplex_sequence_permuted_vertex_translated     =   simplex_sequence_permuted.clone();
        for i in 0..simplex_sequence_permuted_vertex_translated.len() { simplex_sequence_permuted_vertex_translated[i] = compose_f_after_g( & perm_v_o2n, & simplex_sequence_permuted[i]) };
        for i in 0..simplex_sequence_permuted_vertex_translated.len() { simplex_sequence_permuted_vertex_translated[i].sort() };        
        

        println!("sequence_old:          {:?}",     & simplex_sequence_old );
        println!("sequence_old_permuted: {:?}",     & simplex_sequence_permuted );        
        println!("new_sequence:          {:?}",     & simplex_sequence_permuted_vertex_translated );     
        println!("permutation: simplex old -> new {:?}", & perm_s_o2n);           

    }


}    


use crate::polytope::object_def::{Polytope};
use solar::utilities::sequences_and_ordinals::{ordinate_unique_vals, BiMapSequential, reverse_hash_sequential};
use solar::utilities::combinatorics::fixed_sum_sequences;
// use crate::utilities::*;
use solar::utilities::indexing_and_bijection::SuperIndex;
// use ordered_float::OrderedFloat;
// use num::rational::Ratio;
// use std::hash::Hash;
use std::iter::{FromIterator, repeat};
// use std::cmp::Ord;
use itertools::Itertools;
// use itertools::object_defs::{MultiProduct, Combinations};

type Fil = usize;





struct SimprodFaceEnumerator{
    poly:   Polytope,

}

/// Given a sequence of simplices s_1, .., s_n with dimensions d_1, .., d_n,
/// return a vecttor consisting of all sequences (X_1, .., X_n) such that
/// (s_1 - X_1) x ... x (s_n - X_n) is a dimension-`face_dim` face of 
/// s_1 x ... x s_n.
pub fn  vertex_deletion_choices_for_product_of_combinatorial_simplices(
            dims:       &   Vec< usize >,
            face_dim:       usize
        )
        ->
        Vec< Vec< Vec< usize >  > >
        // Flatten< Map< MultiProduct < Map< Combinations < std::ops::Range<usize>  > > > > >
{
    let mut deletion_choices = Vec::new();

    // for each sequence summing to the proper codimension
    for sequence in fixed_sum_sequences( &dims, dims.iter().sum::<usize>() - face_dim ) {

        // enumerate all combinations of vertices to delete
        deletion_choices.extend(
            sequence.iter().enumerate().map( |(y_count, y)| (0..dims[y_count] + 1 ).combinations(*y) ).multi_cartesian_product()
        )
    }

    deletion_choices

    // WHAT I WOULD HAVE LIKED TO RETURN -- BUT FIGURING OUT THE RIGHT TYPE ANNOTATION IS HARD
    // Vec::from_iter(
    //     fixed_sum_sequences( &dims, dims.iter().sum() - face_dim )
    //         .iter()
    //         .map(   |x| 
    //                 x.iter().enumerate().map( |(y_count, y)| (0..dims[y_count + 1]).combinations(*y) ).multi_cartesian_product()
    //         )
    //         .flat_map(|it| it.clone())
    // )
}

// Creates a boolean vector with "true" in slot i if we are supposed to merge
// level set (i-1) into i.  Returns a value of `None` if `verts_to_delete` 
// contains the last vertex of the last simplex.
pub fn  verts_to_delete_per_simp_2_merge_down_flag(
            verts_to_delete_per_simp:   & Vec< Vec< usize > >,
            simplex_dims:               & Vec< usize >,
        )
        ->
        Option< Vec< bool > >
{

    if  verts_to_delete_per_simp.last().unwrap().contains(
            simplex_dims.last().unwrap()
        )
     { return None }

    let num_lsord_old               =   simplex_dims.iter().sum::<usize>() + simplex_dims.len(); // each dim differs from set cardinality by 1, so we must add simplex_dims.len() to compensate

    let mut merge_down_flag         =   Vec::from_iter( repeat(false).take( num_lsord_old ) );

    let mut offset                  =   0;
    for (simplex_count, verts_to_delete) in verts_to_delete_per_simp.iter().enumerate() {
        // we specify vertices to delete for each simplex as a subset of {0, .., n}; however,
        // we really want to think of the polytopes as a product of form
        // {0, .., n_0} x {n_0+1,...,n_1} x ...
        // the "offset" defined below are essentially the n_i's
        offset                      =   simplex_dims.iter().take(simplex_count).sum::<usize>() + simplex_count;
        for vert in verts_to_delete.iter().map(|x| x + offset ) {
            merge_down_flag[ vert + 1 ]   =   true;
        }
    }   
    Some( merge_down_flag )
}



/// Given a boolean vector `v`, interpret `v[i] = true` to mean that level sets
/// `i` and `i - 1` should merge; this function returns a vector X such that 
/// (level set number k    in the original partition) is a subset of 
/// (level set number X[k] in the new, merged partition) (this uniquely defines the vector).
pub fn  merge_down_flag_to_vec_mapping_lsord_old2new(
            merge_down_flag:    &   Vec< bool >
        )
        ->
        Vec< usize >
{
    // Calculate the minimal lev set ord ord with which each lev set merges
    let mut lsord_old2new  =   Vec::from_iter( 0 .. merge_down_flag.len() );
    for lev_set_ord in 0 .. merge_down_flag.len() {
        while   merge_down_flag[  
                    lsord_old2new[
                        lev_set_ord
                    ]
                ] 
        {  
            lsord_old2new[
                lev_set_ord
            ]
            -= 1;
        }
    }

    // If we merged any level sets, then the ordinals of many level sets may have changed
    let ordinal_data                =   ordinate_unique_vals( &lsord_old2new );
    for val in lsord_old2new.iter_mut() { *val = ordinal_data.ord( &val ).unwrap() }

    lsord_old2new    
}


/// Determine whether merging level sets downward as indicated by the merge_down_flag
/// would result in a valid face of the polytope.
/// 
/// Returns false iff at least one of the following conditions holds: (i) all vertices
/// from a single simplex are removed, (ii) level set 0 is supposed to merge down
/// (to say that such a merge were legal would permit pathologies where we think we are
/// reduced the dimension of a simplex but aren't, really)
pub fn  merge_down_flag_is_valid_for_simplex_factor_dims(
            merge_down_flag:        & Vec< bool >,
            simplex_dims:           & Vec< usize >
        )
        ->
        bool
{
    // only the empty vetor is appropriate for an empty product
    if simplex_dims.is_empty() { 
        if merge_down_flag.is_empty() { return true }
        else {return false}
    }

    // can't merge the first level set down
    if merge_down_flag[0] { return false }

    // check that flag vector has length equal to the number of level sets
    let num_lev_sets        =   simplex_dims.iter().sum::<usize>() + simplex_dims.len();
    if merge_down_flag.len() != num_lev_sets {return false }
    
    // check all but the last simplex (it is impossible to remove all vertices from the
    // last simplex by accident)
    let mut offset          =   1;
    for simplex_dim
        in 
        simplex_dims[ 0 .. simplex_dims.len() - 1 ].iter() 
    {        
        if  merge_down_flag[ offset .. offset + simplex_dim + 1 ]
                .iter()
                .all(|&x| x )
        { return false }

        offset  =   offset + simplex_dim + 1
    }

    return true

}

/// Given a vector of "merge-down" flags, determine the vector mapping 
/// level set ordinal to fmin after merging is complete.
/// 
/// NOTE: This only works correctly if you provide a merge flag that's compatible
/// with the given simplex dimensions, as determined by the funciton [`merge_down_flag_is_valid_for_simplex_factor_dims`].
pub fn  merge_down_flag_simplex_dims_to_lsord2fmin(
            merge_down_flag:        & Vec< bool >,
            simplex_dims:           & Vec< usize >
        )
        ->
        Vec< usize >
{
    let num_merged              =   merge_down_flag.iter().filter(|&&x| x).count();
    let num_lev_sets_orig       =   simplex_dims.iter().sum::<usize>() + simplex_dims.len();

    let mut lev_set_ord_to_fmin =   Vec::with_capacity( merge_down_flag.len() - num_merged );

    // Suppose the original lev_set_ord_to_fmin and merge_down_flag vectors look like this
    //      0 0 1 1 1 2
    //      0 0 1 0 0 1
    // then the new lev_set_ord_to_fmin vector should look like
    //      0 * 1 1 * 2
    // in particular we remove 
    let mut offset              =   1;
    for ( simplex_count, simplex_dim ) in simplex_dims.iter().enumerate() {
        let num_removed         =   ( offset .. offset + simplex_dim + 1 )
                                        .filter(
                                            |&x|
                                            merge_down_flag.sindex( x, false ) // we use sindex to handle the edge case where x indexes into the last element of merge_down_flag
                                        )
                                        .count();                             

        lev_set_ord_to_fmin
            .extend( 
                repeat( simplex_count )
                    .take( 1 + simplex_dim -  num_removed ) 
                );

        offset += simplex_dim + 1
    }

    lev_set_ord_to_fmin
}






/// Suppose that we merge level sets as stipulated by `vertices_to_delete_per_simp`.
/// This function returns a vector X such that 
/// (level set number k    in the original partition) is a subset of 
/// (level set number X[k] in the new, merged partition) (this uniquely defines the vector).
/// 
/// Returns a value of `None` if `verts_to_delete` contains the last vertex of the last simplex.
pub fn  vec_mapping_lsord_old2new(            
            verts_to_delete_per_simp:   & Vec< Vec< usize > >,
            simplex_dims:               & Vec< usize >,
        )
        -> 
        Option< Vec< usize > >
{

    if let Some( merge_down_flag ) =    verts_to_delete_per_simp_2_merge_down_flag(
                                            verts_to_delete_per_simp,
                                            simplex_dims,
                                        )
    {
        Some( merge_down_flag_to_vec_mapping_lsord_old2new( & merge_down_flag ) )
    } else {
        None
}

// OLD STUFF -------------------------------
// // This code works by first placing each level set ordinal into an equivalence class.
// // The process works in two steps.  In step 1, each equivalence class is functionally 
// // "labeled" by the smallest level set ordinal that it contains.  After tall the classes 
// // have been computed, however, we relabel them.
// let num_lsord_old               =   simplex_dims.iter().sum::<usize>() + simplex_dims.len(); // each dim differs from set cardinality by 1, so we must add simplex_dims.len() to compensate

// let mut merge_down_flag         =   Vec::from_iter( repeat(false).take( num_lsord_old ) );

// // The following code creates a boolean vector with "true" in slot i if we are supposed to merge
// // level set (i-1) into i.  This essentially builds a directed graph.  The equivalence classes of
// // level sets correspond to the connected components of the graph.
// let mut offset                  =   0;
// for (simplex_count, verts_to_delete) in verts_to_delete_per_simp.iter().enumerate() {
//     // we specify vertices to delete for each simplex as a subset of {0, .., n}; however,
//     // we really want to think of the polytopes as a product of form
//     // {0, .., n_0} x {n_0+1,...,n_1} x ...
//     // the "offset" defined below are essentially the n_i's
//     offset                      =   simplex_dims.iter().take(simplex_count).sum::<usize>() + simplex_count;
//     for vert in verts_to_delete.iter().map(|x| x + offset ) {
//         merge_down_flag[ vert + 1 ]   =   true;
//     }
// } 

// // Calculate the minimal lev set ord ord with which each lev set merges
// let mut lsord_old2new  =   Vec::with_capacity(num_lsord_old ); // level-set-ordinal-old to level-set-ordinal-new    
// for lsord_old in 0 .. num_lsord_old {
//     let mut min_equiv_ord       =   lsord_old.clone();
//     while merge_down_flag[ min_equiv_ord ] {
//         min_equiv_ord -= 1;
//     }
//     lsord_old2new.push( min_equiv_ord );
// }

// // If we merged any level sets, then the ordinals of many level sets may have changed
// let ordinal_data                =   ordinate( &lsord_old2new );
// for val in lsord_old2new.iter_mut() { *val = ordinal_data.ord( &val ).unwrap() }

// lsord_old2new
}


/// Return a face of the underlying polytope, corresponding to deletion of the vertices provided. 
/// 
/// Assumes that user never attempts to delete every vertex from a simplex.
pub fn  poly_face(            
            verts_to_delete_per_simp:   &   Vec< Vec< usize > >,
            poly:                       &   Polytope,
        )
        ->
        Polytope
{
    let mut poly            =   poly.clone();    
    let simplex_dims        =   poly.simplex_factor_dims_cellagnostic();
    let translator          =   vec_mapping_lsord_old2new(            
                                    & verts_to_delete_per_simp,
                                    & simplex_dims,
                                )
                                .unwrap();
    
    // update mapping from cells to level sets                            
    for cell_id in 0 .. poly.num_cells() {
        poly.data_c_to_l[ cell_id ]    =    translator[
                                                poly.data_c_to_l[ cell_id ]
                                            ];
    }

    // update mapping from cells to level sets
    poly.data_l_to_fmin.clear();
    for simplex_count in 0 .. simplex_dims.len() {
        let new_num_lev_sets    =   1 + simplex_dims[ simplex_count] - verts_to_delete_per_simp[ simplex_count ].len();
        poly.data_l_to_fmin.extend(  std::iter::repeat( simplex_count ).take( new_num_lev_sets )  )
    }

    poly
}        

/// Return a face of the underlying polytope, corresponding to deletion of the vertices provided. 
/// 
/// Recall that "vertices" for each simplex correspond to *gaps* between consecutive level sets.
/// 
/// Assumes you do not attempt to delete the last vertex of the last simplex.
/// 
/// We think of polyvv[ k ] as the sequence of level sets with fmin k.
pub fn  polyvvv_face(
            verts_to_delete_per_simp:   &   Vec< Vec< usize > >,
            polyvv:                     &   Vec< Vec< Vec< usize > > >, 
        )   
        ->
        Vec< Vec< Vec< usize > > >
{
    let mut face                    =   polyvv.clone();   
    let num_crit_vals               =   face.len() - 1; // the last bin does not count as a critical value, since topology doesn't change there
    for crit_val_ord in ( 0 .. num_crit_vals ).rev()  
    {

        for gap_ord_to_delete in verts_to_delete_per_simp[ crit_val_ord ].iter().rev() 
        {
            // the level set we will merge into another level set
            let mut lev_set_to_merge    =   face[ crit_val_ord ].remove( gap_ord_to_delete.clone() );   // removing this level sets reduces the length of face[ crit_val_ord ] by 1                     

            // we must now merge our level set with the next level set up

            // if we are deleting the last gap, then we merge into the first level set of the next critical value            
            if *gap_ord_to_delete == face[ crit_val_ord ].len()  {     // note that face[ crit_val_ord ].len() is now shorter by 1 that what it was when we removed lev_set_to_merge
                face[ crit_val_ord + 1 ][ 0 ].append( &mut lev_set_to_merge );
                face[ crit_val_ord + 1 ][ 0 ].sort();
            // otherwise merge it with the next level set with the same critical value                
            } else {                
                face[ crit_val_ord ][ *gap_ord_to_delete ].append( &mut lev_set_to_merge );
                face[ crit_val_ord ][ *gap_ord_to_delete ].sort();                
            }
        }
    }

    face
}     

// TESTING SUPPORT FOR THIS FUNCTION COMES FROM <test_poly_to_face_and_face_enumerate> BELOW
//
//
/// Returns all faces of a polytope of given dimension -- NOT NECESSARILY IN SORTED ORDER.
pub fn  poly_faces( 
            poly: &Polytope, 
            face_dim: usize,
        ) 
        -> 
        Vec< Polytope > 
{
    if face_dim > poly.dim_cellagnostic().unwrap() { return vec![] } // there are no faces of strictly higher dimension

    let simplex_factor_dims         =   poly.simplex_factor_dims_cellagnostic();
    let deletion_choices            =   vertex_deletion_choices_for_product_of_combinatorial_simplices(
                                            & simplex_factor_dims,
                                            face_dim
                                        );    
    Vec::from_iter(
        deletion_choices
            .iter()
            .map(|deletion_choice| poly_face( deletion_choice, poly ))
    )                                       
}

pub fn  poly_faces_by_codim(
            poly: &Polytope, 
            codim_face: usize,
        ) 
        -> 
        Vec< Polytope >     
{
    let dim_poly                    =   poly.dim_cellagnostic().unwrap();
    if codim_face > dim_poly {
        Vec::with_capacity(0)
    } 
    else {
        poly_faces( 
            poly,
            dim_poly - codim_face,
        ) 
    }
}



/// Given a collection of polytopes, enumerate all faces of a given dimension -- IN SORTED ORDER.
pub fn  polys_faces( 
            polys: &Vec< Polytope >, 
            face_dim: usize,
        ) 
        -> 
        Vec< Polytope > 
{
    let mut faces                   =   Vec::new();
    let mut buffer                  =   Vec::new();
    let mut receiver                =   Vec::new();

    for poly in polys.iter() {
        receiver                    =   poly_faces( poly, face_dim );
        receiver.sort();
        buffer.clear();
        buffer.extend(
            faces.iter().cloned().merge( receiver.iter().cloned() ).dedup()
        );
        faces.clear();
        faces.append( &mut buffer );
    }

    faces 
                                     
}


/// Given the facets of a polyhedral complex, (i) enumerate all polyherdra, 
/// (ii) order these polyhedra in ascending order, first by dimension and
/// second by the order derived by Rust on polyhedra, (iii) encode this ordering
/// in an BiMapSequential struct.
pub fn  poly_complex_facets_to_whole_complex_ordinal_data(
            complex_facets: & Vec< Polytope >,       
        )
        ->
        BiMapSequential< Polytope >
{
    let dim_top                 =   complex_facets
                                        .iter()
                                        .map(
                                                |x| 
                                                x.dim_cellagnostic().unwrap() 
                                        )
                                        .max()
                                        .unwrap();
    

    // Collect all faces into a vector                                    
    let mut all_faces           =   Vec::new();
    for poly_dim in 0 .. dim_top + 1 {
        let mut pure_dim_faces  =   polys_faces(
                                        &   complex_facets,
                                            poly_dim,
                                    );
        pure_dim_faces.sort();
        all_faces.append( &mut pure_dim_faces );
    }

    let rev_hash                =   reverse_hash_sequential( & all_faces );

    // build ordinal data
    BiMapSequential{ ord_to_val: all_faces, val_to_ord: rev_hash }
}







#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    use crate::polytope::object_def::{random_polytope, poly_to_polyvvv, polyvvv_to_poly};
    use std::collections::HashSet;

    #[test]
    fn test_vec_mapping_lsord_old2new() {
        
        let simplex_dims                =   vec![1, 3, 0];
        let verts_to_delete_per_simp    =   vec![ vec![1], vec![0], vec![] ];
        
        let translator_test             =   vec_mapping_lsord_old2new(            
                                                & verts_to_delete_per_simp,
                                                & simplex_dims
                                            )
                                            .unwrap();
        
        let translator_true             =   vec![ 0, 1, 1, 1, 2, 3, 4 ];

        assert_eq!( &translator_test, &translator_true );

    }



    #[test]
    fn test_polyvv_to_face() {
        let polyvvv         =   vec![
                                    vec![   vec![ 0, 1 ], vec![ 2, 3], vec![ 4, 5]  ],
                                    vec![   vec![ 8,   ], vec![ 9   ]               ],
                                    vec![   vec![ 10,  ], vec![ 11  ]               ],                                
                                ];
        
        let verts_to_delete_per_simp    
                            =   vec![ vec![ 0 ], vec![ 1 ] ];

        let face_true       =   vec![
                                    vec![   vec![ 0, 1, 2, 3],  vec![ 4, 5]         ],
                                    vec![   vec![ 8,   ],                           ],
                                    vec![   vec![ 9, 10],       vec![ 11  ]         ],                                
                                ];
        let face_test       =   polyvvv_face(
                                    &   verts_to_delete_per_simp,  
                                    &   polyvvv
                                );
                            
        assert_eq!( &face_test, &face_true);
        
    }



    #[test]
    fn test_poly_face() {

        for params 
            in 
            (0..3).map(|_x| 0..3 ).multi_cartesian_product() 
        {            
            let max_num_fmin                =   params[0];
            let max_num_lev_sets_per_fmin   =   params[1];
            let max_num_cells_per_lev_set   =   params[2];
            let attach_dim0_simplex_to_end  =   true;

            let poly                        =   random_polytope(
                                                    max_num_fmin,
                                                    max_num_lev_sets_per_fmin,
                                                    max_num_cells_per_lev_set,
                                                    attach_dim0_simplex_to_end
                                                );                                               
            
            let polyvvv                     =   poly_to_polyvvv( &poly );

            if poly.dim_cellagnostic() == None { continue };
            
            let polytope_dim                =   poly.dim_cellagnostic().unwrap();
            let simplex_factor_dims         =   poly.simplex_factor_dims_cellagnostic();

            assert!( poly.contains( & poly ) ); // sanity check

            for face_dim in 0 .. polytope_dim {
                
                let mut face_pool           =   Vec::new();

                let deletion_choices        =   vertex_deletion_choices_for_product_of_combinatorial_simplices(
                                                    & simplex_factor_dims,
                                                    face_dim
                                                );
                for deletion_choice in deletion_choices.iter() {
                    let face_poly           =   poly_face(
                                                    & deletion_choice,
                                                    & poly,
                                                );
                    let face_polyvvv        =   polyvvv_face(
                                                    & deletion_choice,
                                                    & polyvvv,
                                                );   

                    // check that calculating faces returns equivalent result for both data structures.
                    assert_eq!( &face_poly, &polyvvv_to_poly( &face_polyvvv ));
                    assert_eq!( &poly_to_polyvvv( &face_poly), &face_polyvvv );                    

                    // check that face_poly is indeed a face
                    assert!( poly.contains( & face_poly ) );

                    // check that face_poly has the correct dimension
                    assert_eq!( face_poly.dim_cellagnostic(), Some( face_dim ) );                    

                    // push face_poly to the face pool
                    face_pool.push( face_poly );

                }

                // check that we have not created any duplicate faces
                //
                // note that this only really works if we have assigned at least one 
                // cell to each level set (otherwise the algorithms can return 
                // duplicate faces even when functioning as they are supposed to).
                // therefore we make that a precondition for testing.
                let mut vec_mapping_cell_id_to_lev_set_ord 
                        =   poly.vec_mapping_cell_id_to_level_set_ordinal();
                vec_mapping_cell_id_to_lev_set_ord.sort();
                vec_mapping_cell_id_to_lev_set_ord.dedup();
                if vec_mapping_cell_id_to_lev_set_ord.len() == poly.num_lev_sets() {
                    face_pool.sort();
                    face_pool.dedup();
                    assert_eq!( face_pool.len(), deletion_choices.len() );
                }

            }

        }

    }


    #[test]
    fn test_polys_faces() {

        //  NOTE: THE TOP LEVEL SET ORDINAL IS EMPTY BY DESIGN
        let poly_a      =   Polytope{
                                data_c_to_l:        vec![0,0,0,1,1,1,2,2,3,3],
                                data_l_to_fmin:     vec![0,0,1,1,2]
                            };
        let poly_b      =   Polytope{
                                data_c_to_l:        vec![0,0,1,1,1,2,3,3,4,4],
                                data_l_to_fmin:     vec![0,0,0,1,1,2]
                            };                           
        let poly_c      =   Polytope{
                                data_c_to_l:        vec![0,0,0,1,1,1,2,2,2,2],
                                data_l_to_fmin:     vec![0,1,1,2]
                            };                                
        let poly_d      =   Polytope{
                                data_c_to_l:        vec![0,0,0,0,0,1,1,1,1,1],
                                data_l_to_fmin:     vec![0,1,2]
                            };   
        
        let polys       =   vec![ poly_a, poly_b, poly_c, poly_d ];
                            
        for face_dim in 0 .. 5 {
            let mut faces_test    =   polys_faces( &polys, face_dim );
            let mut faces_true    =   HashSet::new();
            for poly in polys.iter() {
                faces_true.extend( poly_faces( &poly, face_dim ) )
            }
            let mut faces_true : Vec<_>     =   faces_true.iter().cloned().collect();
            
            faces_true.sort();
            faces_test.sort();        

            assert_eq!( &faces_true, &faces_test );

        }                            
    }

    #[test]
    fn  test_merge_down_flag_is_valid_for_simplex_factor_dims() {

        assert!(
            !
            merge_down_flag_is_valid_for_simplex_factor_dims(
                & vec![false, true, false, true, false],
                & vec![0, 2]
            )
        );
        assert!(
            !
            merge_down_flag_is_valid_for_simplex_factor_dims(
                & vec![false, false, true, true, false],
                & vec![0, 1, 1]
            )
        );
        assert!(
            merge_down_flag_is_valid_for_simplex_factor_dims(
                & vec![false, false, true, false, false],
                & vec![0, 1, 1]
            )
        );    
        assert!(
            merge_down_flag_is_valid_for_simplex_factor_dims(
                & vec![false, false, false, true, false],
                & vec![0, 1, 1]
            )
        );  
        assert!(
            merge_down_flag_is_valid_for_simplex_factor_dims(
                & vec![false, false, false, false, true],
                & vec![0, 1, 1]
            )
        );           
        assert!(
            merge_down_flag_is_valid_for_simplex_factor_dims(
                & vec![false, false, true, false],
                & vec![0, 1, 0]
            )
        );  
        assert!(
            merge_down_flag_is_valid_for_simplex_factor_dims(
                & vec![false, false, false, false],
                & vec![0, 1, 0]
            )
        );                                        

    }

    #[test]
    fn test_merge_down_flag_simplex_dims_to_lsord2fmin() {
        assert_eq!(
            vec![0, 1, 2, 2],
            merge_down_flag_simplex_dims_to_lsord2fmin(
                & vec![false, false, true, false, false],
                & vec![0, 1, 1]
            )
        );    
        assert_eq!(
            vec![0, 1, 2, 2],            
            merge_down_flag_simplex_dims_to_lsord2fmin(
                & vec![false, false, false, true, false],
                & vec![0, 1, 1]
            )
        );  
        assert_eq!(
            vec![0, 1, 1, 2],            
            merge_down_flag_simplex_dims_to_lsord2fmin(
                & vec![false, false, false, false, true],
                & vec![0, 1, 1]
            )
        );            
        assert_eq!(
            vec![0, 1, 2],
            merge_down_flag_simplex_dims_to_lsord2fmin(
                & vec![false, false, true, false],
                & vec![0, 1, 0]
            )
        );  
        assert_eq!(
            vec![0, 1, 1, 2],
            merge_down_flag_simplex_dims_to_lsord2fmin(
                & vec![false, false, false, false],
                & vec![0, 1, 0]
            )
        );                 
    }


}    
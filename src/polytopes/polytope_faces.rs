

use crate::polytopes::polytope::{Polytope};
use crate::intervals_and_ordinals::{ordinate};
use solar::utilities::combinatorics::fixed_sum_sequences;
// use crate::utilities::*;
use solar::utilities::index::{histogram};
// use ordered_float::OrderedFloat;
// use num::rational::Ratio;
// use std::collections::{HashMap};
// use std::hash::Hash;
use std::iter::{FromIterator, repeat};
// use std::cmp::Ord;
use itertools::Itertools;
// use itertools::structs::{MultiProduct, Combinations};

type Fil = usize;





struct SimprodFaceEnumerator{
    poly:   Polytope,

}

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

/// Vector mapping the ordinal of each "pre-merge" level set to the ordinal of the
/// corresponding "post-merge" level set (the one into which the original merges).
pub fn  vec_mapping_lsord_old_to_lsord_new(            
            verts_to_delete_per_simp:   & Vec< Vec< usize > >,
            simplex_dims:               & Vec< usize >,
        )
        -> 
        Vec< usize >
{
    // This code works by first placing each level set ordinal into an equivalence class.
    // At first, each equivalence class is functionally "labeled" by the smallest level
    // set ordinal that it contains.  After tall the classes have been computed, however,
    // we relabel them.
    let num_lsord_old               =   simplex_dims.iter().sum::<usize>() + simplex_dims.len(); // each dim differs from set cardinality by 1, so we must add simplex_dims.len() to compensate
    let mut lsord_old_to_lsord_new  =   Vec::with_capacity(num_lsord_old ); // level-set-ordinal-old to level-set-ordinal-new

    let mut merge_down_flag         =   Vec::from_iter( repeat(false).take( num_lsord_old ) );

    let mut offset                  =   0;
    for (simplex_count, verts_to_delete) in verts_to_delete_per_simp.iter().enumerate() {
        offset                      =   simplex_dims.iter().take(simplex_count).sum::<usize>() + simplex_count;
        for vert in verts_to_delete.iter().map(|x| x + offset ) {
            merge_down_flag[ vert + 1 ]   =   true;
        }
    } 

    // Calculate the minimal fmin ord with which each ord merges
    for lsord_old in 0 .. num_lsord_old {
        let mut min_equiv_ord       =   lsord_old.clone();
        while merge_down_flag[ min_equiv_ord ] {
            min_equiv_ord -= 1;
        }
        lsord_old_to_lsord_new.push( min_equiv_ord );
    }

    let ordinal_data                =   ordinate( &lsord_old_to_lsord_new );
    for val in lsord_old_to_lsord_new.iter_mut() { *val = ordinal_data.ord( &val ).unwrap() }

    lsord_old_to_lsord_new




    // lsord_old_to_lsord_new.extend( repeat(0).take( 1 + simplex_dims[0] - verts_to_delete_per_simp[0].len() ) );
    // let mut merged_last             =   verts_to_delete_per_simp[0].contains( &simplex_dims[0] );

    // for simplex_count in 1 .. simplex_dims.len() {
    //     let verts_to_delete         =   & verts_to_delete_per_simp[ simplex_count ]; // the vertices to delete from this simplex, which correspond to **gaps** between level sets
    //     let simplex_dim             =   simplex_dims[ simplex_count ]; // = one less than the number of level sets with this fmin value        
    //     for vertex_count in 0 .. simplex_dim + 1 {
    //         let last_new_ord        =   lsord_old_to_lsord_new.last().unwrap().clone();
    //         match merged_last {
    //             true    =>  lsord_old_to_lsord_new.push( last_new_ord ),
    //             false   =>  lsord_old_to_lsord_new.push( last_new_ord + 1 ),
    //         }
    //         merged_last             =   verts_to_delete.contains( &vertex_count );
    //     }
    // }

    // lsord_old_to_lsord_new


    // let mut num_removed             =   0;
    // let mut global_position         =   0;

    // // run over (# of critical values) = (# of simplices)
    // for simplex_count in 0 .. simplex_dims.len() {
    //     let verts_to_delete         =   & verts_to_delete_per_simp[ simplex_count ]; // the vertices to delete from this simplex, which correspond to **gaps** between level sets
    //     let simplex_dim             =   simplex_dims[ simplex_count ]; // = one less than the number of level sets with this fmin value
        

    //     let mut merged_last         =   false;

    //     for vertex_count in 0 .. simplex_dim  





    //     // remove all the vertices we need to, except the last one, if the last one happens to be the last vetex in the simplex
    //     for vertex_count in 0 .. simplex_dim  
    //     {     
    //         if verts_to_delete.contains( &vertex_count ) { 
    //             num_removed += 1;
    //         }
    //         global_position += 1;
    //         lsord_old_to_lsord_new.push(
    //             global_position - num_removed
    //         );
    //     }
    //     // handle last vertex separately; in this case we merge the next level set up into the last of this simplex's level sets,
    //     // so all cells that used to get an fmin value of x will now get an fmin value of x+1
    //     global_position += 1;
    //     lsord_old_to_lsord_new.push(
    //         global_position - num_removed
    //     );
    //     if verts_to_delete.contains( &simplex_dim ) { 
    //         num_removed += 1;
    //     }        


    //     // if verts_to_delete.contains( &simplex_dims[ simplex_count ] ) { 
    //     //     lsord_old_to_lsord_new.push(
    //     //         global_position - num_removed
    //     //     );            
    //     //     num_removed += 1;
    //     // }        
    // }

    // lsord_old_to_lsord_new
}

/// Return a face of the underlying polytope, corresponding to deletion of the vertices provided. 
pub fn  polytope_face(            
            verts_to_delete_per_simp:   &   Vec< Vec< usize > >,
            poly:                       &   Polytope,
        )
        ->
        Polytope
{
    let mut poly            =   poly.clone();    
    let simplex_dims        =   poly.simplex_factor_dims_cellagnostic();
    let translator          =   vec_mapping_lsord_old_to_lsord_new(            
                                    & verts_to_delete_per_simp,
                                    & simplex_dims,
                                );
    
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
/// Assumes that the last level set corresponds to cells with values fixed at 1 (or whatever the
/// maximum value is).
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







#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    use crate::polytopes::polytope::{random_polytope, poly_to_polyvvv, polyvvv_to_poly};


    #[test]
    fn test_vec_mapping_lsord_old_to_lsord_new() {
        
        let simplex_dims                =   vec![1, 3, 0];
        let verts_to_delete_per_simp    =   vec![ vec![1], vec![0], vec![] ];
        
        let translator_test             =   vec_mapping_lsord_old_to_lsord_new(            
                                                & verts_to_delete_per_simp,
                                                & simplex_dims
                                            );
        
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
    fn test_poly_to_face_and_face_enumerate() {

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

            for face_dim in 0 .. polytope_dim {
                
                let mut face_pool           =   Vec::new();

                let deletion_choices        =   vertex_deletion_choices_for_product_of_combinatorial_simplices(
                                                    & simplex_factor_dims,
                                                    face_dim
                                                );
                for deletion_choice in deletion_choices.iter() {
                    let face_poly           =   polytope_face(
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


}    





pub fn  fixed_sum_sequences(
        caps:   Vec< usizw >,
        target_sum:    usize
        )
        ->
        Vec< Vec< usize> >
{
    let cap_aggregate    =   caps.sum();

    if target_sum > cap_aggregate { Vec::with_capacity(0) }
    else if let Some( last_ind ) = caps.last_index() {

        let trunc_caps   =   Vec::from_iter( caps[0 .. last_ind ] );
        let trunc_cap_agg  =   trunc_caps.sum();
        
        let last_min    =   
            match trunc_cap_agg < sum {
                true    =>  sum - trunc_cap_add,
                false   =>  0
            }
        
        let sequences   =   Vec::new();

        for last_val    in   last_min .. target_sum + 1 {
            
            let trunc_target_sum    =   target_sum - last_val;
            let trunc_sequences     =   fixed_sum_sequences(
                                            trunc_caps,
                                            trunc_target_sum,
                                        );
            for trunc_seq in trunc_sequences { trunc_seq.push!( last_val.clone() ) }  {
            sequences.append( &mut trunc_sequences )
        }

        sequences

    } else { // THIS ONLY HAPPENS IF CAPS IS EMPTY
        vec![vec![]] 
    }

}



struct SimprodFaceEnumerator{
    poly:   polytope,

}

pub fn  face_deletion_iter_from_product_of_combinatorial_simplices(
            dims:       Vec< usize >,
            face_dim:   usize
            fn 
        )
{
    fixed_sum_sequences(dims, dims.sum() - face_dim )
        .map(   |x| 
                x.iter().enumerate().map( |(y_count, y)| (0..dims[y_count + 1]).combinations(y) ).multi_cartesian()
        )
        .chain()
}


pub fn  vec_mapping_lsord_old_to_lsord_new(            
            verts_to_delete_per_simp:   Vec< Vec< usize > >,
            simplex_dims:               Vec< usize >,
        )
{
    let num_lsord_old               =   simplex_dims.len() + simplex_dims.sum(); // dims differ by set cardinality by 1
    let mut lsord_old_to_lsord_new  =   Vec::from_iter( 0 .. num_lsord_old ); // level-set-ordinal-old to level-set-ordinal-new
    let mut num_removed             =   0;
    let mut global_position         =   0;

    for simplex_count in 0 .. simplex_dims.len() {
        let verts_to_delete     =   verts_to_delete_per_simp[ simplex_count ]
        let simplex_dim         =   simplex_dims[ simplex_count ]
        // intentionally excluding last vertex
        for vertex_count in 0 .. simplex_dim  
        {     
            if verts_to_delete.contains( vertex_count ) { 
                num_removed += 1;
            }
            global_position += 1;
            lsord_old_to_lsord_new.push(
                global_position - num_removed
            )
        }
        // handle last vertex separately; in this case we merge the next level set up into the last of this simplex's level sets
        if verts_to_delete.contains( simplex.dims[ simplex_count ] ) { 
            num_removed += 1;
        }        
    }

}

/// Return a face of the underlying polytope, corresponding to deletion of the vertices provided. 
pub fn  polytope_face(            
            verts_to_delete_per_simp:   Vec< Vec< usize > >,
            poly:                       Polytope,
        )
{
    let poly                =   poly.clone();    
    let simplex_dims        =   histogram( poly.data_c_to_l );
    let translator          =   vec_mapping_lsord_old_to_lsord_new(            
                                    verts_to_delete_per_simp,
                                    simplex_dims,
                                )
    
    // update mapping from cells to level sets                            
    for cell_id in 0 .. poly.num_cells() {
        poly.data_c_to_l[ cell_id ]    =    translator[
                                                poly.data_c_to_l[ cell_id ]
                                            ];
    }

    // update mapping from cells to level sets
    poly.data_l_to_fmin.clear();
    for simplex_count in poly.num_lev_sets() {
        new_num_lev_sets    =   1 + simplex_dims[ simplex_count] - verts_to_delete_per_simp[ simplex_count ].len();
        poly.data_l_to_fmin.extend_from(  iter::repeat( simplex_count ).take( new_num_lev_sets )  )
    }

    poly
}        

/// Return a face of the underlying polytope, corresponding to deletion of the vertices provided. 
/// 
pub fn  polytope_setform_face(
            verts_to_delete_per_simp:   Vec< Vec< usize > >,
            poly:                       Vec< Vec< usize > >,
        )    
{
    let face                =   poly.clone();   
    let num_crit_vals       =   face.len() - 1; // the last bin does not count as a critical value, since topology doesn't change there
    for crit_val_ord in ( 0 .. num_crit_vals ).rev()  
    {
        let mut lev_sets_this_val   =   poly[ crit_val ];
        let num_lev_sets_this_val   =   lev_sets_this_val.len();
        let verts_to_delete         =   verts_to_delete_per_simp[ crit_val_ord ];

        for gap_ord_to_delete in verts_to_delete.iter().rev() 
        {
            let cell_ids            =   lev_sets_this_val.remove( gap_ord );   // this reduces length of lev_sets_this_val by 1         
            if gap_ord_to_delete == lev_sets_this_val.len() {
                // if we are deleting the last gap, then merge the last level set into the first level set of the next critical value
                face[ crit_val_ord + 1 ][ 0 ].append( &mut cell_ids );
                face[ crit_val_ord + 1 ][ 0 ].sort();
            } else {
                // otherwise merge it with the next level set within the current critical value
                lev_sets_this_val[ gap_ord_to_delete ].append( &mut cell_ids );
                lev_sets_this_val[ gap_ord_to_delete ].sort();                
            }
        }
    }

    face
}     


pub fn  poly_to_setform_poly( poly: Polytope ) {
    
}
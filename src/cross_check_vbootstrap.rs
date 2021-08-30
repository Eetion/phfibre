



/// Convert a path in the digraph into a polytope
pub fn  digraph_path_to_polytope( 
            path:           &   Vec< usize >,
            edge_grades:    &   Vec< usize >,
            fibre_verts:    &   Vec< Vec< usize > >,
        )
        ->
        Vec< Vec< usize > >
{

    let num_cells_in_filtered_space         =   fibre_verts[0].len();
    let num_crit_vals                       =   fibre_verts[0].iter().max().unwrap();

    // Initialize polytope data
    let mut data_c_to_l                     =   Vec::with_capacity( num_cells_in_filtered_space );
    let mut num_lev_sets_per_crit_val       =   Vec::from_iter( 0 .. num_crit_vals );
    for grade in edge_grades { num_lev_sets_per_crit_val[ grade ] += 1 };  // count: # level sets with critical val k == 1 + # degree-k edges
    let data_l_to_fmin                      =   MAKE LEV SET SIZES FROM num_lev_sets_per_crit_val // we will increase these values as needed
    
    let mut crit_val                        =   num_crit_vals; // this is just a trick, since we can't set it to -1
    let mut within_crit_val_lev_set_ord     =   0;

    FIRST ASSIGN THE MINIMUM POSIBLE (CRITICAL) LEVEL SET ORDINAL TO EACH CELL IN THE FILTERED SPACE; WE WILL THEN "BUMP CELLS UP" IN THE RANKING AS NEEDED

    for path_vertex_count in 1 .. path.len() {

        if edge_grades[ path_vertex_count ] !=  crit_val  {
            crit_val                        =   edge_grades[ path_vertex_count ].clone();
            lev_set_ord_cur                 =   CALCULATE THE LOWEST ORDINAL (THE ACTUALLY CRITICAL ONE) FOR THIS CRITICAL VALUE, THEN ADD 1;
        } else {
            lev_set_ord_cur += 1;
        }

        let fibre_vertex_now                =   fibre_verts[ path_vertex_count ];
        let fibre_vertex_old                =   fibre_verts[ path_vertex_count -1 ];        
        for cell_id in 0 .. num_cells_in_filtered_space {
            if fibre_vertex_now[ cell_id ] > fibre_vertex_old[ cell_id ] {  
                data_c_to_l[ cell_id ]      =   lev_set_ord_cur;
            }
        }
    }

    return  Polytope{   data_c_to_l: data_c_to_l, data_l_to_fmin: data_l_to_fmin }
}


/// Tabulate the max-with-respect-to-inclusion paths in an edge-graded digraph,
///  subject to conidiiton that edge degrees must increase monotonically.
pub fn  maximal_degree_ascending_paths( 
            digraph: &DigraphGraded 
        )
        ->
        Vec< Vec< usize > >
{

}


/// Digraph with edges grouped by degree.
/// 
/// Arcs are stored in "sparse row" format.
pub struct DigraphGraded{
    pub arcs:   Vec< Vec< HashSet< usize > > > 
}


/// The graded digraph G' whose degree-k arcs are the covering relations of the
/// degree-k edges of G.  
/// 
/// Assumes that the degree-k edges of G form a partial order.
fn  digraph_to_atomic_digraph( 
        digraph: & DigraphGraded 
    ) 
    -> 
    DigraphGraded
{
    // RULE: ADD ARC (f,g) IFF THERE EXISTS NO PAIR OF ARCS (f,h),(h,g) in the input digraph
}



/// Given the vertices of the PH fibre, construct the associated digraph.
fn  fibre_vertices_to_digraph( verts: & Vec< Vec< usize > > ) 
    -> 
    BiDigraph
    {

    if verts.is_empty() { return vec![] }

    let num_verts               =   verts.len();
    let num_cells               =   verts[0].len();

    let mut bidigraph           =   BiDigraph::new( num_verts );

    for vert_a in 0 .. num_verts {
        for vert_b in vert_a .. num_verts {
            for cell_id in 0 .. num_cells {
                match   verts[ vert_a ][ cell_id ]
                            .cmp(  
                                verts[ vert_b ][ cell_id ]  
                            )
                    {
                        Equal   =>  {   continue    },
                        Less    =>  {
                                        fibre_vertices_to_digraph_subroutine_check_remainder( 
                                            vert_a,         // vert_ind_mergedn
                                            vert_b,         // vert_ind_mergeup
                                            cell_id,        // cell_id_initial
                                            num_cells,      // num_cells_total
                                            & verts,        // vertices
                                            &mut bidigraph  // bidigraph
                                        )                            
                                    },
                        Greater =>  {
                                        fibre_vertices_to_digraph_subroutine_check_remainder( 
                                            vert_b,         // vert_ind_mergedn
                                            vert_a,         // vert_ind_mergeup
                                            cell_id,        // cell_id_initial
                                            num_cells,      // num_cells_total
                                            & verts,        // vertices
                                            &mut bidigraph  // bidigraph
                                        )                            
                                    }                                    

                    }
                }
            }
        }
    }

    return bidigraph
}


fn  fibre_vertices_to_digraph_subroutine_check_remainder( 
        vert_ind_mergedn:               usize,
        vert_ind_mergeup:               usize,
        cell_id_initial:                usize,
        num_cells_total:                usize,
        vertices:                       & Vec< Vec< usize > >,
        bidigraph:                      & mut BiDigraph,
    )
{
    let vert_dn                     =   vertices[ vert_ind_mergedn ];
    let vert_up                     =   vertices[ vert_ind_mergeup ];  
    let crit_val                    =   vert_dn[ cell_id_initial ];
    let critval_p_1                 =   crit_val + 1;

    let mut push_to_graph           =   false;

    for cell_id in cell_id_initial .. num_cells_total {
        if lower[ cell_id ]     ==  upper[ cell_id ] 
            {   continue } // test continues
        else 
        if  lower[ cell_id ]     !=  crit_val
            || 
            upper[ cell_id ]     !=  crit_val_p_1 
            {   push_to_graph = false;  break } // test fails
        else 
            {   push_to_graph = true    } // test will pass if no future problems are discovered
    }

    if push_to_graph {
        push_edge( digraph, vert_ind_mergedn, vert_ind_mergeup, crit_val );
    }
}
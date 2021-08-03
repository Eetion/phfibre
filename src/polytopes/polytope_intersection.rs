
//  ---------------------------------------------------------------------------  
//  POLYTOPE VERTEX BUFFER
//  --------------------------------------------------------------------------- 


struct PolytopeVertexBuffer<'a> {
    poly:               &'a Polytope,      // polytope to work with
    cell_id_to_fmin:     Vec< Fil >,         // function (cell id) -> min possible filtration value    
    cell_id_to_fmax:     Vec< Fil >,         // function (cell id) -> max possible filtration value
}

fn poly_2_polysolve<'a> ( poly: &'a Polytope ) -> PolytopeVertexBuffer<'a>{

    let num_cells       =   poly.num_cells();
    let num_levels      =   poly.num_lev_sets();    

    let mut cell_id_to_fmin  =   Vec::with_capacity( num_cells );
    let mut cell_id_to_fmax  =   Vec::with_capacity( num_cells );

    // BOUNDARY CASES: #{CELLS} = 0 OR #{LEVEL SET} = 1
    if num_cells == 0 {
        return PolytopeVertexBuffer{ poly:poly, cell_id_to_fmin:cell_id_to_fmin, cell_id_to_fmax:cell_id_to_fmax}        
    } else 
    // this represents the case where there are 1 or more cells, but exactly 1 level set
    // (in this case the level set must be critical, on account of the infinite dim-0 bar)
    if num_levels == 1 {
        for cell_id in 0..num_cells {
            cell_id_to_fmin.push(0);
            cell_id_to_fmax.push(0);
        }
        return PolytopeVertexBuffer{ poly:poly, cell_id_to_fmin:cell_id_to_fmin, cell_id_to_fmax:cell_id_to_fmax}                
    }

    // REMAINING CASE: 
    //      >= 1 cell           AND
    //      >= 2 level sets 

    for i in 0..num_cells { 
        cell_id_to_fmin.push( poly.cell_id_to_fmin( i ).unwrap() );
        cell_id_to_fmax.push( poly.cell_id_to_fmax( i ).unwrap() );        
    }

    PolytopeVertexBuffer{ poly: &poly, cell_id_to_fmin: cell_id_to_fmin, cell_id_to_fmax: cell_id_to_fmax}

}


//  ---------------------------------------------------------------------------  
//  POLYTOPE INTERSECTION FUNCTION
//  --------------------------------------------------------------------------- 


// fn check_intersection(
//     poly_a:     Polytope,
//     poly_b:     Polytope
// )

// {
//     // define solver a
//     let mut 
//     let mut solver_a    =   PolytopeVertexBuffer{   poly: &poly_a, 2 }
// }




#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;


    #[test]
    fn test_cell_id_to_critical_height() { 
        let poly                =   Polytope { 
                                        data_l_to_fmin:     vec![0,1,1,2],
                                        data_c_to_l:        vec![3,2,1,0,1,2,3],
                                    };
        assert_eq!( poly.cell_id_to_critical_height(0), Some(0) );
        assert_eq!( poly.cell_id_to_critical_height(1), Some(1) );        
        assert_eq!( poly.cell_id_to_critical_height(2), Some(0) );                
        assert_eq!( poly.cell_id_to_critical_height(3), Some(0) );                        
        assert_eq!( poly.cell_id_to_critical_height(10), None );                        
    }    

    #[test]
    fn test_making_barocde_from_raw_parts()
    {

        let barcode             =   Barcode::new(
                                        vec![0], // barcode_inf_dim,
                                        to_ordered_float( &vec![0.] ), //barcode_inf_brn,
                                        Vec::new(), //barcode_fin_dim,
                                        Vec::new(), //barcode_fin_brn,
                                        Vec::new(), //barcode_fin_die
                                    );    

        let mut val_to_ord      =   HashMap::new();
        val_to_ord.insert( OrderedFloat(0.0), 0);

        let barcode_true    =   Barcode{ 
                                    inf: vec![ BarInfinite{ dim: 0, birth: 0} ],
                                    fin: vec![],
                                    ordinal: OrdinalData{ 
                                                ord_to_val: vec![OrderedFloat(0.0)],
                                                val_to_ord: val_to_ord    
                                            }
                                };
        
        
        
    }     

    #[test]
    fn test_barcode_from_pairs() {
    
        let pairs   =   vec![ (0,1), (2,3) ];
        let dims    =   vec![ 0,    1,      1,      2,      3   ];
        let births  =   vec![ 0,    1,      1,      1,      2  ];

        let barcode =   pairs_dims_births_to_barcode(
                            & pairs,
                            & dims,
                            & births,
                        );

        println!("{:?}", & barcode);                        

        assert_eq!(     &   barcode.inf, 
                        &   vec![ BarInfinite{ dim: 3, birth: 2 } ] 
        );

        assert_eq!(     &   barcode.fin, 
            &   vec![ 
                    BarFinite{ dim: 0, birth: 0, death: 1 },
                ] 
        );
        
        assert_eq!(     &   barcode.ordinal.ord_to_val, 
                        &   vec![0, 1, 2],
        );

    }

    #[test]
    fn test_is_compatible_with_ordinal_filtration() {

        let polytope                =   Polytope{ 
                                            data_c_to_l: vec![2, 0, 1, 6, 6, 6],
                                            data_l_to_fmin: vec![0, 1, 1],
                                        };

        assert_eq!(     
            true, 
            polytope.min_vertex_is_compatible_with_ordinal_filt( &vec![2, 0, 1, 4, 4, 4] )
        );
        assert_eq!(     
            true, 
            polytope.min_vertex_is_compatible_with_ordinal_filt( &vec![1, 0, 1, 4, 4, 4] )
        );  

        let polytope                =   Polytope{ 
                                            data_c_to_l: vec![2, 0, 1, 6, 6, 6],
                                            data_l_to_fmin: vec![0, 1, 2],
                                        };

        assert_eq!(     
            false, 
            polytope.min_vertex_is_compatible_with_ordinal_filt( &vec![1, 0, 1, 4, 4, 4] )
        );                
            

    }

}    
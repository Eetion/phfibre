
use crate::polytopes::polytope::Polytope;
use crate::polytopes::polytope_faces::{ poly_faces, 
                                        merge_down_flag_to_vec_mapping_lsord_old2new,
                                        merge_down_flag_is_valid_for_simplex_factor_dims,
                                        merge_down_flag_simplex_dims_to_lsord2fmin,
                                    };
// use crate::utilities::*;
use solar::utilities::index::compose_f_after_g;
// use ordered_float::OrderedFloat;
// use num::rational::Ratio;
// use std::collections::{HashMap};
// use std::hash::Hash;
use std::iter::{repeat, FromIterator};
// use std::cmp::Ord;




//  ---------------------------------------------------------------------------  
//  POLYTOPE VERTEX BUFFER
//  --------------------------------------------------------------------------- 


//  UNTESTED
// struct PolytopeVertexBuffer<'a> {
//     poly:               &'a Polytope,      // polytope to work with
//     cell_id_to_fmin:     Vec< Fil >,         // function (cell id) -> min possible filtration value    
//     cell_id_to_fmax:     Vec< Fil >,         // function (cell id) -> max possible filtration value
// }

// impl    < 'a >  
//         PolytopeVertexBuffer
//         < 'a > 
// {


    // UNTESTED
    // fn make_crit_low( &self, cell_id: usize ) -> Option<()> {
    //     let crit_ord        =   self.poly.cell_id_to_fmin( cell_id ).unwrap();
    //     let lev_set_ord     =   self.poly.cell_id_to_lev_set_ord( cell_id ).unwrap();

    //     // for each cell in the complex
    //     for other_cell_id in 0 .. self.poly.num_cells() {
    //         if  crit_ord
    //             == // belonging to the same family of "free" level sets
    //             self.poly.cell_id_to_fmin( &other_cell_id ).unwrap()
    //             &&
    //             lev_set_ord
    //             >= // and appearing at or before the level set that contains this cell
    //             self.poly.cell_id_to_lev_set_ord( &other_cell_id ).unwrap() 
    //         {
    //             // if the other cell has already been pinned to a higher value, then
    //             // report that the inteersection problem is insoluble, and terminate
    //             if  self.cell_id_to_fmin[ other_cell_id ]
    //                 ==
    //                 self.cell_id_to_fmax[ other_cell_id ]
    //                 &&
    //                 self.cell_id_to_fmin[ other_cell_id ]
    //                 >
    //                 crit_ord
    //             {
    //                 return None
    //             }  else {
    //                 // require that this cell take the minimum possible value
    //                 self.cell_id_to_fmax[ other_cell_id ]   =   crit_ord.clone();
    //             }              
    //         }
    //     }
    //     return Some( () )
    // }

    // UNTESTED    
    // fn make_crit_high( &self, cell_id: usize ) -> Option<()> {
    //     let crit_ord        =   self.poly.cell_id_to_fmax( cell_id ).unwrap();
    //     let lev_set_ord     =   self.poly.cell_id_to_lev_set_ord( cell_id ).unwrap();

    //     // for each cell in the complex
    //     for other_cell_id in 0 .. self.poly.num_cells() {
    //         if  crit_ord
    //             == // belonging to the same family of "free" level sets
    //             self.poly.cell_id_to_fmax( &other_cell_id ).unwrap()
    //             &&
    //             lev_set_ord
    //             <= // and appearing at or before the level set that contains this cell
    //             self.poly.cell_id_to_lev_set_ord( &other_cell_id ).unwrap() 
    //         {
    //             // if the other cell has already been pinned to a higher value, then
    //             // report that the inteersection problem is insoluble, and terminate
    //             if  self.cell_id_to_fmin[ other_cell_id ]
    //                 ==
    //                 self.cell_id_to_fmax[ other_cell_id ]
    //                 &&
    //                 self.cell_id_to_fmin[ other_cell_id ]
    //                 <
    //                 crit_ord
    //             {
    //                 return None
    //             }  else {
    //                 // require that this cell take the minimum possible value
    //                 self.cell_id_to_fmin[ other_cell_id ]   =   crit_ord.clone();
    //             }              
    //         }
    //     }
    //     return Some( () )
    // }
// 
// }


//  UNTESTED
// fn poly_2_intersection_solve_buffer<'a> ( poly: &'a Polytope ) -> PolytopeVertexBuffer<'a>{

//     let num_cells       =   poly.num_cells();
//     let num_levels      =   poly.num_lev_sets();    

//     let mut cell_id_to_fmin  =   Vec::with_capacity( num_cells );
//     let mut cell_id_to_fmax  =   Vec::with_capacity( num_cells );

//     // BOUNDARY CASES: #{CELLS} = 0 OR #{LEVEL SET} = 1
//     if num_cells == 0 {
//         return PolytopeVertexBuffer{ poly:poly, cell_id_to_fmin:cell_id_to_fmin, cell_id_to_fmax:cell_id_to_fmax}        
//     } else 
//     // this represents the case where there are 1 or more cells, but exactly 1 level set
//     // (in this case the level set must be critical, on account of the infinite dim-0 bar)
//     if num_levels == 1 {
//         for cell_id in 0..num_cells {
//             cell_id_to_fmin.push(0);
//             cell_id_to_fmax.push(0);
//         }
//         return PolytopeVertexBuffer{ poly:poly, cell_id_to_fmin:cell_id_to_fmin, cell_id_to_fmax:cell_id_to_fmax}                
//     }

//     // REMAINING CASE: 
//     //      >= 1 cell           AND
//     //      >= 2 level sets 

//     for i in 0..num_cells { 
//         cell_id_to_fmin.push( poly.cell_id_to_fmin( i ).unwrap() );
//         cell_id_to_fmax.push( poly.cell_id_to_fmax( i ).unwrap() );        
//     }

//     PolytopeVertexBuffer{ poly: &poly, cell_id_to_fmin: cell_id_to_fmin, cell_id_to_fmax: cell_id_to_fmax}

// }


//  ---------------------------------------------------------------------------  
//  POLYTOPE INTERSECTION FUNCTION
//  --------------------------------------------------------------------------- 



/// If the polytopes have nonempty intersection, return their intersection.  
/// Otherwise return `None`.
pub fn  polytope_intersection(
            poly_a:     &Polytope,
            poly_b:     &Polytope
        )
        ->
        Option< Polytope >
{
    // This function is modeled off of vec_mapping_lsord_old_to_lsord_new; see that
    // code to get a better idea of how this works.

    // PHASE 1: MERGES DETERMINED ONLY BY LEVEL SET ORDINALS

    // Trivial checks
    if  poly_a.num_simplex_factors()    !=  poly_b.num_simplex_factors() 
        ||
        poly_a.num_cells()              !=  poly_b.num_cells()     
    { return None }

    let num_cells                       =   poly_a.num_cells();   
    let num_lev_sets_a                  =   poly_a.num_lev_sets();

    // We will use data about poly_b to decide how to merge level sets in poly_a.
    let mut merge_down_flag         =   Vec::from_iter( 
                                            repeat(false)
                                            .take( 
                                                num_lev_sets_a
                                            ) 
                                        );

    // Each time we see that poly_a allows a cell to take at least one value which
    // poly_b does not, do the only thing we can do: merge classes until that cell
    // joins a critial class, thus restricting its range.
    for cell_id in 0 .. num_cells {
        let fmin_a                  =   poly_a.cell_id_to_fmin( cell_id ).unwrap();
        let fmin_b                  =   poly_b.cell_id_to_fmin( cell_id ).unwrap();            
        if  fmin_a < fmin_b {  
                
            // if the class for poly_a is free to move up to the fmin level of b, then we must make it critical, and move it up to that value
            if fmin_b == poly_a.cell_id_to_fmax( cell_id ).unwrap() {
                let mut dummy_ord       =   poly_a.cell_id_to_lev_set_ord( cell_id ).unwrap();            
                while ! poly_a.lev_set_ord_to_is_critical( dummy_ord ).unwrap() {
                    merge_down_flag[ dummy_ord + 1 ] = true;
                    dummy_ord += 1;
                }
            }
            // otherwise there can be no intersection
            else {
                return None
            }

        }
        let fmax_a                  =   poly_a.cell_id_to_fmax( cell_id ).unwrap();
        let fmax_b                  =   poly_b.cell_id_to_fmax( cell_id ).unwrap();            
        if  fmax_b < fmax_a { 
            // if the class for poly_a is free to move down to the fmax level of b, then we must make it critical, and move it down to that value       
            if  fmax_b == poly_a.cell_id_to_fmin( cell_id ).unwrap() {
                let mut dummy_ord       =   poly_a.cell_id_to_lev_set_ord( cell_id ).unwrap();            
                while ! poly_a.lev_set_ord_to_is_critical( dummy_ord ).unwrap() {
                    merge_down_flag[ dummy_ord ] = true;
                    dummy_ord -= 1;
                }
            }
            // otherwise there can be no intersection            
            else {
                return None
            }
        }

    // Each time we see that to cells belong to the same level set in poly_b, we 
    // will declare that the corresponding two level sets in poly_a (plus all the
    // level sets in between) must merge 
    for cell_id_x in 0 .. num_cells {
        for cell_id_y in 0 .. num_cells {
            if  poly_b.cell_id_to_lev_set_ord( cell_id_x ) 
                <= 
                poly_b.cell_id_to_lev_set_ord( cell_id_y )
            {
                for dummy_ord in 
                        // notice that the order of x and y is reversed here; this is intentional
                        poly_a.cell_id_to_lev_set_ord( cell_id_y ).unwrap() + 1
                        ..
                        poly_a.cell_id_to_lev_set_ord( cell_id_x ).unwrap() + 1
                {
                    merge_down_flag[ dummy_ord ] = true;
                } 
            }
        }
    }

        
    }      

    // We have now decided how level sets of a will merge.  The following test 
    // checkes whether this is a valid merge.
    if  !   
        merge_down_flag_is_valid_for_simplex_factor_dims(
            &   merge_down_flag,
            &   poly_a.simplex_factor_dims_cellagnostic()
        )
    {   return  None    }

    // If we have made it to this point, then (i) we have only performed merges that 
    // absolutely MUST be performed, and (ii) in doing, we have produced a face of poly_a 
    // in which respects the following 
    // (i)  bounds of form filt_val( cell_x ) <= filt_val( cell_ y)  <-- this is acheived by merging according to level set ordinals, i.e. the for-loop that runs over cell_id_x and cell_id_y
    // (ii) bounds of form crit_val_a < filt_val( cell_x ) < crit_val_b
    // I think this suffices to ensure we have obtained a face of poly_b.  
    // !!! THEREFORE WE HAVE COMPUTED THE INTERSECTION OF poly_a AND poly_b, AND THE
    //     INTERSECTION IS NONEMPTY.

    // make new mapping from level set ordinals to critical values
    let data_l_to_fmin_new          =   
        merge_down_flag_simplex_dims_to_lsord2fmin(
            &   merge_down_flag,
            &   poly_a.simplex_factor_dims_cellagnostic()
        );
   
    // make new mapping from cells to level sets
    let translator                  =   merge_down_flag_to_vec_mapping_lsord_old2new( & merge_down_flag );
    let data_c_to_l_new             =   compose_f_after_g(
                                            &   translator,
                                            &   poly_a.data_c_to_l
                                        );

    // let poly_med                    =   Polytope{
    //                                         data_c_to_l:    data_c_to_l_new,
    //                                         data_l_to_fmin: data_l_to_fmin_new
    //                                     };

    // PHASE 2: MERGES THAT INVOLVE FMIN

    // for cell_id_x in 0 .. num_cells {
    //     for cell_id_y in 0 .. num_cells {
    //         if  poly_b.cell_id_to_lev_set_ord( cell_id_x ) 
    //             <= 
    //             poly_b.cell_id_to_lev_set_ord( cell_id_y )
    //         {
    //             for dummy_ord in 
    //                     // notice that the order of x and y is reversed here; this is intentional
    //                     poly_a.cell_id_to_lev_set_ord( cell_id_y ).unwrap() + 1
    //                     ..
    //                     poly_a.cell_id_to_lev_set_ord( cell_id_x ).unwrap() + 1
    //             {
    //                 merge_down_flag[ dummy_ord ] = true;
    //             } 
    //         }
    //     }
    // }   
    
    Some(        
        Polytope{
            data_l_to_fmin:     data_l_to_fmin_new,
            data_c_to_l:        data_c_to_l_new,
        }
    )


}


// OLD VERSION -- COMMENTED TO MAKE SURE I DON'T EDIT IT
// pub fn  polytope_intersection(
//             poly_a:     &Polytope,
//             poly_b:     &Polytope
//         )
//         ->
//         Option< Polytope >

// {

//     let num_cells       =   poly_a.num_cells();

//     // place polytopes into a vector
//     let o_poly          =   vec![ poly_a, poly_b ];

//     // make vectors of length (# cells) to keep track of max and min values

//     let a_fmin          =   poly_a.vec_mapping_cell_id_to_min_filt_ordinal();
//     let a_fmax          =   poly_a.vec_mapping_cell_id_to_max_filt_ordinal();  
    
//     let b_fmin          =   poly_a.vec_mapping_cell_id_to_min_filt_ordinal();
//     let b_fmax          =   poly_a.vec_mapping_cell_id_to_max_filt_ordinal();      

//     // "o" stands for "omnia" (meaning everything)
//     let mut o_fmin      =   [a_fmin, b_fmin];
//     let mut o_fmax      =   [a_fmax, b_fmax];    

//     loop {

//         let mut break_loop     =   true;

//         for cell_id in 0 .. num_cells {

//             // there is nothing to do for this cell in particular if it has the same max/min filtration values in both polytopes
//             if  o_fmin[ 0 ][ cell_id ]  ==  o_fmin[ 1 ][ cell_id ]
//                 &&
//                 o_fmax[ 0 ][ cell_id ]  ==  o_fmax[ 1 ][ cell_id ]
//             { continue }

//             // find the buffer with the (lexicographically) lower pair (cell_id min val, cell_id max val)
//             let i_lower                 =   (0..2)
//                                                 .min_by_key(  
//                                                     |x| 
//                                                     (
//                                                         o_fmin[ *x ][ cell_id ],
//                                                         o_fmax[ *x ][ cell_id ],                                                                
//                                                     )        
//                                                 ).unwrap(); 
//             // get the index for the other buffer                                                        
//             let i_upper                 =   1 - i_lower;

//             // calculate the central value to which all up/downstream cells must now be bound
//             let val_lower               =   o_fmax[ i_lower ][ cell_id ].clone();                
//             let val_upper               =   o_fmin[ i_upper ][ cell_id ].clone();                
            
//             if val_upper != val_lower { return None }
//             let val_center              =   val_lower.clone();                

//             // bind the up/down stream values
//             for upstream_cell_id  in o_poly[ i_lower ].cell_id_to_upstream_cell_ids( cell_id ) {
//                 o_fmin[ i_lower ][ upstream_cell_id ]  =   val_center.clone();
//             }
//             for dnstream_cell_id  in o_poly[ i_upper ].cell_id_to_dnstream_cell_ids( cell_id ) {
//                 o_fmax[ i_upper ][ dnstream_cell_id ]  =   val_center.clone();
//             }   
            
//             break_loop                 =   true;
//         }

//         if break_loop {
//             // in this case we have made it through 
//             return true
//         }
//     }
// }





#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    use crate::polytopes::polytope::{enumerate_poly_test_set_40_total};
    use itertools::Itertools;
    use std::collections::{HashSet};
    use std::iter::FromIterator;

    #[test]
    fn test_polytope_intersection() {

        let polys                   =   enumerate_poly_test_set_40_total();

        for mut poly_pair in ( 0 .. 2 ).map(|x| polys.iter() ).multi_cartesian_product() {
            let poly_a              =   poly_pair.swap_remove( 1 );
            let poly_b              =   poly_pair.swap_remove( 0 );              

            // compute intersectin of vertices
            let d0_faces_a: HashSet<Polytope>   =   poly_faces( &poly_a, 0 ).iter().cloned().collect();
            let d0_faces_b: HashSet<Polytope>   =   poly_faces( &poly_b, 0 ).iter().cloned().collect();
            let cap_of_d0_faces: HashSet<Polytope>  =   d0_faces_a.intersection( & d0_faces_b ).cloned().collect();  

            if let  Some( intersection_polytope ) 
                    =
                    polytope_intersection( &poly_a, &poly_b )
            {
                let d0_faces_of_cap: HashSet<Polytope>     
                                                =   poly_faces( & intersection_polytope, 0 )
                                                        .iter().cloned().collect();
                
                assert_eq!( &d0_faces_of_cap, &cap_of_d0_faces )
            }

            else {              
                assert!( &cap_of_d0_faces.is_empty() );
            }



        }

        // let poly_a                  =   Polytope{ 
        //                                     data_c_to_l: vec![0, 0, 1, 1, 2, 2, 3, 4, 5, 5],
        //                                     data_l_to_fmin: vec![0, 0, 1, 1, 2, 3],
        //                                 };

        // let poly_a                  =   Polytope{ 
        //                                     data_c_to_l: vec![0, 0, 1, 1, 2, 2, 3, 4, 5, 5],
        //                                     data_l_to_fmin: vec![0, 0, 1, 1, 2, 3],
        //                                 };
                                        
                                        
        // MAYBE APPLY BRUTE FORCE -- ENUMERATE ALL FACES OF BOTH POLYTOPES, AND CHECK WHETHER THEY INTERSECT

          
    }  
            


}    
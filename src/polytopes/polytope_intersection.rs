
use crate::polytopes::polytope::Polytope;
// use crate::utilities::*;
// use solar::utilities::index::{SuperVec, EndIndex};
// use ordered_float::OrderedFloat;
// use num::rational::Ratio;
// use std::collections::{HashMap};
// use std::hash::Hash;
// use std::iter::FromIterator;
// use std::cmp::Ord;
// use std::iter;

type Fil = usize;




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


pub fn  intersects(
            poly_a:     &Polytope,
            poly_b:     &Polytope
        )
        ->
        bool

{

    let num_cells       =   poly_a.num_cells();

    // place polytopes into a vector
    let o_poly          =   vec![ poly_a, poly_b ];

    // make vectors of length (# cells) to keep track of max and min values

    let a_fmin          =   poly_a.vec_mapping_cell_id_to_min_filt_ordinal();
    let a_fmax          =   poly_a.vec_mapping_cell_id_to_max_filt_ordinal();  
    
    let b_fmin          =   poly_a.vec_mapping_cell_id_to_min_filt_ordinal();
    let b_fmax          =   poly_a.vec_mapping_cell_id_to_max_filt_ordinal();      

    // "o" stands for "omnia" (meaning everything)
    let mut o_fmin      =   [a_fmin, b_fmin];
    let mut o_fmax      =   [a_fmax, b_fmax];    

    loop {

        let mut break_loop     =   true;

        for cell_id in 0 .. num_cells {

            // there is nothing to do if this cell has the same max/min filtration values in both polytopes
            if  o_fmin[ 0 ][ cell_id ]  ==  o_fmin[ 1 ][ cell_id ]
                &&
                o_fmax[ 0 ][ cell_id ]  ==  o_fmax[ 1 ][ cell_id ]
            { continue }

            // find the buffer with the (lexicographically) lower pair (cell_id min val, cell_id max val)
            let i_lower                 =   (0..2).min_by_key(  |x| 
                                                            (
                                                                o_fmin[ *x ][ cell_id ],
                                                                o_fmax[ *x ][ cell_id ],                                                                
                                                            )        
                                                        ).unwrap(); 
            // get the index for the other buffer                                                        
            let i_upper                 =   1 - i_lower;

            // calculate the central value to which all up/downstream cells must now be bound
            let val_lower               =   o_fmax[ i_lower ][ cell_id ].clone();                
            let val_upper               =   o_fmin[ i_upper ][ cell_id ].clone();                
            
            if val_upper != val_lower { return false }
            let val_center              =   val_lower.clone();                

            // bind the up/down stream values
            for upstream_cell_id  in o_poly[ i_lower ].cell_id_to_upstream_cell_ids( cell_id ) {
                o_fmin[ i_lower ][ upstream_cell_id ]  =   val_center.clone();
            }
            for dnstream_cell_id  in o_poly[ i_upper ].cell_id_to_dnstream_cell_ids( cell_id ) {
                o_fmax[ i_upper ][ dnstream_cell_id ]  =   val_center.clone();
            }   
            
            break_loop                 =   true;
        }

        if break_loop {
            return true
        }
    }
}





#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;


    #[test]
    fn test_polytope_intersection() {

        let poly_a                  =   Polytope{ 
                                            data_c_to_l: vec![0, 0, 1, 1, 2, 2, 3, 4, 5, 5],
                                            data_l_to_fmin: vec![0, 0, 1, 1, 2, 3],
                                        };

        let poly_a                  =   Polytope{ 
                                            data_c_to_l: vec![0, 0, 1, 1, 2, 2, 3, 4, 5, 5],
                                            data_l_to_fmin: vec![0, 0, 1, 1, 2, 3],
                                        };
                                        
                                        
        // MAYBE APPLY BRUTE FORCE -- ENUMERATE ALL FACES OF BOTH POLYTOPES, AND CHECK WHETHER THEY INTERSECT

          
    }  
            


}    
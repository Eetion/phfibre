use crate::utilities::*;
use crate::intervals_and_ordinals::*;
use crate::phfibre::*;
use num::rational::Ratio;
use std::collections::{HashSet, HashMap};
use std::iter::FromIterator;
use std::iter;
use solar;


type Cell = usize;
type Coeff = Ratio< i64 >;
type Ring = solar::rings::ring_native::NativeDivisionRing::< Ratio<i64> >;

#[derive(Clone, Debug)]
pub struct ChainComplexRankNullity{
    data_image:     Vec< usize >,
    data_kernel:    Vec< usize >
}

impl ChainComplexRankNullity {
    fn rank_boundaries( &self, dim: &usize ) -> usize { self.data_image[ *dim ].clone() }
    fn rank_cycles(     &self, dim: &usize ) -> usize { self.data_kernel[ *dim ].clone() }    
    fn rank_homology(   &self, dim: &usize ) -> usize { self.rank_cycles(dim) - self.rank_boundaries(dim) }
    fn rank_chains(     &self, dim: &usize ) -> usize { 
        match *dim == 0 {
            true    =>  self.rank_cycles( dim ),
            false   =>  self.rank_cycles(dim) + self.rank_boundaries( &(dim.clone()-1) ) 
        }        
    }    
}

pub fn  chain_cx_rank_nullity( 
            boundary:       Vec< Vec< (Cell, Coeff) > >,
            ring:           Ring,
            dims:           Vec< usize >
        ) 
        -> 
        ChainComplexRankNullity 
{
    // number of degrees where where compute rank
    let mut num_degrees         =   dims.iter().cloned().max().unwrap() + 1;

    // placeholder for boundary ranks
    let mut dims_b  : Vec< _ >  =   Vec::from_iter( 
                                        iter::repeat(0)
                                        .take( num_degrees )  
                                    );

    // placeholder for cycle ranks
    let mut dims_z              =   dims_b.clone();
    
    // copy of array that will be reduced
    let mut reducindus          =   boundary.clone();

    // reduce the copy
    solar::reduce::vec_of_vec::right_reduce(
        &mut reducindus,
        ring
    );

    // compute ranks
    for (col_count, col) in reducindus.iter().enumerate()  {

        let dim                 =   dims[ col_count ];
        match col.is_empty(){
            true    => dims_z[ dim      ] += 1,
            false   => dims_b[ dim -1   ] += 1
        }

        // // if column is nonempty, add 1 to the rank of boundaries one degree down
        // if let Some( last_entry ) = end_val( col ) {
        //     dims_b[ 
        //         dims[ last_entry.0 ] - 1
        //     ]
        //     +=1;
        // // otherwise increase kernel dimension in this chain degree
        // } else {
        //     dims_z[
        //         dims[ col_count ]
        //     ]
        //     += 1
        // }
    }

    ChainComplexRankNullity { data_image: dims_b, data_kernel: dims_z }
}






//  ---------------------------------------------------------------------------
//  TESTS
//  ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test()
    {

        // ring operator
        let ring                =   Ring::new();

        // boundary matrix of the interval
        let boundary            =   vec![ 
                                        vec![           ],
                                        vec![           ],
                                        vec![  (0, Ratio::new(1,1)),  (1, Ratio::new(-1,1))    ]
                                    ];

        //  DEFINE THE CELL DIMENSIONS

        let dims                =   vec![ 0, 0, 1];
        
        let ranks               =   chain_cx_rank_nullity(
                                        boundary, 
                                        ring, 
                                        dims
                                    );


        println!("{:?}", ranks.clone() );

        assert_eq!( 2, ranks.rank_chains( &0 ));
        assert_eq!( 1, ranks.rank_chains( &1 ));
        assert_eq!( 2, ranks.rank_cycles( &0 ));
        assert_eq!( 0, ranks.rank_cycles( &1 ));
        assert_eq!( 1, ranks.rank_boundaries( &0 ));                                
        assert_eq!( 0, ranks.rank_boundaries( &1 ));                                
        assert_eq!( 1, ranks.rank_homology( &0 ));                                
        assert_eq!( 0, ranks.rank_homology( &1 ));                                
        
    }     

}
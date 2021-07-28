use solar::utilities::index::SuperVec;
use crate::intervals_and_ordinals::*;
use num::rational::Ratio;
use std::iter::FromIterator;
use std::iter;
use solar;
use solar::rings::ring::{Semiring, Ring, DivisionRing};
use std::fmt::Debug;

type Cell = usize;
type Coeff = Ratio< i64 >;
// type Ring = solar::rings::ring_native::NativeDivisionRing::< Ratio<i64> >;

#[derive(Clone, Debug)]
pub struct ChainComplexRankNullity{
    data_image:     SuperVec< usize >,
    data_kernel:    SuperVec< usize >
}

impl ChainComplexRankNullity {
    pub fn rank_boundaries( &self, dim: &usize ) -> usize { self.data_image.val( *dim ).clone() }
    pub fn rank_cycles(     &self, dim: &usize ) -> usize { self.data_kernel.val( *dim ).clone() }    
    pub fn rank_homology(   &self, dim: &usize ) -> usize { self.rank_cycles(dim) - self.rank_boundaries(dim) }
    pub fn rank_chains(     &self, dim: &usize ) -> usize { 
        match *dim == 0 {
            true    =>  self.rank_cycles( dim ),
            false   =>  self.rank_cycles(dim) + self.rank_boundaries( &(dim.clone()-1) ) 
        }        
    }    

    /// Super vector specifying rank of boundary space per degree.        
    pub fn rank_boundaries_supervec( &self ) -> SuperVec< usize > { self.data_image.clone() }

    /// Super vector specifying rank of cycle space per degree.    
    pub fn rank_cycles_supervec( &self )     -> SuperVec< usize > { self.data_kernel.clone() }  
    
    /// Super vector specifying rank of homology space per degree.        
    pub fn rank_homology_supervec( &self )   -> SuperVec< usize > { 
        let vec =       Vec::from_iter(  
                            (0..self.data_image.vec.len())
                                .map( |x| self.rank_homology( &x ) )
                        );
        SuperVec{ vec: vec, val: 0 }
    }

    /// Super vector specifying rank of chain space per degree.
    pub fn rank_vec_chains( &self )   -> SuperVec< usize > { 
        let vec =       Vec::from_iter(  
                            (0..self.data_image.vec.len())
                                .map( |x| self.rank_chains( &x ) )
                        );
        SuperVec{ vec: vec, val: 0 }
    }    

    pub fn top_chain_degree( &self ) -> Option< usize > {
        match self.data_image.vec.is_empty(){
            true    =>  None,
            false   =>  Some( self.data_image.vec.len() - 1 ) // why subtract 1?  consider case where top_chain_degree = 0
        }
    }
}


//  ---------------------------------------------------------------------------
//  CHAIN COMPLEX RANKS
//  ---------------------------------------------------------------------------


pub fn  chain_cx_rank_nullity< RingOp, RingElt >( 
            boundary:       Vec< Vec< (Cell, RingElt) > >,
            ring:           RingOp,
            cell_dims:      & Vec< usize >
        ) 
        -> 
        ChainComplexRankNullity 
    where   RingOp:     Clone + Semiring<RingElt> + Ring<RingElt> + DivisionRing<RingElt>,
            RingElt:    Clone + Debug + PartialOrd
{
    // number of degrees where where compute rank
    let num_degrees         =   cell_dims.iter().cloned().max().unwrap() + 1;

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

        let dim                 =   cell_dims[ col_count ];
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

    ChainComplexRankNullity {   
                                data_image:     SuperVec{ vec: dims_b, val: 0 }, 
                                data_kernel:    SuperVec{ vec: dims_z, val: 0 }
                            }
}


//  ---------------------------------------------------------------------------
//  DEGENERATE BAR QUOTA
//  ---------------------------------------------------------------------------

/// Calculates the exact number of degenerate bars in each dimension, given 
/// the barcode and the dimensions of the boundary spaces.
pub fn num_degenerate_bars_per_degree( 
    ranks:      & ChainComplexRankNullity,
    barcode:    & Barcode,
    ) 
    -> 
    SuperVec< usize > 
{

    // compute dimension of boundary space in each degree
    let mut quota       =   ranks.rank_boundaries_supervec();
    
    // subtract number of bars in each degree
    let mut num_bars_fin_per_dim    =   barcode.num_bars_fin_per_dim().vec;    
    for deg in 0 .. num_bars_fin_per_dim.len() {
        quota.vec[ deg ] -= num_bars_fin_per_dim[ deg ]
    }

    quota // this is a super vec
}


//  ---------------------------------------------------------------------------
//  TESTS
//  ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_chain_complex_rank()
    {

        // ring operator
        let ring                =   solar::rings::ring_native::NativeDivisionRing::< Ratio<i64> >::new();

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
                                        & dims
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
        
        assert_eq!( SuperVec{ vec: vec![2, 1], val: 0 }, ranks.rank_vec_chains() );
        assert_eq!( SuperVec{ vec: vec![2, 0], val: 0 }, ranks.rank_cycles_supervec() );
        assert_eq!( SuperVec{ vec: vec![1, 0], val: 0 }, ranks.rank_boundaries_supervec() );
        assert_eq!( SuperVec{ vec: vec![1, 0], val: 0 }, ranks.rank_homology_supervec() );                                
    }     

}
use solar::utilities::indexing_and_bijection::{SuperIndex};
use crate::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use num::rational::Ratio;
use std::iter::FromIterator;
use std::iter;
use solar;
use solar::rings::ring::{Semiring, Ring, DivisionRing};
use std::fmt::Debug;
use std::hash::Hash;

type Cell = usize;
// type Coeff = Ratio< i64 >;
// type Ring = solar::rings::ring_native::NativeDivisionRing::< Ratio<i64> >;



//  ---------------------------------------------------------------------------
//  COMPUTE PIVOT INDEX PAIRS FROM REDUCED MATRIX
//  ---------------------------------------------------------------------------

pub fn  reduced_mat_to_pivot_index_pairs< T >( reduced_mat: &Vec< Vec< (usize,T)> >) 
    -> 
    Vec< (usize, usize) >  
{

    let mut num_pairs       =   0;
    for col in reduced_mat { if col.len() > 0 { num_pairs += 1} }
    let mut pairs           =   Vec::with_capacity( num_pairs );
    for (col_count, col) in reduced_mat.iter().enumerate() { 
        if col.len() > 0 { 
            pairs.push( ( col.last().unwrap().0, col_count)  ) 
        } 
    }
    pairs
}



//  ---------------------------------------------------------------------------
//  STRUCT: CHAIN COMPLEX RANK
//  ---------------------------------------------------------------------------


#[derive(Clone, Debug)]
pub struct ChainComplexRankNullity{
    data_image:     Vec< usize >,
    data_kernel:    Vec< usize >
}

impl ChainComplexRankNullity {
    pub fn rank_boundaries( &self, dim: &usize ) -> usize { self.data_image.sindex( *dim, 0 ).clone() }
    pub fn rank_cycles(     &self, dim: &usize ) -> usize { self.data_kernel.sindex( *dim, 0 ).clone() }    
    pub fn rank_homology(   &self, dim: &usize ) -> usize { self.rank_cycles(dim) - self.rank_boundaries(dim) }
    pub fn rank_chains(     &self, dim: &usize ) -> usize { 
        match *dim == 0 {
            true    =>  self.rank_cycles( dim ),
            false   =>  self.rank_cycles(dim) + self.rank_boundaries( &(dim.clone()-1) ) 
        }        
    }    

    /// Super vector specifying rank of boundary space per degree.        
    pub fn rank_boundaries_vec( &self ) -> Vec< usize > { self.data_image.clone() }

    /// Super vector specifying rank of cycle space per degree.    
    pub fn rank_cycles_vec( &self )     -> Vec< usize > { self.data_kernel.clone() }  
    
    /// Super vector specifying rank of homology space per degree.        
    pub fn rank_homology_vec( &self )   -> Vec< usize > { 
        let vec =       Vec::from_iter(  
                            (0..self.data_image.len())
                                .map( |x| self.rank_homology( &x ) )
                        );
        vec
    }

    /// Super vector specifying rank of chain space per degree.
    pub fn rank_vec_chains( &self )   -> Vec< usize > { 
        let vec =       Vec::from_iter(  
                            (0..self.data_image.len())
                                .map( |x| self.rank_chains( &x ) )
                        );
        vec
    }    

    pub fn top_chain_degree( &self ) -> Option< usize > {
        match self.data_image.is_empty(){
            true    =>  None,
            false   =>  Some( self.data_image.len() - 1 ) // why subtract 1?  consider case where top_chain_degree = 0
        }
    }
}


//  ---------------------------------------------------------------------------
//  CHAIN COMPLEX RANKS -- FROM BOUNDARY MATRIX
//  ---------------------------------------------------------------------------


pub fn  chain_cx_rank_nullity< RingOp, RingElt >( 
            boundary:       & Vec< Vec< (Cell, RingElt) > >,
            cell_degs:      & Vec< usize >,
            ring:           RingOp,            
        ) 
        -> 
        ChainComplexRankNullity 
    where   RingOp:     Clone + Semiring<RingElt> + Ring<RingElt> + DivisionRing<RingElt>,
            RingElt:    Clone + Debug + PartialOrd
{
    // number of degrees where where compute rank
    let num_degrees         =   cell_degs.iter().cloned().max().unwrap() + 1;

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

        let dim                 =   cell_degs[ col_count ];
        match col.is_empty(){
            true    => dims_z[ dim      ] += 1,
            false   => dims_b[ dim -1   ] += 1
        }

    }

    ChainComplexRankNullity {   
                                data_image:     dims_b,
                                data_kernel:    dims_z,
                            }
}

//  ---------------------------------------------------------------------------
//  CHAIN COMPLEX RANKS -- FROM PAIRS
//  ---------------------------------------------------------------------------

pub fn  pairs_dims_to_chx_rank_nullity( 
    pairs:          &   Vec< (usize, usize) >,
    cell_degs:      &   Vec< usize >
) 
-> 
ChainComplexRankNullity 
{
// number of degrees where where compute rank
let num_degrees         =   cell_degs.iter().cloned().max().unwrap() + 1;

// placeholder for boundary ranks
let mut dims_b  : Vec< _ >  =   Vec::from_iter( 
                                iter::repeat(0)
                                .take( num_degrees )  
                            );

// placeholder for cycle ranks
let mut dims_z              =   dims_b.clone();

// cells form cycles until proven otherwise
for cell_id in 0 .. cell_degs.len() { dims_z[ cell_degs[ cell_id ] ] += 1 }

// now adjust for pivot pairs
for (row, col) in pairs  {
    dims_z[ cell_degs[ *col ] ] -=   1;
    dims_b[ cell_degs[ *row ] ] +=   1;
}

ChainComplexRankNullity {   
                        data_image:     dims_b,
                        data_kernel:    dims_z,
                    }
}


//  ---------------------------------------------------------------------------
//  DEGENERATE BAR QUOTA
//  ---------------------------------------------------------------------------

/// Calculates the exact number of degenerate bars in each dimension, given 
/// the barcode and the dimensions of the boundary spaces.
pub fn num_degenerate_bars_per_degree< FilRaw >( 
    ranks:      & ChainComplexRankNullity,
    barcode:    & Barcode< FilRaw >,
    ) 
    -> 
    Option< Vec< usize > >
    where FilRaw: Ord + Clone + Hash
{

    // compute dimension of boundary space in each degree
    let mut quota       =   ranks.rank_boundaries_vec();
    
    // subtract number of bars in each degree
    let num_bars_fin_per_dim    =   barcode.num_bars_fin_per_dim();    
    for deg in 0 .. num_bars_fin_per_dim.len() {
        if num_bars_fin_per_dim[ deg ] >  quota[ deg ] { 
            println!("\nBARCODE IS INCOMPATIBLE WITH THIS BOUNDARY MATRIX\n");
            return None 
        } else {
            quota[ deg ] -= num_bars_fin_per_dim[ deg ]
        }
    }

    Some( quota ) // this is a super vec
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

        //  define the cell dimensions
        let dims                =   vec![ 0, 0, 1];
        
        //  compute ranks        
        let ranks               =   chain_cx_rank_nullity(
                                        & boundary, 
                                        & dims,
                                        ring.clone(),                                         
                                    );

        //  verify
        println!("{:?}", ranks.clone() );

        assert_eq!( 2, ranks.rank_chains( &0 ));
        assert_eq!( 1, ranks.rank_chains( &1 ));
        assert_eq!( 2, ranks.rank_cycles( &0 ));
        assert_eq!( 0, ranks.rank_cycles( &1 ));
        assert_eq!( 1, ranks.rank_boundaries( &0 ));                                
        assert_eq!( 0, ranks.rank_boundaries( &1 ));                                
        assert_eq!( 1, ranks.rank_homology( &0 ));                                
        assert_eq!( 0, ranks.rank_homology( &1 ));         
        
        assert_eq!( vec![2, 1], ranks.rank_vec_chains() );
        assert_eq!( vec![2, 0], ranks.rank_cycles_vec() );
        assert_eq!( vec![1, 0], ranks.rank_boundaries_vec() );
        assert_eq!( vec![1, 0], ranks.rank_homology_vec() );                                
     

        // ALTERNATE APPROACH !!! 
        // ---------------------- 

        // copy boundary
        let mut reducindus          =   boundary.clone();

        // reduce the copy
        solar::reduce::vec_of_vec::right_reduce(
            &mut reducindus,
            ring.clone(),
        );    

        // extract the pairs + compute ranks
        let pairs                   =   reduced_mat_to_pivot_index_pairs( &reducindus );
        let ranks                   =   pairs_dims_to_chx_rank_nullity( 
                                            & pairs,
                                            & dims
                                        );  
        //  verify
        println!("{:?}", ranks.clone() );

        assert_eq!( 2, ranks.rank_chains( &0 ));
        assert_eq!( 1, ranks.rank_chains( &1 ));
        assert_eq!( 2, ranks.rank_cycles( &0 ));
        assert_eq!( 0, ranks.rank_cycles( &1 ));
        assert_eq!( 1, ranks.rank_boundaries( &0 ));                                
        assert_eq!( 0, ranks.rank_boundaries( &1 ));                                
        assert_eq!( 1, ranks.rank_homology( &0 ));                                
        assert_eq!( 0, ranks.rank_homology( &1 ));         
        
        assert_eq!( vec![2, 1], ranks.rank_vec_chains() );
        assert_eq!( vec![2, 0], ranks.rank_cycles_vec() );
        assert_eq!( vec![1, 0], ranks.rank_boundaries_vec() );
        assert_eq!( vec![1, 0], ranks.rank_homology_vec() );   
    }
}
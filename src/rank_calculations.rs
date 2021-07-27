
use crate::intervals_and_ordinals::*;
use num::rational::Ratio;
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
    pub fn rank_boundaries( &self, dim: &usize ) -> usize { self.data_image[ *dim ].clone() }
    pub fn rank_cycles(     &self, dim: &usize ) -> usize { self.data_kernel[ *dim ].clone() }    
    pub fn rank_homology(   &self, dim: &usize ) -> usize { self.rank_cycles(dim) - self.rank_boundaries(dim) }
    pub fn rank_chains(     &self, dim: &usize ) -> usize { 
        match *dim == 0 {
            true    =>  self.rank_cycles( dim ),
            false   =>  self.rank_cycles(dim) + self.rank_boundaries( &(dim.clone()-1) ) 
        }        
    }    
    pub fn rank_vec_boundaries( &self ) -> Vec< usize > { self.data_image.clone() }
    pub fn rank_vec_cycles( &self )     -> Vec< usize > { self.data_kernel.clone() }    
    pub fn rank_vec_homology( &self )   -> Vec< usize > { 
        Vec::from_iter(  
            (0..self.data_image.len())
                .map( |x| self.rank_homology( &x ) )
        )  
    }
    pub fn rank_vec_chains( &self )   -> Vec< usize > { 
        Vec::from_iter(  
            (0..self.data_image.len())
                .map( |x| self.rank_chains( &x ) )
        )  
    }    

    pub fn top_chain_degree( &self ) -> Option< usize > {
        match self.data_image.is_empty(){
            true    =>  None,
            false   =>  Some( self.data_image.len() )
        }
    }
}


//  ---------------------------------------------------------------------------
//  CHAIN COMPLEX RANKS
//  ---------------------------------------------------------------------------


pub fn  chain_cx_rank_nullity( 
            boundary:       Vec< Vec< (Cell, Coeff) > >,
            ring:           Ring,
            cell_dims:      & Vec< usize >
        ) 
        -> 
        ChainComplexRankNullity 
{
    // number of degrees where where compute rank
    let mut num_degrees         =   cell_dims.iter().cloned().max().unwrap() + 1;

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

    ChainComplexRankNullity { data_image: dims_b, data_kernel: dims_z }
}


//  ---------------------------------------------------------------------------
//  DEGENERATE BAR QUOTA
//  ---------------------------------------------------------------------------

/// Calculates the exact number of degenerate bars in each dimension, given 
/// the barcode and the dimensions of the boundary spaces.
pub fn degenerate_bar_quota( 
    ranks:      & ChainComplexRankNullity,
    barcode:    & Barcode,
    ) 
    -> 
    Vec< usize > 
{
    let top_chain_deg           =   ranks.top_chain_degree().unwrap();

    println!("TOP DIM = {:?}", & top_chain_deg );

    let num_bars_fin_per_dim    =   barcode.num_bars_fin_per_dim();
    let mut quota               =   Vec::with_capacity( top_chain_deg + 1 ); // why +1?  consider the case where top_chain_deg = 0

    for dim in 0..top_chain_deg+1{
        quota
            .push(
                ranks.rank_boundaries( &dim ) - num_bars_fin_per_dim[ dim ]
            )
    }
    quota
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
        
        assert_eq!( vec![2, 1], ranks.rank_vec_chains() );
        assert_eq!( vec![2, 0], ranks.rank_vec_cycles() );
        assert_eq!( vec![1, 0], ranks.rank_vec_boundaries() );
        assert_eq!( vec![1, 0], ranks.rank_vec_homology() );                                
    }     

}
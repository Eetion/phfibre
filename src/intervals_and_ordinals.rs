use crate::utilities::*;
use solar::utilities::indexing_and_bijection::{EndIndex};
use solar::utilities::sequences_and_ordinals::{BiMapSequential, ordinate_unique_vals};
use ordered_float::OrderedFloat;
use num::rational::Ratio;
use std::collections::{HashMap};
use std::hash::Hash;
use std::iter::FromIterator;
use std::cmp::Ord;
use std::iter;
use serde::{Deserialize, Serialize};

// type Cell = usize;
// type Coeff = Ratio< i64 >;
type Fil = usize;
// type FilRaw = OrderedFloat<f64>; // reason for this choice: f64 does not implement the hash trait (cricital for reparametrizing)




//  ---------------------------------------------------------------------------
//  CONVERSION
//  ---------------------------------------------------------------------------


// Produces a vector of ordered floats from a vector of floats
pub fn to_ordered_float( v: & Vec< f64 > ) -> Vec< OrderedFloat< f64> > { v.iter().map(|x| OrderedFloat(x.clone()) ).collect() }





//  ---------------------------------------------------------------------------  
//  BARS + BARCODES
//  ---------------------------------------------------------------------------  


#[derive(Clone, Debug, PartialEq, PartialOrd, Eq, Ord, Serialize, Deserialize)]
pub struct BarFinite 
{ 
    pub dim:        usize, 
    pub birth:      Fil, 
    pub death:      Fil 
}
impl BarFinite{ 
    pub fn dim( &self ) -> usize { self.dim.clone()   }
    pub fn birth( &self ) -> Fil { self.birth.clone() }
    pub fn death( &self ) -> Fil { self.death.clone() }
}

#[derive(Clone, Debug, PartialEq, PartialOrd, Eq, Ord, Serialize, Deserialize)]
pub struct BarInfinite { 
    pub dim:        usize, 
    pub birth:      Fil 
}
impl <'de>
    BarInfinite{ 
    pub fn dim( &self ) -> usize { self.dim.clone()   }
    pub fn birth( &self ) -> Fil { self.birth.clone() }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Barcode< FilRaw > 
    where   FilRaw: Clone + Ord + Hash
{
    pub inf:            Vec< BarInfinite >,     // infinite bars (birth ordinals)
    pub fin:            Vec< BarFinite >,       // finite bars (birth/death ordinals)
    pub ordinal:        BiMapSequential< FilRaw >,  // struct converting endpoints to ordinals and vice versa
}


impl < FilRaw > Barcode< FilRaw >
    where FilRaw: Clone + Ord + Hash
{

    /// The maximum dimension of any bar in the barcode (None if the barcode is empty)
    pub fn top_bar_degree( &self ) -> Option< usize > { 
        self.fin.iter().map(|x| x.dim.clone())
            .chain(
                self.inf.iter().map(|x| x.dim.clone())
            )
            .max()    
    }

    /// The maximum filtration (ordinal) of any bar in the barcode (None if the barcode is empty)
    pub fn top_fil( &self ) -> Option< Fil > { 
        self.fin.iter().map(|x| x.death.clone())
            .chain(
                self.inf.iter().map(|x| x.birth.clone())
            )
            .max()    
    }

    /// Returns a vector `v` such that `v[i] = #{bars of degree i}`
    /// 
    /// The length of `v` is (top_chain_degree + 1), since it's important to be able to 
    pub fn num_bars_fin_per_dim( &self ) -> Vec< usize > { 
  
        if let Some( d ) = self.top_bar_degree() {
            let mut counts      =   Vec::with_capacity( d );
            for _ in 0..d+1 { counts.push(0) }
            for bar in & self.fin { counts[ bar.dim() ] += 1 }  
            return counts      
        } else {
            return Vec::with_capacity(0)
        };
        
        // mut counts          =   Vec::with_capacity( 0 );

        // // if there are any bars in the barcode, count them
        // if let Some( d ) = self.top_bar_degree() {
        //     let mut counts      =   Vec::with_capacity( d );
        //     for _ in 0..d+1 { counts.push(0) }
        //     for bar in & self.fin { counts[ bar.dim() ] += 1 }        
        // }

    }    

    /// Number of finite bars
    pub fn num_bars_fin( &self ) -> usize { self.fin.len() }
    
    /// Number of infinite bars
    pub fn num_bars_inf( &self ) -> usize { self.inf.len() }
    
    /// Get the ith finite bar
    pub fn bar_fin( &self, i: usize ) -> BarFinite { self.fin[ i ].clone() }
    
    /// Get the ith infinite bar
    pub fn bar_inf( &self, i: usize ) -> BarInfinite { self.inf[ i ].clone() }

    /// Create a new barcode from raw parts
    //  COMMENTS:
    //  This function doesn't have to be especially efficient since we only run it once.
    pub fn new( 
        inf_dim:    Vec< usize >,        
        inf_brn:    Vec< FilRaw >,
        fin_dim:    Vec< usize >,        
        fin_brn:    Vec< FilRaw >,
        fin_die:    Vec< FilRaw >,
        ) 
        -> 
        Barcode< FilRaw > 
    {
        // place all endpoints into sequence
        let mut endpoints_raw_unordered   =   Vec::new();
        endpoints_raw_unordered.append( &mut inf_brn.clone() );
        endpoints_raw_unordered.append( &mut fin_brn.clone() );
        endpoints_raw_unordered.append( &mut fin_die.clone() );

        // extract ordinal data about the set of endpoints
        let raw_endpoint_ordinal_data   =   ordinate_unique_vals( &endpoints_raw_unordered );

        // initialize barcode object
        let mut barcode     =   Barcode {
                                inf:        Vec::new(),     // infinite bars (birth ordinals)
                                fin:        Vec::new(),       // finite bars (birth/death ordinals)
                                ordinal:    raw_endpoint_ordinal_data
                            };
                            
        
        for bar_count in 0..fin_brn.len() {
            let bar     =   BarFinite { 
                                dim:        fin_dim[ bar_count ],
                                birth:      barcode.ordinal.ord(
                                                & fin_brn[ bar_count ]
                                            )
                                            .unwrap(), 
                                death:      barcode.ordinal.ord(
                                                & fin_die[ bar_count ]
                                            )
                                            .unwrap(), 
                            };        
            barcode.fin.push( bar );
        }
        
        for bar_count in 0..inf_brn.len() {
            let bar     =   BarInfinite { 
                                dim:        inf_dim[ bar_count ],
                                birth:      barcode.ordinal.ord(
                                    & inf_brn[ bar_count ]
                                )
                                .unwrap(), 
                            };        
            barcode.inf.push( bar );
        }        

        return barcode
    }

    pub fn sort( &mut self ) 
        where FilRaw: Ord 
    {
        self.fin.sort();
        self.inf.sort();
    }
}

//  ---------------------------------------------------------------------------  
//  BARCODE INVERSE
//  ---------------------------------------------------------------------------  


/// Encodes a map sending (ordinal) endpoints to sets of bar id numbers.
#[derive(Clone, Debug)]
pub struct BarcodeInverse {
    pub inf_brn:    Vec< Vec< usize > >,
    pub fin_brn:    Vec< Vec< usize > >,
    pub fin_die:    Vec< Vec< usize > >
}

impl BarcodeInverse{

    /// Returns an object that maps an (ordinal) endpoints back to set of bar id's.
    pub fn from_barcode < FilRaw > (
        barcode: & Barcode< FilRaw >
        ) ->
        BarcodeInverse
        where FilRaw: Ord + Clone + Hash
    {
        let mut endpoint_to_barids = BarcodeInverse {
            inf_brn:    Vec::with_capacity(0),
            fin_brn:    Vec::with_capacity(0),
            fin_die:    Vec::with_capacity(0),
        };    

        if let Some( top_fil )  =  barcode.top_fil() {
            // why have a +1 here?  consider the case where top_fil = 0 (as in, there is exactly one level, equal to 0)
            for _ in 0..(top_fil+1) {
                endpoint_to_barids.inf_brn.push(  Vec::new()  );
                endpoint_to_barids.fin_brn.push(  Vec::new()  );
                endpoint_to_barids.fin_die.push(  Vec::new()  );                        
            }
        } 

        // fill bins with infinite bars
        for (bar_count, bar) in barcode.inf.iter().enumerate() {
            endpoint_to_barids.inf_brn[ bar.birth ].push( bar_count );
        }

        // fill bins with finite bars
        for (bar_count, bar) in barcode.fin.iter().enumerate() {
            endpoint_to_barids.fin_brn[ bar.birth ].push( bar_count );
            endpoint_to_barids.fin_die[ bar.death ].push( bar_count );        
        }

        //  remove excess capacity
        for vec in endpoint_to_barids.inf_brn.iter_mut() { vec.shrink_to_fit() };
        for vec in endpoint_to_barids.fin_brn.iter_mut() { vec.shrink_to_fit() };    
        for vec in endpoint_to_barids.fin_die.iter_mut() { vec.shrink_to_fit() };        

        return endpoint_to_barids
    }
    
    
}


//  ---------------------------------------------------------------------------  
//  COMPUTE BARCODE FROM PAIRS
//  ---------------------------------------------------------------------------  


pub fn  pairs_dims_births_to_barcode< FilRaw >(
    pairs:  & Vec< (usize, usize) >,
    dims:   & Vec< usize >,
    births: & Vec< FilRaw >
) 
-> 
Barcode< FilRaw >
where FilRaw: Ord + Clone + Hash
{

    let ordinal             =   ordinate_unique_vals( &births );

    let num_cells           =   dims.len();
    let mut num_bars_fin    =   0;
    let mut num_bars_inf    =   num_cells;
    let mut essential       =   Vec::from_iter( iter::repeat(true).take(num_cells) );

    // count
    for pair in pairs {
    essential[ pair.0 ]     =   false;
    essential[ pair.1 ]     =   false;
    num_bars_inf            -=  2;
    if births[ pair.0 ] != births[ pair.1 ] { num_bars_fin +=1 }
    }

    let mut inf             =   Vec::with_capacity( num_bars_inf );
    let mut fin             =   Vec::with_capacity( num_bars_fin );

    // push finite bars
    for pair in pairs {
    if births[ pair.0 ] != births[ pair.1 ] { 
        fin.push(
            BarFinite{ 
                dim:    dims[ pair.0 ], 
                birth:  ordinal.ord( &births[ pair.0 ] ).unwrap(),
                death:  ordinal.ord( &births[ pair.1 ] ).unwrap(),
            }
        )
    }
    }

    // push infinite bars
    for cell_id in 0 .. num_cells {
    if essential[ cell_id ] {
        inf.push(
            BarInfinite{
                dim: dims[ cell_id ],
                birth: ordinal.ord( &births[ cell_id ] ).unwrap(),
            }
        )
    }
    }

    Barcode{
    inf:        inf,
    fin:        fin,
    ordinal:    ordinal
}

}     


//  ---------------------------------------------------------------------------  
//  LEVEL SET SIZES
//  ---------------------------------------------------------------------------  


/// A struct that calculates the sizes of a sequence of level sets.
#[derive(Clone, Debug)]
pub struct LevelSetSizes{ pub pointers: Vec< usize > }

impl LevelSetSizes{
    
    /// Size of the kth level set.
    pub fn size( &self, set_index: usize ) -> usize {  
        if set_index == 0 { self.pointers[0] }
        else { self.pointers[ set_index ] - self.pointers[ set_index - 1 ] }
    }

    /// Size of the last level set.
    pub fn size_last( &self ) -> Option< usize > { 
        match self.pointers.is_empty() {
            true    => None,
            false   => Some (   self.size( 
                                    self.pointers.len() - 1
                                ) 
                            )
        }
    }
    
    /// Size of all level sets combined.
    pub fn num_cells_total( &self ) -> usize {
        match self.pointers.last() { Some( val ) => val.clone(), None => 0  }
    }
    
    /// Add a size value for a new (empty) level set at the end.
    pub fn postpend_empty_set( & mut self ) { 
        match self.pointers.last() { 
            Some( val ) => self.pointers.push( val.clone() ) , 
            None => self.pointers.push(0)
        }        
    }

    /// Add one to the size of the last level set.
    pub fn grow_last_set( & mut self ) {
        match self.pointers.end_index() {
            Some( i ) => self.pointers[ i ] += 1 ,
            None => { panic!("There is no set to grow") }
        }
    }

    pub fn last_level_set_ordinal( &self ) -> Option< usize > {
        match self.pointers.is_empty() {
            true    =>  None,
            false   =>  Some( self.pointers.len() -1 )
        }
    }
}





#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
  

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
                                    ordinal: BiMapSequential{ 
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



}    
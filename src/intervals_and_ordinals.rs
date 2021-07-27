use crate::utilities::*;
use ordered_float::OrderedFloat;
use num::rational::Ratio;
use std::collections::{HashMap};
use std::hash::Hash;
use std::iter::FromIterator;

type Cell = usize;
type Coeff = Ratio< i64 >;
type Fil = usize;
type FilRaw = OrderedFloat<f64>; // reason for this choice: f64 does not implement the hash trait (cricital for reparametrizing)




//  ---------------------------------------------------------------------------
//  CONVERSION
//  ---------------------------------------------------------------------------


// Produces a vector of ordered floats from a vector of floats
pub fn to_ordered_float( v: & Vec< f64 > ) -> Vec< OrderedFloat< f64> > { v.iter().map(|x| OrderedFloat(x.clone()) ).collect() }




//  ---------------------------------------------------------------------------
//  PRIMITIVE ORDINALS
//  ---------------------------------------------------------------------------


#[derive(Clone, Debug, PartialEq)]
pub struct OrdinalData < T : Eq + PartialOrd + PartialEq + Hash > {
    ord_to_val:  Vec< T >,
    val_to_ord:  HashMap< T, usize >
}

impl    < T >
        OrdinalData
        < T >
        where T : std::cmp::Eq + PartialOrd + PartialEq + Hash + Clone
{
    /// The ordinal of the raw filtration value
    pub fn ord( &self, a: &T ) -> Option< usize > { 
        self.val_to_ord.get( a ).map(|x| x.clone()) 
    }
    /// The raw filtration value of the ordinal
    pub fn val( &self, a: usize ) -> Option< T > { 
        if a < self.ord_to_val.len() { Some( self.ord_to_val[ a ].clone() ) } else { None }
    }    
}


/// Get the ordinal data for the range of values taken by a vector of ordered floats
pub fn ordinate( v: & Vec< OrderedFloat<f64> > ) -> OrdinalData< OrderedFloat< f64> > {
    let mut a       =   v.clone();
    let mut b       =   HashMap::new();
    a.sort();       // sort entries
    a.dedup();      // remove duplicates

    for (i, t) in a.iter().enumerate() {
        b.insert( t.clone(), i.clone() );
    }

    OrdinalData { ord_to_val: a, val_to_ord: b }
}





//  ---------------------------------------------------------------------------  
//  BARCODES
//  ---------------------------------------------------------------------------  


#[derive(Clone, Debug, PartialEq)]
pub struct BarFinite { 
    dim:        usize, 
    birth:      Fil, 
    death:      Fil 
}
impl BarFinite{ 
    pub fn dim( &self ) -> usize { self.dim.clone()   }
    pub fn birth( &self ) -> Fil { self.birth.clone() }
    pub fn death( &self ) -> Fil { self.death.clone() }
}

#[derive(Clone, Debug, PartialEq)]
pub struct BarInfinite { 
    dim:        usize, 
    birth:      Fil 
}
impl BarInfinite{ 
    pub fn dim( &self ) -> usize { self.dim.clone()   }
    pub fn birth( &self ) -> Fil { self.birth.clone() }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Barcode {
    pub inf:        Vec< BarInfinite >,     // infinite bars (birth ordinals)
    pub fin:        Vec< BarFinite >,       // finite bars (birth/death ordinals)
    pub ordinal:    OrdinalData< FilRaw >,  // struct converting endpoints to ordinals and vice versa
}

impl Barcode{

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
    pub fn num_bars_fin_per_dim( &self ) -> Vec< usize > { 
        // create vector filled with 0's
        let mut counts      =   Vec::from_iter(
                                    std::iter::repeat(0)
                                    .take(
                                        self.top_bar_degree().unwrap() + 1  // why +1?  consider the case where top_bar_degree = 0
                                    ) 
                                ); 
        // with_capacity( self.top_bar_degree().unwrap() +1 ); 
        for bar in & self.fin { counts[ bar.dim() ] += 1 }
        counts
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
        Barcode 
    {
        // place all endpoints into sequence
        let mut endpoints_raw_unordered   =   Vec::new();
        endpoints_raw_unordered.append( &mut inf_brn.clone() );
        endpoints_raw_unordered.append( &mut fin_brn.clone() );
        endpoints_raw_unordered.append( &mut fin_die.clone() );

        // extract ordinal data about the set of endpoints
        let raw_endpoint_ordinal_data   =   ordinate( &endpoints_raw_unordered );

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
}

/// Encodes a map sending (ordinal) endpoints to sets of bar id numbers.
#[derive(Clone, Debug)]
pub struct BarcodeInverse {
    pub inf_brn:    Vec< Vec< usize > >,
    pub fin_brn:    Vec< Vec< usize > >,
    pub fin_die:    Vec< Vec< usize > >
}

impl BarcodeInverse{

    /// Returns an object that maps an (ordinal) endpoints back to set of bar id's.
    pub fn  from_barcode(
        barcode: & Barcode
        ) ->
        BarcodeInverse

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
        match end_val( & self.pointers ) { Some( val ) => val, None => 0  }
    }
    
    /// Add a size value for a new (empty) level set at the end.
    pub fn postpend_empty_set( & mut self ) { 
        match end_val( & self.pointers ) { 
            Some( val ) => self.pointers.push( val ) , 
            None => self.pointers.push(0)
        }        
    }

    /// Add one to the size of the last level set.
    pub fn grow_last_set( & mut self ) {
        match end_index( & self.pointers ) {
            Some( i ) => self.pointers[ i ] += 1 ,
            None => { panic!("There is no set to grow") }
        }
    }
}



//  ---------------------------------------------------------------------------  
//  POLYTOPES (ORDINAL)
//  ---------------------------------------------------------------------------  

/// Represents an (ordinal) polytope, or, equivalently, a family of upper/lower
/// bounds on variable values.
#[derive(Clone, Debug, PartialEq)]
pub struct Polytope {
    pub data_l_to_fmin:       Vec< Fil >,         // function (level set ordinal) -> min possible filtration value
    pub data_c_to_l:          Vec< usize >,       // function (cell id) -> level set ordinal
}

impl Polytope {

    /// Number of cells in the underlying chain complex.
    pub fn  num_cells( &self ) -> usize { self.data_c_to_l.len() }

    /// Number of level sets of the filter function.
    pub fn  num_lev_sets( &self ) -> usize { self.data_l_to_fmin.len() }    

    /// Rightmost finite filtration value
    pub fn  max_filtration_value( &self ) -> Option<Fil> { self.data_l_to_fmin.iter().cloned().last() }

    /// Given an ordinal, determine whether the corresponding level set is critical.
    pub fn  lev_set_ord_to_is_critical( &self, level_set_ordinal: usize ) -> Option<bool> { 
        if level_set_ordinal >= self.num_lev_sets() { return None }
        else if level_set_ordinal == 0 { return Some( true ) }
        else if self.data_l_to_fmin[ level_set_ordinal ] 
                == 
                self.data_l_to_fmin[ level_set_ordinal -1 ] 
        {
            return Some( false )
        }
        else { return Some( true )}
    }

    /// Determine whether the last level set in the polytope is critical.
    pub fn lev_set_last_to_is_critical( &self ) -> Option<bool> {  
        match self.num_lev_sets() == 0 { 
            true => None, 
            false => self.lev_set_ord_to_is_critical( self.num_lev_sets() )
        }
    }

    /// Ensure that the last level set in the polytope is critical 
    /// 
    /// Concretely, this means ensuring that one of the following two conditions hold: 
    /// (1) the level set in question is the first level set,
    /// (2) the level set in question takes a *minimum* filtration value strictly greater 
    ///     than all other level sets. 
    /// This function throws an error if the polytope contains no level sets.
    pub fn ensure_last_lev_set_critical( &mut self ) {
        if let Some( is_crit ) = self.lev_set_last_to_is_critical() {
            if ! is_crit { 
                let i   =   self.num_lev_sets()-1;
                self.data_l_to_fmin[ i ] += 1; // increasing the min filtration value will make the level set critical            
            }
        }
        else { panic!("there are no level sets, so cannot make last critical") } // happens if there *is* no last level set to make critical      
    }

    /// The minimum (ordinal) filtration value associated with the kth level set.
    pub fn  lev_set_ord_to_fmin( &self, level_set_ordinal: usize ) -> Option<Fil> { 
        if level_set_ordinal >= self.num_lev_sets() { return None }
        else {
            return Some ( self.data_l_to_fmin[ level_set_ordinal ].clone() )            
        }
    }

    // NB!!!!  If we assume that barcode endpoints are 0,1,..,N, then we can simplify
    /// The maximum (ordinal) filtration value associated with the kth level set.
    // computation of max filtration values considerably
    pub fn  lev_set_ord_to_fmax( &self, level_set_ordinal: usize ) -> Option<Fil> { 

        // Posit: the ordinal is within legal range
        if let Some( is_critical ) = self.lev_set_ord_to_is_critical( level_set_ordinal.clone() )
        {
            // Posit: the level set ordinal IS critical.  
            // Therefore: max value equals min value.
            if is_critical { return self.lev_set_ord_to_fmin( level_set_ordinal ) }
            // Posit: the level set ordinal IS NOT not critical                
            // Therefore: max value equals min value + 1
            else { return self.lev_set_ord_to_fmin( level_set_ordinal ).map(|x| x+1 ) }


            //  OLD STUFF PROBABLY OK TO DELETE
            //     let fmin                    =   self.lev_set_ord_to_fmin( 
            //                                             level_set_ordinal.clone() 
            //                                         )
            //                                         .unwrap();
            //     if let Some( next_fil_val ) =   & self.data_l_to_fmin[
            //                                         level_set_ordinal.clone()..self.num_lev_sets()
            //                                         ]
            //                                             .iter()
            //                                             .map(|x| self.lev_set_ord_to_fmax(x.clone()).unwrap() )
            //                                             .filter(|&x| x > fmin)
            //                                             .next()
            //     {
            //         return Some( next_fil_val.clone() )
            //     } else {
            //         return Some( self.num_cells() )
            //     }
                                                
            // }
        } 
        else { return None }
    }

    /// Ordinal of this cell's level set.
    pub fn  cell_id_to_lev_set_ord( &self, cell_id: usize ) -> Option<usize> { 
        if cell_id >= self.num_cells() { return None } 
        else {
            return Some( self.data_c_to_l[ cell_id ].clone()  )
        } 
    }

    /// Mainimum ordinal filtration value allowed for this cell.
    pub fn  cell_id_to_fmin( &self, cell_id: usize ) -> Option<Fil> { 
        match self.cell_id_to_lev_set_ord( cell_id ) {
            Some( i )   =>  self.lev_set_ord_to_fmin( i ),
            None        =>  None
        }
    }

    /// Maximum ordinal filtration value allowed for this cell.
    pub fn  cell_id_to_fmax( &self, cell_id: usize ) -> Option<Fil> { 
        match self.cell_id_to_lev_set_ord( cell_id ) {
            Some( i )   =>  self.lev_set_ord_to_fmax( i ),
            None        =>  None
        }
    }    

    /// Post-pend a new (empty) level set.
    pub fn  push_new_lev_set( &mut self )   {
        match self.max_filtration_value() {
            // if the vector of level set filtration values is nonempty, add a repeat of the last element
            Some( fil_max ) =>  self.data_l_to_fmin.push( fil_max ),
            // otherwise the vector is empty, and we should append a 0
            None            =>  self.data_l_to_fmin.push( 0 )
        }
    }

    /// A vector v such that v[ i ] = (the level set ordinal of cell i)
    pub fn  vec_mapping_cell_id_to_filtration_ordinal( &self ) -> Vec< usize > {
        self.data_c_to_l.clone()
    }

}


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

}    
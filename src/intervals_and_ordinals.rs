use crate::utilities::*;
use ordered_float::OrderedFloat;
use num::rational::Ratio;
use std::collections::{HashMap};
use std::hash::Hash;


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


#[derive(Clone, Debug)]
pub struct OrdinalData < T : PartialOrd + PartialEq + Hash > {
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


#[derive(Clone, Debug)]
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

#[derive(Clone, Debug)]
pub struct BarInfinite { 
    dim:        usize, 
    birth:      Fil 
}
impl BarInfinite{ 
    pub fn dim( &self ) -> usize { self.dim.clone()   }
    pub fn birth( &self ) -> Fil { self.birth.clone() }
}

#[derive(Clone, Debug)]
pub struct Barcode {
    pub inf:        Vec< BarInfinite >,     // infinite bars (birth ordinals)
    pub fin:        Vec< BarFinite >,       // finite bars (birth/death ordinals)
    pub ordinal:    OrdinalData< FilRaw >,  // struct converting endpoints to ordinals and vice versa
}

impl Barcode{
    /// The maximum dimension of any bar in the barcode (None if the barcode is empty)
    pub fn top_dim( &self ) -> Option< usize > { 
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
    //  COMMENTS:
    //  This function doesn't have to be especially efficient since we only run it once.
    /// Create a new barcode from raw parts
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
                                    & fin_brn[ bar_count ]
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
pub struct MapEndpoint2BarIDs {
    pub inf_brn:    Vec< Vec< usize > >,
    pub fin_brn:    Vec< Vec< usize > >,
    pub fin_die:    Vec< Vec< usize > >
}

impl MapEndpoint2BarIDs{

    /// Returns an object that maps an (ordinal) endpoints back to set of bar id's.
    pub fn  from_barcode(
        barcode: Barcode
        ) ->
        MapEndpoint2BarIDs

    {
        let mut endpoint__barids = MapEndpoint2BarIDs {
            inf_brn:    Vec::with_capacity(0),
            fin_brn:    Vec::with_capacity(0),
            fin_die:    Vec::with_capacity(0),
        };    

        if let Some( top_fil )  =  barcode.top_fil() {
            for _ in 0..top_fil {
                endpoint__barids.inf_brn.push(  Vec::new()  );
                endpoint__barids.fin_brn.push(  Vec::new()  );
                endpoint__barids.fin_die.push(  Vec::new()  );                        
            }
        } 

        // fill bins with infinite bars
        for (bar_count, bar) in barcode.inf.iter().enumerate() {
            endpoint__barids.inf_brn[ bar.birth ].push( bar_count );
        }

        // fill bins with finite bars
        for (bar_count, bar) in barcode.fin.iter().enumerate() {
            endpoint__barids.fin_brn[ bar.birth ].push( bar_count );
            endpoint__barids.fin_die[ bar.death ].push( bar_count );        
        }

        //  remove excess capacity
        for vec in endpoint__barids.inf_brn.iter_mut() { vec.shrink_to_fit() };
        for vec in endpoint__barids.fin_brn.iter_mut() { vec.shrink_to_fit() };    
        for vec in endpoint__barids.fin_die.iter_mut() { vec.shrink_to_fit() };        

        return endpoint__barids
    }
    
}




//  ---------------------------------------------------------------------------  
//  LEVEL SET SIZES
//  ---------------------------------------------------------------------------  


/// A struct that calculates the sizes of a sequence of level sets.
#[derive(Clone, Debug)]
pub struct LevelSetSizes{ pointers: Vec< usize > }

impl LevelSetSizes{
    
    /// Size of the kth level set.
    pub fn size( &self, set_index: usize ) -> usize {  
        if set_index == 0 { self.pointers[0] }
        else { self.pointers[ set_index ] - self.pointers[ set_index - 1 ] }
    }

    /// Size of the last level set.
    pub fn size_last( &self ) -> usize { self.size( self.pointers.len() - 1) }
    
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
#[derive(Clone, Debug)]
pub struct Polytope {
    data_l__fmin:       Vec< Fil >,         // function (level set ordinal) -> min possible filtration value
    data_c__l:          Vec< usize >,       // function (cell id) -> level set ordinal
}

impl Polytope {

    /// Number of cells in the underlying chain complex.
    fn  num_cells( &self ) -> usize { self.data_c__l.len() }

    /// Number of level sets of the filter function.
    fn  num_lev_sets( &self ) -> usize { self.data_l__fmin.len() }    

    /// Rightmost finite filtration value
    fn  max_filtration_value( &self ) -> Option<Fil> { self.data_l__fmin.iter().cloned().last() }

    fn  cell_id__lev_set_ord( &self, cell_id: usize ) -> Option<usize> { 
        if cell_id >= self.num_cells() { return None } 
        else {
            return Some( self.data_c__l[ cell_id ].clone()  )
        } 
    }

    fn  lev_set_ord__is_critical( &self, level_set_ordinal: usize ) -> Option<bool> { 
        if level_set_ordinal >= self.num_lev_sets() { return None }
        else if level_set_ordinal == 0 { return Some( true ) }
        else if self.data_l__fmin[ level_set_ordinal ] 
                == 
                self.data_l__fmin[ level_set_ordinal -1 ] 
        {
            return Some( false )
        }
        else { return Some( true )}
    }

    pub fn lev_set_last__is_critical( &self ) -> Option<bool> {  
        match self.num_lev_sets() == 0 { 
            true => None, 
            false => self.lev_set_ord__is_critical( self.num_lev_sets() )
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
        if let Some( is_crit ) = self.lev_set_last__is_critical() {
            if ! is_crit { 
                let i   =   self.num_lev_sets()-1;
                self.data_l__fmin[ i ] += 1; // increasing the min filtration value will make the level set critical            
            }
        }
        else { panic!("there are no level sets, so cannot make last critical") } // happens if there *is* no last level set to make critical      
    }

    fn  lev_set_ord__fmin( &self, level_set_ordinal: usize ) -> Option<Fil> { 
        if level_set_ordinal >= self.num_lev_sets() { return None }
        else {
            return Some ( self.data_l__fmin[ level_set_ordinal ].clone() )            
        }
    }

    // NB!!!!  If we assume that barcode endpoints are 0,1,..,N, then we can simplify
    // computation of max filtration values considerably
    fn  lev_set_ord__fmax( &self, level_set_ordinal: usize ) -> Option<Fil> { 

        // Posit: the ordinal is within legal range
        if let Some( is_critical ) = self.lev_set_ord__is_critical( level_set_ordinal.clone() )
        {
            // Posit: the level set ordinal IS critical.  
            // Therefore: max value equals min value.
            if is_critical { return self.lev_set_ord__fmin( level_set_ordinal ) }
            // Posit: the level set ordinal IS NOT not critical                
            // Therefore: max value equals min value + 1
            else { return self.lev_set_ord__fmin( level_set_ordinal ).map(|x| x+1 ) }


            //  OLD STUFF PROBABLY OK TO DELETE
            //     let fmin                    =   self.lev_set_ord__fmin( 
            //                                             level_set_ordinal.clone() 
            //                                         )
            //                                         .unwrap();
            //     if let Some( next_fil_val ) =   & self.data_l__fmin[
            //                                         level_set_ordinal.clone()..self.num_lev_sets()
            //                                         ]
            //                                             .iter()
            //                                             .map(|x| self.lev_set_ord__fmax(x.clone()).unwrap() )
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

    /// Mainimum ordinal filtration value allowed for this cell.
    pub fn  cell_id__fmin( &self, cell_id: usize ) -> Option<Fil> { 
        match self.cell_id__lev_set_ord( cell_id ) {
            Some( i )   =>  self.lev_set_ord__fmin( i ),
            None        =>  None
        }
    }

    /// Maximum ordinal filtration value allowed for this cell.
    fn  cell_id__fmax( &self, cell_id: usize ) -> Option<Fil> { 
        match self.cell_id__lev_set_ord( cell_id ) {
            Some( i )   =>  self.lev_set_ord__fmax( i ),
            None        =>  None
        }
    }    

    pub fn push_new_lev_set( &mut self )   {
        match self.max_filtration_value() {
            // if the vector of level set filtration values is nonempty, add a repeat of the last element
            Some( fil_max ) =>  self.data_l__fmin.push( fil_max ),
            // otherwise the vector is empty, and we should append a 0
            None            =>  self.data_l__fmin.push( 0 )
        }
    }

}


//  ---------------------------------------------------------------------------  
//  POLYTOPE VERTEX BUFFER
//  --------------------------------------------------------------------------- 


struct PolytopeVertexBuffer<'a> {
    poly:               &'a Polytope,      // polytope to work with
    cell_id__fmin:     Vec< Fil >,         // function (cell id) -> min possible filtration value    
    cell_id__fmax:     Vec< Fil >,         // function (cell id) -> max possible filtration value
}

fn poly_2_polysolve<'a> ( poly: &'a Polytope ) -> PolytopeVertexBuffer<'a>{

    let num_cells       =   poly.num_cells();
    let num_levels      =   poly.num_lev_sets();    

    let mut cell_id__fmin  =   Vec::with_capacity( num_cells );
    let mut cell_id__fmax  =   Vec::with_capacity( num_cells );

    // BOUNDARY CASES: #{CELLS} = 0 OR #{LEVEL SET} = 1
    if num_cells == 0 {
        return PolytopeVertexBuffer{ poly:poly, cell_id__fmin:cell_id__fmin, cell_id__fmax:cell_id__fmax}        
    } else 
    // this represents the case where there are 1 or more cells, but exactly 1 level set
    // (in this case the level set must be critical, on account of the infinite dim-0 bar)
    if num_levels == 1 {
        for cell_id in 0..num_cells {
            cell_id__fmin.push(0);
            cell_id__fmax.push(0);
        }
        return PolytopeVertexBuffer{ poly:poly, cell_id__fmin:cell_id__fmin, cell_id__fmax:cell_id__fmax}                
    }

    // REMAINING CASE: 
    //      >= 1 cell           AND
    //      >= 2 level sets 

    for i in 0..num_cells { 
        cell_id__fmin.push( poly.cell_id__fmin( i ).unwrap() );
        cell_id__fmax.push( poly.cell_id__fmax( i ).unwrap() );        
    }

    PolytopeVertexBuffer{ poly: &poly, cell_id__fmin: cell_id__fmin, cell_id__fmax: cell_id__fmax}

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


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
    pub fn  last_of_all_filtration_ordinals( &self ) -> Option<Fil> { self.data_l_to_fmin.iter().cloned().last() }

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
    pub fn lev_set_last_is_critical( &self ) -> Option<bool> {  
        match self.num_lev_sets() == 0 { 
            true => None, 
            false => self.lev_set_ord_to_is_critical( self.num_lev_sets() -1 )
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
        if let Some( is_crit ) = self.lev_set_last_is_critical() {
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

    /// Mainimum ordinal filtration value allowed for this cell.  Returns None if cell has not been assigned to a level set.
    pub fn  cell_id_to_fmin( &self, cell_id: usize ) -> Option<Fil> { 
        match self.cell_id_to_lev_set_ord( cell_id ) {
            Some( i )   =>  self.lev_set_ord_to_fmin( i ),
            None        =>  None
        }
    }

    /// Maximum ordinal filtration value allowed for this cell.  Returns None if cell has not been assigned to a level set.
    pub fn  cell_id_to_fmax( &self, cell_id: usize ) -> Option<Fil> { 
        match self.cell_id_to_lev_set_ord( cell_id ) {
            Some( i )   =>  self.lev_set_ord_to_fmax( i ),
            None        =>  None
        }
    }    

    /// Counts the number of strictly lower level set ordinals mapping to the same fmin value.
    pub fn  cell_id_to_critical_height( &self, cell_id: usize ) -> Option< usize > {
        match self.cell_id_to_lev_set_ord( cell_id ) {
            Some( i )   =>  { 
                let mut height = 0; 
                for k in (0..i).rev() {
                    match self.data_l_to_fmin[ k ] == self.data_l_to_fmin[ i ] {
                        true    => { height += 1; }
                        false   => { break; }
                    }
                }
                Some( height )
            } 
            None        =>  None
        }        
    }

    /// Post-pend a new (empty) level set.
    pub fn  push_new_lev_set( &mut self )   {
        match self.last_of_all_filtration_ordinals() {
            // if the vector of level set filtration values is nonempty, add a repeat of the last element
            Some( fil_max ) =>  self.data_l_to_fmin.push( fil_max ),
            // otherwise the vector is empty, and we should append a 0
            None            =>  self.data_l_to_fmin.push( 0 )
        }
    }

    /// A vector v such that v[ i ] = (the level set ordinal of cell i)
    pub fn  vec_mapping_cell_id_to_level_set_ordinal( &self ) -> Vec< usize > {
        self.data_c_to_l.clone()
    }

    /// A vector v such that v[ i ] = (the level set ordinal of cell i)
    pub fn  vec_mapping_cell_id_to_min_filt_ordinal( &self ) -> Vec< usize > {
        let mut vec         =   Vec::new();
        for cell_id in 0 .. self.num_cells() {
            if let Some( ord ) = self.cell_id_to_fmin( cell_id ) { vec.push( ord ) }
            else { vec.push( self.num_cells() ) }
        }
        vec
    }    

    pub fn  dim( &self ) -> Option< usize > {
        if let Some( max_filt_ord ) = self.data_l_to_fmin.last() {
            Some( self.num_lev_sets() - max_filt_ord - 2 )
        } else {
            None
        }
    }

    /// True or false: if each cell that have been assigned a level set ordinal
    /// takes the minimum possible value it can take, then are the resulting
    /// values equal to those of vector provided?  Cells that have not been assigned
    /// a level set ordinal are ignored.
    pub fn min_vertex_is_compatible_with_ordinal_filt(
                & self,
                ordinal_filtration:  & Vec< Fil >
            )
        -> 
        bool  
    {
        for cell_id in 0 .. self.num_cells() {

            if self.cell_id_to_fmin( cell_id ) == None {
                continue
            } else if   self.lev_set_last_is_critical() == Some( false )
                        &&
                        self.cell_id_to_lev_set_ord( cell_id ) == Some( self.num_lev_sets() -1 )
                        &&
                        self.cell_id_to_fmin( cell_id ).unwrap() + 1 == ordinal_filtration[ cell_id ]
                        {
                continue
            } else if let Some( a ) = self.cell_id_to_fmin( cell_id )  {
                if  a != ordinal_filtration[ cell_id ] 
                {                                      
                    return false 
                }
            }
        }
        return true 
    }

    /// Determine whether `self` contains the other polytope.
    pub fn  contains( &self, other_poly: Polytope ) -> bool {
        if  self.vec_mapping_cell_id_to_min_filt_ordinal()
            !=
            other_poly.vec_mapping_cell_id_to_min_filt_ordinal()
        { return false }
        
        for cell_id_a in 0 .. self.num_cells() {
            for cell_id_b in 0 .. self.num_cells() {
                if        self.cell_id_to_lev_set_ord( cell_id_a ) !=       self.cell_id_to_lev_set_ord( cell_id_b )
                    &&
                    other_poly.cell_id_to_lev_set_ord( cell_id_a ) == other_poly.cell_id_to_lev_set_ord( cell_id_b )
                { return false }
            }
        }

        return true
    }

}
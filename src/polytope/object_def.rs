

// use crate::utilities::*;
// use solar::utilities::indexing_and_bijection::{SuperVec, EndIndex};
// use ordered_float::OrderedFloat;
// use num::rational::Ratio;
// use std::collections::{HashMap};
// use std::hash::Hash;
use std::iter::{FromIterator, repeat};
use itertools::Itertools;
// use std::cmp::Ord;
// use std::iter;

use serde::{Deserialize, Serialize};
type Fil = usize;






//  ---------------------------------------------------------------------------  
//  POLYTOPES (ORDINAL)
//  ---------------------------------------------------------------------------  

/// Represents an (ordinal) polytope, or, equivalently, a family of upper/lower
/// bounds on variable values.
#[derive(Clone, Debug, PartialEq, PartialOrd, Eq, Ord, Hash, Serialize, Deserialize)]
pub struct Polytope {
    pub data_l_to_fmin:       Vec< Fil >,         // function (level set ordinal) -> min possible filtration value
    pub data_c_to_l:          Vec< usize >,       // function (cell id) -> level set ordinal
}

impl Polytope {

    /// Number of cells in the underlying chain complex.
    pub fn  num_cells( &self ) -> usize { self.data_c_to_l.len() }

    /// Number of level sets of the filter function.
    pub fn  num_lev_sets( &self ) -> usize { self.data_l_to_fmin.len() }    

    /// Number of level sets of the filter function.
    pub fn  num_simplex_factors( &self ) -> usize { 
        if self.num_lev_sets() > 0 { 1 + self.data_l_to_fmin.last().unwrap() } 
        else { 0 }
    }     
    
    /// Dimension of the polytope, calculated as (# level sets) - (# fmin values).
    /// 
    /// Here "cellagnostic" refers to the fact that we don't check that every level
    /// set contains at least one cell.
    pub fn  dim_cellagnostic( &self ) -> Option< usize > {
        if let Some( max_filt_ord ) = self.data_l_to_fmin.last() {
            Some( self.num_lev_sets() - max_filt_ord - 1 )
        } else {
            None
        }
    }

    /// Vector of dimensions of the simplices; deach dimension is calculated as 
    /// (# lev sets with this fmin) - 1.
    /// 
    /// Here "cellagnostic" refers to the fact that we don't check that every level
    /// set contains at least one cell.
    pub fn  simplex_factor_dims_cellagnostic( &self ) -> Vec< usize > {
        Vec::from_iter(
            ( 0 .. self.num_simplex_factors() ).map( |x| self.fmin_to_num_lev_sets( x ) - 1 )
        )
    }    
    

    /// Number of level sets with a given `fmin` value.
    pub fn  fmin_to_num_lev_sets( &self, fmin: usize ) -> usize {
        self.data_l_to_fmin.iter().filter(|&&x| x == fmin ).count()
    }

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

    /// Determine whether the given cell has been asigned to a level set.
    pub fn  cell_id_to_has_lev_set( &self, cell_id: usize ) -> Option< bool > { 
        let num_cells   =   self.num_cells();
        if cell_id < num_cells { 
            Some( self.data_c_to_l[ cell_id ] < num_cells ) 
        }  else {
            None
        }
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

    /// Determine whether the level set of the given cell (NOT THE CELL ITSELF) is critical.
    /// 
    /// Returns `None` if cell index is out of range.
    pub fn  cell_id_to_lev_set_is_critical( &self, cell_id: usize ) -> Option< bool > {
        match self.cell_id_to_lev_set_ord( cell_id.clone() ) {
            Some( lev_set_ord )     => self.lev_set_ord_to_is_critical( lev_set_ord ),
            None                    => None
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

    /// A vector v such that v[ i ] = (the minimum filtration ordinal of cell i)
    pub fn  vec_mapping_cell_id_to_min_filt_ordinal( &self ) -> Vec< usize > {
        let mut vec         =   Vec::new();
        for cell_id in 0 .. self.num_cells() {
            if let Some( ord ) = self.cell_id_to_fmin( cell_id ) { vec.push( ord ) }
            else { vec.push( self.num_cells() ) }
        }
        vec
    } 
    
    /// A vector v such that v[ i ] = (the maximum filtration ordinal of cell i)
    pub fn  vec_mapping_cell_id_to_max_filt_ordinal( &self ) -> Vec< usize > {
        let mut vec         =   Vec::new();
        for cell_id in 0 .. self.num_cells() {
            if let Some( ord ) = self.cell_id_to_fmax( cell_id ) { vec.push( ord ) }
            else { vec.push( self.num_cells() ) }
        }
        vec
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


    // TESTED (BUT CHECK THAT TESTS ARE UP TO DATE)
    /// Returns the union of all non-critical level sets with the same fmin value
    /// as the given cell, and which **precede** this cell's level set in order (the
    /// cell's own level set is included).
    /// 
    /// Vectors entries appear in ascending order.
    pub fn  cell_id_to_dnstream_cell_ids( &self, cell_id: usize ) -> Vec< usize > {
        
        let crit_ord        =   self.cell_id_to_fmin( cell_id );
        let lev_set_ord     =   self.cell_id_to_lev_set_ord( cell_id );

        Vec::from_iter(
            ( 0 .. self.num_cells() )
                .filter(    | other_cell_id |
                            crit_ord
                            == // other_cell_id belongs to the same family of "free" level sets
                            self.cell_id_to_fmin( other_cell_id.clone() )
                            &&
                            !  // and does not belong to a critical class
                            self.cell_id_to_lev_set_is_critical( other_cell_id.clone() ).unwrap()
                            &&
                            lev_set_ord
                            >= // and appears at or before the level set that contains this cell
                            self.cell_id_to_lev_set_ord( other_cell_id.clone() )
                )
        )
    }

    // TESTED (BUT CHECK THAT TESTS ARE UP TO DATE)
    /// Returns the union of all non-critical level sets with the same fmin value
    /// as the given cell, and which **succeed** this cell's level set in order (the
    /// cell's own level set is included).
    /// 
    /// Vectors entries appear in ascending order.
    pub fn  cell_id_to_upstream_cell_ids( &self, cell_id: usize ) -> Vec< usize > {
        
        let crit_ord        =   self.cell_id_to_fmin( cell_id );
        let lev_set_ord     =   self.cell_id_to_lev_set_ord( cell_id );

        Vec::from_iter(
            ( 0 .. self.num_cells() )
                .filter(    | other_cell_id |
                            crit_ord
                            == // other_cell_id belongs to the same family of "free" level sets as cell_id
                            self.cell_id_to_fmin( other_cell_id.clone() )
                            &&
                            !  // and does not belong to a critical class
                            self.cell_id_to_lev_set_is_critical( other_cell_id.clone() ).unwrap()
                            &&
                            lev_set_ord
                            <= // and appears at or before the level set that contains this cell
                            self.cell_id_to_lev_set_ord( other_cell_id.clone() )
                )
        )
    }        


    /// Determine whether `self` contains the other polytope.
    pub fn  contains( &self, othr: & Polytope ) -> bool {

        let num_cells       =   self.num_cells();
        
        // polytopes must have the same number of cells and the same number of fmin values
        if  num_cells != othr.num_cells() 
            ||
            self.num_simplex_factors() != othr.num_simplex_factors()
        { return false }

        // the bounds for each cell must be respected
        for cell_id in 0 .. num_cells {
            if  othr.cell_id_to_fmin( cell_id )   <   self.cell_id_to_fmin( cell_id )
                ||
                othr.cell_id_to_fmax( cell_id )   >   self.cell_id_to_fmax( cell_id )
            { return false }                
        }
        
        for cell_id_a in 0 .. self.num_cells() {
            
            for cell_id_b in 0 .. self.num_cells() {
                // let "pull back" refer to the "pull back" preorder on cells induced by the mapping of cells to level set; 
                // this is a binary relation on the set of cells.  the "pull back" for other should be a superset of the 
                // "pull back" for self
                if  (
                        self.cell_id_to_lev_set_ord( cell_id_a )    >=  self.cell_id_to_lev_set_ord( cell_id_b )     
                    )
                    &&
                    ! (
                        othr.cell_id_to_lev_set_ord( cell_id_a )    >=  othr.cell_id_to_lev_set_ord( cell_id_b )
                    )
                    
                { return false }
            }
        }

        return true
    }


    /// Given a partially-built polytope `othr` (meaning that `othr.cell_id_to_lev_set_ord` contains some
    /// values equal to the number of cells in the complex, implying that these `cell_id`s have not been
    /// assigned to a level set), determine whether `self` contains a polytope that the `othr` polytope 
    /// could extend to.
    pub fn  contains_extension( &self, othr: & Polytope ) -> bool {

        let num_cells       =   self.num_cells();
        
        // polytopes must have the same number of cells and the same number of fmin values
        if  num_cells != othr.num_cells() 
        { return false }

        // the bounds for each cell must be respected
        for cell_id in 0 .. num_cells {
            if  othr.cell_id_to_lev_set_ord( cell_id ).unwrap() < num_cells
                &&
                (
                    othr.cell_id_to_fmin( cell_id )   <   self.cell_id_to_fmin( cell_id )
                    ||
                    othr.cell_id_to_fmax( cell_id )   >   self.cell_id_to_fmax( cell_id )
                )

            { return false }                
        }
        
        for cell_id_a in 0 .. self.num_cells() {
            
            for cell_id_b in 0 .. self.num_cells() {
                // let "pull back" refer to the pull back preorder on cells induced by the mapping of cells to level set; 
                // this is a binary relation on the set of cells.  the "pull back" for other should be a superset of the 
                // "pull back" for self
                if  (
                        self.cell_id_to_lev_set_ord( cell_id_a )    >=  self.cell_id_to_lev_set_ord( cell_id_b )     
                    )
                    &&
                    ! (
                        othr.cell_id_to_lev_set_ord( cell_id_a )    >=  othr.cell_id_to_lev_set_ord( cell_id_b )
                    )
                    
                { return false }
            }
        }

        return true
    }  
    
    

    /// Given a partially-built polytope `othr` (meaning that `othr.cell_id_to_lev_set_ord` contains some
    /// values equal to the number of cells in the complex, implying that these `cell_id`s have not been
    /// assigned to a level set), certify that `othr` has been previously explored as a possibility in 
    /// the normal course of the phfibre algorithm
    pub fn  certifies_prior_exploration( &self, othr: & Polytope ) -> bool {

        let num_cells       =   self.num_cells();
        
        // polytopes must have the same number of cells and the same number of fmin values
        if  num_cells != othr.num_cells() 
        { return false }


        // create a mask that makes some elements of the last level set of othr appear to be 
        // unassigned, but returns proper level set assignment for all other cells
        let othr_lev_set_last               =   othr.num_lev_sets() - 1;
        let min_touched_lev_set_self_opt    =   ( 0 .. num_cells )
                                                    .filter( |x| othr.cell_id_to_lev_set_ord( *x ).unwrap() == othr_lev_set_last )
                                                    .map( |x| self.cell_id_to_lev_set_ord( x ).unwrap() )
                                                    .min();
        
        // can't use this trick if the level set currently being built in othr is critical
        //  ???? OR MAYBE YOU CAN???
        // if let Some( min_touched_lev_set_self_inner ) = min_touched_lev_set_self_opt { 
        //     if self.lev_set_ord_to_is_critical( min_touched_lev_set_self_inner ).unwrap() {
        //         return false
        //     } 
        // }       
        // CAN'T GET A CERTIFICATE IF TOO SIMILAR
        // if  min_touched_lev_set_self_opt == Some( othr_lev_set_last ) 
        //     && 
        //     self.data_c_to_l.iter().filter(|x| self.cell_id_to_lev_set_ord(**x) == min_touched_lev_set_self_opt).count()
        //     == 
        //     othr.data_c_to_l.iter().filter(|x| othr.cell_id_to_lev_set_ord(**x) == min_touched_lev_set_self_opt).count()
        //     { return false }
        

        // let min_touched_lev_set_self    =   

        // if let Some( min_touched_lev_set_self_inner ) = min_touched_lev_set_self_opt { 
        //     if self.lev_set_ord_to_is_critical( min_touched_lev_set_self_inner ).unwrap() {
        //         return false
        //     } 
        //     else {
        //         min_touched_lev_set_self_inner                
        //     }
        // };

        let mask = | cell_id: usize | -> usize { 
            if  othr.cell_id_to_lev_set_ord( cell_id ) == Some( othr_lev_set_last )  
                &&
                self.cell_id_to_lev_set_ord( cell_id ) != min_touched_lev_set_self_opt
            { return num_cells.clone() }
            else { return othr.cell_id_to_lev_set_ord( cell_id ).unwrap() }
        };

        // the bounds for each cell must be respected
        for cell_id in 0 .. num_cells {
            if  othr.cell_id_to_lev_set_ord( cell_id ).unwrap() < num_cells
                &&
                (
                    othr.cell_id_to_fmin( cell_id )   <   self.cell_id_to_fmin( cell_id )
                    ||
                    othr.cell_id_to_fmax( cell_id )   >   self.cell_id_to_fmax( cell_id )
                )

            { return false }                
        }
        
        for cell_id_a in 0 .. self.num_cells() {
            
            for cell_id_b in 0 .. self.num_cells() {
                // let "pull back" refer to the pull back preorder on cells induced by the mapping of cells to level set; 
                // this is a binary relation on the set of cells.  the "pull back" for other should be a superset of the 
                // "pull back" for self
                if  (
                        self.cell_id_to_lev_set_ord( cell_id_a )    >=  self.cell_id_to_lev_set_ord( cell_id_b )     
                    )
                    &&
                    ! (
                        mask( cell_id_a )    >=  mask( cell_id_b )
                    )
                    
                { return false }
            }
        }

        return true
    }      

}



//  ---------------------------------------------------------------------------  
//  POLYTOPES (3 VEC FORMAT)
//  --------------------------------------------------------------------------- 


// TESTED (BUT CHECK THAT TESTS ARE UP TO DATE)
/// Convert a Polytope into an equivlant object of type Vec< Vec< Vec< usize >>>
pub fn  poly_to_polyvvv( poly: & Polytope )
        ->
        Vec< Vec< Vec< usize >>> 
    {
    let mut polyvvv    =    Vec::from_iter( 
                                ( 0 .. poly.num_simplex_factors() )
                                    .map( 
                                        |x| 
                                        Vec::from_iter( repeat( vec![] ).take( poly.fmin_to_num_lev_sets(x) ) )
                                    )
                            );                                                    


    // for lev_set_fmin in poly.data_l_to_fmin.iter().enumerate() {
    //     polyvvv[ lev_set_fmin ].push( vec![] )
    // }
    // for crit_count in 0 .. poly.num_simplex_factors() {
    //     polyvvv.push( vec![]  );

    // }

    for cell_id in 0 .. poly.num_cells() {
        polyvvv
            [ poly.cell_id_to_fmin( cell_id ).unwrap() ]
            [ poly.cell_id_to_critical_height( cell_id ).unwrap() ]
            .push( cell_id );
    }
    polyvvv
}

// TESTED (BUT CHECK THAT TESTS ARE UP TO DATE)
/// Convert a polytope in Vec< Vec< Vec< usize >>> format to Polytope format.
pub fn  polyvvv_to_poly( polyvvv: & Vec< Vec< Vec< usize >>>  )
        ->
        Polytope
    {
    let num_cells           =   polyvvv.iter().flatten().flatten().count();
    let num_lev_sets        =   polyvvv.iter().flatten().count();
    
    let mut poly            =   Polytope{   
                                    data_c_to_l:        Vec::from_iter( repeat( num_cells ).take( num_cells) ),
                                    data_l_to_fmin:     Vec::with_capacity( num_lev_sets ),
                                };

    let mut lev_set_ord     =   0;
    for (crit_ord, lev_sets) in polyvvv.iter().enumerate() {
        for lev_set in lev_sets {
            poly.data_l_to_fmin.push( crit_ord );
            for cell_id in lev_set {
                poly.data_c_to_l[ *cell_id ] = lev_set_ord;
            }
            lev_set_ord += 1;
        }
    }

    poly
    
}

//  ---------------------------------------------------------------------------  
//  ENUMERATE ALL POLYTOPES OF GIVEN PARAMS
//  --------------------------------------------------------------------------- 


/// Enumerate all polytopes P such that (i) P has 2 or 3 simplex factors, the last
/// of which is a 0 simplex corresponding to an empty level set (this comports
/// with the output of the main phfibre algorithm), (ii) each simplex has 
/// dimension in {0, 1, 2}, (iii) each level set (except for the last, which is
/// empty) contains exactly one vertex.  For each `data_l_to_fmin`,
/// we generate a polytope for each permutation of `poly.num_cells = poly.fmin_max`.
/// 
/// This produces 40 polytopes total.
pub fn  enumerate_poly_test_set_40_total() -> Vec< Polytope > {

    let mut polytopes   =   Vec::new();
    for num_fmin_vals in 2..4{
        for num_lev_sets_per_fmin in ( 0 .. num_fmin_vals - 1 ).map(|x| 1 .. 3 ).multi_cartesian_product() {
            let mut data_l_to_fmin = Vec::new();
            for (fmin, num_lev_sets) in num_lev_sets_per_fmin.iter().enumerate() { 
                data_l_to_fmin.extend( repeat(fmin).take( * num_lev_sets ) );
            }
            data_l_to_fmin.push( num_fmin_vals - 1 ); // the last fmin val; we give this val a single (empty, critical) level set
            let num_cells   =  data_l_to_fmin.len() - 1; 
            for perm in ( 0 .. num_cells ).permutations(num_cells) {
                polytopes.push(
                    Polytope{
                        data_c_to_l:        perm,
                        data_l_to_fmin:     data_l_to_fmin.clone()
                    }
                )
            }
        }
    }
    return polytopes
}

//  ---------------------------------------------------------------------------  
//  GENERATE RANDOM POLYTOPE
//  --------------------------------------------------------------------------- 

use rand::Rng;

/// Generate a raondom polytope.
pub fn  random_polytope( 
            max_num_fmin: usize, 
            max_num_lev_sets_per_fmin: usize, 
            max_num_cells_per_lev_set: usize,
            attach_dim0_simplex_to_end: bool,
        )
        ->
        Polytope 
    {
    let mut rng = rand::thread_rng();        
    
    let mut data_l_to_fmin      =   Vec::new();
    let mut data_c_to_l         =   Vec::new();

    // fill the data_l_to_fmin vector
    let mut fmin = 0;
    let mut lev_set_ord = 0;
    for _ in 0 .. max_num_fmin {   // for each fmin 
        let num_lev_sets        =   rng.gen_range( 0 .. max_num_lev_sets_per_fmin + 1 ); 
        for k in 0 .. num_lev_sets {     // create random # of level sets (possibly 0, which means deleting this fmin)
            data_l_to_fmin.push( fmin.clone() );  // push a value to data_l_to_fmin for this level set
            let num_cells = rng.gen_range( 0 .. max_num_cells_per_lev_set + 1 ); // choose how many cells to place in this level set
            data_c_to_l.extend( repeat( lev_set_ord ).take( num_cells ) ); // push this number of entries to data_c_to_l
            lev_set_ord += 1; // increment the level set orginal in prepration for the next level set
        }
        if num_lev_sets > 0 { fmin += 1 } // increment fmin value if we are adding a level set as part of a new for-loop
    }

    // if desired, push a dimension-0 simplex to the end of the product (the corresponding level set will be empty)
    if attach_dim0_simplex_to_end {
        data_l_to_fmin.push( fmin );
    }

    Polytope{
        data_c_to_l:        data_c_to_l,
        data_l_to_fmin:     data_l_to_fmin,
    }

}






#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;


    #[test]
    fn test_cell_id_to_critical_height() { 
        let poly                =   Polytope { 
                                        data_l_to_fmin:     vec![0,1,1,2],
                                        data_c_to_l:        vec![3,2,1,0,1,2,3],
                                    };
        assert_eq!( poly.cell_id_to_critical_height(0), Some(0) );
        assert_eq!( poly.cell_id_to_critical_height(1), Some(1) );        
        assert_eq!( poly.cell_id_to_critical_height(2), Some(0) );                
        assert_eq!( poly.cell_id_to_critical_height(3), Some(0) );                        
        assert_eq!( poly.cell_id_to_critical_height(10), None );                        
    }       


    #[test]
    fn test_is_compatible_with_ordinal_filtration() {

        let polytope                =   Polytope{ 
                                            data_c_to_l: vec![2, 0, 1, 6, 6, 6],
                                            data_l_to_fmin: vec![0, 1, 1],
                                        };

        assert_eq!(     
            true, 
            polytope.min_vertex_is_compatible_with_ordinal_filt( &vec![2, 0, 1, 4, 4, 4] )
        );
        assert_eq!(     
            true, 
            polytope.min_vertex_is_compatible_with_ordinal_filt( &vec![1, 0, 1, 4, 4, 4] )
        );  

        let polytope                =   Polytope{ 
                                            data_c_to_l: vec![2, 0, 1, 6, 6, 6],
                                            data_l_to_fmin: vec![0, 1, 2],
                                        };

        assert_eq!(     
            false, 
            polytope.min_vertex_is_compatible_with_ordinal_filt( &vec![1, 0, 1, 4, 4, 4] )
        );                
    }


    #[test]
    fn test_cell_id_to_dnstream_cell_ids() {

        // polytope with the following level sets:
        //  0: crit val 0 { 8 9}   
        //  1: crit val 0 { 6 7 10} 
        //  2: crit val 0 { 4 5 11 12 } 
        //  3: crit val 1 { 2 3 }
        //  4: crit val 1 { 0 1 }
        let poly                =   Polytope{           //   0 1 2 3 4 5 6 7 8 9 101112
                                        data_c_to_l:    vec![4,4,3,3,2,2,1,1,0,0,1,2,2],
                                        data_l_to_fmin: vec![0,0,0,1,1]
                                    };
        assert_eq!(     
            vec![
                vec![ ], // lev set 0
                vec![ ], // lev set 0
                vec![ 6, 7, 10 ], // lev set 1
                vec![ 6, 7, 10 ], // lev set 1
                vec![ 6, 7, 10 ], // lev set 1                                   
                vec![ 4, 5, 6, 7, 10, 11, 12 ], // lev set 2
                vec![ 4, 5, 6, 7, 10, 11, 12 ], // lev set 2
                vec![ 4, 5, 6, 7, 10, 11, 12 ], // lev set 2
                vec![ 4, 5, 6, 7, 10, 11, 12 ], // lev set 2                                                                                            
                vec![ ], // lev set 3
                vec![ ], // lev set 3
                vec![ 0, 1 ], // lev set 4
                vec![ 0, 1 ], // lev set 4
            ],
            Vec::from_iter(
                vec![ 8, 9, 6, 7, 10, 4, 5, 11, 12, 2, 3, 0, 1 ]
                    .iter()
                    .map( |x| poly.cell_id_to_dnstream_cell_ids( x.clone() ) )
            )
        );             
    }

    #[test]
    fn test_cell_id_to_upstream_cell_ids() {

        // polytope with the following level sets:
        //  0: crit val 0 { 8 9}   
        //  1: crit val 0 { 6 7 10} 
        //  2: crit val 0 { 4 5 11 12 } 
        //  3: crit val 1 { 2 3 }
        //  4: crit val 1 { 0 1 }
        let poly                =   Polytope{           //   0 1 2 3 4 5 6 7 8 9 101112
                                        data_c_to_l:    vec![4,4,3,3,2,2,1,1,0,0,1,2,2],
                                        data_l_to_fmin: vec![0,0,0,1,1]
                                    };
        assert_eq!(     
            vec![
                vec![ 4, 5, 6, 7, 10, 11, 12 ], // lev set 0
                vec![ 4, 5, 6, 7, 10, 11, 12 ], // lev set 0
                vec![ 4, 5, 6, 7, 10, 11, 12 ], // lev set 1
                vec![ 4, 5, 6, 7, 10, 11, 12 ], // lev set 1
                vec![ 4, 5, 6, 7, 10, 11, 12 ], // lev set 1                                                                                          
                vec![ 4, 5, 11, 12 ], // lev set 2
                vec![ 4, 5, 11, 12 ], // lev set 2
                vec![ 4, 5, 11, 12 ], // lev set 2
                vec![ 4, 5, 11, 12 ], // lev set 2              
                vec![ 0, 1 ], // lev set 3
                vec![ 0, 1 ], // lev set 3
                vec![ 0, 1 ], // lev set 4
                vec![ 0, 1 ], // lev set 4
            ],
            Vec::from_iter(
                vec![ 8, 9, 6, 7, 10, 4, 5, 11, 12, 2, 3, 0, 1 ]
                    .iter()
                    .map( |x| poly.cell_id_to_upstream_cell_ids( x.clone() ) )
            )
        );           
    }    






    #[test]
    fn test_poly_to_polyvvv_and_vice_versa() { 

        //  NOTE: COULD IMPROVE THIS TEST BY USING SOME RANDOMLY GENERATED POLYTOPES

        let poly                =   Polytope{                  // 0  1  2  3  4  5  6  7  8
                                        data_c_to_l:        vec![ 5, 4, 3, 2, 1, 0, 1, 2, 3],
                                        data_l_to_fmin:     vec![0, 0, 1, 2, 2, 3],
                                    };
        let polyvvv             =   poly_to_polyvvv( &poly );
        let poly_recovered      =   polyvvv_to_poly( & polyvvv );

        assert_eq!( &poly, &poly_recovered );

        assert_eq!( 
            &polyvvv, 
            &vec![  
                vec![ vec![ 5 ], vec![ 4, 6 ] ],
                vec![ vec![ 3, 7 ] ],
                vec![ vec![ 2, 8 ], vec![ 1 ] ],
                vec![ vec![ 0 ] ],                                
            ]
        )
           
    }    
    
    #[test]
    fn test_enumerate_poly_test_set_40_total() {
        let polys               =   enumerate_poly_test_set_40_total();
        for poly in polys.iter().enumerate() { println!("{:?}", poly)}  
    }


    #[test]
    fn test_contains_extension() {
        let poly_a              =   Polytope { 
                                        data_l_to_fmin:     vec![0,1,1,2],
                                        data_c_to_l:        vec![0,1,2,3],
                                    };        
        let poly_b              =   Polytope { 
                                        data_l_to_fmin:     vec![0,1,2],
                                        data_c_to_l:        vec![0,1,1,2],
                                    };                                
        let poly_c              =   Polytope { 
                                        data_l_to_fmin:     vec![0,1],
                                        data_c_to_l:        vec![0,1,1,4],
                                    };    
        let poly_d              =   Polytope { 
                                        data_l_to_fmin:     vec![0,1],
                                        data_c_to_l:        vec![0,0,1,4],
                                    };                                        
                                    
        assert!(    poly_a.contains_extension( &poly_b ) );
        assert!(    poly_a.contains_extension( &poly_c ) );    
        assert!( !  poly_a.contains_extension( &poly_d ) );                
    }


    #[test]
    fn test_certifies_prior_exploration() {
        let poly_a              =   Polytope { 
                                        data_l_to_fmin:     vec![0,1,1,2],
                                        data_c_to_l:        vec![0,1,2,3],
                                    };        
        let poly_b              =   Polytope { 
                                        data_l_to_fmin:     vec![0,1,2],
                                        data_c_to_l:        vec![0,1,1,2],
                                    };                                
        let poly_c              =   Polytope { 
                                        data_l_to_fmin:     vec![0,1],
                                        data_c_to_l:        vec![0,1,1,4],
                                    };    
        let poly_d              =   Polytope { 
                                        data_l_to_fmin:     vec![0,1],
                                        data_c_to_l:        vec![0,0,1,4],
                                    };       
                                    
        
        // This test should return true wherever contains_extension returns true
        assert!(    poly_a.certifies_prior_exploration( &poly_b ) );
        assert!(    poly_a.certifies_prior_exploration( &poly_c ) );    

        // It should return false where there are disagreements about min/max
        // possible filtraiton values
        assert!( !  poly_a.certifies_prior_exploration( &poly_d ) );  
        
        // These tests, however, are unique to <certifies_prior_exploration>
        let poly_aa             =   Polytope{
                                        data_l_to_fmin:     vec![0,0,1,1],
                                        data_c_to_l:        vec![0,0,1,1,2,2,3,3,4,4],
                                    };

        let poly_bb             =   Polytope{
                                        data_l_to_fmin:     vec![0,0,1,1],
                                        data_c_to_l:        vec![0,0,1,1,2,2,2,2,10,10],
                                    };  
                                    
        let poly_cc             =   Polytope{
                                        data_l_to_fmin:     vec![0,0,1,1],
                                        data_c_to_l:        vec![0,0,1,1,2,10,2,2,10,10],
                                    };                                      
        
        assert!(    poly_aa.certifies_prior_exploration( &poly_bb ) );
        assert!( !  poly_aa.certifies_prior_exploration( &poly_cc ) );
     
    }

}    
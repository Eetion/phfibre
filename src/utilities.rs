use ordered_float::OrderedFloat;
use std::collections::{HashMap};
use std::hash::Hash;



//  ---------------------------------------------------------------------------
//  CONVERSION
//  ---------------------------------------------------------------------------


// Produces a vector of ordered floats from a vector of floats
pub fn to_ordered_float( v: & Vec< f64 > ) -> Vec< OrderedFloat< f64> > { v.iter().map(|x| OrderedFloat(x.clone()) ).collect(); }


//  ---------------------------------------------------------------------------
//  WORKING WITH VECTORS
//  ---------------------------------------------------------------------------


/// Last value of a vector.
pub fn end_val< T: Clone >( v: Vec<T> ) -> T  { v[ v.len() - 1].clone() }

/// Mutable reference to last value of a vector.
pub fn end_val_mut< 'a, T >( v: Vec<T> ) -> &'a mut T { &mut v[ v.len() - 1] }

/// Last ordinal for a vector
pub fn last_index< T > ( v: Vec<T> ) -> usize { v.len() - 1 }


//  ---------------------------------------------------------------------------
//  ORDERING
//  ---------------------------------------------------------------------------


#[derive(Clone, Debug)]
pub struct OrdinalData < T : PartialOrd + PartialEq + Hash > {
    ord_to_val:  Vec< T >,
    val_to_ord:  HashMap< T, usize >
}


/// Get the ordinal data for the range of values taken by a vector of ordered floats
pub fn ordinate( v: & Vec< OrderedFloat<f64> > ) -> OrdinalData< OrderedFloat< f64> > {
    // let a   =   to_ordered_float(v);
    let a   =   v.clone();
    let b   =   HashMap::new();
    a.sort();       // sort entries
    a.dedup();      // remove duplicates

    for (i, t) in a.enumerate() {
        b.insert( t.clone(), i.clone() );
    }

    OrdinalData { ord_to_val: a, val_to_ord: b }
}
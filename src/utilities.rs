






//  ---------------------------------------------------------------------------
//  WORKING WITH VECTORS
//  ---------------------------------------------------------------------------


/// Last value of a vector.
pub fn end_val< T: Clone >( v: Vec<T> ) -> T  { v[ v.len() - 1].clone() }

/// Mutable reference to last value of a vector.
pub fn end_val_mut< 'a, T >( v: Vec<T> ) -> &'a mut T { &mut v[ v.len() - 1] }

/// Last ordinal for a vector
pub fn end_index< T > ( v: Vec<T> ) -> Option< usize > { 
    match v.is_empty() { 
        true    =>  None, 
        false   =>  Some(v.len() - 1) 
    }
}


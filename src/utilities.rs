






//  ---------------------------------------------------------------------------
//  WORKING WITH VECTORS
//  ---------------------------------------------------------------------------


/// Last value of a vector.
pub fn end_val< T: Clone >( v: & Vec<T> ) -> Option< T >  {
    match v.is_empty() { 
        true    =>  None, 
        false   =>  Some( v[ v.len() - 1].clone() ) 
    }
}

/// Mutable reference to last value of a vector.
pub fn end_val_mut< 'a, T >( v: &'a mut Vec<T> ) -> Option< &'a mut T > {
    match end_index(v) { 
        None        =>  None, 
        Some(i)   =>  Some( &mut v[i] ) 
    }
}

/// Last ordinal for a vector
pub fn end_index< T > ( v: & Vec<T> ) -> Option< usize > { 
    match v.is_empty() { 
        true    =>  None, 
        false   =>  Some(v.len() - 1) 
    }
}




use solar::utilities::index::{BiMapSequential, compose_f_after_g, sort_perm};
use solar::utilities::ring::{MinusOneToPower};
use solar::rings::ring::{Semiring, Ring};
use std::cmp::Ordering;
use std::collections::{HashSet, HashMap};
use std::hash::Hash;
use itertools::Itertools;
use std::iter::FromIterator;
use std::fmt::Debug;


//  ---------------------------------------------------------------------------
//  COMBINATORIAL SIMPLEX (DEFINITION)
//  ---------------------------------------------------------------------------


// The type we use for vertices
#[derive(Debug, PartialEq, Clone, Eq, Hash)]
pub struct Simplex< Vertex > 
{ 
    pub vertices: Vec< Vertex >     //  vertices should be sorted in ascending order
} 

impl    < Vertex > 
        Simplex
        < Vertex >   
        {
    
    fn num_vertices( &self ) -> usize { self.vertices.len() }
    fn dim( &self ) -> usize { self.vertices.len() - 1 }
}        


impl    < Vertex >           
        PartialOrd for Simplex
        < Vertex >

    where   Vertex: Ord     {

    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl    < Vertex >   
        Ord for Simplex
        < Vertex >

    where Vertex: Ord   {

    fn cmp(&self, other: &Self) -> Ordering {

        // next compare simplex dimensions
        let comp = self.num_vertices().cmp( & other.vertices.len() );
        if comp != Ordering::Equal { return comp }

        // finally, compare simplices lexicographically
        return self.vertices.cmp( & other.vertices )
    }
}

impl    < Vertex >   
        IntoIterator for Simplex
        < Vertex >      {

    type Item = Vertex;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter { self.vertices.into_iter() }
}
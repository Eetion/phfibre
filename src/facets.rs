
use solar::utilities::index::BiMapSequential;
use solar::utilities::ring::{MinusOneToPower};
use solar::rings::ring::{Semiring, Ring};
use std::cmp::Ordering;
use std::collections::{HashSet, HashMap};
use std::hash::Hash;
use itertools::Itertools;
use std::iter::FromIterator;
use std::fmt::Debug;

// type Vertex = u16;



//  ---------------------------------------------------------------------------
//  COMBINATORIAL SIMPLEX (DEFINITION)
//  ---------------------------------------------------------------------------


// The type we use for vertices
#[derive(Debug, PartialEq, Clone, Eq, Hash)]
pub struct Simplex< Vertex > 
{ 
    pub vertices: Vec< Vertex >     //  vertices should be sorted in ascending order
} 


impl < Vertex >   PartialOrd for Simplex< Vertex >
    where Vertex: Ord
{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl < Vertex >   Ord for Simplex< Vertex >
    where Vertex: Ord
{
    fn cmp(&self, other: &Self) -> Ordering {

        // next compare simplex dimensions
        let comp = self.vertices.len().cmp( & other.vertices.len() );
        if comp != Ordering::Equal { return comp }

        // finally, compare simplices lexicographically
        return self.vertices.cmp( & other.vertices )
    }
}

impl < Vertex >   IntoIterator for Simplex< Vertex > 
{
    type Item = Vertex;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter { self.vertices.into_iter() }
}


//  ---------------------------------------------------------------------------
//  FACES FROM FACETS-OF-THE-COMPLEX
//  ---------------------------------------------------------------------------


pub fn  set_of_subsequences< IterFacet, Vertex >( facets: IterFacet ) -> HashSet< Vec< Vertex > > 
    where   IterFacet:      IntoIterator< Item = Vec< Vertex > >,
            Vertex:    Ord + Hash + Clone
{
    let mut faces       =   HashSet::new();
    for facet in facets {
        for seq_length in 1 .. facet.len() {
            for comb in facet.iter().cloned().combinations( seq_length ) {
                faces.insert( comb );
            }
        }
    }
    faces
}


//  NB: THE USE OF SIMPLICES RATHER THAN VECTORS IS IMPORTANT HERE, BECAUSE THE TWO STRUCTS HAVE
//      **DIFFERENT** TOTAL ORDERS
pub fn  ordered_sequence_of_faces< IterFacet, Vertex >( facets: IterFacet ) -> Vec< Simplex< Vertex > > 
    where   IterFacet:      IntoIterator< Item = Vec< Vertex > >,
            Vertex:    Ord + Hash + Clone
{
    let mut faces   =   set_of_subsequences(facets);
    let mut faces   =   Vec::from_iter( faces.drain().map(|x| Simplex{vertices: x}) );
    faces.sort();
    faces
}   

//  ---------------------------------------------------------------------------
//  FACETS-OF-A-SIMPLEX
//  ---------------------------------------------------------------------------

/// Maintains an "inner state" that steps through the facets of a simplex in 
/// ascending lexicographic order; only returns `Some(())` or `None`.
/// 
/// # Examples
/// 
/// ```
/// use phfibre::facets::{Simplex, AscendingFacetIteratorNoReturn};
/// 
/// let mut facet_iterator_noreturn     =   AscendingFacetIteratorNoReturn::new(
/// Simplex{ vertices: vec![0, 1, 2] },
/// None
/// );

/// let answers = vec![
/// AscendingFacetIteratorNoReturn { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![0, 1] }, deleted_vertex_index: Some(2) },
/// AscendingFacetIteratorNoReturn { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![0, 2] }, deleted_vertex_index: Some(1) },
/// AscendingFacetIteratorNoReturn { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![1, 2] }, deleted_vertex_index: Some(0) },
/// AscendingFacetIteratorNoReturn { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![]     }, deleted_vertex_index: None    },            
/// AscendingFacetIteratorNoReturn { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![0, 1] }, deleted_vertex_index: Some(2) },                                                                        
/// ];

/// for i in 0..5 {
/// facet_iterator_noreturn.next();
/// assert_eq!( &facet_iterator_noreturn, &answers[ i ] )    
/// }   
/// ```
#[derive(Clone, Debug, PartialEq)]
pub struct  AscendingFacetIteratorNoReturn< Vertex >
{
    pub simplex: Simplex< Vertex > ,
    pub facet: Simplex< Vertex >,
    pub deleted_vertex_index: Option<usize>
}

impl < Vertex > AscendingFacetIteratorNoReturn < Vertex > 
    where Vertex: Clone
{

    /// Initialize a no-return facet iterator.  
    /// 
    /// Its internal state does not represent a facet; rather, the internal state will
    /// represent the first facet after `next()` is called for the first time.
    pub fn new( simplex: Simplex< Vertex >, buffer: Option< Simplex<Vertex> > ) -> Self {

        let buff = 
            if let Some( vec ) = buffer { vec } 
            else { 
                Simplex{ 
                    vertices: Vec::with_capacity( simplex.vertices.len()-1 )
                }             
            };
        
        AscendingFacetIteratorNoReturn {
            simplex: simplex,
            facet: buff,
            deleted_vertex_index: None
        }
    }
}

impl < Vertex >
    Iterator for 
    AscendingFacetIteratorNoReturn < Vertex >     
    where Vertex : Clone
{
    type Item    =   ();

    fn next( &mut self ) -> Option<()> {

        if let Some( deleted_index ) = self.deleted_vertex_index {

            if deleted_index == 0 {
                // if we start from the facet obtained by deleting vertex 0, then the 
                // next state should **not** represent a facet
                self.deleted_vertex_index   =   None;
                self.facet.vertices.clear();
                return None
                
            } else {
                // if we start from the facet obtained by deleting vertex k > 0, then 
                // the next state should represent the facet obtained by deleting vertex k-1
                let next_deleted_index  =   deleted_index - 1;
                self.facet.vertices[ next_deleted_index ] = self.simplex.vertices[ deleted_index ].clone(); // replace the deleted index and remove the next one
                self.deleted_vertex_index = Some( next_deleted_index );
                return Some( () )
            }
        
        } else {

            // if we start from the state representing no facet, then the next
            // state should represent the first facet (obtained by deleting the
            // last vertex)
            self.facet.vertices.clear();
            for i in 0..self.simplex.vertices.len()-1 { 
                self.facet.vertices.push( self.simplex.vertices[ i ].clone() ) 
            }      
            // set deleted vertex equal to last
            self.deleted_vertex_index = Some( self.simplex.vertices.len() - 1 );                  
            return Some(())            
        }
        // if self.deleted_vertex_index == None {
        //     // reset facet to equal the (lexicographically) first facet
        //     self.facet.vertices.clear();
        //     for i in 0..self.simplex.vertices.len()-1 { 
        //         self.facet.vertices.push( self.simplex.vertices[ i ].clone() ) 
        //     }
        //     // set deleted vertex equal to last
        //     self.deleted_vertex_index = Some( self.simplex.vertices.len() - 1 );
        // } else if self.deleted_vertex_index == 0 {

        // }

        // if self.deleted_vertex_index == 0 { None }
        // else {
        //     self.facet.vertices[ self.deleted_vertex_index -1 ] = self.simplex.vertices[ self.deleted_vertex_index ].clone();
        //     self.deleted_vertex_index -= 1;            
        //     Some(())
        // }
    }
}



//  ---------------------------------------------------------------------------
//  CONSTRUCT BOUNDARY MATRIX FROM FACE SEQUENCE
//  ---------------------------------------------------------------------------



// pub fn  boundary_data_from_facets< Vertex, RingOp, RingElt >( 
//             face_bimap: BiMapSequential< Simplex< Vertex > >,
//             ring:       RingOp
//         ) 
//         ->
//         Vec< Vec < (usize, RingElt) >>
//             where   Vertex:    Ord + Hash + Clone,        
//                     RingOp:     Semiring< RingElt > + Ring< RingElt >,
// {
//     let boundary    =   Vec::with_capacity( face_bimap.ord_to_elt.len() );

//     let mut template_facet  =   Simplex{ vertices: Vec::new() };
//     let mut int_index       =   0;
//     let mut num_vertices    =   0;    

//     for (simplex_count, simplex) in face_bimap.ord_to_elt.iter().cloned().enumerate() {

//         let verts           =   simplex.vertices;

//         if verts.len() == 1 { continue }        

//         template_facet.vertices.clear();

//         // write the first facet (obtained by deleting the last vertex) onto the template
//         for i in 0..verts.len()-1 { template_facet.vertices.push( verts[ i ] ) }

//         // get the corresponding global integer index
//         int_index           =   face_bimap.ord( template_facet );

//         // coefficient 
//         let coeff           =   ring.minus_one_to_power( verts.len() );

//         // push to vector
//         boundary[ simplex_count ].push( (int_index, coeff) )

//         // for every other facet
//         for j in 1 .. verts.len() {

//             // copy the facet to the template
//             template_facet[ j-1 ] = verts[ j ].clone();

//             // get its index
//             int_index           =   face_bimap.elt_to_ord( template_facet );

//             // coefficient 
//             let coeff           =   RingOp::one();

//             // push to vector
//             boundary[ simplex_count ].push( (int_index, coeff) )

//         }

//     }

//     boundary
// }





#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_facet_iterator()
    {

        let mut facet_iterator_noreturn     =   AscendingFacetIteratorNoReturn::new(
                                                    Simplex{ vertices: vec![0, 1, 2] },
                                                    None
                                                );

        let answers = vec![
            AscendingFacetIteratorNoReturn { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![0, 1] }, deleted_vertex_index: Some(2) },
            AscendingFacetIteratorNoReturn { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![0, 2] }, deleted_vertex_index: Some(1) },
            AscendingFacetIteratorNoReturn { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![1, 2] }, deleted_vertex_index: Some(0) },
            AscendingFacetIteratorNoReturn { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![]     }, deleted_vertex_index: None    },            
            AscendingFacetIteratorNoReturn { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![0, 1] }, deleted_vertex_index: Some(2) },                                                                        
        ];

        for i in 0..5 {
            facet_iterator_noreturn.next();
            assert_eq!( &facet_iterator_noreturn, &answers[ i ] )    
        }      
                
    }     

}    
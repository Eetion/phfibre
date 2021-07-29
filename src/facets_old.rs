
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

//  ---------------------------------------------------------------------------
//  FACES FROM FACETS-OF-THE-COMPLEX (NEW)
//  ---------------------------------------------------------------------------

pub fn  subsequences_up_to_card< Vertex >( 
        complex_facets: Vec< Vec< Vertex >>, 
        max_card: usize 
    ) 
    -> 
    Vec< Vec< Vec< Vertex >>> 
    where Vertex: Ord + Clone
{
    let mut seq             =   Vec::with_capacity( max_card );
    seq.push( Vec::with_capacity(0) ); // no empty sequences
    for card in 1 .. max_card + 1 {
        let vec: Vec<_>     =   complex_facets
                                .iter()
                                .map( |x| x.iter().cloned().combinations(card)  )
                                .kmerge()
                                .dedup()
                                .collect();
        seq.push( vec );
    }
    seq
    
}



//  ---------------------------------------------------------------------------
//  FACES FROM FACETS-OF-THE-COMPLEX
//  ---------------------------------------------------------------------------

/// Given something that iterates over vectors (each of which represents a strictly 
/// ascending sequence of vertices), return a HashSet containing all nonempty subsequences.
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

/// Given something that iterates over vectors (each of which represents a strictly 
/// ascending sequence of vertices), return a vector V containing all nonempty ordered
/// subsequences; V is strictly ascending under the order that first compares length of 
/// a sequence, then compares equal-length sequences lexicographically.
/// 
//  NB: THE USE OF SIMPLICES RATHER THAN VECTORS IS IMPORTANT HERE, BECAUSE THE TWO STRUCTS HAVE
//      **DIFFERENT** TOTAL ORDERS
pub fn  ordered_sequence_of_faces< IterFacet, Vertex >( facets: IterFacet ) -> Vec< Simplex< Vertex > > 
    where   IterFacet:  IntoIterator< Item = Vec< Vertex > >,
            Vertex:     Ord + Hash + Clone
{
    println!("THIS FUNCTION COULD BE MADE MUCH MORE EFFICIENT BY FIRST MAKING SEPARATE HASH SETS FOR SEQUENCES OF EACH LENGTH, THEN CONCATENATING");
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
/// use phfibre::facets::{Simplex, FacetIteratorNoReturnAscending};

/// // Create the iterator
/// let mut facet_iterator_noreturn     =   FacetIteratorNoReturnAscending::new(
///                                             Simplex{ vertices: vec![0, 1, 2] },
///                                             None
///                                         );
///
/// // Test it                                                
/// let mut answers = vec![
///     FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![0, 1] }, deleted_vertex_index: Some(2) },
///     FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![0, 2] }, deleted_vertex_index: Some(1) },
///     FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![1, 2] }, deleted_vertex_index: Some(0) },
///     FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![]     }, deleted_vertex_index: None    },            
///     FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![0, 1] }, deleted_vertex_index: Some(2) },                                                                        
/// ];
///
/// for i in 0..5 {
///     facet_iterator_noreturn.next();
///     assert_eq!( &facet_iterator_noreturn, &answers[ i ] )    
/// }      
//
/// // Re-initialize with a new simplex
///
/// facet_iterator_noreturn.reinitialize_with_simplex( Simplex{ vertices: vec![0 ,3]} );
///
/// // Test again        
///
/// answers = vec![
///     FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 3] }, facet: Simplex { vertices: vec![0] }, deleted_vertex_index: Some(1) },
///     FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 3] }, facet: Simplex { vertices: vec![3] }, deleted_vertex_index: Some(0) },
///     FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 3] }, facet: Simplex { vertices: vec![]  }, deleted_vertex_index: None    },
///     FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 3] }, facet: Simplex { vertices: vec![0] }, deleted_vertex_index: Some(1) },                                                                                                      
/// ];    
///
/// for i in 0..4 {
///     facet_iterator_noreturn.next();
///     assert_eq!( &facet_iterator_noreturn, &answers[ i ] )    
/// }   
/// ```
#[derive(Clone, Debug, PartialEq)]
pub struct  FacetIteratorNoReturnAscending< Vertex >
{
    pub simplex: Simplex< Vertex > ,
    pub facet: Simplex< Vertex >,
    pub deleted_vertex_index: Option<usize>
}

impl < Vertex > FacetIteratorNoReturnAscending < Vertex > 
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
                    vertices: Vec::with_capacity( simplex.dim() ) // = 1 less than num_vertices = NUM VERTICES OF A FACET
                }             
            };
        
        FacetIteratorNoReturnAscending {
            simplex: simplex,
            facet: buff,
            deleted_vertex_index: None
        }
    }

    pub fn reinitialize_with_simplex( &mut self, simplex: Simplex< Vertex > ) {

        // if necessary, expand the capacity of the facet vector
        if simplex.dim() > self.facet.vertices.capacity() { 
            self.facet.vertices.reserve_exact(
                simplex.dim() - self.facet.vertices.capacity()
            ) 
        }
        // replace the old simplex with the new
        self.simplex    =   simplex;
        // update the state to indicate that it does not represent a facet
        self.facet.vertices.clear();
        self.deleted_vertex_index = None;
    }
}


impl < Vertex >
    Iterator for 
    FacetIteratorNoReturnAscending < Vertex >     
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

            self.facet.vertices.clear();
            for i in 0..self.simplex.dim() {   // dim = 1 less than num_vertices = INDEX OF LAST VERTEX IN SIMPLEX
                self.facet.vertices.push( self.simplex.vertices[ i ].clone() ) 
            }      
            // set deleted vertex equal to last
            self.deleted_vertex_index = Some( self.simplex.dim() );  // dim = 1 less than num_vertices = INDEX OF LAST VERTEX IN SIMPLEX
            return Some(())            
        }
    }
}



//  ---------------------------------------------------------------------------
//  CONSTRUCT BOUNDARY MATRIX FROM FACE SEQUENCE
//  ---------------------------------------------------------------------------



pub fn  boundary_matrix_from_complex_facets< Vertex, RingOp, RingElt >( 
            simplex_bimap:  BiMapSequential< Simplex< Vertex > >,
            ring:           RingOp
        ) 
        ->
        Vec< Vec < (usize, RingElt) >>

        where   Vertex:    Ord + Hash + Clone + Debug,      
                RingOp:     Semiring< RingElt > + Ring< RingElt >,
{
    if simplex_bimap.ord_to_elt.is_empty() { return vec![] }

    let mut boundary            =   Vec::with_capacity( simplex_bimap.ord_to_elt.len() );  
    
    let mut state_iter          =   FacetIteratorNoReturnAscending{
                                        simplex: Simplex{ vertices: vec![] },
                                        facet: Simplex{ vertices: vec![] },
                                        deleted_vertex_index: None
                                    };

    let mut global_int_index    =   0;
    let mut simplex_dim         =   0;
    let mut simplex_num_verts   =   0;

    for simplex in simplex_bimap.ord_to_elt.iter().cloned() {

        simplex_dim             =   simplex.dim();
        simplex_num_verts       =   simplex.num_vertices();
        state_iter.reinitialize_with_simplex( simplex );

        let mut vec             =   Vec::with_capacity( simplex_num_verts );    // num_vertices = NUMBER OF FACETS
        
        for i in 0 .. simplex_num_verts {
            state_iter.next();
            
            println!("{:?}", &state_iter);
            println!("{:?}", &simplex_bimap);            

            global_int_index    =   simplex_bimap.ord( &state_iter.facet );
            vec.push( 
                (
                    global_int_index.clone(),
                    ring.minus_one_to_power( simplex_dim - i )
                ) 
            )
        }
        boundary.push( vec );
    }

    boundary


    // let mut template_facet  =   Simplex{ vertices: Vec::new() };
    // let mut int_index       =   0;
    // let mut num_vertices    =   0;    

    // for (simplex_count, simplex) in simplex_bimap.ord_to_elt.iter().cloned().enumerate() {

    //     let verts           =   simplex.vertices;

    //     if verts.len() == 1 { continue }        

    //     template_facet.vertices.clear();

    //     // write the first facet (obtained by deleting the last vertex) onto the template
    //     for i in 0..verts.len()-1 { template_facet.vertices.push( verts[ i ] ) }

    //     // get the corresponding global integer index
    //     int_index           =   simplex_bimap.ord( template_facet );

    //     // coefficient 
    //     let coeff           =   ring.minus_one_to_power( verts.len() );

    //     // push to vector
    //     boundary[ simplex_count ].push( (int_index, coeff) )

    //     // for every other facet
    //     for j in 1 .. verts.len() {

    //         // copy the facet to the template
    //         template_facet[ j-1 ] = verts[ j ].clone();

    //         // get its index
    //         int_index           =   simplex_bimap.elt_to_ord( template_facet );

    //         // coefficient 
    //         let coeff           =   RingOp::one();

    //         // push to vector
    //         boundary[ simplex_count ].push( (int_index, coeff) )

    //     }

    // }

    // boundary
}





#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_facet_iterator()
    {

        // Create the iterator
        let mut facet_iterator_noreturn     =   FacetIteratorNoReturnAscending::new(
                                                    Simplex{ vertices: vec![0, 1, 2] },
                                                    None
                                                );

        // Test it                                                
        let mut answers = vec![
            FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![0, 1] }, deleted_vertex_index: Some(2) },
            FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![0, 2] }, deleted_vertex_index: Some(1) },
            FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![1, 2] }, deleted_vertex_index: Some(0) },
            FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![]     }, deleted_vertex_index: None    },            
            FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 1, 2] }, facet: Simplex { vertices: vec![0, 1] }, deleted_vertex_index: Some(2) },                                                                        
        ];

        for i in 0..5 {
            facet_iterator_noreturn.next();
            assert_eq!( &facet_iterator_noreturn, &answers[ i ] )    
        }      

        // Re-initialize with a new simplex

        facet_iterator_noreturn.reinitialize_with_simplex( Simplex{ vertices: vec![0 ,3]} );

        // Test again        

        answers = vec![
            FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 3] }, facet: Simplex { vertices: vec![0] }, deleted_vertex_index: Some(1) },
            FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 3] }, facet: Simplex { vertices: vec![3] }, deleted_vertex_index: Some(0) },
            FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 3] }, facet: Simplex { vertices: vec![]  }, deleted_vertex_index: None    },
            FacetIteratorNoReturnAscending { simplex: Simplex { vertices: vec![0, 3] }, facet: Simplex { vertices: vec![0] }, deleted_vertex_index: Some(1) },                                                                                                      
        ];    

        for i in 0..4 {
            facet_iterator_noreturn.next();
            assert_eq!( &facet_iterator_noreturn, &answers[ i ] )    
        }     
                
    }   
    
    #[test]
    fn test_sequentialize_facets () {

        let ring                    =   solar::rings::ring_native::NativeDivisionRing::< f64 >::new();
        let complex_facets          =   vec![ vec![0,1,2] ];
        let bimap_sequential        =   BiMapSequential::from_vec(
                                            ordered_sequence_of_faces( complex_facets )
                                        );  
        let boundary                =   boundary_matrix_from_complex_facets( bimap_sequential, ring );

        for vec in boundary { println!("{:?}", vec) }
    }

    #[test]
    fn test_subsequences_up_to_card() {

        let facets                  =   vec![ vec![0, 1, 2] ];
        let subsequences            =   subsequences_up_to_card( facets, 3);
        println!("{:?}", & subsequences);

        assert_eq!(     &   subsequences,
                        &   vec![
                                vec![                                           ],
                                vec![  vec![0],     vec![1],    vec![2]         ],                                
                                vec![  vec![0,1],   vec![0,2],  vec![1,2]       ],
                                vec![  vec![0,1,2]                              ]
                            ]
        )
    }


}    
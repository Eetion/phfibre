

use std::iter::FromIterator;
use solar::utilities::indexing_and_bijection::compose_f_after_g;


/// Relabels the vertices of each simplex according to the provided bijection, then sorts
/// the vertices of each simplex.
pub fn  simplices_with_relabeled_vertices( 
            simplices:  &   Vec< Vec< usize >>,
            perm_v_o2n:  &   Vec< usize >,  // map old vertes to new vertex
            ) 
            -> Vec< Vec< usize >> 
{
    let mut permuted_simplices  =   Vec::from_iter(
                                        simplices
                                            .iter()
                                            .map(   |x|
                                                    compose_f_after_g( &perm_v_o2n, x )                            
                                            )
                                    );
    for simplex in permuted_simplices.iter_mut() { simplex.sort() }
    return permuted_simplices
}
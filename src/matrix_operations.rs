use solar::rings::ring_native::NativeDivisionRing;
use solar::rings::ring::{Semiring, Ring, DivisionRing};
use num::rational::Ratio;
use std::iter::IntoIterator;



//let ring : NativeDivisionRing< f64 > =  NativeDivisionRing::<f64>::new();

type Key = usize;
type Val = Ratio< i64 >;


/// Clear an entry of the clearee, using the clearor.
/// 
/// Assumes that all structural nonzero entries are *actually* nonzero
pub fn  clear_if_in< Key, Val, Ring > (
            clearor:        Vec< (Key, Val) >,
            clearee:        Vec< (Key, Val) >,
            buffer:         Vec< (Key, Val) >,
            pivot_entry:         (Key, Val),
            ring:           Ring
        )
    where   Ring: Semiring + Ring + DivisionRing
{
    let entry_to_clear_opt  =   clearee
                                    .iter()
                                    .find( |&x| x.0 == pivot_entry.0 );
    if let Some(entry_to_clear) = entry_to_clear_opt {

        scalar              =   ring.divide( 
                                    ring.negate( pivot_entry.1 ),
                                    entro_to_clear.1
                                );

        merged              =   merge(
                                    clearee.iter().cloned(),
                                    clearor
                                        .iter()
                                        .cloned()
                                        .scale_by( ring, scalar )
                                )
                                .simplify( ring );
        buffer
            .clear()
            .extend( merged );

        clearee
            .clear()
            .append(buffer);
    }

}






pub fn column_pivot<I,R>(   matrix:             Vec< Vec< ( Key, Val ) > >,
                            piv_col:            Vec< ( Key, Val ) >,
                            red_col_indices:    I,
                            ring:               R
                        ) 
    where I : IntoIterator< Item = Key >,
          R : Semiring<Val> + Ring<Val> + DivisionRing<Val>
{

    if piv_col.is_empty() { println!("Error in <column_pivot>: empty pivot column "); return }

    let piv_row_ind     =   piv_col.len() - 1;

    if ring.is_0(    piv_col[ piv_row_ind ].1    )
    {
        println!("Error in <column_pivot>: lowest structural entry is zero"); return
    }
}

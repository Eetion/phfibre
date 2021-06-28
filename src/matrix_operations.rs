use solar::rings::ring_native::NativeDivisionRing;
use solar::rings::ring::{Semiring, Ring, DivisionRing};
use num::rational::Ratio;
use std::iter::IntoIterator;



//let ring : NativeDivisionRing< f64 > =  NativeDivisionRing::<f64>::new();

type Key = usize;
type Val = Ratio< i64 >;



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

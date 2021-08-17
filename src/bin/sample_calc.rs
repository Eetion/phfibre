

use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use phfibre::pipelines::simplex_pipeline;
use solar::utilities::sequences_and_ordinals::{ordinate_unique_vals};
use solar::cell_complexes::simplices_unweighted::facets::{    
    ordered_subsimplices_up_thru_dim_concatenated_vec 
}; 









fn main() {

    //  ----------------------------------------------------------------------------------------------    
    //  B_3_1
    //  ----------------------------------------------------------------------------------------------

    //  CODE
    //  ----

    //  Define the base space, barcode, and ring
    let complex_facets      =   vec![  vec![0,1,2] ];
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 1);    

    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:2} ],
                                    fin: vec![ BarFinite{dim:0,birth:0,death:1}],
                                    ordinal: ordinate_unique_vals( & vec![0, 1, 2] ),
                                };

    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

    simplex_pipeline(
        &   simplex_sequence,
        &   barcode,
        &   ring
    )

    //  RESULTS
    //  -------

    // POLYHEDRAL COMPLEX FACETS
    // number of facets (total): 18
    // number of facets (binned by dimension): [0, 18]
    // POLYHEDRAL COMPLEX CELLS
    // number of polytopes (total): 36
    // number of polytopes (binned by dimension): [18, 18]
    // betti numbers (of polyhedral complex): [1, 1]
    // INTERSECTIONS
    // number of pairs of intersecting facets, binned by the dimension of the intersection polytope: [18, 0]
    // DOWKER NERVE COMPLEX FACETS
    // number of dowker nerve complex facets (total): 36
    // number of dowker nerve complex facets (binned by dimension): [18, 18]
    // DOWKER NERVE COMPLEX CELLS
    // number of nerve dowker complex cells (total): 36
    // number of nerve dowker complex cells (binned by dimension): [18, 18]
    // betti numbers (of dowker nerve): [1, 1]    
    

}
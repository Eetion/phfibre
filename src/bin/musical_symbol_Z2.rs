//  DESCRIPTION

//  A "musical symbol" consisting of a triangle and a dangling edge


//  DEPENDENCIES

use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use phfibre::pipelines::{simplex_pipeline, analyze_fibre};
use phfibre::phfibre::{ConditionNone, ConditionLowerStar, ConditionLowerEdge};
use solar::utilities::sequences_and_ordinals::ordinate_unique_vals;
use solar::utilities::cell_complexes::simplices_unweighted::facets::ordered_subsimplices_up_thru_dim_concatenated_vec; 


//  CODE

fn main() {

    let save_dir_opt        =   None; // we will not save any files    
        

    //  ----------------------------------------------------------------------------------------------    
    //  EMPTY TRIANGLE WITH EDGE ATTACHED
    //  ----------------------------------------------------------------------------------------------


    //  Define the base space, barcode, and ring
    let complex_facets      =   vec![  vec![0,1,2], vec![2, 3] ]; // triangle with edge attached
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 1);    // build all faces up to dimension 1

    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:2} ],
                                    fin: vec![ BarFinite{dim:0,birth:0,death:1} ],
                                    ordinal: ordinate_unique_vals( & vec![0, 1, 2] ),
                                };

    let ring                =   solar::rings::field_prime::GF2{};

    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};      // this struct won't impose any extra conditions on the filtrations we build    
    let analyze_dowker_dual =   true;

    let poly_complex_facets =   simplex_pipeline(
                                    &   simplex_sequence,
                                    &   barcode,
                                    &   ring,
                                    &   precondition_to_make_new_lev_set_lower_none,
                                        analyze_dowker_dual, // do not analyze the dowker dual to the nerve complex
                                        save_dir_opt,
                                );  

    //  ANALYZE    

    let analyze_dowker_dual =   false;
    analyze_fibre( 
        &   poly_complex_facets,
            ring.clone(),
            analyze_dowker_dual,
            save_dir_opt,
    ); 

    //  RESULTS
    //  -------

    // BASE SPACE
    // simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [1, 2], [2, 3]]
    // BARCODE
    // barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 2 }], fin: [BarFinite { dim: 0, birth: 0, death: 1 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2], val_to_ord: {0: 0, 1: 1, 2: 2} } }
    // TIME TO COMPUTE FIBRE FACETS
    // Time elapsed to compute facets of PH fibre: 18.321374ms
    
    // ANALYSIS
    // Each polytope facet has been checked for compatiblity with the given barcode.
    // POLYHEDRAL COMPLEX FACETS
    // number of facets (total): 104
    // number of facets (binned by dimension): [0, 0, 104]
    // number of facets (binned by number of vertices): [0, 0, 104]
    // POLYHEDRAL COMPLEX CELLS
    // "Time elapsed to compute the faces of PH fibre (given the facets): 2.823198ms"
    // number of polytopes (total): 393
    // number of polytopes (binned by dimension): [92, 197, 104]
    // Time elapsed to compute binary-coeff polyhedral betti numbers 224.662µs
    // betti numbers (of polyhedral complex, Z2 coefficients): [1, 2, 0]


}
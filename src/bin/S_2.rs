




use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use phfibre::pipelines::{simplex_pipeline, analyze_fibre};
use phfibre::phfibre::{ConditionNone, ConditionLowerStar, ConditionLowerEdge};
use phfibre::polytope::object_def::Polytope;
use solar::utilities::sequences_and_ordinals::ordinate_unique_vals;
use solar::cell_complexes::simplices_unweighted::facets::ordered_subsimplices_up_thru_dim_concatenated_vec; 
use std::iter::FromIterator;

use std::io::Read;





//  CODE

fn main() {


    //  ----------------------------------------------------------------------------------------------    
    //  2-SKELETON OF 3-SIMPLEX, 1 FINITE BAR IN DIM 0, 1 FINITE BAR IN DIM 1
    //  ----------------------------------------------------------------------------------------------

    //  Define the base space, barcode, and ring

    let complex_facets      =   vec![  vec![0, 1, 2, 3] ];
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 2); 

    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:2,birth:5} ],
                                    fin: vec![ BarFinite{dim:0,birth:1,death:2}, BarFinite{dim:1,birth:3,death:4} ],
                                    ordinal: ordinate_unique_vals( & vec![ 0, 1, 2, 3, 4, 5 ] ),
                                };

    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};      // this struct won't impose any extra conditions on the filtrations we build    

    println!("\n\n2-SKELETON, BARCODE:(0, [0,INF)), (0, [1,2)), (1, [3,4)), (2, [5,INF))");

    let poly_complex_facets =       simplex_pipeline(
                                        &   simplex_sequence,
                                        &   barcode,
                                        &   ring,
                                        &   precondition_to_make_new_lev_set_lower_none,
                                            false, // do not analyze the dowker dual to the nerve complex
                                    );  


    //  SAVE RESULTS
    let json_string         =   serde_json::to_string(&poly_complex_facets ).unwrap();  
    let filepath            =   "/Users/gh10/a/c/pr/xh/pr/phfibre/tmp/s2_bars4.json";
    std::fs::write( &filepath, json_string ).expect("Unable to write file.");  

    let mut data = String::new();
    let mut file = std::fs::File::open( &filepath ).expect("Unable to open file");
    file.read_to_string( &mut data ).unwrap();
    let recovered: Vec< Polytope > = serde_json::from_str( &data ).unwrap();
    assert_eq!( &recovered, &poly_complex_facets );

    //  ANALYZE    
    let analyze_dowker_dual =   true;
    analyze_fibre( 
        &   poly_complex_facets,
            ring.clone(),
            analyze_dowker_dual,
    );    
}



// 2-SKELETON, BARCODE:(0, [0,INF)), (0, [1,2)), (1, [3,4)), (2, [5,INF))
// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 2, birth: 5 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }, BarFinite { dim: 1, birth: 3, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4, 5], val_to_ord: {0: 0, 2: 2, 3: 3, 4: 4, 5: 5, 1: 1} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 5296.645735629s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 162048
// number of facets (binned by dimension): [0, 0, 0, 0, 162048]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 162048]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 43.306967779s"
// number of polytopes (total): 1987032
// number of polytopes (binned by dimension): [87588, 422112, 743892, 571392, 162048]
// Killed: 9


// 2-SKELETON, BARCODE:(0, [0,INF)), (0, [1,2)), (1, [3,4)), (2, [5,INF))
// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 2, birth: 5 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }, BarFinite { dim: 1, birth: 3, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4, 5], val_to_ord: {2: 2, 1: 1, 3: 3, 0: 0, 4: 4, 5: 5} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 5150.821873563s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 162048
// number of facets (binned by dimension): [0, 0, 0, 0, 162048]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 162048]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 40.552123633s"
// number of polytopes (total): 1987032
// number of polytopes (binned by dimension): [87588, 422112, 743892, 571392, 162048]
// Killed: 9


// SUCCESS!!!
// --------------------------------------------------------------------------------------------- 
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 162048
// number of facets (binned by dimension): [0, 0, 0, 0, 162048]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 162048]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 39.787473234s"
// number of polytopes (total): 1987032
// number of polytopes (binned by dimension): [87588, 422112, 743892, 571392, 162048]
// Time elapsed to compute binary-coeff polyhedral betti numbers 50.653351456s
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 2, 27, 2, 0]
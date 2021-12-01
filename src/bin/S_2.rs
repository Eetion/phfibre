




use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use phfibre::pipelines::{simplex_pipeline, analyze_fibre};
use phfibre::phfibre::{ConditionNone, ConditionLowerStar, ConditionLowerEdge};
use phfibre::polytope::object_def::Polytope;
use solar::utilities::sequences_and_ordinals::ordinate_unique_vals;
use solar::utilities::cell_complexes::simplices_unweighted::facets::ordered_subsimplices_up_thru_dim_concatenated_vec; 
use std::iter::FromIterator;

use std::io::Read;





//  CODE

fn main() {

    let save_dir_opt        =   Some("/Users/gh10/a/c/pr/xh/pr/phfibre/tmp/s2_autosaves"); // we will not save any files    


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

    let poly_complex_facets =   simplex_pipeline(
                                    &   simplex_sequence,
                                    &   barcode,
                                    &   ring,
                                    &   precondition_to_make_new_lev_set_lower_none,
                                        false, // do not analyze the dowker dual to the nerve complex
                                        save_dir_opt,
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

    let analyze_dowker_dual =   false;
    analyze_fibre( 
        &   poly_complex_facets,
            ring.clone(),
            analyze_dowker_dual,
            save_dir_opt,
    );   
    
    




    
}


//  SPECIAL CASE (COMPUTED DOWKER NERVE FACET DIMENSIONS BUT MEMORY USE WAS 26 GB BEFORE I STOPPED THE PROCESS FOR COMPUTING BETTI NUMBERS)
// --------------------------------------------------------------------------------------------- 
// BASE SPACE
// simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
// BARCODE
// barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 2, birth: 5 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }, BarFinite { dim: 1, birth: 3, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4, 5], val_to_ord: {0: 0, 4: 4, 3: 3, 5: 5, 2: 2, 1: 1} } }
// TIME TO COMPUTE FIBRE FACETS
// Time elapsed to compute facets of PH fibre: 6510.383176986s

// ANALYSIS
// Each polytope facet has been checked for compatiblity with the given barcode.
// POLYHEDRAL COMPLEX FACETS
// number of facets (total): 162048
// number of facets (binned by dimension): [0, 0, 0, 0, 162048]
// number of facets (binned by number of vertices): [0, 0, 0, 0, 162048]
// POLYHEDRAL COMPLEX CELLS
// "Time elapsed to compute the faces of PH fibre (given the facets): 37.220890372s"
// number of polytopes (total): 1987032
// number of polytopes (binned by dimension): [87588, 422112, 743892, 571392, 162048]
// Time elapsed to compute binary-coeff polyhedral betti numbers 50.033824012s
// betti numbers (of polyhedral complex, Z2 coefficients): [1, 2, 27, 2, 0]
// DOWKER NERVE COMPLEX FACETS
// number of dowker nerve complex facets (total): 162048
// number of dowker nerve complex facets (binned by dimension): [0, 0, 0, 0, 1920, 0, 0, 25344, 35520, 0, 0, 82944, 0, 0, 0, 16320]
// DOWKER NERVE COMPLEX CELLS
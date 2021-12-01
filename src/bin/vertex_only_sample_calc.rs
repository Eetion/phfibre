//  DESCRIPTION

//  A square decomposed into two 2-simplices.


//  DEPENDENCIES

use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use phfibre::hasse_strategy::pipeline_vertex_only::{simplex_pipeline_vertex_only, analyze_fibre_vertex_only};
use phfibre::polytope::object_def::{Polytope};
use phfibre::hasse_strategy::compute_vertex_only::{ConditionNone, ConditionLowerStar, ConditionLowerEdge};
use solar::utilities::sequences_and_ordinals::ordinate_unique_vals;
use solar::utilities::cell_complexes::simplices_unweighted::facets::ordered_subsimplices_up_thru_dim_concatenated_vec; 

use std::io::Read;

fn main() {

    //  ----------------------------------------------------------------------------------------------    
    //  NO FINITE BARS
    //  ----------------------------------------------------------------------------------------------

    //  Define the base space, barcode, and ring

    //  SUBDIVIDED SQUARE
    // let complex_facets      =   vec![  vec![0,1,2], vec![1,2,3] ];
    // let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 2); 
    
    

    //  SOLID 3 SIMPLEX
    let complex_facets      =   vec![  vec![0, 1, 2, 3] ];
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 3); 

    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0} ],
                                    fin: vec![ ],
                                    ordinal: ordinate_unique_vals( & vec![ 0 ] ),
                                };

    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};      // this struct won't impose any extra conditions on the filtrations we build    

    let analyze_dowker_dual =   true;

    let poly_complex_facets =       simplex_pipeline_vertex_only(
                                        &   simplex_sequence,
                                        &   barcode,
                                        &   ring,
                                        &   precondition_to_make_new_lev_set_lower_none,
                                            analyze_dowker_dual
                                    );  
                                    
    let analyze_dowker_dual =   true;
    analyze_fibre_vertex_only( 
        &   poly_complex_facets,
            ring.clone(),
            analyze_dowker_dual,
            // save_dir_opt,
    );  

    let json_string =  serde_json::to_string(&poly_complex_facets ).unwrap();  
    let fp  =   "/Users/gh10/a/c/pr/xh/pr/phfibre/tmp/square_2triangles_stringformat.json";
    std::fs::write( &fp, json_string).expect("Unable to write file.");  

    let mut data = String::new();
    let mut file = std::fs::File::open( &fp ).expect("Unable to open file");
    file.read_to_string( &mut data ).unwrap();
    let recovered: Vec< Polytope > = serde_json::from_str( &data ).unwrap();
    assert_eq!( &recovered, &poly_complex_facets );


    let json_vec =  serde_json::to_vec(&poly_complex_facets ).unwrap();     
    let fp  =   "/Users/gh10/a/c/pr/xh/pr/phfibre/tmp/square_2triangles_u8format.json";
    std::fs::write( &fp, json_vec).expect("Unable to write file.");  

    let mut data = String::new();
    let mut file = std::fs::File::open( &fp ).expect("Unable to open file");
    file.read_to_string( &mut data ).unwrap();
    let recovered: Vec< Polytope > = serde_json::from_str( &data ).unwrap();
    assert_eq!( &recovered, &poly_complex_facets );

    // --------------------------------------------------------------------------------------------- 
    // BASE SPACE
    // simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3], [0, 1, 2, 3]]
    // BARCODE
    // barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0], val_to_ord: {0: 0} } }
    // TIME TO COMPUTE FIBRE FACETS
    // Time elapsed to compute facets of PH fibre: 101.113082685s
    
    // ANALYSIS
    // Each polytope facet has been checked for compatiblity with the given barcode.
    // POLYHEDRAL COMPLEX FACETS
    // number of facets (total): 64
    // number of facets (binned by dimension): [0, 64]
    // number of facets (binned by number of vertices): [0, 64]
    // POLYHEDRAL COMPLEX CELLS
    // "Time elapsed to compute the faces of PH fibre (given the facets): 814.773µs"
    // number of polytopes (total): 129
    // number of polytopes (binned by dimension): [65, 64]
    // Time elapsed to compute binary-coeff polyhedral betti numbers 23.275µs
    // betti numbers (of polyhedral complex, Z2 coefficients): [1, 0]
    // DOWKER NERVE COMPLEX FACETS
    // number of dowker nerve complex facets (total): 64
    // number of dowker nerve complex facets (binned by dimension): [0, 64]
    // DOWKER NERVE COMPLEX CELLS
    // Time elapsed to compute nerve dowker dual betti numbers (user specified coeff) 25.304µs
    // number of nerve dowker complex cells (total): 129
    // number of nerve dowker complex cells (binned by dimension): [65, 64]
    // betti numbers (of dowker nerve, user-specified ring coefficients): [1, 0]
}
//  DESCRIPTION

//  A square decomposed into two 2-simplices.


//  DEPENDENCIES

use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use phfibre::pipelines::{simplex_pipeline, analyze_fibre};
use phfibre::polytope::object_def::{Polytope};
use phfibre::phfibre::{ConditionNone, ConditionLowerStar, ConditionLowerEdge};
use solar::utilities::sequences_and_ordinals::ordinate_unique_vals;
use solar::utilities::cell_complexes::simplices_unweighted::facets::ordered_subsimplices_up_thru_dim_concatenated_vec;

use std::io::Read;

fn main() {

    //  ----------------------------------------------------------------------------------------------    
    //  NO FINITE BARS
    //  ----------------------------------------------------------------------------------------------

    //  Define the base space, barcode, and ring

    let complex_facets      =   vec![  vec![0,1,2], vec![1,2,3] ];

    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 2); 

    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0} ],
                                    fin: vec![ ],
                                    ordinal: ordinate_unique_vals( & vec![ 0 ] ),
                                };

    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};      // this struct won't impose any extra conditions on the filtrations we build    

    let analyze_dowker_dual =   true;

    let save_dir_opt        =   None;
    let poly_complex_facets =   simplex_pipeline(
                                    &   simplex_sequence,
                                    &   barcode,
                                    &   ring,
                                    &   precondition_to_make_new_lev_set_lower_none,
                                        analyze_dowker_dual,
                                        save_dir_opt,
                                );                                      

                                    
    let analyze_dowker_dual =   true;
    analyze_fibre( 
        &   poly_complex_facets,
            ring.clone(),
            analyze_dowker_dual,
            save_dir_opt,
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

    //  RESULTS
    //  -------

    // BASE SPACE
    // simplices of the base space: [[0], [1], [2], [3], [0, 1], [0, 2], [1, 2], [1, 3], [2, 3], [0, 1, 2], [1, 2, 3]]
    // POLYHEDRAL COMPLEX FACETS
    // number of facets (total): 160
    // number of facets (binned by dimension): [0, 0, 0, 0, 0, 160]
    // number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 160]
    // POLYHEDRAL COMPLEX CELLS
    // number of polytopes (total): 2731
    // number of polytopes (binned by dimension): [32, 241, 702, 964, 632, 160]
    // betti numbers (of polyhedral complex): [1, 0, 0, 0, 0, 0]
    // INTERSECTIONS
    // number of pairs of intersecting facets, binned by the dimension of the intersection polytope: [5472, 3872, 2120, 928, 328, 0]
    // DOWKER NERVE COMPLEX FACETS
    // number of dowker nerve complex facets (total): 160
    // number of dowker nerve complex facets (binned by dimension): [0, 0, 0, 0, 0, 160]
    // DOWKER NERVE COMPLEX CELLS
    // number of nerve dowker complex cells (total): 2731
    // number of nerve dowker complex cells (binned by dimension): [32, 241, 702, 964, 632, 160]
    // betti numbers (of dowker nerve): [1, 0, 0, 0, 0, 0]
}
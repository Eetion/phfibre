use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use phfibre::pipelines::{simplex_pipeline, analyze_fibre};
use phfibre::phfibre::{ConditionNone, ConditionLowerStar, ConditionLowerEdge};
use phfibre::polytope::object_def::Polytope;
use solar::utilities::sequences_and_ordinals::ordinate_unique_vals;
use solar::utilities::cell_complexes::simplices_unweighted::facets::ordered_subsimplices_up_thru_dim_concatenated_vec; 
use std::iter::FromIterator;

use std::io::Read;



//  THE PURPOSE OF THIS FILE IS TO EXPOSE THE REASON WHY RUST CONSISTENLY CRASHES WHILE COMPUTING BETTI STATS FOR THIS POLYHEDRAL COMPLEX

fn main() {

    let save_dir_opt        =   None; // we will not save any files    


    //  LOAD FACETS
    let filepath                =   "/Users/gh10/a/c/pr/xh/pr/phfibre/tmp/s2_bars4.json";
    let mut data = String::new();
    let mut file = std::fs::File::open( &filepath ).expect("Unable to open file");
    file.read_to_string( &mut data ).unwrap();
    let poly_complex_facets: Vec< Polytope > = serde_json::from_str( &data ).unwrap();
    
    //  ANALYZE

    let ring                    =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();
    let analyze_dowker_dual     =   false;

    analyze_fibre(
        &   poly_complex_facets,
            ring.clone(),
            analyze_dowker_dual,
            save_dir_opt,
    );     

}
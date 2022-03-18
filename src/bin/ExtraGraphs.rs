use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use phfibre::pipelines::{simplex_pipeline, analyze_fibre};
use phfibre::phfibre::{ConditionNone, ConditionLowerStar, ConditionLowerEdge};
use phfibre::polytope::object_def::Polytope;
use solar::utilities::sequences_and_ordinals::{ordinate_unique_vals, BiMapSequential};
use solar::utilities::cell_complexes::simplices_unweighted::facets::ordered_subsimplices_up_thru_dim_concatenated_vec;
use std::iter::FromIterator;

use std::io::Read;





//  CODE

fn main() {

  let save_dir_opt        = None;

    println!("\n  THIS IS A SUPER FAST TEST : TRIANGLE, 1 FINITE BAR IN DIM 0 ");
    //  ----------------------------------------------------------------------------------------------
    //  TRIANGLE, 1 FINITE BAR IN DIM 0
    //  ----------------------------------------------------------------------------------------------

    //  Define the base space, barcode, and ring


    let complex_facets      =   vec![  vec![0, 1, 2] ];
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 1);

    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:3} ],
                                    fin: vec![ BarFinite{dim:0,birth:1,death:2} ],
                                    ordinal: ordinate_unique_vals( & vec![ 0, 1, 2, 3 ] ),
                                };

    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};      // this struct won't impose any extra conditions on the filtrations we build

    let poly_complex_facets =   simplex_pipeline(
                                    &   simplex_sequence,
                                    &   barcode,
                                    &   ring,
                                    &   precondition_to_make_new_lev_set_lower_none,
                                        false, // do not analyze the dowker dual to the nerve complex
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



        println!("\n THIS IS A SLIGHTLY SLOWER TEST: TWO TRIANGLES JOINED AT A VERTEX  AND NO FINITE BAR");
        println!("===============================================================================================\n");
        //  Define the base space, barcode, and ring

        let complex_facets      =   vec![  vec![0,1,2],vec![2,3,4] ];
        let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 1);    // build all faces up to dimension 1

        let barcode             =   Barcode{
            inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:1},BarInfinite{dim:1,birth:1} ],
            fin: vec![ ],
            ordinal: ordinate_unique_vals( & vec![0, 1] ),
        };

        let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

        let simplex_bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence.clone() );

        //  Define any extra conditions on starting a new level set

        let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};
        let precondition_to_make_new_lev_set_lower_star     =   ConditionLowerStar{ simplex_bimap_sequential: & simplex_bimap_sequential };
        let precondition_to_make_new_lev_set_lower_edge     =   ConditionLowerEdge{ simplex_bimap_sequential: & simplex_bimap_sequential };

        let analyze_dowker_dual                             =   false;


        //  UNCONSTRAINED FILTRATION
        //  -------------------------------------------------------------------------------------------

        println!("\nUNCONSTRAINED FILTRATIONS ");
        let poly_complex_facets =       simplex_pipeline(
                                             &   simplex_sequence,
                                            &   barcode,
                                          &   ring,
                                       &   precondition_to_make_new_lev_set_lower_none,
                                               analyze_dowker_dual,
                                               save_dir_opt
                                         );

        let analyze_dowker_dual =   false;
        analyze_fibre(
           &   poly_complex_facets,
                ring.clone(),
                analyze_dowker_dual,
                save_dir_opt,
        );


         //  LOWER STAR
        //  -------------------------------------------------------------------------------------------

        println!("\n LOWER STARS");
        let poly_complex_facets =       simplex_pipeline(
                                             &   simplex_sequence,
                                            &   barcode,
                                          &   ring,
                                       &   precondition_to_make_new_lev_set_lower_star,
                                               analyze_dowker_dual,
                                               save_dir_opt
                                         );

        let analyze_dowker_dual =   false;
        analyze_fibre(
           &   poly_complex_facets,
                ring.clone(),
                analyze_dowker_dual,
                save_dir_opt,
        );





}

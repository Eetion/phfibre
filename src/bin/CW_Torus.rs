//  DESCRIPTION

//  CW Torus




//  DEPENDENCIES

use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
// use phfibre::pipelines::{simplex_pipeline, analyze_fibre};
use phfibre::pipelines::boundary_matrix_pipeline;
use phfibre::phfibre::{ConditionNone, ConditionLowerStar, ConditionLowerEdge};
use solar::utilities::sequences_and_ordinals::{ordinate_unique_vals, BiMapSequential};
use solar::utilities::cell_complexes::simplices_unweighted::facets::ordered_subsimplices_up_thru_dim_concatenated_vec;
// use std::iter::FromIterator;
use solar::rings::ring::{Ring, Semiring};

/*    OLD DEPENDENCIES
use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use phfibre::pipelines::boundary_matrix_pipeline;
use phfibre::phfibre::ConditionNone;
use solar::utilities::sequences_and_ordinals::ordinate_unique_vals;
use solar::utilities::cell_complexes::simplices_unweighted::facets::ordered_subsimplices_up_thru_dim_concatenated_vec;
use solar::rings::ring::{Ring, Semiring};   */

//  CODE

fn main() {

  let save_dir_opt        = None;


    //  ----------------------------------------------------------------------------------------------
    //  DEFINE COEFFICIENT RING
    //  ----------------------------------------------------------------------------------------------

    //  This is the ring of rationals.
    type CeofficientRing    =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >; // the TYPE of ring object
    let ring                =   CeofficientRing::new(); // the object itself




    //  ----------------------------------------------------------------------------------------------
    //  TWO TORUS -- 6 cells, no finite bars
    //  ----------------------------------------------------------------------------------------------

    //  BOUNDARY MATRIX

    //  Define +1 and -1
    let one_pos             =   CeofficientRing::one();
    let one_neg             =   ring.negate( one_pos.clone() );

    /*

    //  The boundary matrix in vector-of-vector format; each 'internal'vector represents a compressed column of the boundary matrix
    let boundary            =   vec![
                                    vec![                               ],                          // 0-cell
                                    vec![                               ],                          // 0-cell
                                    vec![   ( 0, one_neg.clone() ),      ( 1, one_pos.clone() ) ],  // 1-cell
                                    vec![   ( 0, one_neg.clone() ),      ( 1, one_pos.clone() ) ],  // 1-cell
                                    vec![                               ],                          // 1-cell
                                    vec![                               ]                           // 2-cell
                                ];

    //  DIMENSION OF EACH CELL
    let cell_dims           =   vec![ 0, 0, 1, 1, 1, 2 ];

    //  THE FACES OF EACH CELL (SINCE THE CW COMPLEX IS NOT REGULAR, WE HAVE TO SPECIFY THESE EXPLICITLY)
    let cell_id_to_prereqs  =   vec![
                                    vec![               ],
                                    vec![               ],
                                    vec![ 0, 1,         ],
                                    vec![ 0, 1,         ],
                                    vec![ 0,            ],
                                    vec![ 0, 1, 2, 3, 4 ],
                                ];
    //  BARCODE
    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:1},  BarInfinite{dim:1,birth:2},  BarInfinite{dim:2,birth:3} ],
                                    fin: vec![ ],
                                    ordinal: ordinate_unique_vals( & vec![0, 1, 2, 3] ),
                                };

    //  EXTRA CONDITIONS ON THE FILTRATIONS WE BUILD
    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};

    //  CONFIRM THAT WE WILL ANALYZE THE DUAL DOWKER COMPLEX TO THE NERVE
    let analyze_dowker_dual =   false;

    //  ANALYZE

    println!("");
    println!("T2 (six cells):");
    boundary_matrix_pipeline(
        &   boundary,
        Some( & cell_id_to_prereqs ),
        &   cell_dims,
        &   barcode,
        &   ring,
        &   precondition_to_make_new_lev_set_lower_none,
            analyze_dowker_dual,
            save_dir_opt,
    );

    //  RESULTS
    // ---------------------------------------------------------------------------------------------
    // BOUNDARY MATRIX
    // boundary matrix sparse columns: [[], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [], []]
    // BARCODE
    // barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 1 }, BarInfinite { dim: 1, birth: 2 }, BarInfinite { dim: 2, birth: 3 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3], val_to_ord: {1: 1, 3: 3, 0: 0, 2: 2} } }
    // TIME TO COMPUTE FIBRE FACETS
    // Time elapsed to compute facets of PH fibre: 236.041µs

    // ANALYSIS
    // Each polytope facet has been checked for compatiblity with the given barcode.
    // POLYHEDRAL COMPLEX FACETS
    // number of facets (total): 12
    // number of facets (binned by dimension): [0, 12]
    // number of facets (binned by number of vertices): [0, 12]
    // POLYHEDRAL COMPLEX CELLS
    // "Time elapsed to compute facets of PH fibre: 164.349µs"
    // number of polytopes (total): 24
    // number of polytopes (binned by dimension): [12, 12]
    // betti numbers (of polyhedral complex, Z2 coefficients): [2, 2]
    // DOWKER NERVE COMPLEX FACETS
    // number of dowker nerve complex facets (total): 12
    // number of dowker nerve complex facets (binned by dimension): [0, 12]
    // DOWKER NERVE COMPLEX CELLS
    // number of nerve dowker complex cells (total): 24
    // number of nerve dowker complex cells (binned by dimension): [12, 12]
    // betti numbers (of dowker nerve, user-specified ring coefficients): [2, 2]





    //  ----------------------------------------------------------------------------------------------
    //  TWO TORUS -- 4 cells
    //  ----------------------------------------------------------------------------------------------

    //  Define the boundary matrix

    //  The boundary matrix in vector-of-vector format; each 'internal'vector represents a compressed column of the boundary matrix
    let boundary            =   vec![
                                    vec![                               ],
                                    vec![                               ],
                                    vec![                               ],
                                    vec![                               ]
                                ];

    //  The dimension of each cell.
    let cell_dims           =   vec![ 0, 1, 1, 2 ];

    //  Definition: cell_id_to_prereqs[ m ] = { the cells that must be added to the complex before M }
    let cell_id_to_prereqs  =   vec![
                                    vec![           ],
                                    vec![ 0,        ],
                                    vec![ 0,        ],
                                    vec![ 0, 1, 2,  ]
                                ];

    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:1},  BarInfinite{dim:1,birth:2},  BarInfinite{dim:2,birth:3} ],
                                    fin: vec![ ],
                                    ordinal: ordinate_unique_vals( & vec![0, 1, 2, 3] ),
                                };

    println!("");
    println!("T2 (four cells):");
    boundary_matrix_pipeline(
        &   boundary,
        Some( & cell_id_to_prereqs ),
        &   cell_dims,
        &   barcode,
        &   ring,
        &   precondition_to_make_new_lev_set_lower_none,
            analyze_dowker_dual,
            save_dir_opt,
    );

    //  RESULTS
    // ---------------------------------------------------------------------------------------------
    // BOUNDARY MATRIX
    // boundary matrix sparse columns: [[], [], [], []]
    // BARCODE
    // barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 1 }, BarInfinite { dim: 1, birth: 2 }, BarInfinite { dim: 2, birth: 3 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3], val_to_ord: {0: 0, 2: 2, 3: 3, 1: 1} } }
    // TIME TO COMPUTE FIBRE FACETS
    // Time elapsed to compute facets of PH fibre: 40.27µs

    // ANALYSIS
    // Each polytope facet has been checked for compatiblity with the given barcode.
    // POLYHEDRAL COMPLEX FACETS
    // number of facets (total): 2
    // number of facets (binned by dimension): [2]
    // number of facets (binned by number of vertices): [2]
    // POLYHEDRAL COMPLEX CELLS
    // "Time elapsed to compute facets of PH fibre: 35.272µs"
    // number of polytopes (total): 2
    // number of polytopes (binned by dimension): [2]
    // betti numbers (of polyhedral complex, Z2 coefficients): [2]
    // DOWKER NERVE COMPLEX FACETS
    // number of dowker nerve complex facets (total): 2
    // number of dowker nerve complex facets (binned by dimension): [2]
    // DOWKER NERVE COMPLEX CELLS
    // number of nerve dowker complex cells (total): 2
    // number of nerve dowker complex cells (binned by dimension): [2]
    // betti numbers (of dowker nerve, user-specified ring coefficients): [2]

     //  ----------------------------------------------------------------------------------------------
    //  TWO TORUS -- 4 vertices, 8 edges, 4 faces
     ///////               0---------1---------0
     ///                   |         |         |
     ///                   |         |         |
     ///                   |         |         |
     ///                   2---------3---------2
     ///                   |         |         |
     ///                   |         |         |
     ///                   |         |         |
     ///                   0---------1---------0
     ///

    //  ----------------------------------------------------------------------------------------------
    println!("");
    println!("Torus (CW Structure with four squares, so 16 cells), and barcode with no finite bars:");
    //  Define the boundary matrix
*/
    let analyze_dowker_dual =   false;
    //  The boundary matrix in vector-of-vector format; each 'internal'vector represents a compressed column of the boundary matrix
    let boundary            =   vec![
                                    vec![                               ],                          //0:  0-cell 0
                                    vec![                               ],                          //1:  0-cell 1
                                    vec![                               ],                          //2:  0-cell 2
                                    vec![                               ],                          //3:  0-cell 3
                                    vec![   ( 0, one_neg.clone() ),      ( 1, one_pos.clone() ) ],  //4:  1-cell  0->1
                                    vec![   ( 0, one_pos.clone() ),      ( 1, one_neg.clone() ) ],  //5:  1-cell  1->0
                                    vec![   ( 0, one_neg.clone() ),      ( 2, one_pos.clone() ) ],   //6:  1-cell  0->2
                                    vec![   ( 0, one_pos.clone() ),      ( 2, one_neg.clone() ) ],   //7:  1-cell  2->0
                                    vec![   ( 2, one_neg.clone() ),      ( 3, one_pos.clone() ) ],  //8:  1-cell  2->3
                                    vec![   ( 2, one_pos.clone() ),      ( 3, one_neg.clone() ) ],  //9:  1-cell  3->2
                                    vec![   ( 1, one_neg.clone() ),      ( 3, one_pos.clone() ) ],   //10: 1-cell  1->3
                                    vec![   ( 1, one_pos.clone() ),      ( 3, one_neg.clone() ) ],   //11: 1-cell  3->1
vec![ ( 4, one_pos.clone() ),  ( 6, one_neg.clone()),  ( 8, one_neg.clone() ), ( 10, one_pos.clone())], //12: top left cell
vec![ ( 5, one_pos.clone() ),  ( 6, one_pos.clone()),   ( 9, one_neg.clone() ), ( 10, one_neg.clone())], //13: top right cell                                    vec![   ( 1, one_pos.clone() ),      ( 4, one_neg.clone() ) ],   // 1-cell  3->1
vec![ ( 4, one_neg.clone() ),  ( 7, one_neg.clone()),  ( 8, one_pos.clone() ), ( 11, one_pos.clone())], //14: bottom left cell
vec![ ( 5, one_neg.clone() ),  ( 7, one_pos.clone()),   ( 9, one_pos.clone() ), ( 11, one_neg.clone())] //15: bottom right cell
                                ];

    println!("Built boundary matrix");

    //  The dimension of each cell.
    let cell_dims           =   vec![ 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2 ];

    //  Definition: cell_id_to_prereqs[ m ] = { the cells that must be added to the complex before M }
    let cell_id_to_prereqs  =   vec![
                                    vec![           ], // 0
                                    vec![           ],
                                    vec![           ],
                                    vec![           ], //3
                                    vec![ 0, 1       ],
                                    vec![ 0, 1       ],
                                    vec![ 0, 2       ],
                                    vec![ 0, 2       ], //7
                                    vec![ 2, 3       ],
                                    vec![ 2, 3       ],
                                    vec![ 1, 3       ],
                                    vec![ 1, 3       ], //11
                                    vec![ 4, 6, 8,10  ],
                                    vec![ 5, 6, 9,10  ],
                                    vec![ 4, 7, 8,11  ],
                                    vec![ 5, 7, 9,11  ]
                                ];

   /*
    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:1},  BarInfinite{dim:1,birth:1},  BarInfinite{dim:2,birth:1} ],
                                    fin: vec![ ],
                                    ordinal: ordinate_unique_vals( & vec![0, 1] ),
                                };

     //  ALL Filters
     //  -------------------------------------------------------------------------------------------

     println!("\nAll filters");
     let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};
    boundary_matrix_pipeline(
        &   boundary,
        Some( & cell_id_to_prereqs ),
        &   cell_dims,
        &   barcode,
        &   ring,
        &   precondition_to_make_new_lev_set_lower_none ,
            analyze_dowker_dual,
            save_dir_opt,
    );


    println!("-------------------------------------------------------------------------------------");
    println!("Torus (CW Structure with four squares, so 16 cells), and barcode with no finite bars, first two loops and then cavity:");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:1},  BarInfinite{dim:1,birth:1},  BarInfinite{dim:2,birth:2} ],
        fin: vec![ ],
        ordinal: ordinate_unique_vals( & vec![0, 1, 2] ),
    };

    //  ALL Filters
    //  -------------------------------------------------------------------------------------------

    println!("\nAll filters");
    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};
    boundary_matrix_pipeline(
    &   boundary,
    Some( & cell_id_to_prereqs ),
    &   cell_dims,
    &   barcode,
    &   ring,
    &   precondition_to_make_new_lev_set_lower_none ,
    analyze_dowker_dual,
    save_dir_opt,
    );

    println!("-------------------------------------------------------------------------------------");
    println!("Torus (CW Structure with four squares, so 12 cells), and barcode with no finite bars, first a loop and then a loop and cavity:");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:1},  BarInfinite{dim:1,birth:2},  BarInfinite{dim:2,birth:2} ],
        fin: vec![ ],
        ordinal: ordinate_unique_vals( & vec![0, 1, 2] ),
    };

    //  ALL Filters
    //  -------------------------------------------------------------------------------------------

    println!("\nAll filters");
    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};
    boundary_matrix_pipeline(
    &   boundary,
    Some( & cell_id_to_prereqs ),
    &   cell_dims,
    &   barcode,
    &   ring,
    &   precondition_to_make_new_lev_set_lower_none ,
    analyze_dowker_dual,
    save_dir_opt,
    );


    println!("-------------------------------------------------------------------------------------");
    println!("Torus (CW Structure with four squares, so 16 cells), and barcode with no finite bars, distinct birth times:");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:1},  BarInfinite{dim:1,birth:2},  BarInfinite{dim:2,birth:3} ],
        fin: vec![ ],
        ordinal: ordinate_unique_vals( & vec![0, 1, 2, 3] ),
    };

    //  ALL Filters
    //  -------------------------------------------------------------------------------------------

    println!("\nAll filters");
    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};
    boundary_matrix_pipeline(
    &   boundary,
    Some( & cell_id_to_prereqs ),
    &   cell_dims,
    &   barcode,
    &   ring,
    &   precondition_to_make_new_lev_set_lower_none ,
    analyze_dowker_dual,
    save_dir_opt,
    );

    println!("-------------------------------------------------------------------------------------");
    println!("Torus (CW Structure with four squares, so 16 cells), one dim-0 finite bar, infinite bars same birth times:");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:3},  BarInfinite{dim:1,birth:3},  BarInfinite{dim:2,birth:3} ],
        fin: vec![ BarFinite{dim:0,birth:1,death:2}],
        ordinal: ordinate_unique_vals( & vec![0, 1, 2, 3] ),
    };

    //  ALL Filters
    //  -------------------------------------------------------------------------------------------

    println!("\nAll filters");
    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};
    boundary_matrix_pipeline(
    &   boundary,
    Some( & cell_id_to_prereqs ),
    &   cell_dims,
    &   barcode,
    &   ring,
    &   precondition_to_make_new_lev_set_lower_none ,
    analyze_dowker_dual,
    save_dir_opt,
    );

    println!("-------------------------------------------------------------------------------------");
    println!("Torus (CW Structure with four squares, so 16 cells), one dim-1 finite bar, infinite bars same birth times:");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:3},  BarInfinite{dim:1,birth:3},  BarInfinite{dim:2,birth:3} ],
        fin: vec![ BarFinite{dim:1,birth:1,death:2}],
        ordinal: ordinate_unique_vals( & vec![0, 1, 2, 3] ),
    };

    //  ALL Filters
    //  -------------------------------------------------------------------------------------------

    println!("\nAll filters");
    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};
    boundary_matrix_pipeline(
    &   boundary,
    Some( & cell_id_to_prereqs ),
    &   cell_dims,
    &   barcode,
    &   ring,
    &   precondition_to_make_new_lev_set_lower_none ,
    analyze_dowker_dual,
    save_dir_opt,
    );


    println!("-------------------------------------------------------------------------------------");
    println!("Torus (CW Structure with four squares, so 16 cells), one dim-0 finite bar and then one dim-0 finite bar, infinite bars same birth times:");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:5},  BarInfinite{dim:1,birth:5},  BarInfinite{dim:2,birth:5} ],
        fin: vec![ BarFinite{dim:0,birth:1,death:2}, BarFinite{dim:0,birth:3,death:4}],
        ordinal: ordinate_unique_vals( & vec![0, 1, 2, 3,4,5] ),
    };

    //  ALL Filters
    //  -------------------------------------------------------------------------------------------

    println!("\nAll filters");
    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};
    boundary_matrix_pipeline(
    &   boundary,
    Some( & cell_id_to_prereqs ),
    &   cell_dims,
    &   barcode,
    &   ring,
    &   precondition_to_make_new_lev_set_lower_none ,
    analyze_dowker_dual,
    save_dir_opt,
    );

    println!("-------------------------------------------------------------------------------------");
    println!("Torus (CW Structure with four squares, so 16 cells), one dim-1 finite bar and then one dim-1 finite bar, infinite bars same birth times:");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:5},  BarInfinite{dim:1,birth:5},  BarInfinite{dim:2,birth:5} ],
        fin: vec![ BarFinite{dim:1,birth:1,death:2}, BarFinite{dim:1,birth:3,death:4}],
        ordinal: ordinate_unique_vals( & vec![0, 1, 2, 3,4,5] ),
    };

    //  ALL Filters
    //  -------------------------------------------------------------------------------------------

    println!("\nAll filters");
    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};
    boundary_matrix_pipeline(
    &   boundary,
    Some( & cell_id_to_prereqs ),
    &   cell_dims,
    &   barcode,
    &   ring,
    &   precondition_to_make_new_lev_set_lower_none ,
    analyze_dowker_dual,
    save_dir_opt,
    );
    */
    println!("-------------------------------------------------------------------------------------");
    println!("Torus (CW Structure with four squares, so 16 cells), one dim-0 finite bar and then one dim-1 finite bar, infinite bars same birth times:");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:5},  BarInfinite{dim:1,birth:5},  BarInfinite{dim:2,birth:5} ],
        fin: vec![ BarFinite{dim:0,birth:1,death:2}, BarFinite{dim:1,birth:3,death:4}],
        ordinal: ordinate_unique_vals( & vec![0, 1, 2, 3,4,5] ),
    };

    //  ALL Filters
    //  -------------------------------------------------------------------------------------------

    println!("\nAll filters");
    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};
    boundary_matrix_pipeline(
    &   boundary,
    Some( & cell_id_to_prereqs ),
    &   cell_dims,
    &   barcode,
    &   ring,
    &   precondition_to_make_new_lev_set_lower_none ,
    analyze_dowker_dual,
    save_dir_opt,
    );

    println!("-------------------------------------------------------------------------------------");
    println!("Torus (CW Structure with four squares, so 16 cells), one loop and then one dim-0 finite bar");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:1},  BarInfinite{dim:1,birth:4},  BarInfinite{dim:2,birth:4} ],
        fin: vec![ BarFinite{dim:0,birth:2,death:3}],
        ordinal: ordinate_unique_vals( & vec![0, 1, 2, 3,4] ),
    };

    //  ALL Filters
    //  -------------------------------------------------------------------------------------------

    println!("\nAll filters");
    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};
    boundary_matrix_pipeline(
    &   boundary,
    Some( & cell_id_to_prereqs ),
    &   cell_dims,
    &   barcode,
    &   ring,
    &   precondition_to_make_new_lev_set_lower_none ,
    analyze_dowker_dual,
    save_dir_opt,
    );

    println!("-------------------------------------------------------------------------------------");
    println!("Torus (CW Structure with four squares, so 16 cells), dim-0 finite bar overlapping one loop");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:2},  BarInfinite{dim:1,birth:4},  BarInfinite{dim:2,birth:4} ],
        fin: vec![ BarFinite{dim:0,birth:1,death:3}],
        ordinal: ordinate_unique_vals( & vec![0, 1, 2, 3,4] ),
    };

    //  ALL Filters
    //  -------------------------------------------------------------------------------------------

    println!("\nAll filters");
    let precondition_to_make_new_lev_set_lower_none     =   ConditionNone{};
    boundary_matrix_pipeline(
    &   boundary,
    Some( & cell_id_to_prereqs ),
    &   cell_dims,
    &   barcode,
    &   ring,
    &   precondition_to_make_new_lev_set_lower_none ,
    analyze_dowker_dual,
    save_dir_opt,
    );

}


/*
Torus (CW Structure with four squares, so 16 cells), and barcode with no finite bars
ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 256
number of facets (binned by dimension): [0, 0, 0, 128, 128]
number of facets (binned by number of vertices): [0, 0, 0, 128, 128]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 7.10876ms"
number of polytopes (total): 2048
number of polytopes (binned by dimension): [64, 448, 832, 576, 128]
Time elapsed to compute binary-coeff polyhedral betti numbers 1.176438ms
betti numbers (of polyhedral complex, Z2 coefficients): [1, 2, 1, 0, 0]
DOWKER NERVE COMPLEX FACETS
number of dowker nerve complex facets (total): 256
number of dowker nerve complex facets (binned by dimension): [0, 0, 0, 128, 128]
DOWKER NERVE COMPLEX CELLS
Time elapsed to compute nerve dowker dual betti numbers (user specified coeff) 1.943013ms
number of nerve dowker complex cells (total): 2048
number of nerve dowker complex cells (binned by dimension): [64, 448, 832, 576, 128]
betti numbers (of dowker nerve, user-specified ring coefficients): [1, 2, 1, 0, 0]


Torus (CW Structure with four squares, so 16 cells), and barcode with no finite bars, first two loops and then cavity:

All filters

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 15360
number of facets (binned by dimension): [0, 0, 0, 0, 0, 0, 15360]
number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 0, 15360]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 4.919161391s"
number of polytopes (total): 482432
number of polytopes (binned by dimension): [1824, 19776, 77344, 146688, 146688, 74752, 15360]
Time elapsed to compute binary-coeff polyhedral betti numbers 1.797527276s
betti numbers (of polyhedral complex, Z2 coefficients): [1, 4, 5, 2, 0, 0, 0]
DOWKER NERVE COMPLEX FACETS
number of dowker nerve complex facets (total): 15360
number of dowker nerve complex facets (binned by dimension): [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4096, 11264]
DOWKER NERVE COMPLEX CELLS
/var/spool/slurmd/job1042215/slurm_script: line 21: 1613604 Killed                  cargo run --bin CW_Torus --release
slurmstepd: error: Detected 1 oom-kill event(s) in StepId=1042215.batch cgroup. Some of your processes may have been killed by the cgroup out-of-memory handler.

-------------------------------------------------------------------------------------
Torus (CW Structure with four squares, so 12 cells), and barcode with no finite bars, first a loop and then a loop and cavity:

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: -1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (2, Ratio { numer: -1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: 1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: 1, denom: 1 }), (9, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 1 }, BarInfinite { dim: 1, birth: 2 }, BarInfinite { dim: 2, birth: 2 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0, 1, 2], val_to_ord: {2: 2, 0: 0, 1: 1} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 924.858514925s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 3008
number of facets (binned by dimension): [0, 0, 0, 448, 0, 2560]
number of facets (binned by number of vertices): [0, 0, 0, 448, 0, 2560]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 380.570567ms"
number of polytopes (total): 66688
number of polytopes (binned by dimension): [1472, 8928, 20096, 21856, 11776, 2560]
Time elapsed to compute binary-coeff polyhedral betti numbers 83.905529ms
betti numbers (of polyhedral complex, Z2 coefficients): [6, 8, 2, 0, 0, 0]
-------------------------------------------------------------------------------------

-------------------------------------------------------------------------------------
Torus (CW Structure with four squares, so 16 cells), and barcode with no finite bars, distinct birth times:

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: $
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 1 }, BarInfinite { dim: 1, birth: 2 }, BarInfinite { dim: 2, birth: 3 }], fin: [], ordina$
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 16320.048401729s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 74752
number of facets (binned by dimension): [0, 0, 0, 0, 0, 0, 74752]
number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 0, 74752]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 34.659768234s"
number of polytopes (total): 2770176
number of polytopes (binned by dimension): [18496, 153856, 505280, 852352, 786560, 378880, 74752]
Time elapsed to compute binary-coeff polyhedral betti numbers 17.760308247s
betti numbers (of polyhedral complex, Z2 coefficients): [6, 14, 10, 2, 0, 0, 0]
-------------------------------------------------------------------------------------
Torus (CW Structure with four squares, so 16 cells), one dim-0 finite bar, infinite bars same birth times:

All filters

NEW COMPUTATION

BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: $
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 3 }, BarInfinite { dim: 1, birth: 3 }, BarInfinite { dim: 2, birth: 3 }], fin: [BarFinite$
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 7041.606642676s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 2816
number of facets (binned by dimension): [0, 0, 1408, 1408]
number of facets (binned by number of vertices): [0, 0, 1408, 1408]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 66.64797ms"
number of polytopes (total): 16400
number of polytopes (binned by dimension): [2368, 6800, 5824, 1408]
Time elapsed to compute binary-coeff polyhedral betti numbers 8.677539ms
betti numbers (of polyhedral complex, Z2 coefficients): [1, 17, 0, 0]
-------------------------------------------------------------------------------------
Torus (CW Structure with four squares, so 16 cells), one dim-1 finite bar, infinite bars same birth times:

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: $
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 3 }, BarInfinite { dim: 1, birth: 3 }, BarInfinite { dim: 2, birth: 3 }], fin: [BarFinite$
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 12.678805184s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 128
number of facets (binned by dimension): [0, 0, 0, 128]
number of facets (binned by number of vertices): [0, 0, 0, 128]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 2.779594ms"
number of polytopes (total): 768
number of polytopes (binned by dimension): [64, 256, 320, 128]
Time elapsed to compute binary-coeff polyhedral betti numbers 241.938  s
betti numbers (of polyhedral complex, Z2 coefficients): [4, 4, 0, 0]
-------------------------------------------------------------------------------------

Torus (CW Structure with four squares, so 16 cells), one dim-0 finite bar and then one dim-0 finite bar, infinite bars same birth times:

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: $
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 5 }, BarInfinite { dim: 1, birth: 5 }, BarInfinite { dim: 2, birth: 5 }], fin: [BarFinite$
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 2260.492317076s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 1792
number of facets (binned by dimension): [0, 896, 896]
number of facets (binned by number of vertices): [0, 896, 896]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 22.390401ms"
number of polytopes (total): 6080
number of polytopes (binned by dimension): [2112, 3072, 896]
Time elapsed to compute binary-coeff polyhedral betti numbers 1.554348ms
betti numbers (of polyhedral complex, Z2 coefficients): [2, 66, 0]
-------------------------------------------------------------------------------------
Torus (CW Structure with four squares, so 16 cells), one dim-1 finite bar and then one dim-1 finite bar, infinite bars same birth times:

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: $
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 5 }, BarInfinite { dim: 1, birth: 5 }, BarInfinite { dim: 2, birth: 5 }], fin: [BarFinite$
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 1.576778368s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 0
number of facets (binned by dimension): []
number of facets (binned by number of vertices): []
POLYHEDRAL COMPLEX CELLS
thread 'main' panicked at 'called `Option::unwrap()` on a `None` value', src/polytope/faces.rs:482:42
note: run with `RUST_BACKTRACE=1` environment variable to display a backtrace

Built boundary matrix
-------------------------------------------------------------------------------------
Torus (CW Structure with four squares, so 16 cells), one dim-0 finite bar and then one dim-1 finite bar, infinite bars same birth times:

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: -1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (2, Ratio { numer: -1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: 1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: 1, denom: 1 }), (9, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 5 }, BarInfinite { dim: 1, birth: 5 }, BarInfinite { dim: 2, birth: 5 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }, BarFinite { dim: 1, birth: 3, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4, 5], val_to_ord: {1: 1, 3: 3, 4: 4, 5: 5, 0: 0, 2: 2} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 124.3586993s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 1408
number of facets (binned by dimension): [0, 0, 1408]
number of facets (binned by number of vertices): [0, 0, 1408]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 24.057915ms"
number of polytopes (total): 5120
number of polytopes (binned by dimension): [1152, 2560, 1408]
Time elapsed to compute binary-coeff polyhedral betti numbers 1.490355ms
betti numbers (of polyhedral complex, Z2 coefficients): [8, 8, 0]

Torus (CW Structure with four squares, so 16 cells), one loop and then one dim-0 finite bar

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: -1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (2, Ratio { numer: -1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: 1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: 1, denom: 1 }), (9, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 1 }, BarInfinite { dim: 1, birth: 4 }, BarInfinite { dim: 2, birth: 4 }], fin: [BarFinite { dim: 0, birth: 2, death: 3 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4], val_to_ord: {4: 4, 2: 2, 0: 0, 1: 1, 3: 3} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 128.502780979s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 1856
number of facets (binned by dimension): [0, 0, 320, 0, 1536]
number of facets (binned by number of vertices): [0, 0, 320, 0, 1536]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 159.175797ms"
number of polytopes (total): 28416
number of polytopes (binned by dimension): [2240, 7872, 10432, 6336, 1536]
Time elapsed to compute binary-coeff polyhedral betti numbers 25.492469ms
betti numbers (of polyhedral complex, Z2 coefficients): [8, 16, 8, 0, 0]

Torus (CW Structure with four squares, so 16 cells), dim-0 finite bar overlapping one loop

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: -1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (2, Ratio { numer: -1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: 1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: 1, denom: 1 }), (9, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 2 }, BarInfinite { dim: 1, birth: 4 }, BarInfinite { dim: 2, birth: 4 }], fin: [BarFinite { dim: 0, birth: 1, death: 3 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4], val_to_ord: {0: 0, 1: 1, 3: 3, 2: 2, 4: 4} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 1802.791413147s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 6784
number of facets (binned by dimension): [0, 0, 896, 0, 5888]
number of facets (binned by number of vertices): [0, 0, 896, 0, 5888]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 612.335611ms"
number of polytopes (total): 108288
number of polytopes (binned by dimension): [8320, 29824, 39936, 24320, 5888]
Time elapsed to compute binary-coeff polyhedral betti numbers 149.793679ms
betti numbers (of polyhedral complex, Z2 coefficients): [16, 32, 16, 0, 0]

*/

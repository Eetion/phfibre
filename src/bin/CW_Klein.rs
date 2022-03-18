//  DESCRIPTION

//  CW Klein

// BECAUSE WE TAKE RATIONAL COEFFCIENTS, I THINK HOMOLOGY OF KLEIN IS H_0=Q, H_1=Q, H_n=0 for n>1. SO IN THE ASSOCIATED BARCODE THERE IS JUST 2 INFINITE BARS. If not true we need to modify infinite bars in barcodes below





use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use phfibre::pipelines::boundary_matrix_pipeline;
use phfibre::phfibre::{ConditionNone, ConditionLowerStar, ConditionLowerEdge};
use solar::utilities::sequences_and_ordinals::{ordinate_unique_vals, BiMapSequential};
use solar::utilities::cell_complexes::simplices_unweighted::facets::ordered_subsimplices_up_thru_dim_concatenated_vec;
use solar::rings::ring::{Ring, Semiring};


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
    //  Klein -- 4 vertices, 8 edges, 4 faces
     ///////               0---------1---------0
     ///                   |         |         |
     ///                  |||        |         |
     ///                   |         |         |
     ///                   2---------3---------2
     ///                   |         |         |
     ///                   |         |        |||
     ///                   |         |         |
     ///                   0---------1---------0
     ///



    //  Define +1 and -1
    let one_pos             =   CeofficientRing::one();
    let one_neg             =   ring.negate( one_pos.clone() );

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
vec![ ( 5, one_pos.clone() ),  ( 7, one_neg.clone()),   ( 9, one_neg.clone() ), ( 10, one_neg.clone())], //13: top right cell                                    vec![   ( 1, one_pos.clone() ),      ( 4, one_neg.clone() ) ],   // 1-cell  3->1
vec![ ( 4, one_neg.clone() ),  ( 7, one_neg.clone()),  ( 8, one_pos.clone() ), ( 11, one_pos.clone())], //14: bottom left cell
vec![ ( 5, one_neg.clone() ),  ( 6, one_neg.clone()),   ( 9, one_pos.clone() ), ( 11, one_neg.clone())] //15: bottom right cell
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
                                    vec![ 5, 7, 9,10  ],
                                    vec![ 4, 7, 8,11  ],
                                    vec![ 5, 6, 9,11  ]
                                ];

 /*
    println!("-------------------------------------------------------------------------------------");
    println!("Klein (CW Structure with four squares, so 16 cells), and barcode with no finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:1} ],
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
    println!("Klein (CW Structure with four squares, so 16 cells), and 1 dim-0 finite bar");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:3} ],
        fin: vec![BarFinite{dim:0,birth:1,death:2} ],
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
    println!("Klein (CW Structure with four squares, so 16 cells), and 1 dim-1 finite bar");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:3} ],
        fin: vec![BarFinite{dim:1,birth:1,death:2} ],
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
    println!("Klein (CW Structure with four squares, so 16 cells), and two consecutive dim-0 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:5} ],
        fin: vec![BarFinite{dim:0,birth:1,death:2}, BarFinite{dim:0,birth:3,death:4} ],
        ordinal: ordinate_unique_vals( & vec![0, 1, 2, 3, 4, 5] ),
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
    println!("Klein (CW Structure with four squares, so 16 cells), and two consecutive touching dim-0 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:4} ],
        fin: vec![BarFinite{dim:0,birth:1,death:2}, BarFinite{dim:0,birth:2, death:3} ],
        ordinal: ordinate_unique_vals( & vec![0, 1, 2, 3, 4] ),
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
    println!("Klein (CW Structure with four squares, so 16 cells), and two consecutive dim-1 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:5} ],
        fin: vec![BarFinite{dim:1,birth:1,death:2}, BarFinite{dim:1,birth:3,death:4} ],
        ordinal: ordinate_unique_vals( & vec![0, 1, 2, 3, 4, 5] ),
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
    println!("Klein (CW Structure with four squares, so 16 cells), and two consecutive touching dim-1 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:4} ],
        fin: vec![BarFinite{dim:1,birth:1,death:2}, BarFinite{dim:1,birth:2, death:3} ],
        ordinal: ordinate_unique_vals( & vec![0, 1, 2, 3, 4] ),
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
    println!("Klein (CW Structure with four squares, so 16 cells), dim-0 and then dim-1 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:5} ],
        fin: vec![BarFinite{dim:0,birth:1,death:2}, BarFinite{dim:1,birth:3, death:4} ],
        ordinal: ordinate_unique_vals( & vec![0, 1, 2, 3, 4,5] ),
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
    );   */

    println!("-------------------------------------------------------------------------------------");
    println!("Klein (CW Structure with four squares, so 16 cells), one dim-0 finite bar and then one dim-1 finite bar, infinite bars same birth times:");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:5}],
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
    println!("Klein (CW Structure with four squares, so 16 cells), one loop and then one dim-0 finite bar");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:1} ],
        fin: vec![ BarFinite{dim:0,birth:2,death:3}],
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
    println!("Klein (CW Structure with four squares, so 16 cells), dim-0 finite bar overlapping one loop");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:2}],
        fin: vec![ BarFinite{dim:0,birth:1,death:3}],
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
}

/*

-------------------------------------------------------------------------------------
Klein (CW Structure with four squares, so 16 cells), and barcode with no finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: -1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (2, Ratio { numer: -1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: -1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 1 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0, 1], val_to_ord: {0: 0, 1: 1} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 367.10406117s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 1280
number of facets (binned by dimension): [0, 0, 0, 0, 0, 1280]
number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 1280]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 164.697067ms"
number of polytopes (total): 35584
number of polytopes (binned by dimension): [768, 4768, 10848, 11744, 6176, 1280]
Time elapsed to compute binary-coeff polyhedral betti numbers 43.602218ms
betti numbers (of polyhedral complex, Z2 coefficients): [1, 2, 1, 0, 0, 0]
-------------------------------------------------------------------------------------

Klein (CW Structure with four squares, so 16 cells), and 1 dim-0 finite bar

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: -1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (2, Ratio { numer: -1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: -1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 3 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3], val_to_ord: {2: 2, 1: 1, 3: 3, 0: 0} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 7783.74311249s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 11008
number of facets (binned by dimension): [0, 0, 0, 0, 11008]
number of facets (binned by number of vertices): [0, 0, 0, 0, 11008]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 1.006462121s"
number of polytopes (total): 197856
number of polytopes (binned by dimension): [13512, 52312, 74400, 46624, 11008]
Time elapsed to compute binary-coeff polyhedral betti numbers 570.054592ms
betti numbers (of polyhedral complex, Z2 coefficients): [1, 17, 0, 0, 0]
-------------------------------------------------------------------------------------

Klein (CW Structure with four squares, so 16 cells), and 1 dim-1 finite bar

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: -1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (2, Ratio { numer: -1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: -1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 3 }], fin: [BarFinite { dim: 1, birth: 1, death: 2 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3], val_to_ord: {3: 3, 2: 2, 1: 1, 0: 0} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 124.334913182s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 256
number of facets (binned by dimension): [0, 0, 0, 0, 256]
number of facets (binned by number of vertices): [0, 0, 0, 0, 256]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 16.532442ms"
number of polytopes (total): 3840
number of polytopes (binned by dimension): [192, 896, 1472, 1024, 256]
Time elapsed to compute binary-coeff polyhedral betti numbers 1.996994ms
betti numbers (of polyhedral complex, Z2 coefficients): [4, 4, 0, 0, 0]


Klein (CW Structure with four squares, so 16 cells), and two consecutive dim-0 finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: -1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (2, Ratio { numer: -1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: -1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 5 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }, BarFinite { dim: 0, birth: 3, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4, 5], val_to_ord: {1: 1, 4: 4, 3: 3, 5: 5, 0: 0, 2: 2} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 2141.936502364s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 6528
number of facets (binned by dimension): [0, 0, 0, 6528]
number of facets (binned by number of vertices): [0, 0, 0, 6528]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 346.897707ms"
number of polytopes (total): 66560
number of polytopes (binned by dimension): [10656, 26784, 22592, 6528]
Time elapsed to compute binary-coeff polyhedral betti numbers 93.27435ms
betti numbers (of polyhedral complex, Z2 coefficients): [2, 66, 0, 0]
-------------------------------------------------------------------------------------

Klein (CW Structure with four squares, so 16 cells), and two consecutive touching dim-0 finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: -1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (2, Ratio { numer: -1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: -1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 4 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }, BarFinite { dim: 0, birth: 2, death: 3 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4], val_to_ord: {4: 4, 2: 2, 1: 1, 0: 0, 3: 3} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 1725.560922109s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 5632
number of facets (binned by dimension): [0, 0, 0, 5632]
number of facets (binned by number of vertices): [0, 0, 0, 5632]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 339.328538ms"
number of polytopes (total): 57216
number of polytopes (binned by dimension): [9120, 23008, 19456, 5632]
Time elapsed to compute binary-coeff polyhedral betti numbers 65.748705ms
betti numbers (of polyhedral complex, Z2 coefficients): [2, 66, 0, 0]
-------------------------------------------------------------------------------------

Klein (CW Structure with four squares, so 16 cells), and two consecutive dim-1 finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: -1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (2, Ratio { numer: -1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: -1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 5 }], fin: [BarFinite { dim: 1, birth: 1, death: 2 }, BarFinite { dim: 1, birth: 3, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4, 5], val_to_ord: {3: 3, 4: 4, 1: 1, 5: 5, 0: 0, 2: 2} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 10.121589877s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 0
number of facets (binned by dimension): []
number of facets (binned by number of vertices): []
POLYHEDRAL COMPLEX CELLS
thread 'main' panicked at 'called `Option::unwrap()` on a `None` value', src/polytope/faces.rs:482:42
note: run with `RUST_BACKTRACE=1` environment variable to display a backtrace

-------------------------------------------------------------------------------------
Klein (CW Structure with four squares, so 16 cells), one dim-0 finite bar and then one dim-1 finite bar, infinite bars same birth times:

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: -1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (2, Ratio { numer: -1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: -1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 5 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }, BarFinite { dim: 1, birth: 3, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4, 5], val_to_ord: {1: 1, 2: 2, 4: 4, 5: 5, 0: 0, 3: 3} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 875.796933686s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 2816
number of facets (binned by dimension): [0, 0, 0, 2816]
number of facets (binned by number of vertices): [0, 0, 0, 2816]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 143.007676ms"
number of polytopes (total): 25600
number of polytopes (binned by dimension): [3456, 9984, 9344, 2816]
Time elapsed to compute binary-coeff polyhedral betti numbers 21.858862ms
betti numbers (of polyhedral complex, Z2 coefficients): [8, 8, 0, 0]
-------------------------------------------------------------------------------------

Klein (CW Structure with four squares, so 16 cells), one loop and then one dim-0 finite bar

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: -1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (2, Ratio { numer: -1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: -1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 1 }], fin: [BarFinite { dim: 0, birth: 2, death: 3 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3], val_to_ord: {1: 1, 2: 2, 3: 3, 0: 0} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 64.171548606s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 768
number of facets (binned by dimension): [0, 0, 0, 0, 768]
number of facets (binned by number of vertices): [0, 0, 0, 0, 768]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 64.442938ms"
number of polytopes (total): 14656
number of polytopes (binned by dimension): [1088, 4000, 5472, 3328, 768]
Time elapsed to compute binary-coeff polyhedral betti numbers 12.443578ms
betti numbers (of polyhedral complex, Z2 coefficients): [2, 4, 2, 0, 0]
-------------------------------------------------------------------------------------

Klein (CW Structure with four squares, so 16 cells), dim-0 finite bar overlapping one loop

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: -1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (2, Ratio { numer: -1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: -1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: 1, denom: 1 }), (11, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 2 }], fin: [BarFinite { dim: 0, birth: 1, death: 3 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3], val_to_ord: {3: 3, 2: 2, 1: 1, 0: 0} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 1040.90498307s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 2944
number of facets (binned by dimension): [0, 0, 0, 0, 2944]
number of facets (binned by number of vertices): [0, 0, 0, 0, 2944]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 241.116718ms"
number of polytopes (total): 54720
number of polytopes (binned by dimension): [3936, 14752, 20480, 12608, 2944]
Time elapsed to compute binary-coeff polyhedral betti numbers 71.895071ms
betti numbers (of polyhedral complex, Z2 coefficients): [4, 8, 4, 0, 0]


*/

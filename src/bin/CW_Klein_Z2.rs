//  DESCRIPTION

//  CW Klein Z2

// BECAUSE WE TAKE Z2 COEFFCIENTS, I THINK HOMOLOGY OF KLEIN IS H_0=Z2, H_1=Z2+Z2, H_2=Z2  and H_n=0 for n>2. SO IN THE ASSOCIATED BARCODE THERE IS JUST 2 INFINITE BARS. If not true we need to modify infinite bars in barcodes below





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


    type CeofficientRing    =   solar::rings::field_prime::GF2;
    let ring                =   solar::rings::field_prime::GF2{}; // the object itself

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


    println!("-------------------------------------------------------------------------------------");
    println!("Klein Z2 (CW Structure with four squares, so 16 cells), and barcode with no finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:1}, BarInfinite{dim:1,birth:1},BarInfinite{dim:2,birth:1}  ],
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
    println!("Klein Z2 (CW Structure with four squares, so 16 cells), and 1 dim-0 finite bar");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:3}, BarInfinite{dim:1,birth:3},BarInfinite{dim:2,birth:3}  ],
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
    println!("Klein Z2 (CW Structure with four squares, so 16 cells), and 1 dim-1 finite bar");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:3}, BarInfinite{dim:1,birth:3},BarInfinite{dim:2,birth:3} ],
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
    println!("Klein Z2 (CW Structure with four squares, so 16 cells), and two consecutive dim-0 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:5}, BarInfinite{dim:1,birth:5},BarInfinite{dim:2,birth:5} ],
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
    println!("Klein Z2 (CW Structure with four squares, so 16 cells), and two consecutive touching dim-0 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:4}, BarInfinite{dim:1,birth:4},BarInfinite{dim:2,birth:4}  ],
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

/*
    println!("-------------------------------------------------------------------------------------");
    println!("Klein (CW Structure with four squares, so 16 cells), and two consecutive dim-1 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:5}, BarInfinite{dim:1,birth:5},BarInfinite{dim:2,birth:5}  ],
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
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:4}, BarInfinite{dim:1,birth:4}, BarInfinite{dim:2,birth:4}  ],
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

*/
    println!("-------------------------------------------------------------------------------------");
    println!("Klein Z2 (CW Structure with four squares, so 16 cells), dim-0 and then dim-1 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:5}, BarInfinite{dim:1,birth:5}, BarInfinite{dim:2,birth:5}  ],
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
    );

    println!("-------------------------------------------------------------------------------------");
    println!("Klein Z2 (CW Structure with four squares, so 16 cells), one loop and then one dim-0 finite bar");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:1}, BarInfinite{dim:1,birth:4},BarInfinite{dim:2,birth:4}  ],
        fin: vec![ BarFinite{dim:0,birth:2,death:3}],
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
    println!("Klein Z2 (CW Structure with four squares, so 16 cells), dim-0 finite bar overlapping one loop");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:2}, BarInfinite{dim:1,birth:4},BarInfinite{dim:2,birth:4} ],
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
Built boundary matrix
-------------------------------------------------------------------------------------
Klein Z2 (CW Structure with four squares, so 16 cells), and barcode with no finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, true), (1, true)], [(0, true), (1, true)], [(0, true), (2, true)], [(0, true), (2, true)], [(2, true), (3, true)], [(2, true), (3, true)], [(1, true), (3, true)], [(1, true), (3, true)], [(4, true), (6, true), (8, true), (10, true)], [(5, true), (7, true), (9, true), (10, true)], [(4, true), (7, true), (8, true), (11, true)], [(5, true), (6, true), (9, true), (11, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 1 }, BarInfinite { dim: 1, birth: 1 }, BarInfinite { dim: 2, birth: 1 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0, 1], val_to_ord: {0: 0, 1: 1} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 284.374450779s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 256
number of facets (binned by dimension): [0, 0, 0, 128, 128]
number of facets (binned by number of vertices): [0, 0, 0, 128, 128]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 7.290719ms"
number of polytopes (total): 2048
number of polytopes (binned by dimension): [64, 448, 832, 576, 128]
Time elapsed to compute binary-coeff polyhedral betti numbers 1.162944ms
betti numbers (of polyhedral complex, Z2 coefficients): [1, 2, 1, 0, 0]
-------------------------------------------------------------------------------------

Klein Z2 (CW Structure with four squares, so 16 cells), and 1 dim-0 finite bar

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, true), (1, true)], [(0, true), (1, true)], [(0, true), (2, true)], [(0, true), (2, true)], [(2, true), (3, true)], [(2, true), (3, true)], [(1, true), (3, true)], [(1, true), (3, true)], [(4, true), (6, true), (8, true), (10, true)], [(5, true), (7, true), (9, true), (10, true)], [(4, true), (7, true), (8, true), (11, true)], [(5, true), (6, true), (9, true), (11, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 3 }, BarInfinite { dim: 1, birth: 3 }, BarInfinite { dim: 2, birth: 3 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3], val_to_ord: {3: 3, 1: 1, 0: 0, 2: 2} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 6751.431816555s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 2816
number of facets (binned by dimension): [0, 0, 1408, 1408]
number of facets (binned by number of vertices): [0, 0, 1408, 1408]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 70.044369ms"
number of polytopes (total): 16400
number of polytopes (binned by dimension): [2368, 6800, 5824, 1408]
Time elapsed to compute binary-coeff polyhedral betti numbers 9.192696ms
betti numbers (of polyhedral complex, Z2 coefficients): [1, 17, 0, 0]
-------------------------------------------------------------------------------------
Klein Z2 (CW Structure with four squares, so 16 cells), and 1 dim-1 finite bar

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, true), (1, true)], [(0, true), (1, true)], [(0, true), (2, true)], [(0, true), (2, true)], [(2, true), (3, true)], [(2, true), (3, true)], [(1, true), (3, true)], [(1, true), (3, true)], [(4, true), (6, true), (8, true), (10, true)], [(5, true), (7, true), (9, true), (10, true)], [(4, true), (7, true), (8, true), (11, true)], [(5, true), (6, true), (9, true), (11, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 3 }, BarInfinite { dim: 1, birth: 3 }, BarInfinite { dim: 2, birth: 3 }], fin: [BarFinite { dim: 1, birth: 1, death: 2 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3], val_to_ord: {0: 0, 1: 1, 2: 2, 3: 3} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 11.293032715s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 128
number of facets (binned by dimension): [0, 0, 0, 128]
number of facets (binned by number of vertices): [0, 0, 0, 128]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 2.776045ms"
number of polytopes (total): 768
number of polytopes (binned by dimension): [64, 256, 320, 128]
Time elapsed to compute binary-coeff polyhedral betti numbers 242.924µs
betti numbers (of polyhedral complex, Z2 coefficients): [4, 4, 0, 0]
-------------------------------------------------------------------------------------

Klein Z2 (CW Structure with four squares, so 16 cells), and two consecutive dim-0 finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, true), (1, true)], [(0, true), (1, true)], [(0, true), (2, true)], [(0, true), (2, true)], [(2, true), (3, true)], [(2, true), (3, true)], [(1, true), (3, true)], [(1, true), (3, true)], [(4, true), (6, true), (8, true), (10, true)], [(5, true), (7, true), (9, true), (10, true)], [(4, true), (7, true), (8, true), (11, true)], [(5, true), (6, true), (9, true), (11, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 5 }, BarInfinite { dim: 1, birth: 5 }, BarInfinite { dim: 2, birth: 5 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }, BarFinite { dim: 0, birth: 3, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4, 5], val_to_ord: {0: 0, 1: 1, 4: 4, 5: 5, 3: 3, 2: 2} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 2075.476879377s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 1792
number of facets (binned by dimension): [0, 896, 896]
number of facets (binned by number of vertices): [0, 896, 896]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 22.613276ms"
number of polytopes (total): 6080
number of polytopes (binned by dimension): [2112, 3072, 896]
Time elapsed to compute binary-coeff polyhedral betti numbers 1.563998ms
betti numbers (of polyhedral complex, Z2 coefficients): [2, 66, 0]
-------------------------------------------------------------------------------------
Klein Z2 (CW Structure with four squares, so 16 cells), and two consecutive touching dim-0 finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, true), (1, true)], [(0, true), (1, true)], [(0, true), (2, true)], [(0, true), (2, true)], [(2, true), (3, true)], [(2, true), (3, true)], [(1, true), (3, true)], [(1, true), (3, true)], [(4, true), (6, true), (8, true), (10, true)], [(5, true), (7, true), (9, true), (10, true)], [(4, true), (7, true), (8, true), (11, true)], [(5, true), (6, true), (9, true), (11, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 4 }, BarInfinite { dim: 1, birth: 4 }, BarInfinite { dim: 2, birth: 4 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }, BarFinite { dim: 0, birth: 2, death: 3 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4], val_to_ord: {1: 1, 3: 3, 0: 0, 4: 4, 2: 2} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 1705.101380978s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 1536
number of facets (binned by dimension): [0, 768, 768]
number of facets (binned by number of vertices): [0, 768, 768]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 17.851966ms"
number of polytopes (total): 5184
number of polytopes (binned by dimension): [1792, 2624, 768]
Time elapsed to compute binary-coeff polyhedral betti numbers 1.285774ms
betti numbers (of polyhedral complex, Z2 coefficients): [2, 66, 0]
-------------------------------------------------------------------------------------

Klein Z2 (CW Structure with four squares, so 16 cells), dim-0 and then dim-1 finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, true), (1, true)], [(0, true), (1, true)], [(0, true), (2, true)], [(0, true), (2, true)], [(2, true), (3, true)], [(2, true), (3, true)], [(1, true), (3, true)], [(1, true), (3, true)], [(4, true), (6, true), (8, true), (10, true)], [(5, true), (7, true), (9, true), (10, true)], [(4, true), (7, true), (8, true), (11, true)], [(5, true), (6, true), (9, true), (11, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 5 }, BarInfinite { dim: 1, birth: 5 }, BarInfinite { dim: 2, birth: 5 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }, BarFinite { dim: 1, birth: 3, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4, 5], val_to_ord: {2: 2, 0: 0, 5: 5, 3: 3, 4: 4, 1: 1} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 103.745537482s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 1408
number of facets (binned by dimension): [0, 0, 1408]
number of facets (binned by number of vertices): [0, 0, 1408]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 22.731381ms"
number of polytopes (total): 5120
number of polytopes (binned by dimension): [1152, 2560, 1408]
Time elapsed to compute binary-coeff polyhedral betti numbers 1.381579ms
betti numbers (of polyhedral complex, Z2 coefficients): [8, 8, 0]
-------------------------------------------------------------------------------------
Klein Z2 (CW Structure with four squares, so 16 cells), one loop and then one dim-0 finite bar

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, true), (1, true)], [(0, true), (1, true)], [(0, true), (2, true)], [(0, true), (2, true)], [(2, true), (3, true)], [(2, true), (3, true)], [(1, true), (3, true)], [(1, true), (3, true)], [(4, true), (6, true), (8, true), (10, true)], [(5, true), (7, true), (9, true), (10, true)], [(4, true), (7, true), (8, true), (11, true)], [(5, true), (6, true), (9, true), (11, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 1 }, BarInfinite { dim: 1, birth: 4 }, BarInfinite { dim: 2, birth: 4 }], fin: [BarFinite { dim: 0, birth: 2, death: 3 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4], val_to_ord: {2: 2, 3: 3, 0: 0, 1: 1, 4: 4} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 110.550308113s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 1696
number of facets (binned by dimension): [0, 0, 160, 768, 768]
number of facets (binned by number of vertices): [0, 0, 160, 768, 768]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 115.508981ms"
number of polytopes (total): 21248
number of polytopes (binned by dimension): [2080, 6688, 7776, 3936, 768]
Time elapsed to compute binary-coeff polyhedral betti numbers 16.431488ms
betti numbers (of polyhedral complex, Z2 coefficients): [6, 12, 6, 0, 0]
-------------------------------------------------------------------------------------

Klein Z2 (CW Structure with four squares, so 16 cells), dim-0 finite bar overlapping one loop

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, true), (1, true)], [(0, true), (1, true)], [(0, true), (2, true)], [(0, true), (2, true)], [(2, true), (3, true)], [(2, true), (3, true)], [(1, true), (3, true)], [(1, true), (3, true)], [(4, true), (6, true), (8, true), (10, true)], [(5, true), (7, true), (9, true), (10, true)], [(4, true), (7, true), (8, true), (11, true)], [(5, true), (6, true), (9, true), (11, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 2 }, BarInfinite { dim: 1, birth: 4 }, BarInfinite { dim: 2, birth: 4 }], fin: [BarFinite { dim: 0, birth: 1, death: 3 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4], val_to_ord: {3: 3, 0: 0, 4: 4, 1: 1, 2: 2} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 1707.653436664s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 6336
number of facets (binned by dimension): [0, 0, 448, 2944, 2944]
number of facets (binned by number of vertices): [0, 0, 448, 2944, 2944]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 472.23975ms"
number of polytopes (total): 80384
number of polytopes (binned by dimension): [7616, 25088, 29632, 15104, 2944]
Time elapsed to compute binary-coeff polyhedral betti numbers 98.120642ms
betti numbers (of polyhedral complex, Z2 coefficients): [12, 24, 12, 0, 0]


*/

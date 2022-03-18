//  DESCRIPTION

//  CW Dunce Hat

// BECAUSE Dunce Hat is contractible and WE TAKE RATIONAL COEFFCIENTS, I THINK HOMOLOGY OF MOBIUS IS H_0=Q, H_n=0 for n>0. SO IN THE ASSOCIATED BARCODE THERE IS JUST 1 INFINITE BAR. If not true we need to modify infinite bars in barcodes below





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
    // Dunce Hat 12 cells-- 4 vertices, 6 edges,  2 faces
    //                           0
    ///                        / | \
     ///////                  /  |  \
     ///                     /   |   \
     ///                    /    |    \
     ///                   1-----2-----1
     ///                 //   /  |  \   \\
     ///                /  /     |     \  \
     ///               0----=----1---------0



    //  Define +1 and -1
    let one_pos             =   CeofficientRing::one();
    let one_neg             =   ring.negate( one_pos.clone() );

    let analyze_dowker_dual =   false;
    //  The boundary matrix in vector-of-vector format; each 'internal'vector represents a compressed column of the boundary matrix
    let boundary            =   vec![
                                    vec![                               ],                          //0:  0-cell 0
                                    vec![                               ],                          //1:  0-cell 1
                                    vec![                               ],                          //2:  0-cell 2
                                    vec![   ( 0, one_neg.clone() ),      ( 1, one_pos.clone() ) ],    //3:  1-cell 0->1    (first arrow bottom edge)
                                    vec![   ( 0, one_pos.clone() ),      ( 1, one_neg.clone() ) ],  //4: 1-cell  1->0     (second arrow bottom edge)
                                    vec![   ( 0, one_neg.clone() ),      ( 2, one_pos.clone() ) ],  //5:  1-cell  0->2     (bottom left triangle)
                                    vec![   ( 0, one_pos.clone() ),      ( 2, one_neg.clone() ) ],   //6:  1-cell  2->0     (bottom right triangle)
                                    vec![   ( 0, one_pos.clone() ),      ( 2, one_neg.clone() ) ],  //7:  1-cell  2->0       (top triangle)
                                    vec![   ( 1, one_neg.clone() ),      ( 2, one_pos.clone() ) ],  //8:  1-cell  1->2     (left)
                                    vec![   ( 1, one_pos.clone() ),      ( 2, one_neg.clone() ) ],   //9:  1-cell  2->1     (middle)
                                    vec![   ( 1, one_pos.clone() ),      ( 2, one_neg.clone() ) ],   //10:  1-cell  2->1       (right)
 vec![  ( 3, one_neg.clone() ), ( 5, one_pos.clone() ),  ( 8, one_neg.clone())],  //11:  2-cell  0->2->1 (left most)
 vec![ ( 3, one_pos.clone() ), ( 5, one_neg.clone() ), ( 9, one_neg.clone())  ],  //12: 2-cell  0->1->2 (second left)
vec![ ( 4, one_pos.clone() ),  ( 6, one_neg.clone()),  ( 9, one_pos.clone() )], //13: 2-cell  1->0->2 (right)
vec![ ( 3, one_pos.clone() ),  ( 6, one_pos.clone()),  ( 10, one_neg.clone()) ], //14:2-cell  0->1->2 (right most)
vec![ ( 4, one_neg.clone() ),  ( 7, one_pos.clone()),  ( 8, one_pos.clone()) ], //15: 2-cell  0->1->2 (top left)
vec![ ( 4, one_pos.clone() ),  ( 7, one_neg.clone()),  ( 10, one_pos.clone()) ] //16: 2-cell  0->2->1 (top right)
   ];

    println!("Built boundary matrix");

    //  The dimension of each cell.
    let cell_dims           =   vec![ 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2];

    //  Definition: cell_id_to_prereqs[ m ] = { the cells that must be added to the complex before M }
    let cell_id_to_prereqs  =   vec![
                                    vec![           ], // 0
                                    vec![           ],
                                    vec![           ],  // 2:
                                    vec![  0, 1       ],   //3
                                    vec![ 0, 1       ],
                                    vec![ 0, 2       ],
                                    vec![ 0, 2       ],
                                    vec![ 0, 2       ],  //7
                                    vec![ 1, 2       ],
                                    vec![ 1, 2       ],
                                    vec![ 1, 2       ], //10
                                    vec![ 0, 1, 2  ], //11
                                    vec![ 0, 1, 2  ],
                                    vec![ 0, 1, 2  ],
                                    vec![ 0, 1, 2  ],
                                    vec![ 0, 1, 2  ], //15
                                    vec![ 0, 1, 2  ],
                                ];
    /*

    println!("-------------------------------------------------------------------------------------");
    println!("Dunce Hat (CW Structure with 6 2-simplices, so 17 cells), and barcode with no finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}],
        fin: vec![ ],
        ordinal: ordinate_unique_vals( & vec![0] ),
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
    println!("Dunce Hat (CW Structure with 6 2-simplices, so 17 cells), and barcode with dim-0 finite bar");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0} ],
        fin: vec![ BarFinite{dim:0,birth:1, death:2} ],
        ordinal: ordinate_unique_vals( & vec![0,1,2] ),
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
    println!("Dunce Hat (CW Structure with 6 2-simplices, so 17 cells), and barcode with dim-1 finite bar");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0} ],
        fin: vec![ BarFinite{dim:1,birth:1, death:2} ],
        ordinal: ordinate_unique_vals( & vec![0,1,2] ),
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
    println!("Dunce Hat (CW Structure with 6 2-simplices, so 17 cells), and barcode with two dim-0 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}],
        fin: vec![ BarFinite{dim:0,birth:1, death:2}, BarFinite{dim:0,birth:3, death:4} ],
        ordinal: ordinate_unique_vals( & vec![0,1,2,3,4] ),
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
    println!("Dunce Hat (CW Structure with 6 2-simplices, so 17 cells), and barcode with two overlapping dim-0 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}],
        fin: vec![ BarFinite{dim:0,birth:1, death:3}, BarFinite{dim:0,birth:2, death:4} ],
        ordinal: ordinate_unique_vals( & vec![0,1,2,3,4] ),
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
    println!("Dunce Hat (CW Structure with 6 2-simplices, so 17 cells), and barcode with two dim-1 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}],
        fin: vec![ BarFinite{dim:1,birth:1, death:2}, BarFinite{dim:1,birth:3, death:4} ],
        ordinal: ordinate_unique_vals( & vec![0,1,2,3,4] ),
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
    println!("Dunce Hat (CW Structure with 6 2-simplices, so 17 cells), and barcode with two overlapping dim-1 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}],
        fin: vec![ BarFinite{dim:1,birth:1, death:3}, BarFinite{dim:1,birth:2, death:4} ],
        ordinal: ordinate_unique_vals( & vec![0,1,2,3,4] ),
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
    println!("Dunce Hat (CW Structure with 6 2-simplices, so 17 cells), and barcode with dim-0 and then dim-1 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}],
        fin: vec![ BarFinite{dim:0,birth:1, death:2}, BarFinite{dim:1,birth:3, death:4} ],
        ordinal: ordinate_unique_vals( & vec![0,1,2,3,4] ),
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
    println!("Dunce Hat (CW Structure with 6 2-simplices, so 17 cells), and barcode with two overlapping dim-0,dim-1 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}],
        fin: vec![ BarFinite{dim:0,birth:1, death:3}, BarFinite{dim:1,birth:2, death:4} ],
        ordinal: ordinate_unique_vals( & vec![0,1,2,3,4] ),
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
    println!("Dunce Hat (CW Structure with 6 2-simplices, so 17 cells), and barcode with three dim-1 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}],
        fin: vec![ BarFinite{dim:1,birth:1, death:2}, BarFinite{dim:1,birth:3, death:4}, BarFinite{dim:1,birth:5, death:6} ],
        ordinal: ordinate_unique_vals( & vec![0,1,2,3,4,5,6] ),
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
Dunce Hat (CW Structure with 6 2-simplices, so 17 cells), and barcode with no finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: -1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (2, Ratio { numer: -1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (2, Ratio { numer: -1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: 1, denom: 1 }), (2, Ratio { numer: -1, denom: 1 })], [(1, Ratio { numer: 1, denom: 1 }), (2, Ratio { numer: -1, denom: 1 })], [(3, Ratio { numer: -1, denom: 1 }), (5, Ratio { numer: 1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 })], [(3, Ratio { numer: 1, denom: 1 }), (5, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: 1, denom: 1 })], [(3, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: 1, denom: 1 }), (10, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: 1, denom: 1 }), (8, Ratio { numer: 1, denom: 1 })], [(4, Ratio { numer: 1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: 1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0], val_to_ord: {0: 0} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 1214.895552566s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 92
number of facets (binned by dimension): [0, 0, 0, 20, 24, 48]
number of facets (binned by number of vertices): [0, 0, 0, 20, 24, 48]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 4.544711ms"
number of polytopes (total): 1599
number of polytopes (binned by dimension): [41, 243, 515, 508, 244, 48]
Time elapsed to compute binary-coeff polyhedral betti numbers 586.003µs
betti numbers (of polyhedral complex, Z2 coefficients): [1, 0, 0, 0, 0, 0]
-------------------------------------------------------------------------------------

Dunce Hat (CW Structure with 6 2-simplices, so 17 cells), and barcode with dim-0 finite bar

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: -1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (2, Ratio { numer: -1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (2, Ratio { numer: -1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: 1, denom: 1 }), (2, Ratio { numer: -1, denom: 1 })], [(1, Ratio { numer: 1, denom: 1 }), (2, Ratio { numer: -1, denom: 1 })], [(3, Ratio { numer: -1, denom: 1 }), (5, Ratio { numer: 1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 })], [(3, Ratio { numer: 1, denom: 1 }), (5, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: 1, denom: 1 })], [(3, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: 1, denom: 1 }), (10, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: 1, denom: 1 }), (8, Ratio { numer: 1, denom: 1 })], [(4, Ratio { numer: 1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: 1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2], val_to_ord: {0: 0, 1: 1, 2: 2} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 5439.212181722s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 368
number of facets (binned by dimension): [0, 0, 80, 96, 192]
number of facets (binned by number of vertices): [0, 0, 80, 96, 192]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 15.0577ms"
number of polytopes (total): 4834
number of polytopes (binned by dimension): [500, 1518, 1720, 904, 192]
Time elapsed to compute binary-coeff polyhedral betti numbers 3.553777ms
betti numbers (of polyhedral complex, Z2 coefficients): [1, 11, 0, 0, 0]
-------------------------------------------------------------------------------------


*/

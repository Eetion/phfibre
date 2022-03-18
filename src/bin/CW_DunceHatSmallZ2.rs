//  DESCRIPTION

//  CW Dunce Hat Small triangulation

// BECAUSE Dunce Hat is contractible and WE TAKE RATIONAL COEFFCIENTS, I THINK HOMOLOGY OF MOBIUS IS H_0=Z2, H_n=0 for n>0. SO IN THE ASSOCIATED BARCODE THERE IS JUST 1 INFINITE BAR. If not true we need to modify infinite bars in barcodes below





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
    // Dunce Hat 11 cells-- 3 vertices, 5 edges,  3 faces
    //                           0
    ///                        / | \
     ///////                  /  |  \
     ///                     /   |   \
     ///                    1    |    1        Orientations for edge: up
     ///                   /     2     \
     ///                 //   /    \    \\
     ///                /  /           \  \
     ///               0---------1--------0



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
 vec![ ( 3, one_pos.clone() ), ( 4, one_pos.clone()) ,( 5, one_neg.clone()), ( 6, one_neg.clone())  ],  //8: 2-cell  (bottom)
 vec![  ( 3, one_neg.clone() ), ( 4, one_neg.clone() ),  ( 5, one_pos.clone()),( 7, one_pos.clone()) ],  //9:  2-cell (laft)
vec![ ( 3, one_pos.clone() ), ( 4, one_pos.clone()), ( 6, one_pos.clone()),  ( 7, one_neg.clone())], //10: 2-cell  (right)
   ];

    println!("Built boundary matrix");

    //  The dimension of each cell.
    let cell_dims           =   vec![ 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2];

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
                                    vec![ 0, 1, 2  ], //8
                                    vec![ 0, 1, 2  ],
                                    vec![ 0, 1, 2  ]
                                ];


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

    /*

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
    */

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
}

/*

EMPTY FIBERS FOR
-TWO DIM 1 finite bars
- TWO overlapping DIM 1 finite bars

Dunce Hat (CW Structure with 6 2-simplices, so 17 cells), and barcode with no finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [(0, true), (1, true)], [(0, true), (1, true)], [(0, true), (2, true)], [(0, true), (2, true)], [(0, true), (2, true)], [(3, true), (4, true), (5, true), (6, true)], [(3, true), (4, true), (5, true), (7, true)], [(3, true), (4, true), (6, true), (7, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0], val_to_ord: {0: 0} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 4.731726ms

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 24
number of facets (binned by dimension): [0, 0, 0, 24]
number of facets (binned by number of vertices): [0, 0, 0, 24]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 414.91µs"
number of polytopes (total): 157
number of polytopes (binned by dimension): [15, 54, 64, 24]
Time elapsed to compute binary-coeff polyhedral betti numbers 46.741µs
betti numbers (of polyhedral complex, Z2 coefficients): [1, 0, 0, 0]
-------------------------------------------------------------------------------------
Dunce Hat (CW Structure with 6 2-simplices, so 17 cells), and barcode with dim-0 finite bar

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [(0, true), (1, true)], [(0, true), (1, true)], [(0, true), (2, true)], [(0, true), (2, true)], [(0, true), (2, true)], [(3, true), (4, true), (5, true), (6, true)], [(3, true), (4, true), (5, true), (7, true)], [(3, true), (4, true), (6, true), (7, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2], val_to_ord: {0: 0, 2: 2, 1: 1} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 22.428501ms

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 96
number of facets (binned by dimension): [0, 0, 96]
number of facets (binned by number of vertices): [0, 0, 96]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 1.148938ms"
number of polytopes (total): 436
number of polytopes (binned by dimension): [120, 220, 96]
Time elapsed to compute binary-coeff polyhedral betti numbers 120.859µs
betti numbers (of polyhedral complex, Z2 coefficients): [2, 6, 0]
-------------------------------------------------------------------------------------

Dunce Hat (CW Structure with 6 2-simplices, so 17 cells), and barcode with dim-1 finite bar

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [(0, true), (1, true)], [(0, true), (1, true)], [(0, true), (2, true)], [(0, true), (2, true)], [(0, true), (2, true)], [(3, true), (4, true), (5, true), (6, true)], [(3, true), (4, true), (5, true), (7, true)], [(3, true), (4, true), (6, true), (7, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [BarFinite { dim: 1, birth: 1, death: 2 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2], val_to_ord: {1: 1, 2: 2, 0: 0} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 49.111797ms

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 288
number of facets (binned by dimension): [0, 0, 0, 0, 288]
number of facets (binned by number of vertices): [0, 0, 0, 0, 288]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 17.257736ms"
number of polytopes (total): 3617
number of polytopes (binned by dimension): [163, 790, 1356, 1020, 288]
Time elapsed to compute binary-coeff polyhedral betti numbers 3.242006ms
betti numbers (of polyhedral complex, Z2 coefficients): [1, 4, 0, 0, 0]
-------------------------------------------------------------------------------------
Dunce Hat (CW Structure with 6 2-simplices, so 17 cells), and barcode with two dim-0 finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [(0, true), (1, true)], [(0, true), (1, true)], [(0, true), (2, true)], [(0, true), (2, true)], [(0, true), (2, true)], [(3, true), (4, true), (5, true), (6, true)], [(3, true), (4, true), (5, true), (7, true)], [(3, true), (4, true), (6, true), (7, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }, BarFinite { dim: 0, birth: 3, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4], val_to_ord: {4: 4, 1: 1, 3: 3, 0: 0, 2: 2} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 4.400204ms

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 24
number of facets (binned by dimension): [0, 24]
number of facets (binned by number of vertices): [0, 24]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 172.914µs"
number of polytopes (total): 58
number of polytopes (binned by dimension): [34, 24]
Time elapsed to compute binary-coeff polyhedral betti numbers 5.529µs
betti numbers (of polyhedral complex, Z2 coefficients): [10, 0]

---------------------------

Dunce Hat (CW Structure with 6 2-simplices, so 17 cells), and barcode with two overlapping dim-0 finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [(0, true), (1, true)], [(0, true), (1, true)], [(0, true), (2, true)], [(0, true), (2, true)], [(0, true), (2, true)], [(3, true), (4, true), (5, true), (6, true)], [(3, true), (4, true), (5, true), (7, true)], [(3, true), (4, true), (6, true), (7, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [BarFinite { dim: 0, birth: 1, death: 3 }, BarFinite { dim: 0, birth: 2, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4], val_to_ord: {3: 3, 0: 0, 1: 1, 2: 2, 4: 4} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 3.991743ms

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 24
number of facets (binned by dimension): [0, 24]
number of facets (binned by number of vertices): [0, 24]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 171.451µs"
number of polytopes (total): 58
number of polytopes (binned by dimension): [34, 24]
Time elapsed to compute binary-coeff polyhedral betti numbers 5.096µs
betti numbers (of polyhedral complex, Z2 coefficients): [10, 0]
-------------------------------------------------------------------------------------
Dunce Hat (CW Structure with 6 2-simplices, so 17 cells), and barcode with dim-0 and then dim-1 finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [(0, true), (1, true)], [(0, true), (1, true)], [(0, true), (2, true)], [(0, true), (2, true)], [(0, true), (2, true)], [(3, true), (4, true), (5, true), (6, true)], [(3, true), (4, true), (5, true), (7, true)], [(3, true), (4, true), (6, true), (7, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }, BarFinite { dim: 1, birth: 3, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4], val_to_ord: {3: 3, 1: 1, 2: 2, 4: 4, 0: 0} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 235.196214ms

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 864
number of facets (binned by dimension): [0, 0, 0, 864]
number of facets (binned by number of vertices): [0, 0, 0, 864]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 29.050409ms"
number of polytopes (total): 6466
number of polytopes (binned by dimension): [796, 2370, 2436, 864]
Time elapsed to compute binary-coeff polyhedral betti numbers 5.514777ms
betti numbers (of polyhedral complex, Z2 coefficients): [2, 4, 0, 0]
-------------------------------------------------------------------------------------


Dunce Hat (CW Structure with 6 2-simplices, so 17 cells), and barcode with two overlapping dim-0,dim-1 finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [(0, true), (1, true)], [(0, true), (1, true)], [(0, true), (2, true)], [(0, true), (2, true)], [(0, true), (2, true)], [(3, true), (4, true), (5, true), (6, true)], [(3, true), (4, true), (5, true), (7, true)], [(3, true), (4, true), (6, true), (7, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [BarFinite { dim: 0, birth: 1, death: 3 }, BarFinite { dim: 1, birth: 2, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4], val_to_ord: {3: 3, 1: 1, 2: 2, 4: 4, 0: 0} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 30.573276ms

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 288
number of facets (binned by dimension): [0, 0, 0, 288]
number of facets (binned by number of vertices): [0, 0, 0, 288]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 9.876009ms"
number of polytopes (total): 2232
number of polytopes (binned by dimension): [288, 828, 828, 288]
Time elapsed to compute binary-coeff polyhedral betti numbers 986.527µs
betti numbers (of polyhedral complex, Z2 coefficients): [8, 10, 2, 0]







*/

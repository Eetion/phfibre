//  DESCRIPTION

//  CW RP2 Z2 coeff

// BECAUSE WE TAKE Z2 COEFFCIENTS, I THINK HOMOLOGY OF PROJ PLANE IS H_0=Z2, H_1=Z2, H_2=Z_2  and H_n=0 for n>2. SO IN THE ASSOCIATED BARCODE THERE IS ONE INFINITE BAR. If not true we need to modify infinite bars in barcodes below





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
    type CeofficientRing    =   solar::rings::field_prime::GF2;
    let ring                =   solar::rings::field_prime::GF2{}; // the object itself

     //  ----------------------------------------------------------------------------------------------
    //  ProjPlane -- 4 vertices, 8 edges, 4 faces
     ///////               0---------1---==----4
     ///                   |         |         |
     ///                  |||        |         |
     ///                   |         |         |
     ///                   2---------3---------2
     ///                   |         |         |
     ///                   |         |        |||
     ///                   |         |         |
     ///                   4----==---1---------0
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
                                    vec![                               ],                          //4:  0-cell 4
                                    vec![   ( 0, one_neg.clone() ),      ( 1, one_pos.clone() ) ],  //5:  1-cell  0->1
                                    vec![   ( 1, one_neg.clone() ),      ( 4, one_pos.clone() ) ],  //6:  1-cell  1->4
                                    vec![   ( 0, one_neg.clone() ),      ( 2, one_pos.clone() ) ],   //7:  1-cell  0->2
                                    vec![   ( 2, one_neg.clone() ),      ( 4, one_pos.clone() ) ],   //8:  1-cell  2->4
                                    vec![   ( 2, one_neg.clone() ),      ( 3, one_pos.clone() ) ],  //9:  1-cell  2->3
                                    vec![   ( 2, one_pos.clone() ),      ( 3, one_neg.clone() ) ],  //10:  1-cell  3->2
                                    vec![   ( 1, one_neg.clone() ),      ( 3, one_pos.clone() ) ],   //11: 1-cell  1->3
                                    vec![   ( 1, one_pos.clone() ),      ( 3, one_neg.clone() ) ],   //12: 1-cell  3->1
vec![ ( 5, one_pos.clone() ),  ( 7, one_neg.clone()),  ( 9, one_neg.clone() ), ( 11, one_pos.clone())], //13: top left cell
vec![ ( 6, one_pos.clone() ),  ( 8, one_neg.clone()),   ( 10, one_neg.clone() ), ( 11, one_neg.clone())], //14: top right cell
vec![ ( 6, one_pos.clone() ),  ( 8, one_neg.clone()),  ( 9, one_pos.clone() ), ( 12, one_pos.clone())], //15: bottom left cell
vec![ ( 5, one_pos.clone() ),  ( 7, one_neg.clone()),   ( 10, one_pos.clone() ), ( 12, one_neg.clone())] //16: bottom right cell
   ];

    println!("Built boundary matrix");

    //  The dimension of each cell.
    let cell_dims           =   vec![ 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2 ];

    //  Definition: cell_id_to_prereqs[ m ] = { the cells that must be added to the complex before M }
    let cell_id_to_prereqs  =   vec![
                                    vec![           ], // 0
                                    vec![           ],
                                    vec![           ],
                                    vec![            ],   //3
                                    vec![            ],   //4
                                    vec![ 0, 1       ],
                                    vec![ 1, 4       ],
                                    vec![ 0, 2       ],
                                    vec![ 2, 4       ],     //8
                                    vec![ 2, 3       ],
                                    vec![ 2, 3       ],
                                    vec![ 1, 3       ],
                                    vec![ 1, 3       ],    //12
                                    vec![ 5, 7, 9, 11 ],
                                    vec![ 6, 8, 10, 11 ],
                                    vec![ 6, 8, 9, 12 ],
                                    vec![ 5, 7, 10, 12 ]
                                ];


    println!("-------------------------------------------------------------------------------------");
    println!("Projective Plane Z2 (CW Structure with four squares, so 16 cells), and barcode with no finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:1}, BarInfinite{dim:2,birth:2} ],
        fin: vec![ ],
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
    println!("Projective Plane Z2 (CW Structure with four squares, so 16 cells), and 1 dim-0 finite bar");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:3}, BarInfinite{dim:2,birth:3} ],
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
    &   precondition_to_make_new_lev_set_lower_none,
    analyze_dowker_dual,
    save_dir_opt,
    );


    println!("-------------------------------------------------------------------------------------");
    println!("Projective Plane Z2  (CW Structure with four squares, so 16 cells), and 1 dim-1 finite bar");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:3}, BarInfinite{dim:2,birth:3}],
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
    println!("Projective Plane Z2 (CW Structure with four squares, so 16 cells), and two consecutive dim-0 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:5}, BarInfinite{dim:2,birth:5} ],
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
    println!("Projective Plane Z2 (CW Structure with four squares, so 16 cells), and two consecutive touching dim-0 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:4}, BarInfinite{dim:2,birth:4} ],
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
    println!("Projective Plane Z2 (CW Structure with four squares, so 16 cells), and two consecutive dim-1 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:5}, BarInfinite{dim:2,birth:5} ],
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
    println!("Projective Plane Z2 (CW Structure with four squares, so 16 cells), and two consecutive touching dim-1 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:4}, BarInfinite{dim:2,birth:4} ],
        fin: vec![BarFinite{dim:1,birth:1,death:2}, BarFinite{dim:1,birth:2, death:3} ],
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
    println!("Projective Plane Z2 (CW Structure with four squares, so 16 cells), dim-0 and then dim-1 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:5}, BarInfinite{dim:2,birth:5}  ],
        fin: vec![BarFinite{dim:0,birth:1,death:2}, BarFinite{dim:1,birth:3, death:4} ],
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
    println!("Projective Plane Z2 (CW Structure with four squares, so 16 cells), dim-0 finite bar after loop ");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:1}, BarInfinite{dim:2,birth:4}  ],
        fin: vec![BarFinite{dim:0,birth:2,death:3} ],
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
    println!("Projective Plane Z2 (CW Structure with four squares, so 16 cells), dim-0 finite bar before and one dim-0 finite bar after loop ");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:3}, BarInfinite{dim:2,birth:6}  ],
        fin: vec![BarFinite{dim:0,birth:1,death:2}, BarFinite{dim:0,birth:4,death:5} ],
        ordinal: ordinate_unique_vals( & vec![0, 1, 2, 3, 4, 5,6] ),
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
Projective Plane Z2 (CW Structure with four squares, so 16 cells), and barcode with no finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [], [(0, true), (1, true)], [(1, true), (4, true)], [(0, true), (2, true)], [(2, true), (4, true)], [(2, true), (3, true)], [(2, true), (3, true)], [(1, true), (3, true)], [(1, true), (3, true)], [(5, true), (7, true), (9, true), (11, true)], [(6, true), (8, true), (10, true), (11, true)], [(6, true), (8, true), (9, true), (12, true)], [(5, true), (7, true), (10, true), (12, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 1 }, BarInfinite { dim: 2, birth: 2 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0, 1, 2], val_to_ord: {1: 1, 2: 2, 0: 0} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 93736.575578415s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 64256
number of facets (binned by dimension): [0, 0, 0, 0, 0, 0, 0, 64256]
number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 0, 0, 64256]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 43.570919007s"
number of polytopes (total): 2987264
number of polytopes (binned by dimension): [3312, 46384, 238576, 620656, 908896, 762336, 342848, 64256]
Time elapsed to compute binary-coeff polyhedral betti numbers 18.820972745s
betti numbers (of polyhedral complex, Z2 coefficients): [1, 2, 2, 1, 0, 0, 0, 0]
-------------------------------------------------------------------------------------
Projective Plane Z2 (CW Structure with four squares, so 16 cells), and 1 dim-0 finite bar

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [], [(0, true), (1, true)], [(1, true), (4, true)], [(0, true), (2, true)], [(2, true), (4, true)], [(2, true), (3, true)], [(2, true), (3, true)], [(1, true), (3, true)], [(1, true), (3, true)], [(5, true), (7, true), (9, true), (11, true)], [(6, true), (8, true), (10, true), (11, true)], [(6, true), (8, true), (9, true), (12, true)], [(5, true), (7, true), (10, true), (12, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 3 }, BarInfinite { dim: 2, birth: 3 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3], val_to_ord: {1: 1, 0: 0, 2: 2, 3: 3} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 44068.770210691s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 18880
number of facets (binned by dimension): [0, 0, 0, 0, 18880]
number of facets (binned by number of vertices): [0, 0, 0, 0, 18880]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 1.525315533s"
number of polytopes (total): 242376
number of polytopes (binned by dimension): [10412, 52084, 91896, 69104, 18880]
Time elapsed to compute binary-coeff polyhedral betti numbers 741.018954ms
betti numbers (of polyhedral complex, Z2 coefficients): [1, 3, 2, 0, 0]
-------------------------------------------------------------------------------------

Projective Plane Z2  (CW Structure with four squares, so 16 cells), and 1 dim-1 finite bar

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [], [(0, true), (1, true)], [(1, true), (4, true)], [(0, true), (2, true)], [(2, true), (4, true)], [(2, true), (3, true)], [(2, true), (3, true)], [(1, true), (3, true)], [(1, true), (3, true)], [(5, true), (7, true), (9, true), (11, true)], [(6, true), (8, true), (10, true), (11, true)], [(6, true), (8, true), (9, true), (12, true)], [(5, true), (7, true), (10, true), (12, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 3 }, BarInfinite { dim: 2, birth: 3 }], fin: [BarFinite { dim: 1, birth: 1, death: 2 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3], val_to_ord: {3: 3, 0: 0, 2: 2, 1: 1} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 2022.173455407s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 1216
number of facets (binned by dimension): [0, 0, 0, 0, 1216]
number of facets (binned by number of vertices): [0, 0, 0, 0, 1216]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 64.302196ms"
number of polytopes (total): 12576
number of polytopes (binned by dimension): [408, 2336, 4664, 3952, 1216]
Time elapsed to compute binary-coeff polyhedral betti numbers 7.558906ms
betti numbers (of polyhedral complex, Z2 coefficients): [4, 4, 0, 0, 0]
-------------------------------------------------------------------------------------
Projective Plane Z2 (CW Structure with four squares, so 16 cells), and two consecutive dim-0 finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [], [(0, true), (1, true)], [(1, true), (4, true)], [(0, true), (2, true)], [(2, true), (4, true)], [(2, true), (3, true)], [(2, true), (3, true)], [(1, true), (3, true)], [(1, true), (3, true)], [(5, true), (7, true), (9, true), (11, true)], [(6, true), (8, true), (10, true), (11, true)], [(6, true), (8, true), (9, true), (12, true)], [(5, true), (7, true), (10, true), (12, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 5 }, BarInfinite { dim: 2, birth: 5 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }, BarFinite { dim: 0, birth: 3, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4, 5], val_to_ord: {2: 2, 5: 5, 0: 0, 4: 4, 3: 3, 1: 1} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 35238.053315347s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 26816
number of facets (binned by dimension): [0, 0, 0, 26816]
number of facets (binned by number of vertices): [0, 0, 0, 26816]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 1.409904563s"
number of polytopes (total): 216152
number of polytopes (binned by dimension): [26400, 81288, 81648, 26816]
Time elapsed to compute binary-coeff polyhedral betti numbers 244.347355ms
betti numbers (of polyhedral complex, Z2 coefficients): [1, 57, 0, 0]
-------------------------------------------------------------------------------------

Projective Plane Z2 (CW Structure with four squares, so 16 cells), and two consecutive touching dim-0 finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [], [(0, true), (1, true)], [(1, true), (4, true)], [(0, true), (2, true)], [(2, true), (4, true)], [(2, true), (3, true)], [(2, true), (3, true)], [(1, true), (3, true)], [(1, true), (3, true)], [(5, true), (7, true), (9, true), (11, true)], [(6, true), (8, true), (10, true), (11, true)], [(6, true), (8, true), (9, true), (12, true)], [(5, true), (7, true), (10, true), (12, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 4 }, BarInfinite { dim: 2, birth: 4 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }, BarFinite { dim: 0, birth: 2, death: 3 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4], val_to_ord: {3: 3, 2: 2, 1: 1, 4: 4, 0: 0} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 22408.409787116s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 17856
number of facets (binned by dimension): [0, 0, 0, 17856]
number of facets (binned by number of vertices): [0, 0, 0, 17856]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 825.515603ms"
number of polytopes (total): 143288
number of polytopes (binned by dimension): [17232, 53816, 54384, 17856]
Time elapsed to compute binary-coeff polyhedral betti numbers 139.860551ms
betti numbers (of polyhedral complex, Z2 coefficients): [1, 57, 0, 0]
-------------------------------------------------------------------------------------
Projective Plane Z2 (CW Structure with four squares, so 16 cells), dim-0 and then dim-1 finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [], [(0, true), (1, true)], [(1, true), (4, true)], [(0, true), (2, true)], [(2, true), (4, true)], [(2, true), (3, true)], [(2, true), (3, true)], [(1, true), (3, true)], [(1, true), (3, true)], [(5, true), (7, true), (9, true), (11, true)], [(6, true), (8, true), (10, true), (11, true)], [(6, true), (8, true), (9, true), (12, true)], [(5, true), (7, true), (10, true), (12, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 5 }, BarInfinite { dim: 2, birth: 5 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }, BarFinite { dim: 1, birth: 3, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4, 5], val_to_ord: {4: 4, 2: 2, 5: 5, 0: 0, 3: 3, 1: 1} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 17993.161304558s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 20160
number of facets (binned by dimension): [0, 0, 0, 20160]
number of facets (binned by number of vertices): [0, 0, 0, 20160]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 939.125728ms"
number of polytopes (total): 131968
number of polytopes (binned by dimension): [12904, 45832, 53072, 20160]
Time elapsed to compute binary-coeff polyhedral betti numbers 117.633649ms
betti numbers (of polyhedral complex, Z2 coefficients): [4, 20, 0, 0]
-------------------------------------------------------------------------------------

Projective Plane Z2 (CW Structure with four squares, so 16 cells), dim-0 finite bar after loop

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [], [(0, true), (1, true)], [(1, true), (4, true)], [(0, true), (2, true)], [(2, true), (4, true)], [(2, true), (3, true)], [(2, true), (3, true)], [(1, true), (3, true)], [(1, true), (3, true)], [(5, true), (7, true), (9, true), (11, true)], [(6, true), (8, true), (10, true), (11, true)], [(6, true), (8, true), (9, true), (12, true)], [(5, true), (7, true), (10, true), (12, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 1 }, BarInfinite { dim: 2, birth: 4 }], fin: [BarFinite { dim: 0, birth: 2, death: 3 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4], val_to_ord: {1: 1, 0: 0, 2: 2, 3: 3, 4: 4} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 4514.753282913s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 94208
number of facets (binned by dimension): [0, 0, 0, 0, 0, 0, 94208]
number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 0, 94208]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 60.287462612s"
number of polytopes (total): 3545184
number of polytopes (binned by dimension): [25816, 205296, 656280, 1089408, 996288, 477888, 94208]
Time elapsed to compute binary-coeff polyhedral betti numbers 33.956221525s
betti numbers (of polyhedral complex, Z2 coefficients): [2, 5, 4, 1, 0, 0, 0]
-------------------------------------------------------------------------------------

Projective Plane Z2 (CW Structure with four squares, so 16 cells), dim-0 finite bar before and one dim-0 finite bar after loop

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [], [(0, true), (1, true)], [(1, true), (4, true)], [(0, true), (2, true)], [(2, true), (4, true)], [(2, true), (3, true)], [(2, true), (3, true)], [(1, true), (3, true)], [(1, true), (3, true)], [(5, true), (7, true), (9, true), (11, true)], [(6, true), (8, true), (10, true), (11, true)], [(6, true), (8, true), (9, true), (12, true)], [(5, true), (7, true), (10, true), (12, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 3 }, BarInfinite { dim: 2, birth: 6 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }, BarFinite { dim: 0, birth: 4, death: 5 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4, 5, 6], val_to_ord: {4: 4, 5: 5, 6: 6, 3: 3, 0: 0, 1: 1, 2: 2} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 141036.343159043s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 504832
number of facets (binned by dimension): [0, 0, 0, 0, 0, 504832]
number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 504832]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 185.94028335s"
number of polytopes (total): 11687888
number of polytopes (binned by dimension): [243592, 1443624, 3373600, 3895488, 2226752, 504832]
Time elapsed to compute binary-coeff polyhedral betti numbers 197.70346934s
betti numbers (of polyhedral complex, Z2 coefficients): [3, 30, 27, 0, 0, 0]

*/

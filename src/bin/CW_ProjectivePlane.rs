//  DESCRIPTION

//  CW RP2

// BECAUSE WE TAKE RATIONAL COEFFCIENTS, I THINK HOMOLOGY OF PROJ PLANE IS H_0=Q, H_n=0 for n>0. SO IN THE ASSOCIATED BARCODE THERE IS ONE INFINITE BAR. If not true we need to modify infinite bars in barcodes below





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

 /*
    println!("-------------------------------------------------------------------------------------");
    println!("Projective Plane (CW Structure with four squares, so 16 cells), and barcode with no finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0} ],
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
    println!("Projective Plane (CW Structure with four squares, so 16 cells), and 1 dim-0 finite bar");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0} ],
        fin: vec![BarFinite{dim:0,birth:1,death:2} ],
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
    println!("Projective Plane  (CW Structure with four squares, so 16 cells), and 1 dim-1 finite bar");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}],
        fin: vec![BarFinite{dim:1,birth:1,death:2} ],
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
    println!("Projective Plane (CW Structure with four squares, so 16 cells), and two consecutive dim-0 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0} ],
        fin: vec![BarFinite{dim:0,birth:1,death:2}, BarFinite{dim:0,birth:3,death:4} ],
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
    println!("Projective Plane (CW Structure with four squares, so 16 cells), and two consecutive touching dim-0 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0} ],
        fin: vec![BarFinite{dim:0,birth:1,death:2}, BarFinite{dim:0,birth:2, death:3} ],
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
    println!("Projective Plane (CW Structure with four squares, so 16 cells), and two consecutive dim-1 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0}],
        fin: vec![BarFinite{dim:1,birth:1,death:2}, BarFinite{dim:1,birth:3,death:4} ],
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
    println!("Projective Plane (CW Structure with four squares, so 16 cells), and two consecutive touching dim-1 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0} ],
        fin: vec![BarFinite{dim:1,birth:1,death:2}, BarFinite{dim:1,birth:2, death:3} ],
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

*/

    println!("-------------------------------------------------------------------------------------");
    println!("Projective Plane (CW Structure with four squares, so 16 cells), dim-0 and then dim-1 finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0} ],
        fin: vec![BarFinite{dim:0,birth:1,death:2}, BarFinite{dim:1,birth:3, death:4} ],
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

}


/*

Built boundary matrix
-------------------------------------------------------------------------------------
Projective Plane (CW Structure with four squares, so 16 cells), and barcode with no finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (4, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (4, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 }), (11, Ratio { numer: 1, denom: 1 })], [(6, Ratio { numer: 1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: -1, denom: 1 }), (11, Ratio { numer: -1, denom: 1 })], [(6, Ratio { numer: 1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: 1, denom: 1 }), (12, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: 1, denom: 1 }), (12, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0], val_to_ord: {0: 0} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 151.656677095s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 960
number of facets (binned by dimension): [0, 0, 0, 0, 0, 0, 960]
number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 0, 960]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 136.96533ms"
number of polytopes (total): 30835
number of polytopes (binned by dimension): [106, 1261, 5000, 9436, 9352, 4720, 960]
Time elapsed to compute binary-coeff polyhedral betti numbers 16.383181ms
betti numbers (of polyhedral complex, Z2 coefficients): [1, 0, 0, 0, 0, 0, 0]
-------------------------------------------------------------------------------------

Projective Plane (CW Structure with four squares, so 16 cells), and 1 dim-0 finite bar

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (4, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (4, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 }), (11, Ratio { numer: 1, denom: 1 })], [(6, Ratio { numer: 1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: -1, denom: 1 }), (11, Ratio { numer: -1, denom: 1 })], [(6, Ratio { numer: 1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: 1, denom: 1 }), (12, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: 1, denom: 1 }), (12, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2], val_to_ord: {0: 0, 1: 1, 2: 2} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 10923.25197433s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 18880
number of facets (binned by dimension): [0, 0, 0, 0, 0, 18880]
number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 18880]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 3.557397999s"
number of polytopes (total): 503072
number of polytopes (binned by dimension): [11348, 67224, 151436, 165432, 88752, 18880]
Time elapsed to compute binary-coeff polyhedral betti numbers 2.196099379s
betti numbers (of polyhedral complex, Z2 coefficients): [1, 3, 2, 0, 0, 0]
-------------------------------------------------------------------------------------

Projective Plane  (CW Structure with four squares, so 16 cells), and 1 dim-1 finite bar

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (4, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (4, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 }), (11, Ratio { numer: 1, denom: 1 })], [(6, Ratio { numer: 1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: -1, denom: 1 }), (11, Ratio { numer: -1, denom: 1 })], [(6, Ratio { numer: 1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: 1, denom: 1 }), (12, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: 1, denom: 1 }), (12, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [BarFinite { dim: 1, birth: 1, death: 2 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2], val_to_ord: {1: 1, 0: 0, 2: 2} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 130329.381744557s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 66432
number of facets (binned by dimension): [0, 0, 0, 0, 0, 2176, 0, 64256]
number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 2176, 0, 64256]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 43.255722094s"
number of polytopes (total): 3031424
number of polytopes (binned by dimension): [4000, 51104, 250784, 635840, 918080, 764512, 342848, 64256]
Time elapsed to compute binary-coeff polyhedral betti numbers 18.662131308s
betti numbers (of polyhedral complex, Z2 coefficients): [1, 2, 2, 1, 0, 0, 0, 0]
-------------------------------------------------------------------------------------
Projective Plane (CW Structure with four squares, so 16 cells), and two consecutive dim-0 finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (4, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (4, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 }), (11, Ratio { numer: 1, denom: 1 })], [(6, Ratio { numer: 1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: -1, denom: 1 }), (11, Ratio { numer: -1, denom: 1 })], [(6, Ratio { numer: 1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: 1, denom: 1 }), (12, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: 1, denom: 1 }), (12, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }, BarFinite { dim: 0, birth: 3, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4], val_to_ord: {1: 1, 3: 3, 0: 0, 2: 2, 4: 4} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 13058.841511906s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 26816
number of facets (binned by dimension): [0, 0, 0, 0, 26816]
number of facets (binned by number of vertices): [0, 0, 0, 0, 26816]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 3.149973297s"
number of polytopes (total): 469208
number of polytopes (binned by dimension): [31732, 122700, 176040, 111920, 26816]
Time elapsed to compute binary-coeff polyhedral betti numbers 1.69312394s
betti numbers (of polyhedral complex, Z2 coefficients): [1, 33, 0, 0, 0]
-------------------------------------------------------------------------------------

Projective Plane (CW Structure with four squares, so 16 cells), and two consecutive touching dim-0 finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (4, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (4, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 }), (11, Ratio { numer: 1, denom: 1 })], [(6, Ratio { numer: 1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: -1, denom: 1 }), (11, Ratio { numer: -1, denom: 1 })], [(6, Ratio { numer: 1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: 1, denom: 1 }), (12, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: 1, denom: 1 }), (12, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }, BarFinite { dim: 0, birth: 2, death: 3 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3], val_to_ord: {0: 0, 2: 2, 3: 3, 1: 1} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 8333.25060812s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 17856
number of facets (binned by dimension): [0, 0, 0, 0, 17856]
number of facets (binned by number of vertices): [0, 0, 0, 0, 17856]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 1.504254619s"
number of polytopes (total): 305192
number of polytopes (binned by dimension): [20092, 78964, 114632, 73648, 17856]
Time elapsed to compute binary-coeff polyhedral betti numbers 866.011194ms
betti numbers (of polyhedral complex, Z2 coefficients): [1, 33, 0, 0, 0]
-------------------------------------------------------------------------------------

Projective Plane (CW Structure with four squares, so 16 cells), and two consecutive dim-1 finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (4, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (4, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 }), (11, Ratio { numer: 1, denom: 1 })], [(6, Ratio { numer: 1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: -1, denom: 1 }), (11, Ratio { numer: -1, denom: 1 })], [(6, Ratio { numer: 1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: 1, denom: 1 }), (12, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: 1, denom: 1 }), (12, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [BarFinite { dim: 1, birth: 1, death: 2 }, BarFinite { dim: 1, birth: 3, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4], val_to_ord: {0: 0, 2: 2, 3: 3, 4: 4, 1: 1} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 31976.494813164s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 17664
number of facets (binned by dimension): [0, 0, 0, 0, 0, 0, 17664]
number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 0, 17664]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 9.210763794s"
number of polytopes (total): 617280
number of polytopes (binned by dimension): [4056, 33304, 109880, 188040, 177040, 87296, 17664]
Time elapsed to compute binary-coeff polyhedral betti numbers 1.40438765s
betti numbers (of polyhedral complex, Z2 coefficients): [4, 8, 4, 0, 0, 0, 0]

---------------------

Projective Plane (CW Structure with four squares, so 16 cells), and two consecutive touching dim-1 finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (4, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (4, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 }), (11, Ratio { numer: 1, denom: 1 })], [(6, Ratio { numer: 1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: -1, denom: 1 }), (11, Ratio { numer: -1, denom: 1 })], [(6, Ratio { numer: 1, denom: 1 }), (8, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: 1, denom: 1 }), (12, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (10, Ratio { numer: 1, denom: 1 }), (12, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }], fin: [BarFinite { dim: 1, birth: 1, death: 2 }, BarFinite { dim: 1, birth: 2, death: 3 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3], val_to_ord: {0: 0, 1: 1, 3: 3, 2: 2} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 62495.61269044s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 14592
number of facets (binned by dimension): [0, 0, 0, 0, 0, 0, 14592]
number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 0, 14592]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 5.773946455s"
number of polytopes (total): 488256
number of polytopes (binned by dimension): [2904, 24856, 84536, 148104, 142096, 71168, 14592]
Time elapsed to compute binary-coeff polyhedral betti numbers 1.045319527s
betti numbers (of polyhedral complex, Z2 coefficients): [4, 8, 4, 0, 0, 0, 0]
-------------------------------------------------------------------------------------




*/

//  DESCRIPTION

//  CW Cylinder

// BECAUSE WE TAKE RATIONAL COEFFCIENTS, I THINK HOMOLOGY OF MOBIUS IS H_0=Q, H_1=Q and H_n=0 for n>0. SO IN THE ASSOCIATED BARCODE THERE IS JUST 2 INFINITE BAR. If not true we need to modify infinite bars in barcodes below





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
    // Cylinder 12 cells-- 4 vertices, 6 edges,  2 faces
     ///////               0---------1---------0
     ///                   |         |         |
     ///                   |         |         |
     ///                   |         |         |
     ///                   2---------3---------2
     ///
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
                                    vec![       ( 0, one_pos.clone() ), ( 1, one_neg.clone() ) ],  //5:  1-cell  1->0
                                    vec![   ( 0, one_neg.clone() ),      ( 2, one_pos.clone() ) ],   //6:  1-cell  0->2
                                    vec![   ( 1, one_neg.clone() ),      ( 3, one_pos.clone() ) ],   //7:  1-cell  1->3
                                    vec![   ( 2, one_neg.clone() ),      ( 3, one_pos.clone() ) ],  //8:  1-cell  2->3
                                    vec![   ( 2, one_pos.clone() ),      ( 3, one_neg.clone() ) ],  //9:  1-cell  3->2
vec![ ( 4, one_pos.clone() ),  ( 6, one_neg.clone()),  ( 7, one_neg.clone() ), ( 8, one_pos.clone())], //10: left cell
vec![ ( 5, one_pos.clone() ),  ( 6, one_neg.clone()),  ( 7, one_neg.clone()), ( 9, one_neg.clone() ) ], //11: right cell
   ];

    println!("Built boundary matrix");

    //  The dimension of each cell.
    let cell_dims           =   vec![ 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2];

    //  Definition: cell_id_to_prereqs[ m ] = { the cells that must be added to the complex before M }
    let cell_id_to_prereqs  =   vec![
                                    vec![           ], // 0
                                    vec![           ],
                                    vec![           ],
                                    vec![            ],   //3
                                    vec![ 0, 1       ],
                                    vec![ 0, 1       ],
                                    vec![ 0, 2       ],
                                    vec![ 1, 3       ],     //8
                                    vec![ 2, 3       ],
                                    vec![ 2, 3       ],
                                    vec![ 4, 6, 7, 8 ],
                                    vec![ 5, 6, 7,9 ] //11
                                ];


    println!("-------------------------------------------------------------------------------------");
    println!("Cylinder (CW Structure with four squares, so 16 cells), and barcode with no finite bars");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0},  BarInfinite{dim:1,birth:1} ],
        fin: vec![ ],
        ordinal: ordinate_unique_vals( & vec![0,1] ),
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
    println!("Cylinder (CW Structure with four squares, so 16 cells), and barcode with dim-0 finite bar");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0},  BarInfinite{dim:1,birth:3} ],
        fin: vec![ BarFinite{dim:0,birth:1, death:2} ],
        ordinal: ordinate_unique_vals( & vec![0,1,2,3] ),
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
    println!("Cylinder (CW Structure with four squares, so 16 cells), and barcode with dim-1 finite bar");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0},  BarInfinite{dim:1,birth:3} ],
        fin: vec![ BarFinite{dim:1,birth:1, death:2} ],
        ordinal: ordinate_unique_vals( & vec![0,1,2,3] ),
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
    println!("Cylinder (CW Structure with four squares, so 16 cells), and barcode with two dim-0 finite bars before infinite");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0},  BarInfinite{dim:1,birth:5} ],
        fin: vec![ BarFinite{dim:0,birth:1, death:2}, BarFinite{dim:0,birth:3, death:4} ],
        ordinal: ordinate_unique_vals( & vec![0,1,2,3,4,5] ),
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
    println!("Cylinder (CW Structure with four squares, so 16 cells), and barcode with two overlapping dim-0 finite bars before infinite");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0},  BarInfinite{dim:1,birth:5} ],
        fin: vec![ BarFinite{dim:0,birth:1, death:3}, BarFinite{dim:0,birth:2, death:4} ],
        ordinal: ordinate_unique_vals( & vec![0,1,2,3,4,5] ),
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
    println!("Cylinder (CW Structure with four squares, so 16 cells), and barcode with dim-0 finite bar overlapping infinite");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0},  BarInfinite{dim:1,birth:2} ],
        fin: vec![ BarFinite{dim:0,birth:1, death:3} ],
        ordinal: ordinate_unique_vals( & vec![0,1,2,3] ),
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
    println!("Cylinder (CW Structure with four squares, so 16 cells), and barcode with dim-0 finite bar after infinite");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0},  BarInfinite{dim:1,birth:1} ],
        fin: vec![ BarFinite{dim:0,birth:2, death:3} ],
        ordinal: ordinate_unique_vals( & vec![0,1,2,3] ),
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
    println!("Cylinder (CW Structure with four squares, so 16 cells), and barcode with one dim-0 finite bar before and one after infinite loop");

    let barcode             =   Barcode{
        inf: vec![ BarInfinite{dim:0,birth:0},  BarInfinite{dim:1,birth:3} ],
        fin: vec![ BarFinite{dim:0,birth:1, death:2}, BarFinite{dim:0,birth:4, death:5} ],
        ordinal: ordinate_unique_vals( & vec![0,1,2,3,4,5] ),
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
Cylinder (CW Structure with four squares, so 16 cells), and barcode with no finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: -1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 1 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0, 1], val_to_ord: {0: 0, 1: 1} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 600.155348ms

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 640
number of facets (binned by dimension): [0, 0, 0, 0, 0, 640]
number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 640]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 77.368962ms"
number of polytopes (total): 16064
number of polytopes (binned by dimension): [328, 2040, 4760, 5352, 2944, 640]
Time elapsed to compute binary-coeff polyhedral betti numbers 19.920961ms
betti numbers (of polyhedral complex, Z2 coefficients): [1, 1, 0, 0, 0, 0]
-------------------------------------------------------------------------------------

Cylinder (CW Structure with four squares, so 16 cells), and barcode with dim-0 finite bar

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: -1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 3 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3], val_to_ord: {1: 1, 3: 3, 2: 2, 0: 0} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 6.299268189s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 5504
number of facets (binned by dimension): [0, 0, 0, 0, 5504]
number of facets (binned by number of vertices): [0, 0, 0, 0, 5504]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 472.504047ms"
number of polytopes (total): 89276
number of polytopes (binned by dimension): [5612, 22528, 33520, 22112, 5504]
Time elapsed to compute binary-coeff polyhedral betti numbers 280.799323ms
betti numbers (of polyhedral complex, Z2 coefficients): [1, 5, 0, 0, 0]
-------------------------------------------------------------------------------------

Cylinder (CW Structure with four squares, so 16 cells), and barcode with dim-1 finite bar

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: -1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 3 }], fin: [BarFinite { dim: 1, birth: 1, death: 2 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3], val_to_ord: {3: 3, 0: 0, 2: 2, 1: 1} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 204.107526ms

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 128
number of facets (binned by dimension): [0, 0, 0, 0, 128]
number of facets (binned by number of vertices): [0, 0, 0, 0, 128]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 8.015674ms"
number of polytopes (total): 1920
number of polytopes (binned by dimension): [96, 448, 736, 512, 128]
Time elapsed to compute binary-coeff polyhedral betti numbers 980.713µs
betti numbers (of polyhedral complex, Z2 coefficients): [2, 2, 0, 0, 0]
-------------------------------------------------------------------------------------

Cylinder (CW Structure with four squares, so 16 cells), and barcode with two dim-0 finite bars before infinite

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: -1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 5 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }, BarFinite { dim: 0, birth: 3, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4, 5], val_to_ord: {0: 0, 1: 1, 2: 2, 5: 5, 3: 3, 4: 4} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 1.506240248s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 3264
number of facets (binned by dimension): [0, 0, 0, 3264]
number of facets (binned by number of vertices): [0, 0, 0, 3264]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 145.398782ms"
number of polytopes (total): 30064
number of polytopes (binned by dimension): [4432, 11776, 10592, 3264]
Time elapsed to compute binary-coeff polyhedral betti numbers 42.358796ms
betti numbers (of polyhedral complex, Z2 coefficients): [2, 18, 0, 0]
-------------------------------------------------------------------------------------

Cylinder (CW Structure with four squares, so 16 cells), and barcode with two overlapping dim-0 finite bars before infinite

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: -1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 5 }], fin: [BarFinite { dim: 0, birth: 1, death: 3 }, BarFinite { dim: 0, birth: 2, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4, 5], val_to_ord: {3: 3, 5: 5, 0: 0, 4: 4, 1: 1, 2: 2} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 2.456192718s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 4160
number of facets (binned by dimension): [0, 0, 0, 4160]
number of facets (binned by number of vertices): [0, 0, 0, 4160]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 189.33587ms"
number of polytopes (total): 38512
number of polytopes (binned by dimension): [5712, 15104, 13536, 4160]
Time elapsed to compute binary-coeff polyhedral betti numbers 57.850606ms
betti numbers (of polyhedral complex, Z2 coefficients): [2, 18, 0, 0]

Cylinder (CW Structure with four squares, so 16 cells), and barcode with dim-0 finite bar overlapping infinite

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: -1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 2 }], fin: [BarFinite { dim: 0, birth: 1, death: 3 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3], val_to_ord: {2: 2, 1: 1, 3: 3, 0: 0} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 335.243621ms

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 1472
number of facets (binned by dimension): [0, 0, 0, 0, 1472]
number of facets (binned by number of vertices): [0, 0, 0, 0, 1472]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 114.107624ms"
number of polytopes (total): 26272
number of polytopes (binned by dimension): [1904, 7056, 9760, 6080, 1472]
Time elapsed to compute binary-coeff polyhedral betti numbers 33.007125ms
betti numbers (of polyhedral complex, Z2 coefficients): [4, 8, 4, 0, 0]
-------------------------------------------------------------------------------------
Cylinder (CW Structure with four squares, so 16 cells), and barcode with dim-0 finite bar after infinite

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: -1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 1 }], fin: [BarFinite { dim: 0, birth: 2, death: 3 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3], val_to_ord: {0: 0, 1: 1, 3: 3, 2: 2} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 24.916746ms

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 384
number of facets (binned by dimension): [0, 0, 0, 0, 384]
number of facets (binned by number of vertices): [0, 0, 0, 0, 384]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 29.062924ms"
number of polytopes (total): 6816
number of polytopes (binned by dimension): [496, 1824, 2528, 1584, 384]
Time elapsed to compute binary-coeff polyhedral betti numbers 5.692229ms
betti numbers (of polyhedral complex, Z2 coefficients): [2, 4, 2, 0, 0]
-------------------------------------------------------------------------------------

Cylinder (CW Structure with four squares, so 16 cells), and barcode with one dim-0 finite bar before and one after infinite loop

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, Ratio { numer: -1, denom: 1 }), (1, Ratio { numer: 1, denom: 1 })], [(0, Ratio { numer: 1, denom: 1 }), (1, Ratio { numer: -1, denom: 1 })], [(0, Ratio { numer: -1, denom: 1 }), (2, Ratio { numer: 1, denom: 1 })], [(1, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: -1, denom: 1 }), (3, Ratio { numer: 1, denom: 1 })], [(2, Ratio { numer: 1, denom: 1 }), (3, Ratio { numer: -1, denom: 1 })], [(4, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (8, Ratio { numer: 1, denom: 1 })], [(5, Ratio { numer: 1, denom: 1 }), (6, Ratio { numer: -1, denom: 1 }), (7, Ratio { numer: -1, denom: 1 }), (9, Ratio { numer: -1, denom: 1 })]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 3 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }, BarFinite { dim: 0, birth: 4, death: 5 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4, 5], val_to_ord: {2: 2, 0: 0, 1: 1, 4: 4, 3: 3, 5: 5} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 63.782188ms

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 768
number of facets (binned by dimension): [0, 0, 0, 768]
number of facets (binned by number of vertices): [0, 0, 0, 768]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 30.139869ms"
number of polytopes (total): 7208
number of polytopes (binned by dimension): [1120, 2840, 2480, 768]
Time elapsed to compute binary-coeff polyhedral betti numbers 5.498229ms
betti numbers (of polyhedral complex, Z2 coefficients): [4, 12, 0, 0]
*/

//  DESCRIPTION

//  CW Mobius Z2

// BECAUSE WE TAKE Z2 COEFFCIENTS, I THINK HOMOLOGY OF MOBIUS IS H_0=Z2, H_1=Z2 and H_n=0 for n>0. SO IN THE ASSOCIATED BARCODE THERE IS JUST 2 INFINITE BAR. If not true we need to modify infinite bars in barcodes below




/*
use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use phfibre::pipelines::boundary_matrix_pipeline;
use phfibre::phfibre::{ConditionNone, ConditionLowerStar, ConditionLowerEdge};
use solar::utilities::sequences_and_ordinals::{ordinate_unique_vals, BiMapSequential};
use solar::utilities::cell_complexes::simplices_unweighted::facets::ordered_subsimplices_up_thru_dim_concatenated_vec;
use solar::rings::ring::{Ring, Semiring};
*/

use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use phfibre::pipelines::boundary_matrix_pipeline;
use phfibre::phfibre::ConditionNone;
use solar::utilities::sequences_and_ordinals::ordinate_unique_vals;
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
    //type CeofficientRing    =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >; // the TYPE of ring object
    //let ring                =   CeofficientRing::new(); // the object itself

     //  ----------------------------------------------------------------------------------------------
    // Mobius 12 cells-- 4 vertices, 6 edges,  2 faces
     ///////               0---------1---------2
     ///                   |         |         |
     ///                   |         |         |
     ///                   |         |         |
     ///                   2---------3---------0
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
                                    vec![   ( 1, one_neg.clone() ),      ( 2, one_pos.clone() ) ],  //5:  1-cell  1->2
                                    vec![   ( 0, one_neg.clone() ),      ( 2, one_pos.clone() ) ],   //6:  1-cell  0->2
                                    vec![   ( 1, one_neg.clone() ),      ( 3, one_pos.clone() ) ],   //7:  1-cell  1->3
                                    vec![   ( 2, one_neg.clone() ),      ( 3, one_pos.clone() ) ],  //8:  1-cell  2->3
                                    vec![   ( 0, one_pos.clone() ),      ( 3, one_neg.clone() ) ],  //9:  1-cell  3->0
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
                                    vec![ 1, 2       ],
                                    vec![ 0, 2       ],
                                    vec![ 1, 3       ],  //7
                                    vec![ 2, 3       ],
                                    vec![ 0, 3       ],
                                    vec![ 4, 6, 7, 8 ],
                                    vec![ 5, 6, 7,9 ] //11
                                ];


    println!("-------------------------------------------------------------------------------------");
    println!("Mobius Z2 (CW Structure with four squares, so 16 cells), and barcode with no finite bars");

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
    println!("Mobius Z2 (CW Structure with four squares, so 16 cells), and barcode with dim-0 finite bar");

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
    println!("Mobius Z2 (CW Structure with four squares, so 16 cells), and barcode with dim-1 finite bar");

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
    println!("Mobius Z2 (CW Structure with four squares, so 16 cells), and barcode with two dim-0 finite bars before infinite");

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
    println!("Mobius Z2 (CW Structure with four squares, so 16 cells), and barcode with two overlapping dim-0 finite bars before infinite");

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
    println!("Mobius Z2 (CW Structure with four squares, so 16 cells), and barcode with dim-0 finite bar overlapping infinite");

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
    println!("Mobius Z2 (CW Structure with four squares, so 16 cells), and barcode with dim-0 finite bar after infinite");

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
    println!("Mobius Z2 (CW Structure with four squares, so 16 cells), and barcode with one dim-0 finite bar before and one after infinite loop");

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
Mobius Z2 (CW Structure with four squares, so 16 cells), and barcode with no finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, true), (1, true)], [(1, true), (2, true)], [(0, true), (2, true)], [(1, true), (3, true)], [(2, true), (3, true)], [(0, true), (3, true)], [(4, true), (6, true), (7, true), (8, true)], [(5, true), (6, true), (7, true), (9, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 1 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0, 1], val_to_ord: {0: 0, 1: 1} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 374.638215ms

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 768
number of facets (binned by dimension): [0, 0, 0, 0, 0, 768]
number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 768]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 96.715065ms"
number of polytopes (total): 19112
number of polytopes (binned by dimension): [372, 2396, 5664, 6392, 3520, 768]
Time elapsed to compute binary-coeff polyhedral betti numbers 27.149008ms
betti numbers (of polyhedral complex, Z2 coefficients): [1, 1, 0, 0, 0, 0]
-------------------------------------------------------------------------------------
Mobius Z2 (CW Structure with four squares, so 16 cells), and barcode with dim-0 finite bar

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, true), (1, true)], [(1, true), (2, true)], [(0, true), (2, true)], [(1, true), (3, true)], [(2, true), (3, true)], [(0, true), (3, true)], [(4, true), (6, true), (7, true), (8, true)], [(5, true), (6, true), (7, true), (9, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 3 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3], val_to_ord: {2: 2, 0: 0, 1: 1, 3: 3} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 8.224600326s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 6528
number of facets (binned by dimension): [0, 0, 0, 0, 6528]
number of facets (binned by number of vertices): [0, 0, 0, 0, 6528]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 593.092938ms"
number of polytopes (total): 108412
number of polytopes (binned by dimension): [7028, 27764, 40644, 26448, 6528]
Time elapsed to compute binary-coeff polyhedral betti numbers 684.080811ms
betti numbers (of polyhedral complex, Z2 coefficients): [1, 13, 0, 0, 0]
-------------------------------------------------------------------------------------

Mobius Z2 (CW Structure with four squares, so 16 cells), and barcode with dim-1 finite bar

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, true), (1, true)], [(1, true), (2, true)], [(0, true), (2, true)], [(1, true), (3, true)], [(2, true), (3, true)], [(0, true), (3, true)], [(4, true), (6, true), (7, true), (8, true)], [(5, true), (6, true), (7, true), (9, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 3 }], fin: [BarFinite { dim: 1, birth: 1, death: 2 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3], val_to_ord: {1: 1, 0: 0, 2: 2, 3: 3} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 165.361671ms

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 128
number of facets (binned by dimension): [0, 0, 0, 0, 128]
number of facets (binned by number of vertices): [0, 0, 0, 0, 128]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 8.080527ms"
number of polytopes (total): 1920
number of polytopes (binned by dimension): [96, 448, 736, 512, 128]
Time elapsed to compute binary-coeff polyhedral betti numbers 1.07372ms
betti numbers (of polyhedral complex, Z2 coefficients): [2, 2, 0, 0, 0]
-------------------------------------------------------------------------------------
Mobius Z2 (CW Structure with four squares, so 16 cells), and barcode with two dim-0 finite bars before infinite

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, true), (1, true)], [(1, true), (2, true)], [(0, true), (2, true)], [(1, true), (3, true)], [(2, true), (3, true)], [(0, true), (3, true)], [(4, true), (6, true), (7, true), (8, true)], [(5, true), (6, true), (7, true), (9, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 5 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }, BarFinite { dim: 0, birth: 3, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4, 5], val_to_ord: {5: 5, 2: 2, 1: 1, 0: 0, 3: 3, 4: 4} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 2.368035082s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 4032
number of facets (binned by dimension): [0, 0, 0, 4032]
number of facets (binned by number of vertices): [0, 0, 0, 4032]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 184.829648ms"
number of polytopes (total): 37848
number of polytopes (binned by dimension): [5744, 14904, 13168, 4032]
Time elapsed to compute binary-coeff polyhedral betti numbers 58.099583ms
betti numbers (of polyhedral complex, Z2 coefficients): [1, 25, 0, 0]
-------------------------------------------------------------------------------------

Mobius Z2 (CW Structure with four squares, so 16 cells), and barcode with two overlapping dim-0 finite bars before infinite

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, true), (1, true)], [(1, true), (2, true)], [(0, true), (2, true)], [(1, true), (3, true)], [(2, true), (3, true)], [(0, true), (3, true)], [(4, true), (6, true), (7, true), (8, true)], [(5, true), (6, true), (7, true), (9, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 5 }], fin: [BarFinite { dim: 0, birth: 1, death: 3 }, BarFinite { dim: 0, birth: 2, death: 4 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4, 5], val_to_ord: {1: 1, 0: 0, 5: 5, 4: 4, 3: 3, 2: 2} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 3.616433793s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 4992
number of facets (binned by dimension): [0, 0, 0, 4992]
number of facets (binned by number of vertices): [0, 0, 0, 4992]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 231.303478ms"
number of polytopes (total): 47096
number of polytopes (binned by dimension): [7184, 18568, 16352, 4992]
Time elapsed to compute binary-coeff polyhedral betti numbers 75.878038ms
betti numbers (of polyhedral complex, Z2 coefficients): [1, 25, 0, 0]
-------------------------------------------------------------------------------------
Mobius Z2 (CW Structure with four squares, so 16 cells), and barcode with dim-0 finite bar overlapping infinite

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, true), (1, true)], [(1, true), (2, true)], [(0, true), (2, true)], [(1, true), (3, true)], [(2, true), (3, true)], [(0, true), (3, true)], [(4, true), (6, true), (7, true), (8, true)], [(5, true), (6, true), (7, true), (9, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 2 }], fin: [BarFinite { dim: 0, birth: 1, death: 3 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3], val_to_ord: {1: 1, 2: 2, 0: 0, 3: 3} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 96.319111ms

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 768
number of facets (binned by dimension): [0, 0, 0, 0, 768]
number of facets (binned by number of vertices): [0, 0, 0, 0, 768]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 56.777681ms"
number of polytopes (total): 13224
number of polytopes (binned by dimension): [936, 3492, 4908, 3120, 768]
Time elapsed to compute binary-coeff polyhedral betti numbers 10.969264ms
betti numbers (of polyhedral complex, Z2 coefficients): [8, 8, 0, 0, 0]
-------------------------------------------------------------------------------------

Mobius Z2 (CW Structure with four squares, so 16 cells), and barcode with dim-0 finite bar after infinite

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, true), (1, true)], [(1, true), (2, true)], [(0, true), (2, true)], [(1, true), (3, true)], [(2, true), (3, true)], [(0, true), (3, true)], [(4, true), (6, true), (7, true), (8, true)], [(5, true), (6, true), (7, true), (9, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 1 }], fin: [BarFinite { dim: 0, birth: 2, death: 3 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3], val_to_ord: {3: 3, 2: 2, 0: 0, 1: 1} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 9.533502ms

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 192
number of facets (binned by dimension): [0, 0, 0, 0, 192]
number of facets (binned by number of vertices): [0, 0, 0, 0, 192]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 12.727491ms"
number of polytopes (total): 3192
number of polytopes (binned by dimension): [216, 828, 1188, 768, 192]
Time elapsed to compute binary-coeff polyhedral betti numbers 1.820713ms
betti numbers (of polyhedral complex, Z2 coefficients): [4, 4, 0, 0, 0]
-------------------------------------------------------------------------------------
Mobius Z2 (CW Structure with four squares, so 16 cells), and barcode with one dim-0 finite bar before and one after infinite loop

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [(0, true), (1, true)], [(1, true), (2, true)], [(0, true), (2, true)], [(1, true), (3, true)], [(2, true), (3, true)], [(0, true), (3, true)], [(4, true), (6, true), (7, true), (8, true)], [(5, true), (6, true), (7, true), (9, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 3 }], fin: [BarFinite { dim: 0, birth: 1, death: 2 }, BarFinite { dim: 0, birth: 4, death: 5 }], ordinal: BiMapSequential { ord_to_val: [0, 1, 2, 3, 4, 5], val_to_ord: {1: 1, 2: 2, 5: 5, 4: 4, 0: 0, 3: 3} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 58.22658ms

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 768
number of facets (binned by dimension): [0, 0, 0, 768]
number of facets (binned by number of vertices): [0, 0, 0, 768]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 31.260592ms"
number of polytopes (total): 7296
number of polytopes (binned by dimension): [1152, 2880, 2496, 768]
Time elapsed to compute binary-coeff polyhedral betti numbers 4.443247ms
betti numbers (of polyhedral complex, Z2 coefficients): [8, 8, 0, 0]
*/

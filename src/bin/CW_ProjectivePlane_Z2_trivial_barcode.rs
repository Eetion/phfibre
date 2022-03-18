//  DESCRIPTION

//  CW RP2 Z2 coeff just trivial barcode

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
        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:1}, BarInfinite{dim:2,birth:1} ],
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

}

/*
Built boundary matrix
-------------------------------------------------------------------------------------
Projective Plane Z2 (CW Structure with four squares, so 16 cells), and barcode with no finite bars

All filters

NEW COMPUTATION
---------------------------------------------------------------------------------------------
BOUNDARY MATRIX
boundary matrix sparse columns: [[], [], [], [], [], [(0, true), (1, true)], [(1, true), (4, true)], [(0, true), (2, true)], [(2, true), (4, true)], [(2, true), (3, true)], [(2, true), (3, true)], [(1, true), (3, true)], [(1, true), (3, true)], [(5, true), (7, true), (9, true), (11, true)], [(6, true), (8, true), (10, true), (11, true)], [(6, true), (8, true), (9, true), (12, true)], [(5, true), (7, true), (10, true), (12, true)]]
BARCODE
barcode: Barcode { inf: [BarInfinite { dim: 0, birth: 0 }, BarInfinite { dim: 1, birth: 1 }, BarInfinite { dim: 2, birth: 1 }], fin: [], ordinal: BiMapSequential { ord_to_val: [0, 1], val_to_ord: {1: 1, 0: 0} } }
TIME TO COMPUTE FIBRE FACETS
Time elapsed to compute facets of PH fibre: 2045.028814351s

ANALYSIS
Each polytope facet has been checked for compatiblity with the given barcode.
POLYHEDRAL COMPLEX FACETS
number of facets (total): 960
number of facets (binned by dimension): [0, 0, 0, 0, 0, 960]
number of facets (binned by number of vertices): [0, 0, 0, 0, 0, 960]
POLYHEDRAL COMPLEX CELLS
"Time elapsed to compute the faces of PH fibre (given the facets): 74.402006ms"
number of polytopes (total): 15417
number of polytopes (binned by dimension): [105, 1156, 3844, 5592, 3760, 960]
Time elapsed to compute binary-coeff polyhedral betti numbers 10.393109ms
betti numbers (of polyhedral complex, Z2 coefficients): [1, 1, 1, 0, 0, 0]
*/

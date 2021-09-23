//  DESCRIPTION

//  These examples reproduce some results from Leygonie and Tillmann, "The Fier of Persistent Homology for Simplicial Complexes"


//  DEPENDENCIES

use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use phfibre::pipelines::boundary_matrix_pipeline;
use phfibre::phfibre::ConditionNone;
use solar::utilities::sequences_and_ordinals::ordinate_unique_vals;
use solar::cell_complexes::simplices_unweighted::facets::ordered_subsimplices_up_thru_dim_concatenated_vec; 
use solar::rings::ring::{Ring, Semiring};

//  CODE

fn main() {

    let save_dir_opt        =   None;  // we will not save any files


    //  ----------------------------------------------------------------------------------------------    
    //  DEFINE COEFFICIENT RING
    //  ----------------------------------------------------------------------------------------------    

    //  This is the ring of rationals.
    type CeofficientRing    =   solar::rings::field_prime::GF2;
    let ring                =   solar::rings::field_prime::GF2{}; // the object itself
    
    


    //  ----------------------------------------------------------------------------------------------    
    //  TWO TORUS -- 6 cells, no finite bars
    //  ----------------------------------------------------------------------------------------------

    //  BOUNDARY MATRIX

    //  Define +1 and -1
    let one_pos             =   CeofficientRing::one();
    let one_neg             =   ring.negate( one_pos.clone() );

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
    let analyze_dowker_dual =   true;

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



  
}
 

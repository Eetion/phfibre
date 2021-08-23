//  DESCRIPTION

//  These examples reproduce some results from Leygonie and Tillmann, "The Fier of Persistent Homology for Simplicial Complexes"


//  DEPENDENCIES

use phfibre::intervals_and_ordinals::{Barcode, BarFinite, BarInfinite};
use phfibre::pipelines::simplex_pipeline;
use phfibre::phfibre::ConditionNone;
use solar::utilities::sequences_and_ordinals::ordinate_unique_vals;
use solar::cell_complexes::simplices_unweighted::facets::ordered_subsimplices_up_thru_dim_concatenated_vec; 


//  CODE

fn main() {

    //  ----------------------------------------------------------------------------------------------    
    //  FILTRATION RESTRICTIONS (LOWER-STAR, LOWER-EDGE)
    //  ----------------------------------------------------------------------------------------------

    //  By using the following struct, we impose not "lower-X" conditions
    let precondition_to_make_new_lev_set     =   ConditionNone{};


    //  ----------------------------------------------------------------------------------------------    
    //  B_3^1
    //  ----------------------------------------------------------------------------------------------


    //  Define the base space, barcode, and ring
    let complex_facets      =   vec![  vec![0,1,2] ]; // this represents a list containing exactly one simplex, whose vertices are {0, 1, 2} 
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 1);    // build all faces up to dimension 1

    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:2} ],
                                    fin: vec![ BarFinite{dim:0,birth:0,death:1} ],
                                    ordinal: ordinate_unique_vals( & vec![0, 1, 2] ),
                                };

    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

    println!("");    
    println!("B_3^1:");        
    simplex_pipeline(
        &   simplex_sequence,
        &   barcode,
        &   ring,
        &   precondition_to_make_new_lev_set,
    );

    //  RESULTS
    //  -------

    // POLYHEDRAL COMPLEX FACETS
    // number of facets (total): 18
    // number of facets (binned by dimension): [0, 18]
    // POLYHEDRAL COMPLEX CELLS
    // number of polytopes (total): 36
    // number of polytopes (binned by dimension): [18, 18]
    // betti numbers (of polyhedral complex): [1, 1]
    // INTERSECTIONS
    // number of pairs of intersecting facets, binned by the dimension of the intersection polytope: [18, 0]
    // DOWKER NERVE COMPLEX FACETS
    // number of dowker nerve complex facets (total): 36
    // number of dowker nerve complex facets (binned by dimension): [18, 18]
    // DOWKER NERVE COMPLEX CELLS
    // number of nerve dowker complex cells (total): 36
    // number of nerve dowker complex cells (binned by dimension): [18, 18]
    // betti numbers (of dowker nerve): [1, 1]   

    
    //  ----------------------------------------------------------------------------------------------    
    //  B_3^2
    //  ----------------------------------------------------------------------------------------------


    //  Define the base space, barcode, and ring
    let complex_facets      =   vec![  vec![0,1,2] ];
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 1);    

    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:2} ],
                                    fin: vec![ BarFinite{dim:0,birth:1,death:2}],
                                    ordinal: ordinate_unique_vals( & vec![0, 1, 2] ),
                                };

    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

    println!("");    
    println!("B_3^2:");        
    simplex_pipeline(
        &   simplex_sequence,
        &   barcode,
        &   ring,
        &   precondition_to_make_new_lev_set,
    );

    //  RESULTS
    //  -------

    // POLYHEDRAL COMPLEX FACETS
    // number of facets (total): 18
    // number of facets (binned by dimension): [0, 18]
    // POLYHEDRAL COMPLEX CELLS
    // number of polytopes (total): 36
    // number of polytopes (binned by dimension): [18, 18]
    // betti numbers (of polyhedral complex): [1, 1]
    // INTERSECTIONS
    // number of pairs of intersecting facets, binned by the dimension of the intersection polytope: [18, 0]
    // DOWKER NERVE COMPLEX FACETS
    // number of dowker nerve complex facets (total): 36
    // number of dowker nerve complex facets (binned by dimension): [18, 18]
    // DOWKER NERVE COMPLEX CELLS
    // number of nerve dowker complex cells (total): 36
    // number of nerve dowker complex cells (binned by dimension): [18, 18]
    // betti numbers (of dowker nerve): [1, 1]      
    
    
    //  ----------------------------------------------------------------------------------------------    
    //  B_3^3
    //  ----------------------------------------------------------------------------------------------


    //  Define the base space, barcode, and ring
    let complex_facets      =   vec![  vec![0,1,2] ];
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 1);    

    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:2} ],
                                    fin: vec![ BarFinite{dim:0,birth:1,death:2}, BarFinite{dim:0,birth:0,death:1}],
                                    ordinal: ordinate_unique_vals( & vec![0, 1, 2] ),
                                };

    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

    println!("");    
    println!("B_3^3:");
    simplex_pipeline(
        &   simplex_sequence,
        &   barcode,
        &   ring,
        &   precondition_to_make_new_lev_set,
    );

    //  RESULTS
    //  -------

    // POLYHEDRAL COMPLEX FACETS
    // number of facets (total): 3
    // number of facets (binned by dimension): [3]
    // POLYHEDRAL COMPLEX CELLS
    // number of polytopes (total): 3
    // number of polytopes (binned by dimension): [3]
    // betti numbers (of polyhedral complex): [3]
    // INTERSECTIONS
    // number of pairs of intersecting facets, binned by the dimension of the intersection polytope: [0]
    // DOWKER NERVE COMPLEX FACETS
    // number of dowker nerve complex facets (total): 3
    // number of dowker nerve complex facets (binned by dimension): [3]
    // DOWKER NERVE COMPLEX CELLS
    // number of nerve dowker complex cells (total): 3
    // number of nerve dowker complex cells (binned by dimension): [3]
    // betti numbers (of dowker nerve): [3]    
    

    //  ----------------------------------------------------------------------------------------------    
    //  B_2^1
    //  ----------------------------------------------------------------------------------------------


    //  Define the base space, barcode, and ring
    let complex_facets      =   vec![  vec![0,1,2] ];
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 1);    

    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:3} ],
                                    fin: vec![ BarFinite{dim:0,birth:1,death:2} ],
                                    ordinal: ordinate_unique_vals( & vec![0, 1, 2, 3] ),
                                };

    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new();

    println!("");    
    println!("B_2^1:");    
    simplex_pipeline(
        &   simplex_sequence,
        &   barcode,
        &   ring,
        &   precondition_to_make_new_lev_set,
    );

    //  RESULTS
    //  -------
       
        

    //  ----------------------------------------------------------------------------------------------    
    //  B_4^1
    //  ----------------------------------------------------------------------------------------------


    //  Define the base space, barcode, and ring
    let complex_facets      =   vec![  vec![0,1,2] ]; // triangle with edge attached
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 1);    // build all faces up to dimension 1

    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:1} ],
                                    fin: vec![ ],
                                    ordinal: ordinate_unique_vals( & vec![ 0, 1 ] ),
                                };

    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new(); // this just defines the rational numbers

    println!("");    
    println!("B_4^1:");        
    simplex_pipeline(
        &   simplex_sequence,
        &   barcode,
        &   ring,
        &   precondition_to_make_new_lev_set,
    );

    //  RESULTS
    //  -------

    // BASE SPACE
    // simplices of the base space: [[0], [1], [2], [0, 1], [0, 2], [1, 2]]
    // POLYHEDRAL COMPLEX FACETS
    // number of facets (total): 12
    // number of facets (binned by dimension): [0, 0, 12]
    // number of facets (binned by number of vertices): [0, 0, 12]
    // POLYHEDRAL COMPLEX CELLS
    // number of polytopes (total): 42
    // number of polytopes (binned by dimension): [9, 21, 12]
    // betti numbers (of polyhedral complex): [1, 1, 0]
    // INTERSECTIONS
    // number of pairs of intersecting facets, binned by the dimension of the intersection polytope: [24, 15, 0]
    // DOWKER NERVE COMPLEX FACETS
    // number of dowker nerve complex facets (total): 12
    // number of dowker nerve complex facets (binned by dimension): [0, 0, 12]
    // DOWKER NERVE COMPLEX CELLS
    // number of nerve dowker complex cells (total): 42
    // number of nerve dowker complex cells (binned by dimension): [9, 21, 12]
    // betti numbers (of dowker nerve): [1, 1, 0]


    //  ----------------------------------------------------------------------------------------------    
    //  B_5
    //  ----------------------------------------------------------------------------------------------


    //  Define the base space, barcode, and ring
    let complex_facets      =   vec![  vec![0,1,2] ]; // triangle with edge attached
    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 1);    // build all faces up to dimension 1

    let barcode             =   Barcode{
                                    inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:0} ],
                                    fin: vec![ ],
                                    ordinal: ordinate_unique_vals( & vec![ 0 ] ),
                                };

    let ring                =   solar::rings::ring_native::NativeDivisionRing::< num::rational::Ratio<i64> >::new(); // this just defines the rational numbers

    println!("");    
    println!("B_5:");        
    simplex_pipeline(
        &   simplex_sequence,
        &   barcode,
        &   ring,
        &   precondition_to_make_new_lev_set,
    );

    //  RESULTS
    //  -------    

    // BASE SPACE
    // simplices of the base space: [[0], [1], [2], [0, 1], [0, 2], [1, 2]]
    // POLYHEDRAL COMPLEX FACETS
    // number of facets (total): 1
    // number of facets (binned by dimension): [1]
    // number of facets (binned by number of vertices): [1]
    // POLYHEDRAL COMPLEX CELLS
    // number of polytopes (total): 1
    // number of polytopes (binned by dimension): [1]
    // betti numbers (of polyhedral complex): [1]
    // INTERSECTIONS
    // number of pairs of intersecting facets, binned by the dimension of the intersection polytope: [0]
    // DOWKER NERVE COMPLEX FACETS
    // number of dowker nerve complex facets (total): 1
    // number of dowker nerve complex facets (binned by dimension): [1]
    // DOWKER NERVE COMPLEX CELLS
    // number of nerve dowker complex cells (total): 1
    // number of nerve dowker complex cells (binned by dimension): [1]
    // betti numbers (of dowker nerve): [1]

}
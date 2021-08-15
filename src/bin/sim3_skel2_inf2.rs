
use solar::utilities::sequences_and_ordinals::BiMapSequential;
use phfibre::phfibre::{Node, explore, verify_that_barcode_is_compatible};
use phfibre::intervals_and_ordinals::{Barcode, BarcodeInverse, to_ordered_float};
use phfibre::polytope::faces::{polys_faces};
use solar::cell_complexes::simplices_unweighted::facets::{    
    ordered_subsimplices_up_thru_dim_concatenated_vec     }; 
use solar::cell_complexes::simplices_unweighted::boundary_matrices::{    
    boundary_matrix_from_complex_facets     }; 
use num::rational::Ratio;
use std::iter::FromIterator;


type RingEltRational = Ratio<i64>;
type RingOpRational = solar::rings::ring_native::NativeDivisionRing<Ratio<i64>>;




fn main() { 

    //  DEFINE COEFFICIENT FIELD
    //  --------------------------------------------
    let ring                =   RingOpRational::new();

    //  DEFINE THE CELL DIMENSIONS + BOUNDARY MATRIX
    //  --------------------------------------------

    let complex_facets      =   vec![  vec![0, 1, 2, 3] ];

    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 2);    
    let cell_dims: Vec<_>   =   simplex_sequence.iter().map(|x| x.len()-1 ).collect();

    let bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence );
    let boundary            =   boundary_matrix_from_complex_facets(bimap_sequential, ring.clone());

    println!("{:?}", & cell_dims);
    println!("{:?}", & boundary);    


    //  DEFINE THE BARCODE + INVERSE BARCODE
    //  ------------------------------------

    let barcode_inf_dim     =   vec![0, 2];    
    let barcode_inf_brn     =   vec![0, 1];

    let barcode_fin_dim     =   Vec::new();    
    let barcode_fin_brn     =   Vec::new();
    let barcode_fin_die     =   Vec::new();


    let barcode             =   Barcode::new(
                                    barcode_inf_dim,
                                    barcode_inf_brn,
                                    barcode_fin_dim,
                                    barcode_fin_brn,
                                    barcode_fin_die
                                );   
    
    let barcode_inverse     =   BarcodeInverse::from_barcode( & barcode );
            

    //  ALLOW FOR DEGENERATE CELLS AFTER THE LAST FINITE BARCODE ENDPOINT
    //  -----------------------------------------------------------------

    // let last_must_be_crit   =   false;


    //  DEFINE THE ROOT NODE + RESULTS VECTOR
    //  -------------------------------------

    let root =  Node::make_root(
                        boundary.clone(),
                    &   barcode,
                    &   barcode_inverse,
                    &   cell_dims,  
                        RingOpRational::new(),
                        // last_must_be_crit,
                );

    let mut results         =   Vec::new();                                

    //  GATHER RESULTS
    //  --------------

    explore( & root, &mut results );

    //  DOUBLE CHECK RESULTS
    //  --------------------

    for (result_count, result) in results.iter().cloned().enumerate() {
        verify_that_barcode_is_compatible( 
            & root,
            & result
        );
    }    

    //  DISPLAY RESULTS
    //  --------------------    

    println!("\nRESULTS\n");

    println!("number of facets: {:?}", results.len() ); 

    let dim_top             =   results.iter().map(|x| x.dim_cellagnostic().unwrap() ).max().unwrap();
    let mut dim_to_polycount   =   Vec::from_iter( std::iter::repeat(0).take(dim_top + 1) );
    for dim in 0 .. dim_top + 1 {
        dim_to_polycount[ dim ]    =   polys_faces( &results, dim ).len();
    }
    println!("number of polytopes total: {:?}", dim_to_polycount.iter().sum::<usize>() );
    println!("number of polytopes by dimension (total): {:?}", dim_to_polycount );

    let mut dim_to_facetcount   =   Vec::from_iter( std::iter::repeat(0).take(dim_top + 1) );
    for facet in results.iter() {
        dim_to_facetcount[ facet.dim_cellagnostic().unwrap() ] += 1;
    }
    println!("number of facets total: {:?}", dim_to_facetcount.iter().sum::<usize>() );    
    println!("number of facets by dimension (total): {:?}", dim_to_facetcount );           



}  


// number of facets: 1920
// number of polytopes total: 39566
// number of polytopes by dimension (total): [64, 858, 4480, 10860, 13320, 8064, 1920]
// number of facets total: 1920
// number of facets by dimension (total): [0, 0, 0, 0, 0, 0, 1920]
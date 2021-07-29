
use solar::utilities::index::BiMapSequential;
use phfibre::phfibre::{Node, explore};
use phfibre::intervals_and_ordinals::{Barcode, BarcodeInverse, to_ordered_float};
use phfibre::facets::{ordered_subsimplices_up_to_dim_concatenated, boundary_matrix_from_complex_facets};
use num::rational::Ratio;



type RingEltRational = Ratio<i64>;
type RingOpRational = solar::rings::ring_native::NativeDivisionRing<Ratio<i64>>;



fn main() { 

    //  DEFINE COEFFICIENT FIELD
    //  --------------------------------------------
    let ring                =   RingOpRational::new();

    //  DEFINE THE CELL DIMENSIONS + BOUNDARY MATRIX
    //  --------------------------------------------

    let complex_facets      =   vec![  vec![0,1,2] ];

    let simplex_sequence    =   ordered_subsimplices_up_to_dim_concatenated( &complex_facets, 1);    
    let cell_dims: Vec<_>   =   simplex_sequence.iter().map(|x| x.len()-1 ).collect();

    let bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence );
    let boundary            =   boundary_matrix_from_complex_facets(bimap_sequential, ring.clone());


    //  DEFINE THE BARCODE + INVERSE BARCODE
    //  ------------------------------------

    let barcode_inf_dim     =   vec![0];    
    let barcode_inf_brn     =   to_ordered_float( & vec![0.] );

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
    
    let last_must_be_crit   =   false;    


    //  DEFINE THE ROOT NODE + RESULTS VECTOR
    //  -------------------------------------

    let root =  Node::make_root(
                        boundary.clone(),
                    &   barcode,
                    &   barcode_inverse,
                    &   cell_dims,  
                        RingOpRational::new(),
                        last_must_be_crit,
                );

    let mut results         =   Vec::new();                                

    //  GATHER RESULTS
    //  --------------

    // println!("{:?}", &root );

    explore( & root, &mut results );

    println!("\nRESULTS\n");
    
    for result in results.iter().cloned().enumerate() {
        println!("{:?}", result)
    }    


}  
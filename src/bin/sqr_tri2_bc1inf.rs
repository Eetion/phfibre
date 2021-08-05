
use solar::utilities::index::{BiMapSequential, histogram};
use phfibre::phfibre::{Node, explore, verify_that_barcode_is_compatible};
use phfibre::intervals_and_ordinals::{Barcode, BarcodeInverse, to_ordered_float};
use solar::cell_complexes::simplices_unweighted::maximal_cliques::{    
    ordered_subsimplices_up_thru_dim_concatenated_vec 
}; 
use solar::cell_complexes::simplices_unweighted::boundary_matrices::{    
    boundary_matrix_from_complex_facets 
};   
use num::rational::Ratio;
use ordered_float::OrderedFloat;


type RingEltRational = OrderedFloat<f64>;
type RingOpRational = solar::rings::ring_native::NativeDivisionRing<RingEltRational>;



fn main() { 

    //  DEFINE COEFFICIENT FIELD
    //  --------------------------------------------
    let ring                =   RingOpRational::new();

    //  DEFINE THE CELL DIMENSIONS + BOUNDARY MATRIX
    //  --------------------------------------------

    let complex_facets      =   vec![  vec![0,1,2], vec![1,2,3] ];

    let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 2);    
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

    //  REPORT
    //  ------

    println!("number of results: {:?}", results.len() ); 

    let hist = histogram( results.iter().map(|x| x.dim_cellagnostic().unwrap() ) );
    println!("number of polytopes by dimension: {:?}", &hist );

    // number of results: 1366
    // number of polytopes by dimension: [1, 31, 210, 492, 472, 160]    

    //  CHECK
    //  ------    

    

}  
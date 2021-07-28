
use crate::phfibre::{Node, explore};
use crate::intervals_and_ordinals::{Barcode, BarcodeInverse, to_ordered_float};
use num::rational::Ratio;


type RingEltRational = Ratio<i64>;
type RingOpRational = solar::rings::ring_native::NativeDivisionRing<Ratio<i64>>;


pub fn interval_d1_v2_barcode_trivial< 'a > () {

    //  DEFINE THE CELL DIMENSIONS
    //  --------------------------

    let cell_dims           =   vec![ 0, 0, 1];      

    //  DEFINE THE RING OPERATOR + BOUNDARY MATRIX
    //  ------------------------------------------

    // boundary matrix of the interval
    let boundary            =   vec![ 
                                    vec![           ],
                                    vec![           ],
                                    vec![  (0, Ratio::new(1,1)),  (1, Ratio::new(-1,1))    ]
                                ];

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
            

    //  DEFINE THE ROOT NODE + RESULTS VECTOR
    //  -------------------------------------

    let root =  Node::make_root(
                        boundary.clone(),
                    &   barcode,
                    &   barcode_inverse,
                    &   cell_dims,  
                        RingOpRational::new()                                  
                );

    let mut results         =   Vec::new();                                

    //  GATHER RESULTS
    //  --------------

    // println!("{:?}", &root );

    explore( & root, &mut results );

    println!("\nRESULTS\n");
    
    for result in results {
        println!("{:?}", result)
    }    


}                                
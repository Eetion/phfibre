

use phfibre::utilities::*;
use phfibre::intervals_and_ordinals::{Barcode, BarcodeInverse, Polytope, to_ordered_float};
use phfibre::phfibre::{explore, Node};
use num::rational::Ratio;


type Cell = usize;
type Coeff = Ratio< i64 >;
type Fil = usize;
type Ring = solar::rings::ring_native::NativeDivisionRing::< Ratio<i64> >;



fn main() {
    println!("Hello, world!");


    //  DEFINE THE CELL DIMENSIONS
    //  --------------------------

    let cell_dims           =   vec![ 0, 0, 1];      

    //  DEFINE THE RING OPERATOR + BOUNDARY MATRIX
    //  ------------------------------------------

    // ring operator
    let ring                =   Ring::new();

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

    let root                =   Node::make_root(
                                        boundary.clone(),
                                    &   barcode,
                                    &   barcode_inverse,
                                    &   cell_dims,                                    
                                );

    let mut results         =   Vec::new();                                

    //  GATHER RESULTS
    //  --------------

    println!("{:?}", &root );

    explore( & root, &mut results );

    println!("{:?}", &results)
    
}





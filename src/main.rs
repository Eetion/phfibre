

use phfibre::utilities::*;
use phfibre::intervals_and_ordinals::*;
use phfibre::phfibre::*;
use num::rational::Ratio;
use std::collections::{HashSet, HashMap};
use std::iter::FromIterator;
use std::iter;
use ordered_float::OrderedFloat;
use solar;


type Cell = usize;
type Coeff = Ratio< i64 >;
type Fil = usize;
type FilRaw = OrderedFloat<f64>; // reason for this choice: f64 does not implement the hash trait (cricital for reparametrizing)
type Ring = solar::rings::ring_native::NativeDivisionRing::< Ratio<i64> >;



fn main() {
    println!("Hello, world!");


    //  DEFINE THE BARCODE
    let barcode_inf_brn     =   to_ordered_float( & vec![0.] );
    let barcode_fin_brn     =   Vec::new();
    let barcode_fin_die     =   Vec::new();

    let barcode_inf_dim     =   vec![0];
    let barcode_fin_dim     =   Vec::new();    

    let barcode             =   Barcode::new(
                                    barcode_inf_dim,
                                    barcode_inf_brn,
                                    barcode_fin_dim,
                                    barcode_fin_brn,
                                    barcode_fin_die
                                );

    //  DEFINE THE RING OPERATOR + BOUNDARY MATRIX

    // ring operator
    let ring                =   Ring::new();

    // boundary matrix of the interval
    let boundary            =   vec![ 
                                    vec![           ],
                                    vec![           ],
                                    vec![  (0, Ratio::new(1,1)),  (1, Ratio::new(-1,1))    ]
                                ];

    //  DEFINE THE CELL DIMENSIONS

    let dims                =   vec![ 0, 0, 1];
    

}






use phfibre::phfibre::{Node, explore, verify_that_barcode_is_compatible};
use phfibre::intervals_and_ordinals::{Barcode, BarcodeInverse, to_ordered_float};
use phfibre::polytopes::polytope_faces::{polys_faces, poly_complex_facets_to_whole_complex_ordinal_data};
use phfibre::polytopes::polytope_intersection::{polytope_intersection};
use phfibre::polytopes::polytope_differential::polyhedral_boundary_matrix_binary_coeff;
use phfibre::rank_calculations::chain_cx_rank_nullity;
use solar::utilities::index::{BiMapSequential, histogram};
use solar::cell_complexes::simplices_unweighted::maximal_cliques::{    
    ordered_subsimplices_up_thru_dim_concatenated_vec 
}; 
use solar::cell_complexes::simplices_unweighted::boundary_matrices::{    
    boundary_matrix_from_complex_facets 
};   
use solar::rings::field_prime::GF2;
use std::iter::FromIterator;
use ordered_float::OrderedFloat;
use num::rational::Ratio;


use exhact::matrix::{SmOracle, RingSpec, RingMetadata, MajorDimension};
use exhact::cubical::{Cube, CubicalComplex};
use exhact::chx::{ChainComplex, factor_chain_complex, ChxTransformKind};
use exhact::clique::Simplex;




type RingEltFixed = OrderedFloat<f64>;
type RingOpFixed = solar::rings::ring_native::NativeDivisionRing<RingEltFixed>;



fn main() { 

    //  DEFINE COEFFICIENT FIELD
    //  --------------------------------------------
    let ring                =   RingOpFixed::new();

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

    // let last_must_be_crit   =   false;


    //  DEFINE THE ROOT NODE + RESULTS VECTOR
    //  -------------------------------------

    let root =  Node::make_root(
                        boundary.clone(),
                    &   barcode,
                    &   barcode_inverse,
                    &   cell_dims,  
                        RingOpFixed::new(),
                        // last_must_be_crit,
                );

    let mut poly_complex_facets         =   Vec::new();                                

    //  GATHER RESULTS
    //  --------------

    // println!("{:?}", &root );

    explore( & root, &mut poly_complex_facets );

    //  REPORT
    //  ------

    println!("number of facets (total): {:?}", poly_complex_facets.len() ); 
    println!("number of facets (binned by dimension): {:?}", histogram( poly_complex_facets.iter().map(|x| x.dim_cellagnostic().unwrap() ) ) );  

    let poly_complex_ordinal_data   =   poly_complex_facets_to_whole_complex_ordinal_data( & poly_complex_facets );
    let poly_complex_dims           =   Vec::from_iter( poly_complex_ordinal_data.ord_to_val.iter().map(|x| x.dim_cellagnostic().unwrap() ) );    
    let poly_complex_dim_top        =   poly_complex_dims.iter().max().unwrap();    
    println!("number of polytopes (total): {:?}", poly_complex_ordinal_data.ord_to_val.len() ); 
    println!("number of polytopes (binned by dimension): {:?}", histogram( poly_complex_dims.iter().cloned() ) );      




    // let mut dim_to_polycount   =   Vec::from_iter( std::iter::repeat(0).take(poly_complex_dim_top + 1) );
    // for dim in 0 .. poly_complex_dim_top + 1 {
    //     dim_to_polycount[ dim ]    =   polys_faces( &poly_complex_facets, dim ).len();
    // }
    // println!("number of polytopes total: {:?}", dim_to_polycount.iter().sum::<usize>() );
    // println!("number of polytopes by dimension (total): {:?}", dim_to_polycount );
    
    // let mut dim_to_facetcount   =   Vec::from_iter( std::iter::repeat(0).take(poly_complex_dim_top + 1) );
    // for facet in poly_complex_facets.iter() {
    //     dim_to_facetcount[ facet.dim_cellagnostic().unwrap() ] += 1;
    // }
    // println!("number of facets total: {:?}", dim_to_facetcount.iter().sum::<usize>() );    
    // println!("number of facets by dimension (total): {:?}", dim_to_facetcount );   

    // // number of facets: 160
    // // number of polytopes total: 2731
    // // number of polytopes by dimension (total): [32, 241, 702, 964, 632, 160]
    // // number of facets total: 160
    // // number of facets by dimension (total): [0, 0, 0, 0, 0, 160] 

    let num_facets              =   poly_complex_facets.len();
    let mut dismat              =   Vec::from_iter(
                                        std::iter::repeat(
                                            Vec::from_iter(
                                                std::iter::repeat( OrderedFloat(1.) )
                                                .take( num_facets)                                                
                                            )
                                        )
                                        .take(num_facets)
                                    );
    for facet_id in 0 .. num_facets{ dismat[facet_id][facet_id] = OrderedFloat(0.) }

    let mut intersection_dim_bins  =   Vec::from_iter( std::iter::repeat(0).take( poly_complex_dim_top + 1) );

    for count_a in 0 .. num_facets{
        for count_b in count_a + 1 .. num_facets {
            if let Some( intersection_poly ) = polytope_intersection( &poly_complex_facets[count_a], &poly_complex_facets[count_b] ){
                dismat[ count_a ][ count_b ] = OrderedFloat( 0.1 );
                dismat[ count_b ][ count_a ] = OrderedFloat( 0.1 );  
                let dim                 =   intersection_poly.dim_cellagnostic().unwrap();                
                intersection_dim_bins[ dim ] += 1              
            }
        }
    }

    println!("number of pairs of intersecting facets, binned by the dimension of the intersection polytope: {:?}", &intersection_dim_bins );

    // ================================================================================== 
    //      THIS CODE IS ON PROBATION -- I DON'T THINK THE CLIQUE COMPLEX IS A NERVE COMPLEX
    //
    // ----------------------------------------------------------------------------------
    // Define an object to represent the ring Z/3Z
    // ----------------------------------------------------------------------------------
    let ringmetadata = exhact::matrix::RingMetadata{
    	ringspec: RingSpec::Rational,
    	identity_additive: Ratio::from_integer(0),
    	identity_multiplicative: Ratio::from_integer(1),
    };


    // ----------------------------------------------------------------------------------
    // Construct the corresponding filtered clique complex
    // ----------------------------------------------------------------------------------
    let chx = exhact::clique::CliqueComplex {
        // the distance/dissimilarity matrix
        dissimilarity_matrix: dismat,
        // threshold to stop the filtration
        dissimilarity_value_max: OrderedFloat(0.5),
        // sets "safeguards" on dimension; we'll get warnings if we try to
        // get boundary matrices in dimension higher than dim+1
        safe_homology_degrees_to_build_boundaries: vec![0],
        // set the default major dimension (for sparse matrices) to be row
        major_dimension: MajorDimension::Row,
        // indicates we want Z/3Z coefficients
        ringmetadata: ringmetadata,
        // don't worry about this
        simplex_count: Vec::new()
    };    

    // ----------------------------------------------------------------------------------
    // Count cliques in each dimension
    // ----------------------------------------------------------------------------------    

    let mut cliques_bin     =   Vec::new();
    for dim in 0 .. 4 {
        cliques_bin.push(
            chx.keys_unordered_itr( dim ).count()
        )
    }

    println!("number of cliques per dimension: {:?}", &cliques_bin );

// ==================================================================================

    let poly_complex_differential           =   polyhedral_boundary_matrix_binary_coeff( & poly_complex_ordinal_data ); 
    let poly_complex_rank_nullity           =   chain_cx_rank_nullity(
                                                    & poly_complex_differential,
                                                    & poly_complex_dims,
                                                    GF2{}
                                                );                                                   
    let poly_complex_betti_vec              =   poly_complex_rank_nullity.rank_homology_vec();         
    println!("betti numbers: {:?}", &poly_complex_betti_vec)




}  


// number of facets (total): 160
// number of facets (binned by dimension): [0, 0, 0, 0, 0, 160]
// number of polytopes (total): 2731
// number of polytopes (binned by dimension): [32, 241, 702, 964, 632, 160]
// number of pairs of intersecting facets, binned by the dimension of the intersection polytope: [5472, 3872, 2120, 928, 328, 0]
// number of cliques per dimension: [160, 12720, 669920, 26294360]
// betti numbers: [1, 0, 0, 0, 0, 0]
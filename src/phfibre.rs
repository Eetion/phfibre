
use crate::utilities::*;
use crate::intervals_and_ordinals::{Barcode, BarcodeInverse, BarFinite, BarInfinite, LevelSetSizes, pairs_dims_births_to_barcode};
use crate::polytope::object_def::{Polytope};
use crate::rank_calculations::{reduced_mat_to_pivot_index_pairs, chain_cx_rank_nullity, num_degenerate_bars_per_degree};
use solar::reduce::vec_of_vec::{clear_cols};
use solar::utilities::indexing_and_bijection::{SuperVec, SuperIndex, sort_perm, inverse_perm};
use solar::utilities::sequences_and_ordinals::{BiMapSequential, ordinate_unique_vals};
use solar::utilities::statistics::histogram;
use solar::rings::ring::{Semiring, Ring, DivisionRing};
use solar::cell_complexes::simplices_unweighted::facets::{    
    ordered_subsimplices_up_thru_dim_concatenated_vec     }; 
use solar::cell_complexes::simplices_unweighted::boundary_matrices::{    
    boundary_matrix_from_complex_facets     }; 
use num::rational::Ratio;
use std;
use std::collections::{HashSet};
use std::hash::Hash;
use std::iter::FromIterator;
use std::fmt::Debug;
use ordered_float::OrderedFloat;
use auto_impl::auto_impl;
use itertools::Itertools;
use serde::{Serialize, Deserialize};

type Cell = usize;
type Coeff = Ratio< i64 >;
type Fil = usize;
// type FilRaw = OrderedFloat<f64>; // reason for this choice: f64 does not implement the hash trait (cricital for reparametrizing)
// type Ring = solar::rings::ring_native::NativeDivisionRing::< Ratio<i64> >;

// use indicatif::{HumanDuration, ProgressBar, ProgressStyle};






// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  TO DO 








//  ---------------------------------------------------------------------------  
//  TREE NODE
//  ---------------------------------------------------------------------------  


/// Represents an "entry" for a sincle cell; raw data about cells is held in a
/// central repository (a vector of these entries) within each node
#[derive(Clone, Debug)]
pub struct CellEntry{
    birth_ordinal:      usize,
    bounding_cell_id:   usize,
    dim:                usize,
}

/// Represents a node of the search tree.
#[derive(Clone, Debug)]
pub struct Node < 'a, FilRaw, RingOp, RingElt, PreconditionToMakeNewLevSet> 
    where   RingOp:             Clone + Semiring<RingElt> + Ring<RingElt> + DivisionRing<RingElt>,
            FilRaw:             Ord + Clone + Hash + Debug,
            PreconditionToMakeNewLevSet:    ExtraConditionToStartNewLevSet + Clone + Debug,
{
    boundary:               Vec< Vec< ( Cell, RingElt ) > >,  // boundary matrix reprsented as a vector of vectors
    barcode:                &'a Barcode< FilRaw >,          // the target barcode (contains useful ordinal data)
    barcode_inverse:        &'a BarcodeInverse,             // pointer to a central register that maps bar endpoints to bar ids     
    cell_id_to_prereqs:     Option< &'a Vec< Vec< usize > > >,        //  cell_id_to_prereqs[m] = cells that must be in K_{in} before cell m can be added to K_{in}

    ring:                   RingOp,
    boundary_buffer:        Vec< ( Cell, Coeff) >,          // a "holding space" for matrix entries, when these need to be moved around

    bars_degn_quota:        Option< Vec< usize > >,         // number of degenerate cells anticipated in each dimension    

    lev_set_sizes:          LevelSetSizes,                  // Kth value = # cells in Kth level set
    polytope:               Polytope,    
    
    cells_all:              Vec< CellEntry >,               // all cells
    cell_ids_out:           Vec< Vec< usize > >,            // cells not yet assigned a birth,
    
    cell_ids_pos_crit:      HashSet< usize >,               // unmatched positive critical cells, grouped by (birth_time, dimension)
    cell_ids_pos_degn:      HashSet< usize >,               // unmatched positive degenerate cells
                                                            // NB: doesn't need to be hash; we know birth time
    
    bar_ids_dun_fin:        Vec< usize >,                   // nonempty bars with all endpoints accounted for (finite)
    bar_ids_dun_inf:        Vec< usize >,                   // nonempty bars with all endpoints accounted for (infinite)

    bar_ids_now_inf_brn:    Vec< usize >,                   // bars still to match in this batch
    bar_ids_now_fin_brn:    Vec< usize >,                   // bars to be born with this level set
    bar_ids_now_fin_die:    Vec< usize >,                   // bars to be bounded with this level set    

    precondition_to_make_new_lev_set: &'a PreconditionToMakeNewLevSet,

    // last_must_be_crit:      bool,                           // a flag to indicate whether the last level set must be critical

}

impl    < 'a, RingOp, RingElt, FilRaw, PreconditionToMakeNewLevSet > 
        Node 
        < 'a, FilRaw, RingOp, RingElt, PreconditionToMakeNewLevSet > 
    where   RingOp:     Clone + Semiring<RingElt> + Ring<RingElt> + DivisionRing<RingElt>,
            RingElt:    Clone + Debug + PartialOrd,
            FilRaw:     Ord + Clone + Hash + Debug,
            PreconditionToMakeNewLevSet:    ExtraConditionToStartNewLevSet + Clone + Debug,                    
{
    /// Generate a root node from boundary data
    pub fn make_root(   
        boundary:               Vec< Vec< (Cell, RingElt)>>,
        barcode:            &'a Barcode< FilRaw >,
        barcode_inverse:    &'a BarcodeInverse,    
        cell_id_to_prereqs:     Option< &'a Vec< Vec< usize > > >,
        cell_dims:          &   Vec< usize >,   
        ring:                   RingOp,  
        precondition_to_make_new_lev_set: &'a PreconditionToMakeNewLevSet,
        // last_must_be_crit:      bool,        
    ) 
    -> 
    Node< 'a, FilRaw, RingOp, RingElt, PreconditionToMakeNewLevSet >    
    {

        let num_cells               =   cell_dims.len();       
        let boundary_buffer         =   Vec::new();

        // compute degenrate bar quotas
        let ranks                   =   chain_cx_rank_nullity(
                                            & boundary, 
                                            & cell_dims,
                                            ring.clone(),                                            
                                        );
        let bars_degn_quota         =   num_degenerate_bars_per_degree(
                                            & ranks,
                                            & barcode
                                        );
        // println!("initial quota: {:?}", &bars_degn_quota);
        // println!("rank vec: {:?}", ranks.rank_boundaries_vec());
        // println!("barcode.num_bars_fin_per_dim: {:?}", barcode.num_bars_fin_per_dim());

        // level set sizes
        let lev_set_sizes           =   LevelSetSizes{ pointers: vec![0] };

        // polytope
        let polytope                =   Polytope{ 
                                            data_l_to_fmin:   vec![0],   // why initialize this way?  because, when recursing, we often initialize a new empty set; indeed, this is **the way** we start new level sets; this is consistent with that
                                            data_c_to_l:      Vec::from_iter( 
                                                                std::iter::repeat( num_cells )
                                                                .take( num_cells )
                                                            )
                                        };              


        // cell registry
        let mut cells_all           =   Vec::with_capacity( num_cells );
        for dim in cell_dims {
            cells_all.push( CellEntry{ 
                                birth_ordinal: num_cells, 
                                bounding_cell_id: num_cells, 
                                dim: *dim } );
        }

        // cell_ids_out
        let chain_space_dim_per_deg =   ranks.rank_vec_chains();    // precompute sizes of component vectors
        let mut cell_ids_out        =   Vec::new();                 // initialize empty sequence of vectors
        for dim in chain_space_dim_per_deg {                                  // initialize constituent vectors
            cell_ids_out.push( Vec::with_capacity(dim) )
        }
        for cell_id in 0..num_cells {                               // populate vectors
            cell_ids_out[ cell_dims[ cell_id ] ].push( cell_id );
        }



        let cell_ids_pos_crit       =   HashSet::new();
        let cell_ids_pos_degn       =   HashSet::new();        
        let bar_ids_dun_fin         =   Vec::new();
        let bar_ids_dun_inf         =   Vec::new();


        let bc_endpoint_now         =   polytope.last_of_all_filtration_ordinals().unwrap(); // equals 0
        let bar_ids_now_inf_brn     =   barcode_inverse
                                            .inf_brn
                                            [ bc_endpoint_now ]
                                            .clone();

        let bar_ids_now_fin_brn     =   barcode_inverse
                                            .fin_brn
                                            [ bc_endpoint_now ]
                                            .clone();            

        let bar_ids_now_fin_die     =   Vec::new(); // no bars die at the first step in the filtration
        
        Node { 
            boundary:               boundary.clone(),
            barcode:            &   barcode,
            barcode_inverse:    &   barcode_inverse,
            cell_id_to_prereqs:     cell_id_to_prereqs,
            ring:                   ring.clone(),
            boundary_buffer:        boundary_buffer,
            bars_degn_quota:        bars_degn_quota,
            lev_set_sizes:          lev_set_sizes,
            polytope:               polytope,
            cells_all:              cells_all,
            cell_ids_out:           cell_ids_out,
            cell_ids_pos_crit:      cell_ids_pos_crit,
            cell_ids_pos_degn:      cell_ids_pos_degn,                            
            bar_ids_dun_fin:        bar_ids_dun_fin,
            bar_ids_dun_inf:        bar_ids_dun_inf,
            bar_ids_now_inf_brn:    bar_ids_now_inf_brn,
            bar_ids_now_fin_brn:    bar_ids_now_fin_brn,
            bar_ids_now_fin_die:    bar_ids_now_fin_die,
            precondition_to_make_new_lev_set: precondition_to_make_new_lev_set,
            // last_must_be_crit:      last_must_be_crit,
        }     
    }



}

//  ---------------------------------------------------------------------------  
//  EXTRA CONDITIONS TO IMPOSE ON STARTING A NEW LEVEL SET
//  ---------------------------------------------------------------------------  

#[auto_impl(&, Box)]
pub trait ExtraConditionToStartNewLevSet {
    fn extra_condition_to_start_new_lev_set( &self, poly: & Polytope ) -> bool;
}




//  NO RESTRICTION

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ConditionNone { }

impl 
    ExtraConditionToStartNewLevSet 
    for 
    ConditionNone
{
    fn extra_condition_to_start_new_lev_set( &self, poly: & Polytope) -> bool { true }
}


//  LOWER STAR FILTRATION

/// Returns false if every cell in the o-skeleton of a simplex S has been assigned a level set, 
/// but S has not yet been assigned to a level set.
#[derive(Clone, Debug)]
pub struct ConditionLowerStar<'a> { pub simplex_bimap_sequential: &'a BiMapSequential< Vec<usize> > }

impl < 'a >
    ExtraConditionToStartNewLevSet 
    for 
    ConditionLowerStar
    < 'a > 
{
    fn extra_condition_to_start_new_lev_set( &self, poly: & Polytope) -> bool {
        for ( simplex_index, simplex ) in self.simplex_bimap_sequential.ord_to_val.iter().enumerate() {
            if  ( 
                    ! poly.cell_id_to_has_lev_set( simplex_index ).unwrap()   // simplex has not been assinged to a level set
                )
                &&
                simplex                                             // but all its vertices have
                    .iter()
                    .cloned()                                       
                    .combinations( 1 )
                    .all(   
                        |x| 
                        {
                            poly.cell_id_to_has_lev_set( 
                                self.simplex_bimap_sequential.ord( &x ).unwrap()                                
                            )
                            .unwrap() 
                        }
                    ) 
            {
                return false
            }
        }
        true 
    }
}

//  LOWER EDGE FILTRATION


/// Returns false if every cell in the 1-skeleton of a simplex S has been assigned a level set, 
/// but S has not yet been assigned to a level set.
#[derive(Clone, Debug)]
pub struct ConditionLowerEdge<'a> { pub simplex_bimap_sequential: &'a BiMapSequential< Vec<usize> > }

impl < 'a >
    ExtraConditionToStartNewLevSet 
    for 
    ConditionLowerEdge
    < 'a > 
{
    fn extra_condition_to_start_new_lev_set( &self, poly: & Polytope) -> bool {
        for ( simplex_index, simplex ) in self.simplex_bimap_sequential.ord_to_val.iter().enumerate() {        

            if  (
                    ! poly.cell_id_to_has_lev_set( simplex_index ).unwrap()  // simplex has not been assinged to a level set
                )
                &&
                simplex                                         // but all its vertices have
                    .iter() 
                    .cloned()                                         
                    .combinations( 1 )
                    .all(   
                        |x| 
                        {
                            poly.cell_id_to_has_lev_set( 
                                self.simplex_bimap_sequential.ord( &x ).unwrap()                                
                            )
                            .unwrap()
                        }
                    ) 
                &&
                simplex                                         // but all its edges have, also
                    .iter()
                    .cloned()
                    .combinations( 2 )
                    .all(   
                        |x| 
                        {
                            poly.cell_id_to_has_lev_set( 
                                self.simplex_bimap_sequential.ord( &x ).unwrap()                                
                            )
                            .unwrap()
                        }
                    ) 
            {
                return false
            }
        }
        true  // return true if no violations have been found
    }
}





//  ---------------------------------------------------------------------------  
//  SEED DATA
//  ---------------------------------------------------------------------------  

// pub struct SeedOfNode< FilRaw, RingOp, RingElt > 
//     where   RingOp:     Clone + Semiring<RingElt> + Ring<RingElt> + DivisionRing<RingElt>,
//             RingElt:    Clone + Debug + PartialOrd,
//             FilRaw:     Ord + Clone + Hash + Debug  
// {
//     boundary:           Vec< Vec< (Cell, RingElt) > >,
//     cell_dims:          Vec< usize >
    
// }    


    // pub fn seed_from_base_space_facets(   
    //     complex_facets:         Vec< Vec< usize > >,
    //     barcode:            &'a Barcode< FilRaw >,
    //     ring:                   RingOp,   
    // )      
    // -> 
    // Node< 'a, FilRaw, RingOp, RingElt >      
    
    // {
    //     let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 1);    
    //     let cell_dims: Vec<_>   =   simplex_sequence.iter().map(|x| x.len()-1 ).collect();
    
    //     let bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence );
    //     let boundary            =   boundary_matrix_from_complex_facets(bimap_sequential, ring.clone());

    //     let barcode_inverse     =   BarcodeInverse::from_barcode( & barcode );        

    //     Node::make_root(
    //             boundary.clone(),
    //         &   barcode,
    //         &   barcode_inverse,
    //         &   cell_dims,  
    //             ring.clone(),
    //             // last_must_be_crit,
    //     )        
    // }


//  ---------------------------------------------------------------------------  
//  EXPLORE TREE
//  ---------------------------------------------------------------------------  


pub fn  explore< FilRaw, RingOp, RingElt, PreconditionToMakeNewLevSet >( 
        node:       &       Node < FilRaw, RingOp, RingElt, PreconditionToMakeNewLevSet >, 
        results:    &mut    Vec < Polytope > 
    )
    where   RingOp:                     Clone + Semiring<RingElt> + Ring<RingElt> + DivisionRing<RingElt>,
            RingElt:                    Clone + Debug + PartialOrd,
            FilRaw:                     Clone + Debug + Ord + Hash,
            PreconditionToMakeNewLevSet:    ExtraConditionToStartNewLevSet + Clone + Debug,
{ 

    //  TERMINATE IMMEDIATELY IF BARCODE IS INCOMPATIBLE

    if node.bars_degn_quota.is_none() { 
        println!("\nFUNCTIN <explore> CAN ADD NOTHING TO RESULTS BECAUSE BARCDE IS INCOMPATIBLE WITH THIS BOUNDARY MATRIX\n");
        return 
    }

    //  RECORD SOME USEFUL VARIABLES
    let node_lev_set_last_is_crit   =   node.polytope.lev_set_last_is_critical().unwrap();
    let node_fmin_last              =   node.polytope.last_of_all_filtration_ordinals().unwrap();    
    let node_fmin_last_intersects_a_finite_bar    =   
                                            // evaluatest to true iff every bar that passes through this level set lasts forever (in which case, even if the class is critical and could be built up to a larger critical class in the future, the inclusion of the smaller into the larger will induce monomophrism, hence isomorphism)
                                            node.barcode.fin.iter().any(|x| x.birth <= node_fmin_last && node_fmin_last < x.death );    


    //  SHORT CIRCUIT IF WE HAVE ALREADY BEEN DOWN THIS ROAD
        // BENCHMARKS ON THE SQUARE SUBDIVIDED INTO TWO TRIANGLES (NO FINITE BARS)
        // with    this short circuit: 12.926543132
        // without this short circuit: 44.080082316s  
        
    if node.lev_set_sizes.num_cells_total() < node.cells_all.len() {
    // original short circuit:
        if node.lev_set_sizes.size_last()          ==  Some( 0 ) {
            for poly in results.iter() {
                if poly.contains_extension( & node.polytope ) { return }
            }
        }
        // experimental new short circuit:
        else if node.lev_set_sizes.size_last() > Some( 1 )  
        {
            for poly in results.iter() {
                if  poly.certifies_prior_exploration( & node.polytope, node_fmin_last_intersects_a_finite_bar ) { return }
            } 
        }   

        // for poly in results.iter() {
        //     if  poly.certifies_prior_exploration( & node.polytope ) { return }
        // }        
    }


    //  PUSH LEAF NODE TO RESULTS

    if  node.lev_set_sizes.size_last()          ==  Some( 0 ) &&
        node.lev_set_sizes.num_cells_total()    ==  node.cells_all.len() &&
        node.bar_ids_dun_fin.len()              ==  node.barcode.num_bars_fin() &&
        node.bar_ids_dun_inf.len()              ==  node.barcode.num_bars_inf() //&&
        // (            
        //     node.polytope.lev_set_last_is_critical() == Some( true )
        //     ||
        //     ! node.last_must_be_crit
        // )           
        
        {

        let mut poly                            =   node.polytope.clone();

        // the last level set is non-critical and empty; switch it to critical.  this empty, newly-made-critical
        // level set will represent the last, dimension-0 simplex factor of the polytope
        *poly.data_l_to_fmin.last_mut().unwrap() += 1;

        // we may have already constructed an equivalent filtration (just with a different order
        // on cells within a level set).  if we haven't then push to results.        
        let mut push                            =   true;
        for result_count in ( 0 .. results.len() ).rev()  {
            // if the current result is a face of a prior one, then don't bother adding to results
            if results[ result_count ].contains( & poly ) { 
                push                            =   false; 
                break 
            }
            // if conversely the current polytope contains a prexisting polytope, then remove the smaller
            else if poly.contains( &results[ result_count ] ) { 
                results.remove( result_count );
            }
        }
        if push {
            // println!("num results: {:?}",  histogram( results.iter().map(|x| x.dim_cellagnostic().unwrap() ) )   );
            results.push( poly ); 
        }

    }

    else {    

        //  IF FEASIBLE, INITIALIZE NEW LEVEL SET 
        //  ---------------------------------------------------------------------
     
        if node.lev_set_sizes.size_last() != Some( 0 )      &&  // current level set is nonempty 
                node.precondition_to_make_new_lev_set.extra_condition_to_start_new_lev_set( & node.polytope ) && 
                node.cell_ids_pos_degn.is_empty()           &&  // C-rule (degenerate) there are no unmatched "degenerate" positive cells
                (                                               // C-rule (critical) 
                    (
                    !node.polytope                                  // the level set contains no "critical" cell
                            .lev_set_last_is_critical()
                            .unwrap() 
                    )                          
                    ||                                              // OR
                    (                                               // I_cur is empty
                        node.bar_ids_now_inf_brn.is_empty() &&     
                        node.bar_ids_now_fin_brn.is_empty() &&    
                        node.bar_ids_now_fin_die.is_empty()    
                    )
                )                                       

                //  THIS CODE CAUSED PROBLEMS AND I THINK IT'S OK TO DELETE
                // &&
                // node.barcode                                    // the current level set does not map to 1.0
                //     .ordinal
                //     .val( node.bc_endpoint_now )
                //     != Some( OrderedFloat( 1.0 ) )

            {

            // CREATE CHILD NODE    
        
            let mut child                   =   node.clone();

            // APPEND NEW COUNTER FOR LEVEL SET SIZE

            child.lev_set_sizes.postpend_empty_set();       // post-pend new entry (= 0) to vector of level set sizes
            child.polytope.push_new_lev_set();              // post-pend new level set to polytope object (note this object doesn't have a well formed notion of an empty level set)


            // IF PARENT NODE IS CRITICAL, THEN (1) UPDATE BARCODE ENDPOINT, AND (2) UPDATE SET OF BARS TO ACCOUNT FOR

            // NB: we don't check that bc_endpoint_now < maximum_endpoint; nothing bad happens
            // if this is indeed the case; one simply obtains three empty sets (moreover, this doesn't
            // produce an infinite loop, because we only reach this point after creating a new
            // nonempty level set.

            if  node_lev_set_last_is_crit  { 
                
                // update barcode endpoint

                // child.polytope.ensure_last_lev_set_critical();
                
                let bc_endpoint_now         =   node_fmin_last + 1;


                // update set of bars to account for
                child.bar_ids_now_inf_brn   =   child.barcode_inverse
                                                    .inf_brn
                                                    .sindex( bc_endpoint_now, vec![] )
                                                    .clone();

                child.bar_ids_now_fin_brn   =   child.barcode_inverse
                                                    .fin_brn
                                                    .sindex( bc_endpoint_now, vec![] )
                                                    .clone();            
                
                child.bar_ids_now_fin_die   =   child.barcode_inverse
                                                    .fin_die
                                                    .sindex( bc_endpoint_now, vec![] )
                                                    .clone();                        

            }

            // RUN EXPLORE ON CHILD

            // println!("BOUNDARY:");
            // for col in child.boundary.iter() {
            //     println!("{:?}", &col)
            // }
            // println!("POLYTOPE: {:?}", & child.polytope.data_c_to_l);                    
            // println!("CELL IDS OUT: {:?}", & child.cell_ids_out );                       
            // println!("BIRTH ORDINALS: {:?}", Vec::from_iter( child.cells_all.iter().cloned().map(|x| x.birth_ordinal) ));        
            
            // println!("CALL 1: INITIALIZED NEW LEV SET");
            explore( & child, results );

            // I DON'T THINK WE CAN USE THIS SHORT CIRCUIT AT ALL ANY MORE
            // we can only use this short-circuit for degenerate sets; for critical sets it's possible we may need to keep going
            if  (   // either the level set is critical (so any future additions could be simulated by merging this class upward)
                    ! 
                    node_lev_set_last_is_crit
                )
                ||
                (   // or every bar that passes through this level set lasts forever (thus, even if the class is critical and could be built up to a larger critical class in the future, the inclusion of the smaller into the larger will induce monomophrism, hence isomorphism)
                    !
                    node_fmin_last_intersects_a_finite_bar
                )                


            { return } 
            
        }

        //  ENLARGE CURRENT LEVEL SET 
        //  ---------------------------------------------------------------------


        //  ADD DEATH CELL FOR A *DEGENERATE* BAR 
        //  -------------------------------------

        // loop over all dimensions
        for dim in 1..node.cell_ids_out.len() {

            // (error - deprecated) short circuit if the current level set is critical
            // if node.lev_set_sizes.size_last().unwrap() > 0  && node.polytope.num_lev_sets() > 1 && ( node.polytope.lev_set_last_is_critical().unwrap() ) { continue };                  

            // loop over chains of given dimension with nonzero boundary, skipping empty ones but maintaining a global count of chains checked
            // (this is why filter has to come last)
            for (neg_id_count, neg_id) in   node.cell_ids_out
                                                [ dim ]
                                                .iter()
                                                .cloned()
                                                .enumerate()
                                                .filter(|(_x_count, x)|  ! node.boundary[ *x ].is_empty() )                                                
            {


                let low_entry   =   node.boundary
                                        [ neg_id ]
                                        .iter()
                                        .cloned()
                                        .max_by(
                                            |(xind, xcoeff),  (yind, ycoeff)|
                                            node.cells_all[ *xind ].birth_ordinal.cmp(
                                                & node.cells_all[ *yind ].birth_ordinal
                                            )
                                        )
                                        .unwrap();  
                
                let low_cell_id     =   low_entry.0.clone();
                let low_birth_ord   =   node.cells_all[ low_cell_id ].birth_ordinal.clone();            


                //   IF  every cell in the support of this column has been assigned a birth ordinal 
                //  AND  low_id corresponds to a positive degenerate cell
                // THEN  add a negative degenerate cell
                if  low_birth_ord < node.cells_all.len() &&  // this first condition helps to short circuit / avoid expensive look-ups in the hash set
                    node.cell_ids_pos_degn.contains( & low_cell_id ) 
                {

                    // update: 
                    //  lev_set_sizes 
                    //  cells_all
                    //  cell_ids_out 
                    //  bar_ids_dun 
                    //  bar_ids_now
                    
                    // create child node
                    let mut child   =   node.clone();                                           // this must come before swap_remove

                    // move cells
                    let _           =   child.cell_ids_out[ dim ].swap_remove( neg_id_count );         // remove negative cell from "out list"
                    child.cell_ids_pos_degn.remove( & low_cell_id );                     // remove positive cell from set of unmatched positive cells

                    // record the birth ordinal of the added cell in the central registry
                    child
                        .cells_all
                        [ neg_id ]
                        .birth_ordinal
                                    =   child.lev_set_sizes.num_cells_total();

                    // record the birth **level set ordinal** of the added cell in the polytope
                    child
                        .polytope
                        .data_c_to_l    
                        [ neg_id ]
                                    =   child.lev_set_sizes.last_level_set_ordinal().unwrap();                                       
                   
                    // update number of cells born
                    child.lev_set_sizes.grow_last_set(); 

                    // record the pairing in the central registry
                    child.cells_all
                        [ low_cell_id ]
                        .bounding_cell_id
                                    =   neg_id.clone();

                    // UPDATE THE BOUNDARY MATRIX

                    // define the ingredients for clearing
                    let clearor     =   child.boundary[ neg_id ].clone();
                    let pivot_entry =   low_entry.clone();
                    let ring        =   child.ring.clone();

                    // deallocate the negative column (it's no longer needed)
                    child.boundary[ neg_id ].clear();
                    child.boundary[ neg_id ].shrink_to_fit();    

                    // clear bottom entries from columns                    
                    clear_cols(
                        &       clearor,
                        & mut   child.boundary,
                                child.cell_ids_out
                                    [ dim ]
                                    .clone(),
                        &       pivot_entry,
                                ring
                    );

                    // RUN EXPLORE ON CHILD

                    // println!("BOUNDARY:");
                    // for col in child.boundary.iter() {
                    //     println!("{:?}", &col)
                    // }
                    // println!("POLYTOPE: {:?}", & child.polytope.data_c_to_l);                    
                    // println!("CELL IDS OUT: {:?}", & child.cell_ids_out );                       
                    // println!("BIRTH ORDINALS: {:?}", Vec::from_iter( child.cells_all.iter().cloned().map(|x| x.birth_ordinal) ));
                    // println!("node.cell_ids_pos_degn: {:?}", &node.cell_ids_pos_degn);     
                    // println!("child.cell_ids_pos_degn: {:?}", &child.cell_ids_pos_degn);      
                                        
                    // println!("CALL 6: ADDED DEATH CELL FOR NON-CRICIAL BAR OF DIM {:?}", dim.clone()-1 );                                                                      
                    explore( & child, results );
                }
            }
        }     
        
        //  ADD DEATH CELL FOR A *FINITE* BAR 
        //  ------------------------------------

        // NB:  PROBABLY AN OPORTUNITY TO IMPROVE PERFORMANCE HERE: rather than loop bar by bar, just create a 
        //      HASHMAP : (BIRTH, DEATH) -> MULTIPLICITY; loop once over all columns and for each column check 
        //      WHETHER THE BAR IT WOULD CREATE WOULD LIE IN THE DESIRED SET

        for (bar_id_count, bar_id) in node.bar_ids_now_fin_die.iter().cloned().enumerate() {

            // (error - deprecated) short circuit if the current level set is non-empty and non-critical
            // if node.lev_set_sizes.size_last().unwrap() > 0  && node.polytope.num_lev_sets() > 1 && ( ! node.polytope.lev_set_last_is_critical().unwrap() ) { continue };            

            // clone info about bar to be added
            let bar_new         =   node.barcode
                                        .bar_fin( bar_id );

            // loop over all dimension (dim+1) chains, skipping empty ones but maintaining a global count of chains checked
            // (this is why filter has to come last)
            for (neg_id_count, neg_id) in   node.cell_ids_out
                                                [ bar_new.dim() + 1 ]
                                                .iter()
                                                .cloned()
                                                .enumerate()
                                                .filter(|(_x_count, x)|  ! node.boundary[ *x ].is_empty() )
            {


                // find the lowest entry of the vector, according to birth order
                let low_entry   =   node.boundary
                                        [ neg_id ]
                                        .iter()
                                        .cloned()
                                        .max_by(
                                            |(xind, xcoeff),  (yind, ycoeff)|
                                            node.cells_all[ *xind ].birth_ordinal.cmp(
                                                & node.cells_all[ *yind ].birth_ordinal
                                            )
                                        )
                                        .unwrap();  
                
                let low_cell_id     =   low_entry.0.clone();
                let low_birth_ord   =   node.cells_all[ low_cell_id ].birth_ordinal.clone();
                                    
                //   IF  every cell in the support of this column has been assigned a birth ordinal 
                //  AND  low_id corresponds to a positive degenerate cell
                // THEN  add a negative degenerate cell
                if  low_birth_ord < node.cells_all.len()  // this first condition helps to short circuit / avoid expensive look-ups in the hash set
                    &&  
                    node.polytope.cell_id_to_fmin( low_cell_id.clone() ).unwrap() == bar_new.birth()
                    &&
                    node.cell_ids_pos_crit.contains( & low_cell_id )

                {

                    // update: 
                    //  lev_set_sizes 
                    //  cells_all
                    //  cell_ids_out 
                    //  bar_ids_dun 
                    //  bar_ids_now
                    
                    // create child node
                    let mut child   =   node.clone();                                           // this must come before swap_remove

                    // mark child's level set as critical (posisbly redundantly)
                    child.polytope.ensure_last_lev_set_critical();

                    // move cells
                    let _           =   child.cell_ids_out[ bar_new.dim() + 1 ].swap_remove(neg_id_count);       // remove negative cell from "out list"
                    child.cell_ids_pos_crit.remove( & low_cell_id );                                   // remove positive cell from set of unmatched positive cells

                    // move bar
                    let _           =   child.bar_ids_now_fin_die.swap_remove(bar_id_count);    // remove bar from "unadded list"
                    child.bar_ids_dun_fin.push( bar_id.clone() );                               // add bar to "done" list (consumes variable)

                   
                    // record the birth ordinal of the added cell in the central registry
                    child
                        .cells_all
                        [ neg_id ]
                        .birth_ordinal
                                    =   child.lev_set_sizes.num_cells_total();

                    // record the birth *level set ordinal** of the added cell in the polytope
                    child
                        .polytope
                        .data_c_to_l    
                        [ neg_id ]
                                        =   child.lev_set_sizes.last_level_set_ordinal().unwrap();                                      
                   
                    // update number of cells born
                    child.lev_set_sizes.grow_last_set(); 

                    // record the pairing in the central registry
                    child.cells_all
                        [ low_cell_id ]
                        .bounding_cell_id
                                    =   neg_id.clone();


                    // UPDATE THE BOUNDARY MATRIX

                    // define the ingredients for clearing
                    let clearor     =   child.boundary[ neg_id ].clone();
                    let pivot_entry =   low_entry.clone();
                    let ring        =   child.ring.clone();

                    // deallocate the negative column (it's no longer needed)
                    child.boundary[ neg_id ].clear();
                    child.boundary[ neg_id ].shrink_to_fit();    

                    // clear bottom entries from columns                    
                    clear_cols(
                        &       clearor,
                        & mut   child.boundary,
                                child.cell_ids_out
                                    [ bar_new.dim() +1 ]  // only columns indexed by cells of degree dim+1 could be changed
                                    .clone(),
                        &       pivot_entry,
                                ring
                    );
                    
                    


                    // RUN EXPLORE ON CHILD

                    // println!("BOUNDARY:");
                    // for col in child.boundary.iter() {
                    //     println!("{:?}", &col)
                    // }
                    // println!("POLYTOPE: {:?}", & child.polytope.data_c_to_l);                    
                    // println!("CELL IDS OUT: {:?}", & child.cell_ids_out );                       
                    // println!("BIRTH ORDINALS: {:?}", Vec::from_iter( child.cells_all.iter().cloned().map(|x| x.birth_ordinal) ));                       
                    
                    // println!("CALL 4: ADDED DEATH CELL FOR FIN BAR OF DIM {:?}", bar_new.dim());                                                    
                    explore( & child, results );
                }
            }
        }        

        //  ADD BIRTH CELL FOR AN *INFINITE* BAR 
        //  ------------------------------------
        
        for (bar_id_count, bar_id) in node.bar_ids_now_inf_brn.iter().cloned().enumerate() { 

            // (error - deprecated) short circuit if the current level set is non-empty and non-critical
            // if node.lev_set_sizes.size_last().unwrap() > 0  && node.polytope.num_lev_sets() > 1 && ( ! node.polytope.lev_set_last_is_critical().unwrap() ) { continue };


            // clone info about bar to be added
            let bar_new     =   node
                                    .barcode
                                    .bar_inf( bar_id );  

            // loop over all cycles of appropriate dimension 
            for (pos_id_out_count, pos_id_out) in node.cell_ids_out
                                                    [ bar_new.dim()  ]
                                                    .iter()
                                                    .cloned()
                                                    .enumerate()
                                                    .filter(    |x| 
                                                                node.boundary[x.1.clone()]
                                                                    .is_empty()
                    )  
            {   

                // update: 
                //  lev_set_sizes 
                //  cells_all
                //  cell_ids_out 
                //  bar_ids_dun 
                //  bar_ids_now

                // create child node
                let mut child   =   node.clone();                                   // this must come before swap_remove

                // mark child's current level set as critical (posisbly redundantly)
                child.polytope.ensure_last_lev_set_critical();

                // move bar and positive cell
                let _    =   child.cell_ids_out[ bar_new.dim() ].swap_remove(pos_id_out_count);      // remove cell from "out list"
                let _    =   child.bar_ids_now_inf_brn.swap_remove(bar_id_count);   // remove bar from "unadded list"
                child.bar_ids_dun_inf.push( bar_id );                               // add bar to "done" list (consumes variable)
               
                // record the birth ordinal of the added cell in the central registry
                child
                    .cells_all
                    [ pos_id_out ]
                    .birth_ordinal
                                        =   child.lev_set_sizes.num_cells_total();

                // record the birth *level set ordinal** of the added cell in the polytope
                child
                    .polytope
                    .data_c_to_l    
                    [ pos_id_out ]
                                        =   child.lev_set_sizes.last_level_set_ordinal().unwrap();                                        
               
                // update number of cells born
                child.lev_set_sizes.grow_last_set();

                // RUN EXPLORE ON CHILD

                // println!("BOUNDARY:");
                // for col in child.boundary.iter().skip(2) {
                //     println!("{:?}", &col);
                // }
                // println!("POLYTOPE: {:?}", & child.polytope.data_c_to_l);                    
                // println!("CELL IDS OUT: {:?}", & child.cell_ids_out );                       
                // println!("BIRTH ORDINALS: {:?}", Vec::from_iter( child.cells_all.iter().cloned().map(|x| x.birth_ordinal) ));   
                // println!("child.bar_ids_now_inf_brn: {:?}", &child.bar_ids_now_inf_brn);   
                // println!("child.cell_ids_pos_degn: {:?}", &child.cell_ids_pos_degn);                   
                
                
                // println!("CALL 2: ADDED BIRTH CELL FOR INF BAR OF DIM {:?}", bar_new.dim());                
                explore( & child, results );
            }
        } 
        
        //  ADD BIRTH CELL FOR A *FINITE* BAR 
        //  ------------------------------------

        for (bar_id_count, bar_id) in node.bar_ids_now_fin_brn.iter().cloned().enumerate() {

            // update: 
            //  lev_set_sizes 
            //  cells_all
            //  cell_ids_out 
            //  bar_ids_dun 
            //  bar_ids_now

            // (error, deprecated) short circuit if the current level set is non-empty and non-critical
            // if node.lev_set_sizes.size_last().unwrap() > 0  && node.polytope.num_lev_sets() > 1 &&  ( ! node.polytope.lev_set_last_is_critical().unwrap() ) { continue };            

            // clone info about bar to be added
            let bar_new     =   node.barcode
                                    .bar_fin( bar_id );

            // loop over all cycles of appropriate dimension 
            for (pos_id_out_count, pos_id_out) in node.cell_ids_out
                                                    [ bar_new.dim() ]
                                                    .iter()
                                                    .cloned()
                                                    .enumerate()
                                                    .filter(    |x| 
                                                                node.boundary[ x.1 ]
                                                                    .is_empty()
                                                    )  
            {   

                // create child node
                let mut child   =   node.clone();                                   // this must come before swap_remove

                // mark child's level set as critical (posisbly redundantly)
                child.polytope.ensure_last_lev_set_critical();
                
                // move bar_id
                let _    =   child.bar_ids_now_fin_brn.swap_remove(bar_id_count);   // remove bar from list of bars that have unaccounted endpoints here

                // move pos_id_out
                let _    =   child.cell_ids_out[ bar_new.dim() ].swap_remove(pos_id_out_count);      // remove cell from "out list"
                
                child.cell_ids_pos_crit.insert(
                    pos_id_out 
                );
               
                // record the birth ordinal of the added cell in the central registry
                child
                    .cells_all
                    [ pos_id_out ]
                    .birth_ordinal
                                =   child.lev_set_sizes.num_cells_total();

                // record the birth *level set ordinal** of the added cell in the polytope
                child
                    .polytope
                    .data_c_to_l    
                    [ pos_id_out ]
                                        =   child.lev_set_sizes.last_level_set_ordinal().unwrap();                                      

                // update number of cells born
                child.lev_set_sizes.grow_last_set(); 
                

                // RUN EXPLORE ON CHILD

                // println!("BOUNDARY:");
                // for col in child.boundary.iter() {
                //     println!("{:?}", &col)
                // }
                // println!("POLYTOPE: {:?}", & child.polytope.data_c_to_l);                    
                // println!("CELL IDS OUT: {:?}", & child.cell_ids_out );                       
                // println!("BIRTH ORDINALS: {:?}", Vec::from_iter( child.cells_all.iter().cloned().map(|x| x.birth_ordinal) ));                
                
                // println!("CALL 3: ADDED BIRTH CELL FOR FIN BAR OF DIM {:?}", bar_new.dim());                
                explore( & child, results );

            }
        } 

        //  ADD BIRTH CELL FOR A *DEGENERATE* BAR 
        //  -------------------------------------

        // if node.polytope.min_vertex_is_compatible_with_ordinal_filt( &vec![2, 0, 1, 2, 2, 3] ) {   
        //     println!("WILL TRY TO ADD SOME DEGENERATE CELLS!");
        //     println!("NODE WAS COMPATIBLE FROM THE START: {:?}", first_pass_was_true)
        // }          

        // loop over all dimensions
        for dim in 0..node.cell_ids_out.len() {

            // (error - deprecated) short circuit if the current level set is critical
            // if node.lev_set_sizes.size_last().unwrap() > 0  && node.polytope.num_lev_sets() > 1 && ( node.polytope.lev_set_last_is_critical().unwrap() ) { continue };            

            // short-circuit if we have already met the quota for degenerate bars in this dimension
            if node.bars_degn_quota.as_ref().map(|x| x.sindex(dim,0) ) == Some(0) { 
                // println!("skipping: quota {:?}, dim {:?}", &node.bars_degn_quota, &dim); 
                continue             
            }

            // loop over all cycles of given dimension 
            for (pos_id_out_count, pos_id_out) in node.cell_ids_out
                                                    [ dim ]
                                                    .iter()
                                                    .cloned()
                                                    .enumerate()
                                                    .filter(    |x| 
                                                                node.boundary[ x.1 ]
                                                                    .is_empty()
                                                    )  
            {   
                
             
                // update: 
                //  lev_set_sizes 
                //  cells_all
                //  cell_ids_out 
                //  bar_ids_dun 
                //  bar_ids_now

                // create child node
                let mut child   =   node.clone();                                   // this must come before swap_remove


                // move pos_id_out
                let _    =   child.cell_ids_out[ dim ].swap_remove(pos_id_out_count);      // remove cell from "out list"
                child.cell_ids_pos_degn.insert( pos_id_out );                         // add cell to list of positive degenerate cells
               
                // record the birth ordinal of the added cell in the central registry
                child
                    .cells_all
                    [ pos_id_out ]
                    .birth_ordinal
                                =   child.lev_set_sizes.num_cells_total();

                // record the birth **level set ordinal** of the added cell in the polytope
                child
                    .polytope
                    .data_c_to_l    
                    [ pos_id_out ]
                                        =   child.lev_set_sizes.last_level_set_ordinal().unwrap();                                  
               
                // update counters
                child.lev_set_sizes.grow_last_set();    // number of cells born

                // drop the quota for degenerate bars by 1
                match child.bars_degn_quota.as_mut() {
                    None => { return },
                    Some(quota_vec) => { quota_vec[dim] -= 1; }         
                }

                // RUN EXPLORE ON CHILD

                // if node.polytope.min_vertex_is_compatible_with_ordinal_filt( &vec![2, 0, 1, 2, 2, 3] ) {   
                //     println!("PARENT polytope:  {:?} {:?}", & node.polytope.vec_mapping_cell_id_to_min_filt_ordinal(), & node.polytope.data_c_to_l );                     
                //     println!("CHILD  polytope:  {:?} {:?}", & child.polytope.vec_mapping_cell_id_to_min_filt_ordinal(),  & child.polytope.data_c_to_l );
                //     println!("child full polytope:  {:?}", & child.polytope );                                                           
                //     println!("{:?}", child.polytope.min_vertex_is_compatible_with_ordinal_filt( &vec![2, 0, 1, 2, 2, 3] ));
                // }                                
                
                // println!("BOUNDARY:");
                // for col in child.boundary.iter() {
                //     println!("{:?}", &col)
                // }
                // println!("POLYTOPE: {:?}", & child.polytope.data_c_to_l);                    
                // println!("CELL IDS OUT: {:?}", & child.cell_ids_out );                       
                // println!("BIRTH ORDINALS: {:?}", Vec::from_iter( child.cells_all.iter().cloned().map(|x| x.birth_ordinal) )); 

                // println!("CALL 5: ADDED BIRTH CELL FOR NON-CRICIAL BAR OF DIM {:?}", dim.clone() );                                    
                explore( & child, results );
            }
        }
    }

    // spinner.finish_and_clear();
}

//  ---------------------------------------------------------------------------  
//  CHECK BARCODE
//  ---------------------------------------------------------------------------  


pub fn  is_barcode_compatible< FilRaw, RingOp, RingElt, PreconditionToMakeNewLevSet >( 
            root_node:  &   Node< FilRaw, RingOp, RingElt, PreconditionToMakeNewLevSet >, 
            result:     &   Polytope
        )
        ->
        bool
    where   RingOp:                     Clone + Semiring<RingElt> + Ring<RingElt> + DivisionRing<RingElt>,
            RingElt:                    Clone + Debug + PartialOrd + Ord,
            FilRaw:                     Clone + Debug + Ord + Hash,
            PreconditionToMakeNewLevSet:    ExtraConditionToStartNewLevSet + Clone + Debug,
{

    let num_cells               =   root_node.cells_all.len();

    // compute the sorting permutation
    let new_to_old              =   sort_perm( &result.data_c_to_l  ); 
    let old_to_new              =   inverse_perm( &new_to_old );

    // make copy of boundary matrix with permuted columns -- note we must convert FROM new indices TO old indices
    let mut reducindus_permuted =   Vec::from_iter( new_to_old.iter().map(|&x| root_node.boundary[x].clone()) );

    // permute row indices and sort rows
    for col_index in 0 .. num_cells {
        let mut vec             =   reducindus_permuted[ col_index ].clone();
        for row_pointer in 0 .. vec.len() {
            let row_index           =   vec[ row_pointer ].0.clone();
            vec[ row_pointer ].0    =   old_to_new[ row_index ];       // note we must convert FROM old indices TO new indices     
        }
        vec.sort();
        reducindus_permuted[ col_index ]    =   vec;
    }

    // reduce
    solar::reduce::vec_of_vec::right_reduce(
        &mut reducindus_permuted,
        root_node.ring.clone(),
    );    

    // extract the pairs + compute the barcode
    let pairs                   =   reduced_mat_to_pivot_index_pairs( &reducindus_permuted );

    // obtain permuted dimension vector (WE MUST ACCOUNT FOR THE PERMUTATION)
    let cell_dims               =   Vec::from_iter(
                                        new_to_old.iter().map(|&x| root_node.cells_all[ x ].dim )
                                    );    

    // calculate cell births
    // obtain permuted dimension vector (WE MUST ACCOUNT FOR THE PERMUTATION)
    let births                  =   Vec::from_iter(
                                        new_to_old
                                            .iter()
                                            .map(   |&x| 
                                                    result.cell_id_to_fmin( x ).unwrap()  // recall that this returns ordinals of barcode endpoints
                                            )
                                            .map(   |x|
                                                    root_node.barcode.ordinal.val( x ).unwrap() // in light of preceding comment, we must convert back to raw filtration values
                                            )
                                    );        

    // compute the proposed barcode (with sorted entries)
    let mut barcode_proposed    =   pairs_dims_births_to_barcode(
                                        & pairs,
                                        & cell_dims,
                                        & births
                                    );

    barcode_proposed.inf.sort();
    barcode_proposed.fin.sort();

    // compute the objective barcode (with sorted entries)

    let mut barcode_true        =   root_node.barcode.clone();
    barcode_true.inf.sort();
    barcode_true.fin.sort();


    // print out any discrepencies

    if &barcode_proposed.inf != &barcode_true.inf {
        println!("barcode_proposed.inf: {:?}", &barcode_proposed.inf);
        println!("barcode_true.inf: {:?}",     &barcode_true.inf);        
    }

    if &barcode_proposed.fin != &barcode_true.fin {
        println!("barcode_proposed.fin: {:?}", &barcode_proposed.fin);
        println!("barcode_true.fin: {:?}",     &barcode_true.fin);        
    }  
    
    if &barcode_proposed.ordinal.ord_to_val != &barcode_true.ordinal.ord_to_val {
        println!("barcode_proposed.ordinal.ord_to_val: {:?}", &barcode_proposed.ordinal.ord_to_val);
        println!("barcode_true.ordinal.ord_to_val: {:?}",     &barcode_true.ordinal.ord_to_val);        
    }       

    // return true iff there are no discrepencies

    return  &barcode_proposed.inf   ==  &barcode_true.inf
            &&
            &barcode_proposed.fin   ==  &barcode_true.fin
            &&
            &barcode_proposed.ordinal.ord_to_val    ==  &barcode_true.ordinal.ord_to_val

    // assert_eq!( &barcode_proposed.inf, &barcode_true.inf );
    // assert_eq!( &barcode_proposed.fin, &barcode_true.fin );    
    // assert_eq!( &barcode_proposed.ordinal.ord_to_val, &barcode_true.ordinal.ord_to_val );        

}   


pub fn  verify_that_barcode_is_compatible< FilRaw, RingOp, RingElt, PreconditionToMakeNewLevSet >( 
            root_node:  &   Node< FilRaw, RingOp, RingElt, PreconditionToMakeNewLevSet >, 
            result:     &   Polytope
        )
    where   RingOp:     Clone + Semiring<RingElt> + Ring<RingElt> + DivisionRing<RingElt>,
            RingElt:    Clone + Debug + PartialOrd + Ord,
            FilRaw:     Clone + Debug + Ord + Hash,
            PreconditionToMakeNewLevSet:    ExtraConditionToStartNewLevSet + Clone + Debug,

{
    assert!( is_barcode_compatible( root_node, result ) );
}






#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;


    #[test]
    fn test_verify_that_barcode_is_compatible() {

        type RingEltFixed       =   OrderedFloat<f64>;
        type RingOpFixed        =   solar::rings::ring_native::NativeDivisionRing<RingEltFixed>;
        let ring                =   RingOpFixed::new();

        let complex_facets      =   vec![ vec![0, 1], vec![1, 2], vec![0, 2] ];
        
        let simplex_sequence    =   ordered_subsimplices_up_thru_dim_concatenated_vec( &complex_facets, 1);    
        let cell_dims: Vec<_>   =   simplex_sequence.iter().map(|x| x.len()-1 ).collect();
    
        let bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence );
        let boundary            =   boundary_matrix_from_complex_facets( &bimap_sequential, ring.clone());
        
        let barcode             =   Barcode{
                                        inf: vec![ BarInfinite{dim:0,birth:0}, BarInfinite{dim:1,birth:3} ],
                                        fin: vec![ BarFinite{dim:0,birth:1,death:2} ],
                                        ordinal: ordinate_unique_vals( & vec![0,1,2,3] )
                                    };
        
        let barcode_inverse     =   BarcodeInverse::from_barcode( & barcode );    
        
        let cell_id_to_prereqs  =   None; // no prerequisites for any cell

        let ok_to_start_new_lev_set     =   ConditionNone{};

        let root_node           =  Node::make_root(
                                            boundary.clone(),
                                        &   barcode,
                                        &   barcode_inverse,
                                            cell_id_to_prereqs,
                                        &   cell_dims,  
                                            ring.clone(),
                                        &   ok_to_start_new_lev_set,
                                            // last_must_be_crit,
                                    );

        let poly_correct_a      =   Polytope{
                                        data_c_to_l:    vec![ 0, 0, 1, 0, 2, 3],
                                        data_l_to_fmin: vec![ 0, 1, 2, 3, 4]
                                    };

        let poly_correct_b      =   Polytope{
                                        data_c_to_l:    vec![ 0, 0, 1, 0, 3, 2],
                                        data_l_to_fmin: vec![ 0, 1, 2, 3, 4]
                                    };  
        let poly_incorrect      =   Polytope{
                                        data_c_to_l:    vec![ 0, 0, 1, 3, 2, 0],
                                        data_l_to_fmin: vec![ 0, 1, 2, 3, 4]
                                    };                                      
                                    
        assert!(        is_barcode_compatible( &root_node, &poly_correct_a    ) );
        assert!(        is_barcode_compatible( &root_node, &poly_correct_b    ) );        
        assert!(    !   is_barcode_compatible( &root_node, &poly_incorrect    ) );              
    }  
    

    #[test]
    fn test_extra_conditions_to_start_new_level_set() {

        let simplex_sequence    =   vec![ vec![0], vec![1], vec![2], vec![0, 1], vec![0, 2], vec![1, 2], vec![ 0, 1, 2 ]];
        let simplex_bimap_sequential    =   BiMapSequential::from_vec( simplex_sequence );

        let condition_none      =   ConditionNone{};
        let condition_star      =   ConditionLowerStar{ simplex_bimap_sequential: & simplex_bimap_sequential };
        let condition_edge      =   ConditionLowerEdge{ simplex_bimap_sequential: & simplex_bimap_sequential };        

        let poly_star_n_edge_n  =   Polytope{
                                        data_c_to_l:    vec![ 0, 0, 0, 1, 1, 1, 7 ],
                                        data_l_to_fmin: vec![ 0, 1 ]
                                    };

        let poly_star_n_edge_y  =   Polytope{
                                        data_c_to_l:    vec![ 0, 1, 2, 7, 7, 7, 7 ],
                                        data_l_to_fmin: vec![ 0, 1, 2 ]
                                    };

        let poly_star_y_edge_y  =   Polytope{
                                        data_c_to_l:    vec![ 0, 1, 7, 1, 7, 7, 7 ],
                                        data_l_to_fmin: vec![ 0, 1, ]
                                    };    

        println!(
            "{:?} \n{:?}",
            poly_star_n_edge_n.cell_id_to_has_lev_set( 6 ), 
            vec![0, 1, 2]                                             // but all its vertices have
                .iter()
                .cloned()                                       
                .combinations( 1 )
                .all(   |x| 
                        poly_star_n_edge_n.cell_id_to_has_lev_set( 
                            simplex_bimap_sequential.ord( &x ).unwrap()                                
                        )
                        .unwrap()
                ) 
        );
                                    
        assert!(        condition_none.extra_condition_to_start_new_lev_set( &poly_star_n_edge_n)   );
        assert!(        condition_none.extra_condition_to_start_new_lev_set( &poly_star_n_edge_y)   );
        assert!(        condition_none.extra_condition_to_start_new_lev_set( &poly_star_y_edge_y)   );                

        assert!(    !   condition_star.extra_condition_to_start_new_lev_set( &poly_star_n_edge_n)   );
        assert!(    !   condition_star.extra_condition_to_start_new_lev_set( &poly_star_n_edge_y)   );
        assert!(        condition_star.extra_condition_to_start_new_lev_set( &poly_star_y_edge_y)   );                        

        assert!(    !   condition_edge.extra_condition_to_start_new_lev_set( &poly_star_n_edge_n)   );
        assert!(        condition_edge.extra_condition_to_start_new_lev_set( &poly_star_n_edge_y)   );
        assert!(        condition_edge.extra_condition_to_start_new_lev_set( &poly_star_y_edge_y)   );                                
    }

  
}
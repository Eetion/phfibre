
use crate::utilities::*;
use crate::intervals_and_ordinals::*;
use crate::rank_calculations::*;
use solar::reduce::vec_of_vec::{clear_cols};
use solar::utilities::index::{SuperVec, SuperIndex};
use solar::rings::ring::{Semiring, Ring, DivisionRing};
use num::rational::Ratio;
use std::collections::{HashSet, HashMap};
use std::iter::FromIterator;
use std;
use std::fmt::Debug;
use ordered_float::OrderedFloat;

type Cell = usize;
type Coeff = Ratio< i64 >;
type Fil = usize;
type FilRaw = OrderedFloat<f64>; // reason for this choice: f64 does not implement the hash trait (cricital for reparametrizing)
// type Ring = solar::rings::ring_native::NativeDivisionRing::< Ratio<i64> >;








// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  TO DO 

// PRE-CHECK THAT BARCODE AND COMPLEX ARE COMPATIBLE
// CHECK INTERSECTION 
// CALCULATE STATISTICS ABOUT POLYTOPES
// FUNCTION TO CHECK THAT BARCODE HAS BEEN COMPUTED CORRECTLY (IE COMPUTE BARCODE AFTER THE FACT, USING THE POLYHEDRON OBJECT)
// (NOW, UPON REFLECTION, I THINK WE HAVE TO BE CAREFUL ABOUT WHETHER THIS IS TRULY WELL-FOUNDED, MATHEMATICALLY) ADD A SCREENER FOR WHETHER A NEW POSTIVE CELL SHOULD HAVE FINITE OR INFINITE LIFE
// IMPLEMENT THE "CROSS-CHECK" STRATEGY AND COMPARE WITH ORIGINAL IMPLEMENTATION

    //  PRECOMPUTE
    //  compute: bars that stop and start at each endpoint
   
    //  FEASIBILITY CHECK
    //  compute: ranks of boundary operators (concomitantly, betti numbers)
    //  compute: bars_degn_quota[ dim ] = rank(boundary_(dim+1)) - #(finte bars of dimension dim)    
    //  check: #(inf bars of dimension dim) == betti_dim( num_cells_total space )
    //  check: #(fin bars of dimension dim) <= rank( num_cells_total space )

    //  INITIALIZE PARTIAL DATA

    //  EXPLORE THE TREE
    //  remove duplicate polyhedra as we go





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
pub struct Node < 'a, RingOp, RingElt > 
    where RingOp: Clone + Semiring<RingElt> + Ring<RingElt> + DivisionRing<RingElt>
{
    boundary:               Vec< Vec< ( Cell, RingElt ) > >,  // boundary matrix reprsented as a vector of vectors
    barcode:                &'a Barcode,                    // the target barcode (contains useful ordinal data)
    barcode_inverse:        &'a BarcodeInverse,             // pointer to a central register that maps bar endpoints to bar ids       

    ring:                   RingOp,
    boundary_buffer:        Vec< ( Cell, Coeff) >,          // a "holding space" for matrix entries, when these need to be moved around
    bc_endpoint_now:        usize,                          // the (ordinal) barcode endpoint of concern for this (or some future) level set

    bars_degn_quota:        SuperVec< usize >,              // number of degenerate cells anticipated in each dimension    

    lev_set_sizes:          LevelSetSizes,                  // Kth value = # cells in Kth level set
    polytope:               Polytope,    

    lev_set_is_crit:        bool,                           // true iff current level set contains a "critical cell"
    
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

}

impl < 'a, RingOp, RingElt > Node < 'a, RingOp, RingElt > 
    where   RingOp:     Clone + Semiring<RingElt> + Ring<RingElt> + DivisionRing<RingElt>,
            RingElt:    Clone + Debug + PartialOrd
{

    pub fn make_root(   
        boundary:               Vec< Vec< (Cell, RingElt)>>,
        barcode:            &'a Barcode,
        barcode_inverse:    &'a BarcodeInverse,        
        cell_dims:          &   Vec< usize >,   
        ring:                   RingOp,          
    ) 
    -> 
    Node< 'a, RingOp, RingElt >    
    {

        let num_cells               =   cell_dims.len();       
        let boundary_buffer         =   Vec::new();
        let bc_endpoint_now         =   0;

        // compute degenrate bar quotas
        let ranks                   =   chain_cx_rank_nullity(
                                            boundary.clone(), 
                                            ring.clone(),
                                            & cell_dims
                                        );
        let bars_degn_quota         =   num_degenerate_bars_per_degree(
                                            & ranks,
                                            & barcode
                                        );
        // level set sizes
        let lev_set_sizes           =   LevelSetSizes{ pointers: vec![0] };

        // polytope
        let polytope                =   Polytope{ 
                                            data_l_to_fmin:   Vec::new(), 
                                            data_c_to_l:      Vec::from_iter( 
                                                                std::iter::repeat( num_cells )
                                                                .take( num_cells )
                                                            )
                                        };
        
        // all empty level sets are degenerate (i.e. not critical)
        let lev_set_is_crit         =   false;

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
        for dim in chain_space_dim_per_deg.vec {                                  // initialize constituent vectors
            cell_ids_out.push( Vec::with_capacity(dim) )
        }
        for cell_id in 0..num_cells {                               // populate vectors
            cell_ids_out[ cell_dims[ cell_id ] ].push( cell_id );
        }



        let cell_ids_pos_crit       =   HashSet::new();
        let cell_ids_pos_degn       =   HashSet::new();        
        let bar_ids_dun_fin         =   Vec::new();
        let bar_ids_dun_inf         =   Vec::new();

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
            ring:                   ring.clone(),
            boundary_buffer:        boundary_buffer,
            bc_endpoint_now:        bc_endpoint_now,
            bars_degn_quota:        bars_degn_quota,
            lev_set_sizes:          lev_set_sizes,
            polytope:               polytope,
            lev_set_is_crit:        lev_set_is_crit,
            cells_all:              cells_all,
            cell_ids_out:           cell_ids_out,
            cell_ids_pos_crit:      cell_ids_pos_crit,
            cell_ids_pos_degn:      cell_ids_pos_degn,                            
            bar_ids_dun_fin:        bar_ids_dun_fin,
            bar_ids_dun_inf:        bar_ids_dun_inf,
            bar_ids_now_inf_brn:    bar_ids_now_inf_brn,
            bar_ids_now_fin_brn:    bar_ids_now_fin_brn,
            bar_ids_now_fin_die:    bar_ids_now_fin_die
        }     
    }

}


//  ---------------------------------------------------------------------------  
//  EXPLORE TREE
//  ---------------------------------------------------------------------------  


pub fn explore< RingOp, RingElt >( node: & Node< RingOp, RingElt >, results: &mut Vec< Polytope > )
    where   RingOp:     Clone + Semiring<RingElt> + Ring<RingElt> + DivisionRing<RingElt>,
            RingElt:    Clone + Debug + PartialOrd
{ 

    //  PUSH LEAF NODE TO RESULTS

    if  node.lev_set_sizes.size_last()          ==  Some( 0 ) &&
        node.lev_set_sizes.num_cells_total()    ==  node.cells_all.len() &&
        node.bar_ids_dun_fin.len()              ==  node.barcode.num_bars_fin() &&
        node.bar_ids_dun_inf.len()              ==  node.barcode.num_bars_inf()
    {
        // we may have alread constructed an equivalent filtration (just with a different order
        // on cells within a level set).  if we haven't then push to results.
        if ! results.contains( & node.polytope ) {
            results.push( node.polytope.clone() ); 
        }
    }

    //  TERMINATE CONSTRUCTION OF PRESENT LEVEL SET; INITIALIZE NEW LEVEL SET 
    //  ---------------------------------------------------------------------
     
    else if node.lev_set_sizes.size_last() != Some( 0 ) &&  // current level set is nonempty 
            node.cell_ids_pos_degn.is_empty()       &&      // C-rule (degenerate) there are no unmatched "degenerate" positive cells
            (                                               // C-rule (critical) 
                ! node.lev_set_is_crit                          // the level set contains no "critical" cell
                ||                                              // OR
                (                                               // I_cur is empty
                node.bar_ids_now_inf_brn.is_empty() &&     
                node.bar_ids_now_fin_brn.is_empty() &&    
                node.bar_ids_now_fin_die.is_empty()    
                )
            )                                       &&
            node.barcode                                    // the current level set does not map to 1.0
                .ordinal
                .val( node.bc_endpoint_now )
                != Some( OrderedFloat( 1.0 ) )

        {

        // CREATE CHILD NODE    
    
        let mut child                   =   node.clone();

        // APPEND NEW COUNTER FOR LEVEL SET SIZE

        child.lev_set_sizes.postpend_empty_set();       // post-pend new entry (= 0) to vector of level set sizes
        child.polytope.push_new_lev_set();              // post-pend new level set to polytope object (note this object doesn't have a well formed notion of an empty level set)

        // IF PARENT NODE IS CRITICAL, THEN (1) UPDATE BARCODE ENDPOINT, AND (2) UPDATE SET OF BARS TO ACCOUNT FOR

        // NB: we don't check that child.bc_endpoint_now < maximum_endpoint; nothing bad happens
        // if this is indeed the case; one simply obtains three empty sets (moreover, this doesn't
        // produce an infinite loop, because we only reach this point after creating a new
        // nonempty level set.

        if  node.lev_set_is_crit  { 

            // since we are creating a new level set, place child.lev_set_is_crit = false
            // this value remains false as long as the new level set remains empty 
            child.lev_set_is_crit       =   false;          
            
            // update barcode endpoint
            child.bc_endpoint_now       +=  1;

            // update set of bars to account for
            child.bar_ids_now_inf_brn   =   child.barcode_inverse
                                                .inf_brn
                                                .sindex( child.bc_endpoint_now, vec![] )
                                                .clone();

            child.bar_ids_now_fin_brn   =   child.barcode_inverse
                                                .fin_brn
                                                .sindex( child.bc_endpoint_now, vec![] )
                                                .clone();            
            
            child.bar_ids_now_fin_die   =   child.barcode_inverse
                                                .fin_die
                                                .sindex( child.bc_endpoint_now, vec![] )
                                                .clone();                        

        }

        // RUN EXPLORE ON CHILD
        
        println!("CALL 1");
        explore( & child, results );
        
    }

    //  ENLARGE LEVEL SET 
    //  ---------------------------------------------------------------------
    else {

        //  ADD BIRTH CELL FOR AN *INFINITE* BAR 
        //  ------------------------------------
        
        for (bar_id_count, bar_id) in node.bar_ids_now_inf_brn.iter().cloned().enumerate() { 

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

                println!("ADDING INFINITE POSITIVE CELL");

                // update: 
                //  lev_set_sizes 
                //  cells_all
                //  cell_ids_out 
                //  bar_ids_dun 
                //  bar_ids_now

                // create child node
                let mut child   =   node.clone();                                   // this must come before swap_remove

                // mark child's level set as critical (posisbly redundantly)
                child.lev_set_is_crit   =   true;

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
                
                println!("CALL 2");                
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
                println!("ADDING FINITE POSITIVE CELL");

                // create child node
                let mut child   =   node.clone();                                   // this must come before swap_remove

                // mark child's level set as critical (posisbly redundantly)
                child.lev_set_is_crit   =   true;
                
                // move bar_id
                let _    =   child.bar_ids_now_inf_brn.swap_remove(bar_id_count);   // remove bar from list of bars that have unaccounted endpoints here
                child.bar_ids_dun_inf.push( bar_id );                               // add bar to "done" list (consumes variable)

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
                
                println!("CALL 3");                
                explore( & child, results );

            }
        } 
        
        //  ADD DEATH CELL FOR A *FINITE* BAR 
        //  ------------------------------------

        // NB:  PROBABLY AN OPORTUNITY TO IMPROVE PERFORMANCE HERE: rather than loop bar by bar, just create a 
        //      HASHMAP : (BIRTH, DEATH) -> MULTIPLICITY; loop once over all columns and for each column check 
        //      WHETHER THE BAR IT WOULD CREATE WOULD LIE IN THE DESIRED SET

        for (bar_id_count, bar_id) in node.bar_ids_now_fin_die.iter().cloned().enumerate() {

            // clone info about bar to be added
            let bar_new         =   node.barcode
                                        .bar_fin( bar_id );

            // loop over all dimension (dim+1) chains
            for (neg_id_count, neg_id) in   node.cell_ids_out
                                                [ bar_new.dim() + 1 ]
                                                .iter()
                                                .cloned()
                                                .enumerate()
            {

                // find the lowest nonzero entry in this column
                let low_id_opt        =    node.boundary
                                                [ neg_id ]  
                                                .iter()
                                                .map( |x| x.0 )
                                                .max()
                                                .map( |x| x.clone() );

                // short-circuit if the column is empty                                                 
                if low_id_opt == None { continue }

                // unwrap the bottom index
                let low_id              =   low_id_opt.unwrap();

                //   IF  every cell in the support of this column has been assigned a birth ordinal 
                //  AND  low_id corresponds to a positive cell with the correct birth time
                // THEN  add a positive critical cell
                if  low_id < node.cells_all.len()   // this first condition helps to short circuit / avoid expensive look-ups in the hash set
                    &&  
                    node.polytope.cell_id_to_fmin( low_id.clone() ).unwrap() == bar_new.birth()
                    &&
                    node.cell_ids_pos_crit
                        .contains( & low_id )
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
                    child.lev_set_is_crit   =   true;

                    // move cells
                    let _           =   child.cell_ids_out[ bar_new.dim() ].swap_remove(neg_id_count);           // remove negative cell from "out list"
                    child.cell_ids_pos_crit.remove( & low_id );                                   // remove positive cell from set of unmatched positive cells

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
                        [ low_id ]
                        .bounding_cell_id
                                    =   neg_id.clone();


                    // UPDATE THE BOUNDARY MATRIX

                    // define the ingredients for clearing
                    let clearor     =   child.boundary[ neg_id ].clone();
                    let pivot_entry =   clearor.last().unwrap();
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
                    
                    println!("CALL 4");                    
                    explore( & child, results );
                }
            }
        }

        //  ADD BIRTH CELL FOR A *DEGENERATE* BAR 
        //  -------------------------------------

        // loop over all dimensions
        for dim in 0..node.cell_ids_out.len() {

            // short-circuit if we have already met the quota for degenerate bars in this dimension
            if node.bars_degn_quota.val( dim ) == 0 { continue }

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
                child.bars_degn_quota.vec[dim] -= 1;         // drop the quota for degenerate bars by 1

                // RUN EXPLORE ON CHILD
                
                println!("CALL 5");                
                explore( & child, results );
            }
        }
        
        //  ADD DEATH CELL FOR A *DEGENERATE* BAR 
        //  -------------------------------------

        // loop over all dimensions
        for dim in 1..node.cell_ids_out.len() {

            // loop over chains of given dimension with nonzero boundary
            for (neg_id_count, neg_id) in   node.cell_ids_out
                                                [ dim ]
                                                .iter()
                                                .cloned()
                                                .enumerate()
            {

                // find the lowest nonzero entry in this column
                let low_birth_and_id_opt    =   node.boundary
                                                    [ neg_id ]  
                                                    .iter()
                                                    .map(   |x| 
                                                            (   
                                                                node.cells_all[ x.0 ].birth_ordinal.clone(),   
                                                                x.0.clone())                            
                                                            )
                                                    .max();                

                // short-circuit if the column is empty                                                 
                if low_birth_and_id_opt == None { continue }

                // unwrap the bottom index
                let low_birth_and_id            =   low_birth_and_id_opt.unwrap();


                //   IF  every cell in the support of this column has been assigned a birth ordinal 
                //  AND  low_id corresponds to a positive degenerate cell
                // THEN  add a negative degenerate cell
                if  low_birth_and_id.0 < node.cells_all.len() &&  // this first condition helps to short circuit / avoid expensive look-ups in the hash set
                    node.cell_ids_pos_degn.contains( & low_birth_and_id.1 ) 
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
                    child.cell_ids_pos_degn.remove( & low_birth_and_id.1 );                     // remove positive cell from set of unmatched positive cells

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
                        [ low_birth_and_id.1 ]
                        .bounding_cell_id
                                    =   neg_id.clone();

                    // UPDATE THE BOUNDARY MATRIX

                    // define the ingredients for clearing
                    let clearor     =   child.boundary[ neg_id ].clone();
                    let pivot_entry =   clearor.last().unwrap();
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
                    
                    println!("CALL 6");                    
                    explore( & child, results );
                }
            }
        }
    }
}


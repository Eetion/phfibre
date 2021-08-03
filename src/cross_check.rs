
use crate::utilities::*;
use crate::intervals_and_ordinals::*;
use crate::rank_calculations::{reduced_mat_to_pivot_index_pairs, chain_cx_rank_nullity, num_degenerate_bars_per_degree};
use solar::reduce::vec_of_vec::{clear_cols};
use solar::utilities::index::{SuperVec, SuperIndex, sort_perm, inverse_perm, histogram};
use solar::rings::ring::{Semiring, Ring, DivisionRing};
use num::rational::Ratio;
use std::collections::{HashSet, HashMap};
use std::hash::Hash;
use std::iter::FromIterator;
use std;
use std::fmt::Debug;
use ordered_float::OrderedFloat;

type Cell = usize;
type Coeff = Ratio< i64 >;
type Fil = usize;
type FilRaw = OrderedFloat<f64>; // reason for this choice: f64 does not implement the hash trait (cricital for reparametrizing)
// type Ring = solar::rings::ring_native::NativeDivisionRing::< Ratio<i64> >;



//  possibly push result (involves a check of whether the current level set is critical and if so contains all necessary cells)
//  possibly start new level set
//  for s in cell_ids_out
//      if boundary(s) = 0
//          if there's space in pair_quota + crit_quota
//              add s to complex and mark ask positive
//      else if max supp( boundary(s) ) < #(K_in)
//          if there space in pair_quota + degn_quota
//              add s to complex and clear boundary matrix







//  let N = # simplices in total space
//  let E = # barcode endpoints
//  for each perm in  
//      [permutations that respect face relations]
//      reduce boundary matrix
//          for each [partitions of N] x [N choose E] 
//             compute full barcode
//             add to set of feasible data if barcode correct 



fn pos_quota< FilRaw, RingOp, RingElt >(  node: & Node< FilRaw, RingOp, RingElt >, dim ) -> usize {
    
}




//  ---------------------------------------------------------------------------  
//  EXPLORE TREE
//  ---------------------------------------------------------------------------  


pub fn explore< FilRaw, RingOp, RingElt >( node: & Node< FilRaw, RingOp, RingElt >, results: &mut Vec< Polytope > )
    where   RingOp:     Clone + Semiring<RingElt> + Ring<RingElt> + DivisionRing<RingElt>,
            RingElt:    Clone + Debug + PartialOrd,
            FilRaw:     Clone + Debug + Ord + Hash
{ 

    if node.partial_coface_exists_in(& results ) {
        return
    } elseif node.is_feasible_leaf() {
        results.push( node.polytope )
        return
    }

    if node.can_start_new_lev_set() {
        node.start_new_lev_set();
        return
    }

    for cell_id in node.cell_ids_out {
        if node.boundary[ cell ].is_empty{
            if node.quota_pos_cell( node.cell_dim( cell_id ) ) > 0 {  // quota encorporates both crit and degn
                node.add_pos_cell()
            }
        } else {
            calc birth + death time
            if node.quota_neg_cell ( birth + death time + dim) // quota encorporates both crit and degn  !!! AND EXCLUDES DEGEN CELLS OF UNHELPFUL DEGREES FROM CRITICAL CLASSES
                node.add_neg_cell()
        }
    }

      

    //  PUSH LEAF NODE TO RESULTS

    if  node.lev_set_sizes.size_last()          ==  Some( 0 ) &&
        node.lev_set_sizes.num_cells_total()    ==  node.cells_all.len() &&
        node.bar_ids_dun_fin.len()              ==  node.barcode.num_bars_fin() &&
        node.bar_ids_dun_inf.len()              ==  node.barcode.num_bars_inf() &&
        (            
            node.polytope.lev_set_last_is_critical() == Some( true )
            ||
            ! node.last_must_be_crit
        )           
        
        {

        // we may have alread constructed an equivalent filtration (just with a different order
        // on cells within a level set).  if we haven't then push to results.
        if ! results.contains( & node.polytope ) {
            println!("num results: {:?}",  histogram( results.iter().map(|x| x.dim().unwrap() ) )   );
            results.push( node.polytope.clone() ); 
        }

    }

    else {    

        //  IF FEASIBLE, INITIALIZE NEW LEVEL SET 
        //  ---------------------------------------------------------------------
     
        if node.lev_set_sizes.size_last() != Some( 0 ) &&  // current level set is nonempty 
                node.cell_ids_pos_degn.is_empty()       &&      // C-rule (degenerate) there are no unmatched "degenerate" positive cells
                (                                               // C-rule (critical) 
                    (
                    !node.polytope                           // the level set contains no "critical" cell
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

            if  node.polytope.lev_set_last_is_critical().unwrap()  { 
                
                // update barcode endpoint

                // child.polytope.ensure_last_lev_set_critical();
                
                let bc_endpoint_now         =   child.polytope.last_of_all_filtration_ordinals().unwrap() + 1;


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
            return
            
        }

        //  ENLARGE CURRENT LEVEL SET 
        //  ---------------------------------------------------------------------

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
        
        //  ADD DEATH CELL FOR A *FINITE* BAR 
        //  ------------------------------------

        // NB:  PROBABLY AN OPORTUNITY TO IMPROVE PERFORMANCE HERE: rather than loop bar by bar, just create a 
        //      HASHMAP : (BIRTH, DEATH) -> MULTIPLICITY; loop once over all columns and for each column check 
        //      WHETHER THE BAR IT WOULD CREATE WOULD LIE IN THE DESIRED SET

        for (bar_id_count, bar_id) in node.bar_ids_now_fin_die.iter().cloned().enumerate() {

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
  
//-------------------------------------------------------------------------------------------------                
                
                // if node.polytope.min_vertex_is_compatible_with_ordinal_filt( &vec![2, 0, 1, 2, 2, 3] ) {  
                    
                //     println!("!!!!!!!!!! ATTEMPTING TO ADD NEGATIVE CELL PARENT polytope:  {:?} {:?}", & node.polytope.vec_mapping_cell_id_to_min_filt_ordinal(), & node.polytope.data_c_to_l );                     
                //     println!("{:?} {:?} {:?} {:?}",  
                //                                 low_birth_ord < node.cells_all.len(),
                //                                 node.polytope.cell_id_to_fmin( low_cell_id.clone() ) ==  Some( bar_new.birth() ),
                //                                 node.cell_ids_pos_crit.contains( & low_cell_id ),
                //                                 & low_cell_id 
                //                             )
                // }                   


//-------------------------------------------------------------------------------------------------                                    
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

        //  ADD BIRTH CELL FOR A *DEGENERATE* BAR 
        //  -------------------------------------

        // if node.polytope.min_vertex_is_compatible_with_ordinal_filt( &vec![2, 0, 1, 2, 2, 3] ) {   
        //     println!("WILL TRY TO ADD SOME DEGENERATE CELLS!");
        //     println!("NODE WAS COMPATIBLE FROM THE START: {:?}", first_pass_was_true)
        // }          

        // loop over all dimensions
        for dim in 0..node.cell_ids_out.len() {

            // short-circuit if we have already met the quota for degenerate bars in this dimension
            if node.bars_degn_quota.sindex( dim, 0 ) == 0 { 
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
                child.bars_degn_quota[dim] -= 1;         // drop the quota for degenerate bars by 1

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
        
        //  ADD DEATH CELL FOR A *DEGENERATE* BAR 
        //  -------------------------------------

        // loop over all dimensions
        for dim in 1..node.cell_ids_out.len() {

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
    }
}
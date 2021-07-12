


type Cell = usize;
type Val = Ratio< i64 >;
type Fil = f64;


//  ---------------------------------------------------------------------------  
//  UTILITIES 
//  ---------------------------------------------------------------------------  


fn end_val< T >( v: Vec<T> )  { (v[ v.len() - 1]).clone() }
fn end_val_mut< T >( v: Vec<T> )  { (v[ v.len() - 1]) }

fn lev_set_size< T >( v: Vec<T>, ind: usize ) {
    if ind == 0 { return v[0] }
    else { return v[ind] - v[ind-1] }
}


//  ---------------------------------------------------------------------------  
//  ENTRIES 
//  ---------------------------------------------------------------------------  


struct EBarcodeTotal{
    birth_time:     Fil,
    death_time:     Fil,
    birth_cell:     Option< Cell >,
    death_cell:     Option< Cell >,
    dim:            usize
}


struct ETotalComplex{
    birth_order:    usize,
    birth_class:    usize,
    birth_bar:      Option< usize >,
    death_bar:      Option< usize >,
    dim:            usize
}


// struct ETotalBarcode{
//     birth_index:    usize,
//     death_index:    Option< usize >,
//     dim:            usize,
// }
// 
// 
// struct EPositiveCell{
//     cell:           Cell,
//     death_cells:    Vec< Cell >,
//     interval:       Option< usize >, // uuid of the corresponding interval 
// }


struct EResults{
    birth_class:    usize,
    birth_order:    usize,
    birth_itv:      Option< usize >,
    death_itv:      Option< usize >
}


//  ---------------------------------------------------------------------------  
//  CONTAINERS
//  ---------------------------------------------------------------------------  



// struct Node
// {
//     boundary:           Vec< Vec< ( Cell, Val ) > >,
//     num_cells_born    usize,
//     bc_endpoint_now         usize,
//     
//     cells_all:          Vec< ETotalComplex >,   // all cells
//     cell_ids_out:          Vec< usize >            // cells not yet assigned a birth,
//     cell_ids_pos_crit:     Vec< EPositiveCell >,   // unmatched positive critical cells
//     cell_ids_pos_degn:     Vec< EPositiveCell >,   // unmatched positive degenerate cells
// 
//     bar_ids_dun:           Vec< usize >,           // bars with all endpoints accounted for
//     bar_ids_und:           Vec< usize >,           // bars without all endpoints accounted for 
//     bars_all:           Vec< ETotalBarcode>,    // all bars
//     bar_ids_now:           Vec< usize >,           // bars still to match in this batch
// 
//     class_id:           usize,
//     class_size:         usize,
//     class_crit:         bool 
// }


// struct Node
// {
//     boundary:           Vec< Vec< ( Cell, Val ) > >,
//     num_cells_born      usize,
//     bc_endpoint_now        usize,
// 
//     cells_all:          Vec< ETotalComplex >,   // all cells
//     cell_ids_out:          Vec< usize >            // cells not yet assigned a birth,
//   
//     bar_ids_crit_all:          Vec< BarPair >,     // all bars in desired barcode; 
//     bar_ids_crit_needbirth:    Vec< usize >,       // bars not yet assigned a birth cell
//     bar_ids_crit_needdeath:    Vec< usize >,       // bars assigned a birth cell but no death cell
// 
//     bar_ids_degn_all:          Vec< BarPair >,     // all degenerate bars 
//     bar_ids_degn_needdeath:    Vec< usize >,       // bars assigned a birth cell but no death cell
// 
// 
//     bar_ids_fut_crit:      Vec< usize >,           // bars without all endpoints accounted for 
//     bar_ids_fut_degn:      Vec< usize >,           // bars without all endpoints accounted for 
//     
//     bar_ids_now_crit:      Vec< usize >,           // bars still to match in this batch
//     bar_ids_now_degn:      Vec< usize >,           // bars still to match in this batch
// 
//     class_id:           usize,
//     class_size:         usize,
//     class_crit:         bool 
// }


struct Node
{

    boundary:           Vec< Vec< ( Cell, Val ) > >,
    bc_endpoint_now:    usize,
    
    lev_set_sizes:      Vec< usize >,           // Kth value = # cells in Kth level set
    lev_set_is_crit:    bool,                   // true iff current level set contains a "critical cell"
    
    cells_all:          Vec< ETotalComplex >,   // all cells
    cell_ids_out:          Vec< usize >            // cells not yet assigned a birth,
    
    cell_ids_pos_crit:     Hash< (Fil, usize), Vec< usize > >,   // unmatched positive critical cells, grouped by (birth_time, dimension)
    cell_ids_pos_degn:     Vec< usize >,           // unmatched positive degenerate cells
                                                // NB: doesn't need to be hash; we know birth time
    bars_all_inf:       Vec< Fil >,
    bars_all_fin:       Vec< (Fil, Fil) >,
    
    bar_ids_dun_fin:       Vec< usize >,           // nonempty bars with all endpoints accounted for (finite)
    bar_ids_dun_inf:       Vec< usize >,           // nonempty bars with all endpoints accounted for (infinite)
    
    bar_ids_und_fin:       Vec< usize >,           // nonempty bars with 1 endpoint left to account for 

    bar_ids_now_inf:       Vec< usize >,           // bars still to match in this batch
    bar_ids_now_fin:       Vec< usize >,           // bars still to match in this batch

}




//  ---------------------------------------------------------------------------  
//  FUNCTIONS 
//  ---------------------------------------------------------------------------  


fn format_result( node: Node ) -> EResults
{
    
}

fn explore( node: &mut Node, results: &mut Vec< Vec< ETotalComplex> > )
{
   
    //  PRECOMPUTE SOME SIZES
    //  ** ADD A CHECK TO MAKE SURE WE HAVE THE RIGHT NUMBER OF CELLS IN ALL DIMENSIONS
    //          - CAN CALCULATE THIS RECURSIVELY STARTING WITH DIM 0
    //  ** COULD ALSO DO AN INITIAL CHECK TO MAKE SURE HOMOLOGY OF LAST SPACE IN SEQUENCE IS
    //  COMPATIBLE WITH THE BARCODE PROVIDED

    //  PUSH LEAF NODE TO RESULTS

    if  node.lev_set_sizes.sum()  ==  node.all_cells.len() &&
        node.bar_ids_dun_fin.len() == node.bars_all_fin.len() &&
        node.bar_ids_dun_inf.len() == node.bar_ids_dun_inf.len() &&
        node.class_size         == 0
        {
        results.push( node.cells_all.clone() ); 
    }

    //  TERMINATE CONSTRUCTION OF PRESENT LEVEL SET; INITIALIZE NEW LEVEL SET 
     
    else if end_val(node.lev_set_sizes) != 0 &&     // current class is nonempty 
            node.cell_ids_pos_degn.is_empty()   &&     // C-rule (degenerate) there are no unmatched "degenerate" positive cells
            (                                       // C-rule (critical) 
                ! node.lev_set_is_crit                  // the level set contains no "critical" cell
                ||                                      // OR
                (                                       // I_cur is empty
                node.bar_ids_now_inf.is_empty()     &&     
                node.bar_ids_now_fin_brn.is_empty() &&    
                node.bar_ids_now_fin_die.is_empty() &&   
                )
            )
        {
        
        // CREATE CHILD NODE    
    
        let mut child           =   node.clone();
        child.lev_set_sizes.push(0);                // post-pend new entry (= 0) to vector of level set sizes
        chile.lev_set_is_crit   =   false;          // this value remains false as long as the new level set remains empty

        // IF PARENT NODE IS CRITICAL, THEN (1) UPDATE BARCODE ENDPOINT, AND (2) UPDATE SET OF BARS TO ACCOUNT FOR

        // NB: we don't check that child.bc_endpoint_now < maximum_endpoint; nothing bad happens
        // if this is indeed the case; one simply obtains three empty sets (moreover, this doesn't
        // produce an infinite loop, because we only reach this point after creating a new
        // nonempty level set.

        if  child.lev_set_is_crit  { 

            // update barcode endpoint
            child.bc_endpoint_now       +=  1;

            // update set of bars to account for
            for bar in node.bars_all_inf.filter(|x| x.birth == child.bc_endpoint_now) {
                child.bar_ids_now_inf.push( bar );
            }
            for bar in node.bars_all_fin.filter(|x| x.birth == child.bc_endpoint_now) {
                child.bar_ids_now_brn.push( bar );
            }
            for bar in node.bars_all_fin.filter(|x| x.death == child.bc_endpoint_now) {
                child.bar_ids_now_die.push( bar );
            }

        }

        results                 =   explore( child, results );

    
        // // fork 1 for child: new critical class
        // if ! child.bar_ids_und.is_empty() {            // if there are bars yet to match
        //     child.class_crit    =   true;           // flag that we are starting a critical class
        //     child.bar_ids_now      =   bar_ids_und        // calculate bars with this endpoint
        //                                 .iter()
        //                                 .map(|x| child.bars_all[ x ] )
        //                                 .filter(|x|     x.birth_time == child.bc_endpoint_now
        //                                                 ||
        //                                                 x.death_time == Some(child.bc_endpoint_now)
        //                                         )
        //                                 .collect();
        //     results             =   explore( child.clone(), results );
        // }

        // // fork 2 for child: new degenerate class
        // child.class_crit        =   false;      // flag that we're creating a degenerate class
        // child.bar_ids_now.clear();                 // remove all bars we may have added to the  "docket"
        // results                 =   explore( child, results ); 

        // }
        
    }

    //  ENLARGE LEVEL SET 
    else {

        //  ADD BIRTH CELL FOR AN *INFINITE* BAR 
        
        for ind_bar_new in node.bar_ids_now_inf {

            // update: 
            //  lev_set_sizes 
            //  cells_all
            //  cell_ids_out 
            //  bar_ids_dun 
            //  bar_ids_now

            // clone info about bar to be added
            let mut bar_new     =   node
                                        .bars_all_inf
                                        [ind_bar_new]
                                        .clone();  

            // loop over all cycles of appropriate dimension 
            for ind_out in 0..node.cell_ids_out[ bar_new.dim ].len()
                .filter(    |x| 
                            node.boundary
                                [  
                                    node.cell_ids_out
                                        [ bar_new.dim ]
                                        [ ind_out ]
                                ]
                            .is_empty()
                    )  
            {   
                
                // create child node
                let mut child   =   node.clone();                               // this must come before swap_remove

                // mark child's level set as critical (posisbly redundantly)
                child.lev_set_is_crit   =   true;

                // move bar and positive cell
                let mut cell    =   child.cell_ids_out.swap_remove(ind_out);       // remove cell from "out list"
                let mut _       =   child.bar_ids_now_inf.swap_remove(ind_out);    // remove bar from "unadded list"
                child.bar_ids_dun_inf.push( bar_new );                             // add bar to "done" list (consumes variable)
               
                // record the birth ordinal of the added cell in the central registry
                child
                    .cells_all
                        [ cell ]
                    .birth
                                =   end_val( child.lev_set_sizes );
               
                // update number of cells born
                end_val_mut( child.lev_set_sizes ) += 1; 

            }
        } 
        
        //  ADD BIRTH CELL FOR A *FINITE* BAR 

        for ind_bar_new in node.bar_ids_now_fin_brn {

            // update: 
            //  lev_set_sizes 
            //  cells_all
            //  cell_ids_out 
            //  bar_ids_dun 
            //  bar_ids_now

            // clone info about bar to be added
            let mut bar_new     =   node.bars_all_fin
                                        [ind_bar_new]
                                        .clone();  

            // loop over all cycles of appropriate dimension 
            for ind_out in 0..node.cell_ids_out[ bar_new.dim ].len()
                .filter(    |x| 
                            node.boundary
                                [  
                                    node.cell_ids_out
                                        [ bar_new.dim ]
                                        [ ind_out ]
                                ]
                            .is_empty()
                    )  
            {   
                
                // create child node
                let mut child   =   node.clone();                               // this must come before swap_remove

                // mark child's level set as critical (posisbly redundantly)
                child.lev_set_is_crit   =   true;

                // move the positive cell
                let mut cell    =   child.cell_ids_out.swap_remove(ind_out);       // remove cell from "out list"
                child.cell_ids_pos_crit
                    .get_mut(
                                ( child.bc_endpoint_now.clone(), bar_new.dim.clone() )
                        )
                    .push( cell );

                // remove bar from list of bars that have unaccounted endpoints here
                let mut _       =   child.bar_ids_now_fin.swap_remove(ind_bar_new);    // remove bar from "unadded list"
               
                // record the birth ordinal of the added cell in the central registry
                child
                    .cells_all
                        [ cell ]
                    .birth
                                =   end_val( child.lev_set_sizes );
               
                // update number of cells born
                end_val_mut( child.lev_set_sizes ) += 1; 

            }
        } 
        
        //  ADD DEATH  CELL FOR A *FINITE* BAR 

        for ind_bar_new in node.bar_ids_now_fin_die {

            // update: 
            //  lev_set_sizes 
            //  cells_all
            //  cell_ids_out 
            //  bar_ids_dun 
            //  bar_ids_now

            // clone info about bar to be added
            let mut bar_new     =   node.bars_all_fin
                                        [ind_bar_new]
                                        .clone();  

            // loop over all cycles of appropriate dimension 
            for ind_out in 0..node.cell_ids_out[ bar_new.dim ].len() {

                snz_ind_max     =   
                                    node.boundary
                                        [  
                                            node.cell_ids_out
                                                [ bar_new.dim ]
                                                [ ind_out ]
                                        ]
                                    .iter()
                                    .enumerate()
                                    .max_by( 
                                                |(x,y)| 
                                                node.cells_all[ y.0 ].birth
                                        )



                                        node.boundary
                                            [  
                                                node.cell_ids_out
                                                    [ bar_new.dim ]
                                                    [ ind_out ]
                                            ]


                .filter(    |x| 
                            node.boundary
                                [  
                                    node.cell_ids_out
                                        [ bar_new.dim ]
                                        [ ind_out ]
                                ]
                            .iter()
                            .map( |x| 
                                  node.cells_all
                                    [x.0]
                                    .birth
                                )
                            .empty_max_equals( bar_new.birth )
                    )

            {   
                
                // create child node
                let mut child   =   node.clone();                               // this must come before swap_remove

                // mark child's level set as critical (posisbly redundantly)
                child.lev_set_is_crit   =   true;

                // move the positive cell
                let mut cell    =   child.cell_ids_out.swap_remove(ind_out);       // remove cell from "out list"
                child.cell_ids_pos_crit
                    .get_mut(
                                ( child.bc_endpoint_now.clone(), bar_new.dim.clone() )
                        )
                    .push( cell );

                // remove bar from list of bars that have unaccounted endpoints here
                let mut _       =   child.bar_ids_now_fin.swap_remove(ind_bar_new);    // remove bar from "unadded list"
               
                // record the birth ordinal of the added cell in the central registry
                child
                    .cells_all
                        [ cell ]
                    .birth
                                =   end_val( child.lev_set_sizes );
               
                // update number of cells born
                end_val_mut( child.lev_set_sizes ) += 1; 

            }
        } 
// -----------------------------------------------------------------------------------------
        for ind_bar_new in node.bar_ids_now_fin_brn {

            let mut bar_new     =   node
                                        .bars_all
                                        [ ind_bar_new ]
                                        .clone();   
            let mut pos_key     =   ( 
                                        bar_new.left, 
                                        bar_new.dim 
                                    )

            for ind_out in 0..node.cell_ids_out.len()
                .filter(    |x| 
                            node.boundary[  node.cell_ids_out[ ind_out ] ].is_empty() &&
                            node.all_cells[ node.cell_ids_out[ ind_out ] ].dim == bar_new.dim     )

                {
                    let mut child   =   node.clone(); // this must come before swap_remove
                    let cell        =   child.cell_ids_out.swap_remove(out_ind);

                    if child.cell_ids_pos_crit.has_key(    ( 
                                                        bar_new.left, 
                                                        bar_new.dim 
                                                        ) 
                                                    ) 
                    {
                        chold.cell_ids_pos_crit.get_mut(
                    }
                    {
                        child.cell_ids_pos_crit.push( cell.clone() ); // DO WE REALLY WANT TO DO THIS?
                        // update: 
                        //  lev_set_sizes
                        //  cells_all
                        //  cell_ids_out (already done)
                        //  cell_ids_pos_crit (already done -- but do we want this?)
                        //  bar_ids_dun (if bar never dies)
                        //  bar_ids_und (if bar never dies)
                        //  bars_all
                        //  bar_ids_now
                        //  class_size
                    }
                }
            } 
        }

        for bar in node.bar_ids_now_end {
            for cell in node.cell_ids_pos_crit.get(    (
                                                    node.bc_endpoint_now,
                                                    node.bars_all[ bar ]
                                                        .dim()

        }






        for bar in node.bar_ids_now {
            
            //  ADD A BAR_DEATH CELL
            if bar.birth_time  ==  node.bc_endpoint_now{
                for out_ind in 0..node.cell_ids_out.len()
                    .filter(    |x| 
                                node.boundary[  node.cell_ids_out[ind] ].is_empty() &&
                                node.all_cells[ node.cell_ids_out[ind] ].dim == bar.dim     )

                    {
                    let mut child   =   node.clone(); // this must come before swap_remove
                    let cell    =   child.cell_ids_out.swap_remove(out_ind);
                    child.cell_ids_pos_crit.push( cell.clone() ); // DO WE REALLY WANT TO DO THIS?
                    // update: 
                    //  lev_set_sizes
                    //  cells_all
                    //  cell_ids_out (already done)
                    //  cell_ids_pos_crit (already done -- but do we want this?)
                    //  bar_ids_dun (if bar never dies)
                    //  bar_ids_und (if bar never dies)
                    //  bars_all
                    //  bar_ids_now
                    //  class_size
                }
            } 

            //  ADD A BAR_BITH CELL
            else if bar.death_time == Some( node.bc_endpoint_now ) {
                
                for cell_birth in node.cell_pos_crit { // ASK JACOB IF THIS FOR-LOOP IS OK
                    if  {
                            
                    }
                }
                // STUFF
            }


        }

        for pos_degen in node.cell_ids_pos_degn {
            // STUFF
        }

        for cell in node.cell_ids_out
                        .iter()
                        .cloned()
                        .filter(|x| node.boundary[ x ].is_empty() )
            {
            // STUFF
        }
    }


}

// the added class is born in interval (-inf, birth_of_optimized_class]
// AND
// the added class dies in    interval (birth_of_optimized_class, death_of_optimized class)
let condition_a =   chx.key_2_filtration(&index_birth)    <=  chx.key_2_filtration(&birth) &&
                    chx.key_2_filtration(&index_death)    >   chx.key_2_filtration(&birth) &&
                    chx.key_2_filtration(&index_death)    <   chx.key_2_filtration(&death) 

// the added class is born strictly before the birth simplex (in lexicographic order) 
// AND
// the added class dies when the optimized class dies
let condition_b =   index_birth                           <   birth                        && 
                    chx.key_2_filtration(&index_death)    ==  chx.key_2_filtration(&death) 



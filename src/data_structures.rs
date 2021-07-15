


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

fn last_index< T > ( v: Vec<T> ) -> usize { v.len() - 1 }

//  ---------------------------------------------------------------------------  
//  SPARSE POINTERS 
//  ---------------------------------------------------------------------------  

#[derive(Clone, Debug)]
struct LevelSetSizes{ pointers: Vec< usize > }

impl LevelSetSizes{
    
    fn size( &self, set_index: usize ) -> usize {  
        if set_index == 0 { self.pointers[0] }
        else { self.pointers[ set_index ] - self.pointers[ set_index - 1 ] }
    }

    fn size_last( &self ) -> usize { self.set_size( self.pointers.len() - 1) }
    
    fn total( &self ) -> usize {
        if self.pointers.is_empty() { 0 }
        else { end_val( self.pointers ).clone() }
    }
    
    fn postpend_empty_set( &self ) { 
        if self.pointers.is_empty() { self.pointers.push(0) }
        else { self.pointers.push( end_val(self.pointers) ) }
    }

    fn grow_last_set( &self ) {
        self.pointers
            [ last_index(self.pointers) ]  
        =
        self.pointers
            [ last_index(self.pointers) ]  
        + 1 
    }
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

#[derive(Clone, Debug)]
struct Node
{

    boundary:           Vec< Vec< ( Cell, Val ) > >,
    bc_endpoint_now:    usize,
    

    bars_degn_quota:    Vec< usize >,

    lev_set_sizes:      LevelSetSizes,           // Kth value = # cells in Kth level set

    lev_set_values:     Vec< Fil >,             // the function phi
    lev_set_is_crit:    bool,                   // true iff current level set contains a "critical cell"
    
    cells_all:              Vec< ETotalComplex >,   // all cells
    cell_ids_out:           Vec< Vec< usize > >           // cells not yet assigned a birth,
    
    cell_ids_pos_crit:      HashSet< (Fil, usize, usize) >,   // unmatched positive critical cells, grouped by (birth_time, dimension)
    cell_ids_pos_degn:      Vec< usize >,           // unmatched positive degenerate cells
                                                // NB: doesn't need to be hash; we know birth time
    bars_all_inf:           Vec< Fil >,
    bars_all_fin:           Vec< (Fil, Fil) >,
    
    bar_ids_dun_fin:        Vec< usize >,           // nonempty bars with all endpoints accounted for (finite)
    bar_ids_dun_inf:        Vec< usize >,           // nonempty bars with all endpoints accounted for (infinite)
    
    bar_ids_und_fin:        Vec< usize >,           // nonempty bars with 1 endpoint left to account for 

    bar_ids_now_inf:        Vec< usize >,           // bars still to match in this batch
    bar_ids_now_fin:        Vec< usize >,           // bars still to match in this batch

}




//  ---------------------------------------------------------------------------  
//  FUNCTIONS 
//  ---------------------------------------------------------------------------  


fn format_result( node: Node ) -> EResults
{
    
}

fn enumerate_compatible_filtrations( boundary, cells ) {
   
    //  PRECOMPUTE
    //  compute: ranks of boundary operators (concomitantly, betti numbers)
    //  compute: bars_degn_quota[ dim ] = rank(boundary_(dim+1)) - #(finte bars of dimension dim)

    //  FEASIBILITY CHECK
    //  check: #(inf bars of dimension dim) == betti_dim( total space )
    //  check: #(fin bars of dimension dim) <= rank( total space )

    //  INITIALIZE PARTIAL DATA
    //  RUN THE EXPLORE ALGORITHM
}

fn explore( node: &mut Node, results: &mut Vec< Vec< ETotalComplex> > )
{
   
    //  PUSH LEAF NODE TO RESULTS

    if  node.lev_set_sizes.size_last()  ==  0
        node.lev_set_sizes.total()      ==  node.all_cells.len() &&
        node.bar_ids_dun_fin.len()      ==  node.bars_all_fin.len() &&
        node.bar_ids_dun_inf.len()      ==  node.bars_all_inf.len() &&
        {
        results.push( node.cells_all.clone() ); 
    }

    //  TERMINATE CONSTRUCTION OF PRESENT LEVEL SET; INITIALIZE NEW LEVEL SET 
    //  ---------------------------------------------------------------------
     
    else if node.lev_set_sizes.size_last() != 0     &&      // current class is nonempty 
            node.cell_ids_pos_degn.is_empty()       &&      // C-rule (degenerate) there are no unmatched "degenerate" positive cells
            (                                           // C-rule (critical) 
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
    
        let mut child                   =   node.clone();

        // APPEND NEW COUNTER FOR LEVEL SET SIZE

        child.lev_set_sizes.postpend_empty();        // post-pend new entry (= 0) to vector of level set sizes

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

        // UPDATE FUNCTION f
        
        // note that we must ask whether the node level set is critical, since we have already
        // declared the child level set to be degenerate
        if      node.lev_set_is_crit  { child.f.push_crit() }
        else    { child.f.push_degn() }

        // RUN EXPLORE ON CHILD
        
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
    //  ---------------------------------------------------------------------
    else {

        //  ADD BIRTH CELL FOR AN *INFINITE* BAR 
        //  ------------------------------------
        
        for (bar_id_count, bar_id) in node.bar_ids_now_inf.iter().enumerate() { 

            // clone info about bar to be added
            let mut bar_new     =   node
                                        .bars_all_inf
                                        [bar_id]
                                        .clone();  

            // loop over all cycles of appropriate dimension 
            for (pos_id_out_count, pos_id_out) in node.cell_ids_out
                                                    [ bar_new.dim  ]
                                                    .iter()
                                                    .enumerate()
                .filter(    |x| 
                            node.boundary[x.1]
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

                // mark child's level set as critical (posisbly redundantly)
                child.lev_set_is_crit   =   true;

                // move bar and positive cell
                let _    =   child.cell_ids_out.swap_remove(pos_id_out_count);      // remove cell from "out list"
                let _    =   child.bar_ids_now_inf.swap_remove(bar_id_count);       // remove bar from "unadded list"
                child.bar_ids_dun_inf.push( bar_id );                               // add bar to "done" list (consumes variable)
               
                // record the birth ordinal of the added cell in the central registry
                child
                    .cells_all
                    [ pos_id_out ]
                    .birth
                                =   end_val( child.lev_set_sizes );
               
                // update number of cells born
                child.lev_set_sizes.grow_last_set();

                // RUN EXPLORE ON CHILD
                
                results                 =   explore( child, results );
            }
        } 
        
        //  ADD BIRTH CELL FOR A *FINITE* BAR 
        //  ------------------------------------

        for (bar_id_count, bar_id) in node.bar_ids_now_fin_brn {

            // update: 
            //  lev_set_sizes 
            //  cells_all
            //  cell_ids_out 
            //  bar_ids_dun 
            //  bar_ids_now

            // clone info about bar to be added
            let mut bar_new     =   node.bars_all_fin
                                        [bar_id]
                                        .clone();  

            // loop over all cycles of appropriate dimension 
            for (pos_id_out_count, pos_id_out) in node.cell_ids_out
                                                    [ bar_new.dim ]
                                                    .iter()
                                                    .enumerate()
                                                    .filter(    |x| 
                                                                node.boundary[x.1]
                                                                    .is_empty()
                                                    )  
            {   
                
                // create child node
                let mut child   =   node.clone();                               // this must come before swap_remove

                // mark child's level set as critical (posisbly redundantly)
                child.lev_set_is_crit   =   true;
                
                // move bar_id
                let _    =   child.bar_ids_now_inf.swap_remove(bar_id_count);       // remove bar from list of bars that have unaccounted endpoints here
                child.bar_ids_dun_inf.push( bar_id );                               // add bar to "done" list (consumes variable)

                // move pos_id_out
                let _    =   child.cell_ids_out.swap_remove(pos_id_out_count);      // remove cell from "out list"
                
                child.cell_ids_pos_crit.insert(
                    pow_id_out 
                    // ( 
                    // bar_new.birth.clone(), 
                    // bar_new.dim.clone(),
                    // pos_id_out
                    // )
                )
               
                // record the birth ordinal of the added cell in the central registry
                child
                    .cells_all
                    [ pos_id_out ]
                    .birth
                                =   child.lev_set_sizes.total();
               
                // update number of cells born
                child.lev_set_sizes.grow_last_set(); 

                // RUN EXPLORE ON CHILD
                
                results                 =   explore( child, results );

            }
        } 
        
        //  ADD DEATH CELL FOR A *FINITE* BAR 
        //  ------------------------------------

        for (bar_id_count, bar_id) in node.bar_ids_now_fin_die {

            // clone info about bar to be added
            let mut bar_new     =   node.bars_all_fin
                                        [bar_id]
                                        .clone();  

            // loop over all dimension (n+1) chains with nonzero boundary
            for (neg_id_count, neg_id) in   node.cell_ids_out
                                                [ bar_new.dim + 1 ]
                                                .iter()
                                                .enumerate()
                                                .filter( |(i, x)| !node.boundary[x].is_empty() )
            {

                // find the lowest nonzero entry in this column
                let mut low_id      =    node.boundary
                                                [ neg_id ]  
                                                .iter()
                                                .map( |x| x.0 )
                                                .max()
                                                .unwrap()   // we are guaranteed to get well-defined max, since we filtered out empty columns
                                                .clone();

                //   IF  every cell in the support of this column has been assigned a birth ordinal 
                //  AND  low_id corresponds to a positive cell with the correct birth time
                // THEN  add a positive critical cell
                if  low_id < node.all_cells.len()   // this first condition helps to short circuit / avoid expensive look-ups in the hash set
                    &&  
                    node.cells_all[ low_id ].birth == Some(bar_new.birth)
                    &&
                    node.cell_ids_pos_crit
                        .contains( low_id )
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
                    let _           =   child.cell_ids_out.swap_remove(neg_id_count);           // remove negative cell from "out list"
                    child.cell_ids_pos_crit.remove( low_id );                                   // remove positive cell from set of unmatched positive cells

                    // move bar
                    let _           =   child.bar_ids_now_fin_die.swap_remove(bar_id_count);    // remove bar from "unadded list"
                    child.bar_ids_dun_fin.push( bar_id.clone() );                               // add bar to "done" list (consumes variable)

                   
                    // record the birth ordinal of the added cell in the central registry
                    child
                        .cells_all
                        [ neg_id ]
                        .birth
                                    =   child.lev_set_sizes.total();
                   
                    // update number of cells born
                    child.lev_set_sizes.grow_last_set(); 

                    // record the pairing in the central registry
                    child.cells_all
                        [ low_id ]
                        .bounding_cell_id
                                    =   neg_id.clone();


                    // UPDATE THE BOUNDARY MATRIX
                    
                    // clear bottom entries from columns
                    //      !!! must insert code here !!! 

                    // deallocate the negative column (it's no longer needed)
                    child.boundary
                        [neg_id]
                        .clear()
                        .shrink_to_fit()

                    // RUN EXPLORE ON CHILD
                    
                    results                 =   explore( child, results );
                }
            }
        }

        //  ADD BIRTH CELL FOR A *DEGENERATE* BAR 
        //  -------------------------------------

        // loop over all dimensions
        for dim in 0..node.cell_ids_out.len() {

            // short-circuit if we have already met the quota for degenerate bars in this dimension
            if node.bars_degn_quota[ dim ] == 0 { continue }

            // loop over all cycles of given dimension 
            for (pos_id_out_count, pos_id_out) in node.cell_ids_out
                                                    [ dim ]
                                                    .iter()
                                                    .enumerate()
                                                    .filter(    |x| 
                                                                node.boundary[x.1]
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
                let _    =   child.cell_ids_out.swap_remove(pos_id_out_count);      // remove cell from "out list"
                child.cell_ids_pos_degn.push( pos_id_out );                         // add cell to list of positive degenerate cells
               
                // record the birth ordinal of the added cell in the central registry
                child
                    .cells_all
                    [ pos_id_out ]
                    .birth
                                =   child.lev_set_sizes.total();
               
                // update counters
                child.lev_set_sizes.grow_last_set();    // number of cells born
                child.bar_degn_quota[dim] -= 1;         // drop the quota for degenerate bars by 1

                // RUN EXPLORE ON CHILD
                
                results         =   explore( child, results );
            }
        }
        
        //  ADD DEATH CELL FOR A *DEGENERATE* BAR 
        //  -------------------------------------

        // loop over all dimensions
        for dim in 1..node.cell_ids_out.len() {

            // short-circuit if there exists no positive cell of appropriate dimension
            if node.cell_ids_pos[ dim-1 ].is_empty() { continue }

            // loop over chains of given dimension with nonzero boundary
            for (neg_id_count, neg_id) in   node.cell_ids_out
                                                [ dim ]
                                                .iter()
                                                .enumerate()
                                                .filter( |(i, x)| !node.boundary[x].is_empty() )
            {

                // find the lowest nonzero entry in this column
                let mut low_id      =   node.boundary
                                            [ neg_id ]  
                                            .iter()
                                            .map( |x| x.0 )
                                            .max()
                                            .unwrap()   // we are guaranteed to get well-defined max, since we filtered out empty columns
                                            .clone();

                //   IF  every cell in the support of this column has been assigned a birth ordinal 
                //  AND  low_id corresponds to a positive degenerate cell
                // THEN  add a negative degenerate cell
                if  low_id < node.all_cells.len() &&  // this first condition helps to short circuit / avoid expensive look-ups in the hash set
                    node.cell_ids_pos_degn.contains( low_id ) 
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
                    let _           =   child.cell_ids_out.swap_remove(neg_id_count);           // remove negative cell from "out list"
                    child.cell_ids_pos_degn.remove( low_id );                                   // remove positive cell from set of unmatched positive cells

                    // record the birth ordinal of the added cell in the central registry
                    child
                        .cells_all
                        [ neg_id ]
                        .birth
                                    =   child.lev_set_sizes.total();
                   
                    // update number of cells born
                    child.lev_set_sizes.grow_last_set(); 

                    // record the pairing in the central registry
                    child.cells_all
                        [ low_id ]
                        .bounding_cell_id
                                    =   neg_id.clone();

                    // UPDATE THE BOUNDARY MATRIX
                    
                    // clear bottom entries from columns
                    //      !!! must insert code here !!! 

                    // deallocate the negative column (it's no longer needed)
                    child.boundary
                        [neg_id]
                        .clear()
                        .shrink_to_fit()

                    // RUN EXPLORE ON CHILD
                    
                    results                 =   explore( child, results );
                }
            }
        }
    }

    return results
}


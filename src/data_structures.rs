


type Cell = usize;
type Val = Ratio< i64 >;
type Fil = f64;




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
//     endpoint_now         usize,
//     
//     cells_all:          Vec< ETotalComplex >,   // all cells
//     cells_out:          Vec< usize >            // cells not yet assigned a birth,
//     cells_pos_crit:     Vec< EPositiveCell >,   // unmatched positive critical cells
//     cells_pos_degn:     Vec< EPositiveCell >,   // unmatched positive degenerate cells
// 
//     bars_dun:           Vec< usize >,           // bars with all endpoints accounted for
//     bars_und:           Vec< usize >,           // bars without all endpoints accounted for 
//     bars_all:           Vec< ETotalBarcode>,    // all bars
//     bars_now:           Vec< usize >,           // bars still to match in this batch
// 
//     class_id:           usize,
//     class_size:         usize,
//     class_crit:         bool 
// }


// struct Node
// {
//     boundary:           Vec< Vec< ( Cell, Val ) > >,
//     num_cells_born      usize,
//     endpoint_now        usize,
// 
//     cells_all:          Vec< ETotalComplex >,   // all cells
//     cells_out:          Vec< usize >            // cells not yet assigned a birth,
//   
//     bars_crit_all:          Vec< BarPair >,     // all bars in desired barcode; 
//     bars_crit_needbirth:    Vec< usize >,       // bars not yet assigned a birth cell
//     bars_crit_needdeath:    Vec< usize >,       // bars assigned a birth cell but no death cell
// 
//     bars_degn_all:          Vec< BarPair >,     // all degenerate bars 
//     bars_degn_needdeath:    Vec< usize >,       // bars assigned a birth cell but no death cell
// 
// 
//     bars_fut_crit:      Vec< usize >,           // bars without all endpoints accounted for 
//     bars_fut_degn:      Vec< usize >,           // bars without all endpoints accounted for 
//     
//     bars_now_crit:      Vec< usize >,           // bars still to match in this batch
//     bars_now_degn:      Vec< usize >,           // bars still to match in this batch
// 
//     class_id:           usize,
//     class_size:         usize,
//     class_crit:         bool 
// }


struct Node
{
    boundary:           Vec< Vec< ( Cell, Val ) > >,
    num_cells_born      usize,
    endpoint_now        usize,
    
    cells_all:          Vec< ETotalComplex >,   // all cells
    cells_out:          Vec< usize >            // cells not yet assigned a birth,
    cells_pos_crit:     Hash< FParam, Vec< usize > >,   // unmatched positive critical cells
    cells_pos_degn:     Vec< usize >,           // unmatched positive degenerate cells
                                                // NB: doesn't need to be hash; we know birth time

    bars_all_inf:       Vec< FParam >,
    bars_all_fin:       Vec< (Fparam, Fparam) >,
    
    bars_und:           Vec< usize >,           // nonempty bars with 1 endpoint left to account for 
    
    bars_dun_fin:       Vec< usize >,           // nonempty bars with all endpoints accounted for (finite)
    bars_dun_inf:       Vec< usize >,           // nonempty bars with all endpoints accounted for (infinite)

    bars_now_inf:       Vec< usize >,           // bars still to match in this batch
    bars_now_fin:       Vec< usize >,           // bars still to match in this batch

    class_id:           usize,
    class_size:         usize,
    class_crit:         bool 
}




//  ---------------------------------------------------------------------------  
//  FUNCTIONS 
//  ---------------------------------------------------------------------------  


fn format_result( node: Node ) -> EResults
{
    
}

fn explore( node: &mut Node, results: &mut Vec< Vec< ETotalComplex> > )
{
    
    //  PUSH TO RESULTS
    if  node.num_cells_born == node.all_cells.len() &&
        node.bars_und.is_empty() &&
        node.class_size == 0
        {
        results.push( node.cells_all.clone() ); 
    }

    //  BEGIN NEW CLASS
    else if node.class_size != 0 &&         // current class is nonempty 
            node.bars_now.is_empty() &      // I_cur is empty
            node.cells_pos_degn.is_empty(); // \bar C_cur contains no unmatched elements
        {
        // create child node
        let mut child = node.clone();
        child.class_id          +=  1;      // increase class #
        child.class_size        =   0;      // set size to 0 
        if child.class_crit 
           child.endpoint_now    +=  1;      // if parent was critical, increase bc endpoint index 
        }
        
        // fork 1 for child: new critical class
        if ! child.bars_und.is_empty() {    // if there are bars yet to match
            child.class_crit    =   true;   // flag that we are starting a critical class
            child.bars_now      =   bars_und    // calculate bars with this endpoint
                                        .iter()
                                        .map(|x| child.bars_all[ x ] )
                                        .filter(|x|     x.birth_time == child.endpoint_now
                                                        ||
                                                        x.death_time == Some(child.endpoint_now)
                                                )
                                        .collect();
            results             =   explore( child.clone(), results );
        }

        // fork 2 for child: new degenerate class
        child.class_crit        =   false;      // flag that we're creating a degenerate class
        child.bars_now.clear();                 // remove all bars we may have added to the  "docket"
        results                 =   explore( child, results ); 
    }

    //  ADD TO CLASS
    else {
        
        for bar in node.bars_now {
            
            //  ADD A BAR_DEATH CELL
            if bar.birth_time  ==  node.endpoint_now{
                for out_ind in 0..node.cells_out.len()
                    .filter(    |x| 
                                node.boundary[  node.cells_out[ind] ].is_empty() &&
                                node.all_cells[ node.cells_out[ind] ].dim == bar.dim     )

                    {
                    let mut child   =   node.clone(); // this must come before swap_remove
                    let cell    =   child.cells_out.swap_remove(out_ind);
                    child.cells_pos_crit.push( cell.clone() ); // DO WE REALLY WANT TO DO THIS?
                    // update: 
                    //  num_cells_born
                    //  cells_all
                    //  cells_out (already done)
                    //  cells_pos_crit (already done -- but do we want this?)
                    //  bars_dun (if bar never dies)
                    //  bars_und (if bar never dies)
                    //  bars_all
                    //  bars_now
                    //  class_size
                }
            } 

            //  ADD A BAR_BITH CELL
            else if bar.death_time == Some( node.endpoint_now ) {
                
                for cell_birth in node.cell_pos_crit { // ASK JACOB IF THIS FOR-LOOP IS OK
                    if  {
                            
                    }
                }
                // STUFF
            }


        }

        for pos_degen in node.cells_pos_degn {
            // STUFF
        }

        for cell in node.cells_out
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



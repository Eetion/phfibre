


type Cell = usize;
type Val = Ratio< i64 >;
type Fil = f64;




//  ---------------------------------------------------------------------------  
//  ENTRIES 
//  ---------------------------------------------------------------------------  


struct EIntervalsCur{
    birth_index:    Fil,
    death_index:    Fil,
    birth_cell:     Option< Cell >,
    death_cell:     Option< Cell >,
    dim:            usize
}


struct ETotalComplex{
    birth_order:    usize,
    birth_class:    usize,
    birth_bc_ind:   Option< usize >,
    dim:            usize
}


struct ETotalBarcode{
    birth_index:    usize,
    death_index:    Option< usize >,
    dim:            usize,
}


struct EPositiveCell{
    cell:           Cell,
    death_cells:    Vec< Cell >,
    interval:       Option< usize >, // uuid of the corresponding interval 
}


struct EResults{
    birth_class:    usize,
    birth_order:    usize,
    birth_itv:      Option< usize >,
    death_itv:      Option< usize >
}


//  ---------------------------------------------------------------------------  
//  CONTAINERS
//  ---------------------------------------------------------------------------  



struct Node
{
    global_boundary:    Vec< Vec< ( Cell, Val ) > >,
    
    cells_all:          Vec< ETotalComplex >,   // all cells
    cells_out:          Vec< Cell >             // cells not yet assigned a birth,
    cells_pos_crit:     Vec< EPositiveCell >,   // unmatched positive critical cells
    cells_pos_degn:     Vec< EPositiveCell >,   // unmatched positive degenerate cells

    bars_all:           Vec< ETotalBarcode>,    // all bars
    bars_out:           Vec< usize >            // bars never placed in bars_now
    bars_now:           Vec< usize >            // bars to match in this batch

    cur_class_id:       usize,
    cur_num_cells:      usize,
    cur_bc_endpoint:    usize,
}



//  ---------------------------------------------------------------------------  
//  FUNCTIONS 
//  ---------------------------------------------------------------------------  


fn explore( node: &mut Node, results: &mut Vec< EResults > )
{
    if  node.cur_cells_num == node.global_cells.len() &&
        node.bars_out.is_empty() &&

}



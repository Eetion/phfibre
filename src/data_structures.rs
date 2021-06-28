


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
    birth_bar:      Option< usize >,
    death_bar:      Option< usize >,
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
    boundary:           Vec< Vec< ( Cell, Val ) > >,
    total_cells_born    usize,
    bc_endpoint         usize,
    
    cells_all:          Vec< ETotalComplex >,   // all cells
    cells_out:          Vec< usize >            // cells not yet assigned a birth,
    cells_pos_crit:     Vec< EPositiveCell >,   // unmatched positive critical cells
    cells_pos_degn:     Vec< EPositiveCell >,   // unmatched positive degenerate cells

    bars_dun:           Vec< usize >,           // bars with all endpoints accounted for
    bars_und:           Vec< usize >,           // bars without all endpoints accounted for 
    bars_all:           Vec< ETotalBarcode>,    // all bars
    bars_now:           Vec< usize >,           // bars still to match in this batch

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
    if  node.total_cells_born == node.all_cells.len() &&
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
           child.bc_endpoint    +=  1;      // if parent was critical, increase bc endpoint index 
        }
        
        // fork 1 for child: new critical class
        if ! child.bars_und.is_empty() {    // if there are bars yet to match
            child.class_crit    =   true;   // flag that we are starting a critical class
            child.bars_now      =   bars_und    // calculate bars with this endpoint
                                        .iter()
                                        .map(|x| child.bars_all[ x ] )
                                        .filter(|x|     x.birth_index == child.bc_endpoint
                                                        ||
                                                        x.death_index == Some(child.bc_endpoint)
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
            
            //  ADD A CRITICAL BIRTH CELL
            if bar.birth_index  ==  node.bc_endpoint{
                 
            } 

            //  ADD A CRITICAL DEATH CELL
            else if bar.death_index == Some( node.bc_endpoint ) {
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



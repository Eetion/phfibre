
use crate::utilities::*;
use num::rational::Ratio;
use std::collections::{HashSet, HashMap};
use std::iter::FromIterator;
use ordered_float::OrderedFloat;

type Cell = usize;
type Coeff = Ratio< i64 >;
type Fil = usize;
type FilRaw = OrderedFloat<f64>; // reason for this choice: f64 does not implement the hash trait (cricital for reparametrizing)




//  ---------------------------------------------------------------------------  
//  BARCODES
//  ---------------------------------------------------------------------------  


#[derive(Clone, Debug)]
struct FiniteBar { 
    dim:        usize, 
    birth:      Fil, 
    death:      Fil 
}

#[derive(Clone, Debug)]
struct InfiniteBar { 
    dim:        usize, 
    birth:      Fil 
}

#[derive(Clone, Debug)]
struct Barcode {
    inf:        Vec< InfiniteBar >,     // infinite bars (birth ordinals)
    fin:        Vec< FiniteBar >,       // finite bars (birth/death ordinals)
    ordinal:    OrdinalData< FilRaw >,            // struct converting endpoints to ordinals and vice versa
}

impl Barcode{
    /// The maximum dimension of any bar in the barcode (None if the barcode is empty)
    fn top_dim( &self ) -> Option< usize > { 
        self.fin.iter().map(|&x| x.dim.clone())
            .chain(
                self.inf.iter().map(|&x| x.dim.clone())
            )
            .max()    
    }
    /// The maximum filtration (ordinal) of any bar in the barcode (None if the barcode is empty)
    fn top_fil( &self ) -> Option< Fil > { 
        self.fin.iter().map(|&x| x.death.clone())
            .chain(
                self.inf.iter().map(|&x| x.birth.clone())
            )
            .max()    
    }    
    /// Create a new barcode from raw parts
    fn new( 
        inf_brn:    Vec< FilRaw >,
        inf_dim:    Vec< usize >,
        fin_brn:    Vec< FilRaw >,
        fin_die:    Vec< FilRaw >,
        fin_dim:    Vec< usize >
        ) 
        -> 
        Barcode 
    {
        let endpoints_raw_unordered   =   Vec::new();
        endpoints_raw_unordered.append( &mut inf_brn.clone() );
        endpoints_raw_unordered.append( &mut fin_brn.clone() );
        endpoints_raw_unordered.append( &mut fin_die.clone() );

        let raw_endpoint_ordinal_data   =   ordinate( &endpoints_raw_unordered );

        let barcode     =   Barcode {
                                inf:        Vec::new(),     // infinite bars (birth ordinals)
                                fin:        Vec::new(),       // finite bars (birth/death ordinals)
                                ordinal:    raw_endpoint_ordinal_data
                            };
                            
        for bar_count in 0..fin_brn.len() {
            let bar     =   FiniteBar { 
                                dim:        fin_dim[ bar_count ],
                                birth:      fin_brn[ bar_count ], 
                                death:      fin_die[ bar_count ] 
                            };        
            barcode.fin.push( bar );
        }
        
        for bar_count in 0..inf_brn.len() {
            let bar     =   InfiniteBar { 
                                dim:        inf_dim[ bar_count ],
                                birth:      inf_brn[ bar_count ]
                            };        
            barcode.inf.push( bar );
        }        

        return barcode
    }
}

/// Encodes a map sending (ordinal) endpoints to sets of bar id numbers.
#[derive(Clone, Debug)]
struct MapEndpoint2BarIDs {
    inf:        Vec< Vec< usize > >,
    fin_brn:    Vec< Vec< usize > >,
    fin_die:    Vec< Vec< usize > >
}

/// Returns an object that maps an (ordinal) endpoints back to set of bar id's.
fn  make_endpoint_to_barids_map(
    barcode: Barcode
    ) ->
    MapEndpoint2BarIDs

{
    let mut endpoint__barids = MapEndpoint2BarIDs {
        inf:        Vec::with_capacity(0),
        fin_brn:    Vec::with_capacity(0),
        fin_die:    Vec::with_capacity(0),
    };    

    if let Some( top_fil )  =  barcode.top_fil() {
        for _ in 0..top_fil {
            endpoint__barids.inf.push(      Vec::new()  );
            endpoint__barids.fin_brn.push(  Vec::new()  );
            endpoint__barids.fin_die.push(  Vec::new()  );                        
        }
    } 

    // fill bins with infinite bars
    for (bar_count, bar) in barcode.inf.iter().enumerate() {
        endpoint__barids.inf[ bar.birth ].push( bar_count );
    }

    // fill bins with finite bars
    for (bar_count, bar) in barcode.fin.iter().enumerate() {
        endpoint__barids.fin_brn[ bar.birth ].push( bar_count );
        endpoint__barids.fin_die[ bar.death ].push( bar_count );        
    }

    return endpoint__barids
}


//  ---------------------------------------------------------------------------  
//  LEVEL SET SIZES
//  ---------------------------------------------------------------------------  


/// A struct representing the sizes of a sequence of level sets.
#[derive(Clone, Debug)]
struct LevelSetSizes{ pointers: Vec< usize > }

impl LevelSetSizes{
    
    /// Size of the kth level set.
    fn size( &self, set_index: usize ) -> usize {  
        if set_index == 0 { self.pointers[0] }
        else { self.pointers[ set_index ] - self.pointers[ set_index - 1 ] }
    }

    /// Size of the last level set.
    fn size_last( &self ) -> usize { self.size( self.pointers.len() - 1) }
    
    /// Size of all level sets combined.
    fn total( &self ) -> usize {
        if self.pointers.is_empty() { 0 }
        else { end_val( self.pointers ).clone() }
    }
    
    /// Add a size value for a new (empty) level set at the end.
    fn postpend_empty_set( &self ) { 
        if self.pointers.is_empty() { self.pointers.push(0) }
        else { self.pointers.push( end_val(self.pointers) ) }
    }

    /// Add one to the size of the last level set.
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
//  SEARCH TREE
//  ---------------------------------------------------------------------------  


/// Represents an "entry" for a sincle cell within a central repository of 
/// information about individual cells.
#[derive(Clone, Debug)]
struct CellEntry{
    birth_ordinal:  usize,
    birth_class:    usize,
    birth_bar:      Option< usize >,
    death_bar:      Option< usize >,
    dim:            usize
}

/// Represents a node of the search tree.
#[derive(Clone, Debug)]
struct Node<'a >
{
    boundary:               Vec< Vec< ( Cell, Coeff ) > >,  // boundary matrix reprsented as a vector of vectors
    boundary_buffer:        Vec< ( Cell, Coeff) >,          // a "holding space" for matrix entries, when these need to be moved around
    bc_endpoint_now:        usize,                          // the (ordinal) barcode endpoint of concern for this (or some future) level set
    
    bars_degn_quota:        Vec< usize >,                   // number of degenerate cells anticipated in each dimension

    lev_set_sizes:          LevelSetSizes,                  // Kth value = # cells in Kth level set

    lev_set_values:         Vec< Fil >,                     // the function phi
    lev_set_is_crit:        bool,                           // true iff current level set contains a "critical cell"
    
    cells_all:              Vec< CellEntry >,               // all cells
    cell_ids_out:           Vec< Vec< usize > >,            // cells not yet assigned a birth,
    
    cell_ids_pos_crit:      HashSet< (Fil, usize, usize) >, // unmatched positive critical cells, grouped by (birth_time, dimension)
    cell_ids_pos_degn:      Vec< usize >,                   // unmatched positive degenerate cells
                                                            // NB: doesn't need to be hash; we know birth time
    bars_all_inf:           Vec< Fil >,                     
    bars_all_fin:           Vec< (Fil, Fil) >,
    
    bar_ids_dun_fin:        Vec< usize >,                   // nonempty bars with all endpoints accounted for (finite)
    bar_ids_dun_inf:        Vec< usize >,                   // nonempty bars with all endpoints accounted for (infinite)

    bar_ids_now_inf:        &'a Vec<Vec< usize >>,          // bars still to match in this batch
    bar_ids_now_fin_brn:    &'a Vec<Vec< usize >>,          // bars to be born with this level set
    bar_ids_now_fin_die:    &'a Vec<Vec< usize >>,          // bars to be bounded with this level set    

}


//  ---------------------------------------------------------------------------  
//  ENUMERATE POLYHEDRA
//  ---------------------------------------------------------------------------  


fn  enumerate_compatible_filtrations( 
        boundary: Vec<Vec<Coeff>>, 
        cell_registry: Vec<Cell>,
        barcode_inf_births: Vec< Fil >,
        barcode_fin_births: Vec< Fil >,
        barcode_fin_deaths: Vec< Fil >                
     ) 
{

    //  PRECOMPUTE
    //  compute: bars that stop and start at each endpoint
   
    //  FEASIBILITY CHECK
    //  compute: ranks of boundary operators (concomitantly, betti numbers)
    //  compute: bars_degn_quota[ dim ] = rank(boundary_(dim+1)) - #(finte bars of dimension dim)    
    //  check: #(inf bars of dimension dim) == betti_dim( total space )
    //  check: #(fin bars of dimension dim) <= rank( total space )

    //  INITIALIZE PARTIAL DATA

    //  EXPLORE THE TREE
    //  remove duplicate polyhedra as we go

}


//  ---------------------------------------------------------------------------  
//  EXPLORE TREE
//  ---------------------------------------------------------------------------  


fn explore( node: &mut Node, results: &mut Vec< Vec< CellEntry> > )
{
   
    //  PUSH LEAF NODE TO RESULTS

    if  node.lev_set_sizes.size_last()  ==  0 &&
        node.lev_set_sizes.total()      ==  node.cells_all.len() &&
        node.bar_ids_dun_fin.len()      ==  node.bars_all_fin.len() &&
        node.bar_ids_dun_inf.len()      ==  node.bars_all_inf.len()  
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
                node.bar_ids_now_fin_die.is_empty()    
                )
            )
        {

        // CREATE CHILD NODE    
    
        let mut child                   =   node.clone();

        // APPEND NEW COUNTER FOR LEVEL SET SIZE

        child.lev_set_sizes.postpend_empty_set();        // post-pend new entry (= 0) to vector of level set sizes

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
            for bar in node.bars_all_inf.iter().filter(|&&x| x == child.bc_endpoint_now) {
                child.bar_ids_now_inf.push( bar.clone() );
            }
            for bar in node.bars_all_fin.iter().filter(|&&x| x.birth == child.bc_endpoint_now) {
                child.bar_ids_now_brn.push( bar.clone() );
            }
            for bar in node.bars_all_fin..iter().filter(|&&x| x.death == child.bc_endpoint_now) {
                child.bar_ids_now_die.push( bar.clone() );
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
                    pos_id_out 
                    // ( 
                    // bar_new.birth.clone(), 
                    // bar_new.dim.clone(),
                    // pos_id_out
                    // )
                );
               
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
                        .shrink_to_fit();

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
                        .shrink_to_fit();

                    // RUN EXPLORE ON CHILD
                    
                    results                 =   explore( child, results );
                }
            }
        }
    }

    return results
}


//  ---------------------------------------------------------------------------  
//  CHECK INTERSECTION OF POLYTOPES
//  ---------------------------------------------------------------------------  


struct Polytope {
    data_l__fmin:    Vec< Fil >,    // function (level set ordinal) -> min possible filtration value
    data_c__l:       Vec< usize >,     // function (cell id) -> level set ordinal
}

impl Polytope {

    /// Number of cells in the underlying chain complex.
    fn num_cells( &self ) -> usize { self.data_c__l.len() }

    /// Number of level sets of the filter function.
    fn num_lev_sets( &self ) -> usize { self.data_l__fmin.len() }    

    /// Rightmost finite filtration value
    fn max_filtration_value( &self ) -> Option<Fil> { self.data_l__fmin.iter().cloned().last() }

    fn  c__l( &self, cell_id: usize ) -> Option<usize> { 
        if cell_id >= self.num_cells() { return None } 
        else {
            return Some( self.data_c__l[ cell_id ].clone()  )
        } 
    }

    fn  l__is_critical( &self, level_set_ordinal: usize ) -> Option<bool> { 
        if level_set_ordinal >= self.num_lev_sets() { return None }
        else if level_set_ordinal == 0 { return Some( true ) }
        else if self.data_l__fmin[ level_set_ordinal ] 
                == 
                self.data_l__fmin[ level_set_ordinal -1 ] 
        {
            return Some( false )
        }
        else { return Some( true )}
    }

    fn  l__fmin( &self, level_set_ordinal: usize ) -> Option<Fil> { 
        if level_set_ordinal >= self.num_lev_sets() { return None }
        else {
            return Some ( self.data_l__fmin[ level_set_ordinal ].clone() )            
        }
    }

    // NB!!!!  If we assume that barcode endpoints are 0,1,..,N, then we can simplify
    // computation of max filtration values considerably
    fn  l__fmax( &self, level_set_ordinal: usize ) -> Option<Fil> { 

        // Posit: the ordinal is within legal range
        if let Some( is_critical ) = self.l__is_critical( level_set_ordinal.clone() )
        {
            // Posit: the level set ordinal IS critical.  Then max value equals min value.
            if is_critical { return self.l__fmin( level_set_ordinal ) }
            // Posit: the level set ordinal IS NOT not critical                
            else {
                let fmin                    =   self.l__fmin( 
                                                        level_set_ordinal.clone() 
                                                    )
                                                    .unwrap();
                if let Some( next_fil_val ) =   & self.data_l__fmin[
                                                    level_set_ordinal.clone()..self.num_lev_sets()
                                                    ]
                                                        .iter()
                                                        .map(|x| self.l__fmax(x.clone()).unwrap() )
                                                        .filter(|&x| x > fmin)
                                                        .next()
                {
                    return Some( next_fil_val.clone() )
                } else {
                    return Some( self.num_cells() )
                }
                                                    
                // // Posit: this is the last legal level set ordinal
                // // This implies that  not critical, and lies in open interval (0, max_legal_ordinal]
                // if level_set_ordinal == self.num_lev_sets()-1 { 
                //     return Some (self.max_filtration_value().unwrap() + 1 ) 
                // } 
                // // Posit: this is NOT the last legal ordinal
                // // This implies that there is one above:
                // else                 
                // {
                //     // SEARCH FORWARD FOR FIRST DISTINCT VALUE
                //     return Some( self.data_l__fmin[ level_set_ordinal +1 ])
                // }   
            }
        } 
        else { return None }
    }

}


struct PolytopeVertexSolver<'a> {
    poly:           &'a Polytope,       // polytope to work with
    c__fmin:     Vec< Fil >,         // function (cell id) -> min possible filtration value    
    c__fmax:     Vec< Fil >,         // function (cell id) -> max possible filtration value
}

fn poly_2_polysolve<'a> ( poly: &'a Polytope ) -> PolytopeVertexSolver<'a>{

    let num_cells       =   poly.c__l.len();
    let num_levels      =   poly.l__fmin.len();    

    let mut c__fmin  =   Vec::with_capacity( num_cells );
    let mut c__fmax  =   Vec::with_capacity( num_cells );

    // DEGENERATE CASES: 0 CELLS OR EXACTLY 1 LEVEL SET
    if num_cells == 0 {
        return PolytopeVertexSolver{ poly:poly, c__fmin:c__fmin, c__fmax:c__fmax}        
    } else 
    // this represents the case where there are 1 or more cells, but exactly 1 level set
    // (in this case the level set must be critical, on account of the infinite dim-0 bar)
    if num_levels == 1 {
        for cell_id in 0..num_cells {
            c__fmin.push(0);
            c__fmax.push(0);
        }
        return PolytopeVertexSolver{ poly:poly, c__fmin:c__fmin, c__fmax:c__fmax}                
    }

    // REMAINING CASES: 
    //      >= 1 cell           AND
    //      >= 2 level sets 

    // min filtration value for each cell is the min filtration value for its level set
    c__fmin.extend(       
        poly.c__l
            .iter()
            .map(
                |&x| 
                poly.l__fmin[ x ].clone()
            )
    );

    // max filtration value for each cell is either (a) the min filtration value of the
    // succeeding level set, or (b) infinity (realizes as # cells, which is strictly 
    // greater than the index of any cell)

    c__fmax.extend(  
        poly.c__l[0..(num_cells-1)]
            .iter()
            .map( |&i|
                poly.l__fmin[ i + 1 ]
                .clone()
            )     
    );
    

    for cell_id in 0..num_cells {


    }
}


fn check_intersection(
    poly_a:     Polytope,
    poly_b:     Polytope
)

{
    // define solver a
    let mut 
    let mut solver_a    =   PolytopeVertexSolver{   poly: &poly_a, 2 }
}
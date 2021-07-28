




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

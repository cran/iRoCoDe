# Generates number of treatments in a Design
get_treat <- function(design){
  d = design
  blks = nrow(d)
  blk_size = ncol(d)
  trt_list = c()
  for (i in 1:blk_size){
    for (j in 1:blks){
      if (length(trt_list)==0){
        trt_list = append(trt_list, d[j,i])
      }
      else {
        check_t = d[j,i] %in% trt_list
        if(check_t == FALSE){
          trt_list = append(trt_list, d[j,i])
        }
      }
    }
  }
  trt_list = sort(trt_list)
  return(trt_list)
}

# Generates number of replications in a Design
get_rep <- function(design){
  d = design
  blks = nrow(d)
  blk_size = ncol(d)
  trt_list = get_treat(d)
  repli_list = c()
  for (t in (1:length(trt_list))){
    flag_r =0
    for (i in 1:blks){
      for (j in 1:blk_size){
        check_r= trt_list[t] %in% d[i,j]
        if (check_r ==TRUE){
          flag_r= flag_r+1
        }
      }
    }
    repli_list = append(repli_list, flag_r)
  }
  for (i in (1: (length(repli_list)-1))){
    if(repli_list[i]!= repli_list[i+1]){
      return('error')
    }
  }
  rep  = repli_list[1]
  return(rep)
}

# Generates the Incidence Matrix (N_matrix) of a Design
n_matrix <- function(design, treatments, blocks){
  d = design
  v = treatments
  b = blocks
  cols = ncol(d)
  rows = nrow(d)
  n_mat = matrix(,nrow = v, ncol = b)
  for (i in (1:v)){
    pos =c()
    for (k in (1:rows)){
      for (l in (1:cols)){
        if (d[k,l]==i){
          pos = append(pos,k)
          break
        }
      }
    }
    for (j in (1:length(pos))){
      n_mat[i,pos[j]] = 1
    }
  }
  n_mat[is.na(n_mat)]=0
  return(n_mat)
}


#  Generates the A_matrix of a design
a_matrix <- function(design2, n_matrix){
  n_mat = n_matrix
  d2 = design2
  cols = ncol(n_mat)
  rows = nrow(n_mat)
  a_mat = matrix(,nrow = rows, ncol = cols)
  for (i in (1:cols)){
    k = 1
    for (j in (1:rows)){
      if(n_mat[j,i]!=0){
        a_mat[j,i] = d2[i, k]
        k=k+1
      }
    }
  }
  a_mat[is.na(a_mat)]=0
  #print(a_mat)
  return(a_mat)
}

# Rearranges the a_matrix
rearrange_mat <- function(a_matrix, treatments){
  a_mat = a_matrix
  v = treatments
  rows = nrow(a_mat)
  cols = ncol(a_mat)

  ### are colummns unique???
  col_unique = TRUE
  for (i in (1:cols)){
    for (j in (1:rows)){
      if(j != rows){
        if(a_mat[j,i]!=0){
          first = a_mat[j,i]
          for (k in ((j+1):rows)){
            if (a_mat[k,i]!=0){
              second = a_mat[k,i]
              if(first==second){
                col_unique =FALSE
                break
              }
            }
          }
        }
      }
    }
  }
  if (col_unique == TRUE){
    r_unique = row_unique(a_mat, v, rows, cols)
    if(r_unique==TRUE){
      return(a_mat)
    }
    else {
      trt_list = c(1:v)
      flag = 0
      for (i in (1:rows)){
        ### creations of the necessary empty lists
        exist_list = c()
        miss_list = c()
        dup_list = c()
        for (j in (1:cols)){
          if(a_mat[i,j]!=0){
            val= a_mat[i,j] %in% exist_list
            if (val == TRUE){
              dup_list = append(dup_list, a_mat[i,j])
              exist_list = append(exist_list, a_mat[i,j])
            } else {
              exist_list = append(exist_list, a_mat[i,j])
            }
          }
        }
        for (t in (1:v)) {
          e= trt_list[t] %in% exist_list
          if(e==FALSE){
            miss_list= append(miss_list, trt_list[t])
          }
        }
        ### re-arranging columns row-wise
        spcl_col = c()
        flag_2 =0
        for (j in (1:cols)){
          if (a_mat[i,j]!=0){
            if (flag==0){
              first_mat_element = a_mat[i,j]
              flag=1
            } else {
              val_check_1 = a_mat[i,j] %in% dup_list
              if(val_check_1== TRUE){
                flag_2 =0
                curr_val_1 = a_mat[i,j]
                for(x in ((i+1):rows)){
                  if(x<=rows){
                    if(a_mat[x,j]!=0){
                      val_check_2 = a_mat[x,j] %in% miss_list
                      if(val_check_2==TRUE){
                        next_val_1 = a_mat[x,j]
                        a_mat[i,j] = next_val_1
                        a_mat[x,j] = curr_val_1
                        miss_list = miss_list[miss_list !=next_val_1]
                        exist_list = rmv_list(exist_list, curr_val_1)
                        exist_list = append(exist_list, next_val_1)
                        flag_3 =0
                        for (ij in (1:length(exist_list))){
                          if (curr_val_1 == exist_list[ij]){
                            flag_3= flag_3+1
                          }
                        }
                        if (flag_3==1){
                          dup_list = dup_list[dup_list != curr_val_1]
                        }
                        flag_2 = 1
                        break
                      }
                    }
                  }
                }
                if (flag_2!=1){
                  spcl_col = append(spcl_col, j)
                }
              }
            }
          }
          if(length(miss_list) ==0){
            break
          }
        } ## end of second loop (j)
        while(length(miss_list) !=0){
          for (j in (1:cols)){
            if (a_mat[i,j]!=0){
              curr_val_3 = a_mat[i,j]
              for (y in (i+1:rows)){
                if (y<=rows){
                  if(a_mat[y,j]!=0){
                    val_check_4 = a_mat[y,j] %in% miss_list
                    if (val_check_4== TRUE){
                      next_val_3 = a_mat[y,j]
                      a_mat[y,j] = curr_val_3
                      a_mat[i,j] = next_val_3
                      miss_list = miss_list[miss_list !=next_val_3]
                      val_check_5 = curr_val_3 %in% dup_list
                      if (val_check_5 == TRUE){
                        exist_list = rmv_list(exist_list, curr_val_3)
                        exist_list = append(exist_list, next_val_3)
                        flag_3 =0
                        for (ij in (1:length(exist_list))){
                          if (curr_val_1 == exist_list[ij]){
                            flag_3= flag_3+1
                          }
                        }
                        if (flag_3==1){
                          dup_list = dup_list[dup_list != curr_val_1]
                        }
                      }else {
                        miss_list = append(miss_list, curr_val_3)
                        exist_list = rmv_list(exist_list, curr_val_3)
                        exist_list = append(exist_list, next_val_3)
                      }
                      ### check spcl columns now
                      if(length(miss_list) !=0 & length(spcl_col)!=0){
                        for (k in (1: cols)){
                          val_check_3 = k %in% spcl_col
                          if (val_check_3==TRUE){
                            curr_val_2 = a_mat[i,k]
                            for (l in ((i+1): rows)){
                              if (l<=rows){
                                if (a_mat[l,k]!=0){
                                  val_check_7 = a_mat[l,k] %in% miss_list
                                  if(val_check_7 == TRUE){
                                    next_val_2 = a_mat[l,k]
                                    a_mat[l,k] = curr_val_2
                                    a_mat[i,k] = next_val_2
                                    miss_list = miss_list[miss_list !=next_val_2]
                                    exist_list = append(exist_list, next_val_2)
                                    exist_list = rmv_list(exist_list, curr_val_2)
                                    flag_3 =0
                                    for (ij in (1:length(exist_list))){
                                      if (curr_val_1 == exist_list[ij]){
                                        flag_3= flag_3+1
                                      }
                                      if (flag_3<=1){
                                        dup_list = dup_list[dup_list != curr_val_1]
                                      }
                                    }
                                    spcl_col = spcl_col[spcl_col!=k]
                                    break
                                  }
                                }
                              }
                            }
                          }
                          if(length(miss_list) ==0){
                            break
                          }
                        }
                      }
                      break
                    }
                  }
                }
              }
            }
          }
        }
      }# print(a_mat)
    } ## end of 1st for loop (i)
  } else {
    return('error')
  }
  return(a_mat)
}

# Deletes an item from a List
rmv_list <- function(list, val){
  list_1 = list
  for (i in (1: length(list_1))){
    if (list_1[i]==val){
      list_1 = list_1[-i]
      break
    }
  }
  return(list_1)
}

# checks the row uniqueness of matrix
row_unique <- function(a_matrix, treatments, rows, cols){
  a_mat = a_matrix
  v = treatments
  rows = rows
  cols = cols
  con = FALSE
  trt_list = c(1:v)
  for (a in (1:rows)) {
    c = 1
    row_con = c()
    for (b in (1:cols)){
      if(a_mat[a,b]!=0){
        row_con[c] = a_mat[a,b]
        c = c+1
      }
    }
    if(length(trt_list)!= length(row_con)){
      con = FALSE
      return (con)
    } else{
      row_con = sort(row_con)
      trt_list = sort(trt_list)
      for (x in (1:length(trt_list))){
        con = TRUE
        if(row_con[x]!=trt_list[x]){
          con = FALSE
          return(con)
        }
      }
    }
  }
  return(con)
}

# Generates a row-column Design with incomplete Rows and Columns
iRoCoDe <- function(design1, design2){
  # Input designs
  d1 = design1
  d2 = design2

  # treatments
  trt_list1 = get_treat(d1)
  trt_list2 = get_treat(d2)
  v1 = length(trt_list1)
  v2 = length(trt_list2)

  # number of blocks
  b1 = nrow(d1)
  b2 = nrow(d2)

  # replications
  r1 = get_rep(d1)
  r2 = get_rep(d2)
  if (r1 == 'error' || r2 =='error'){
    return('replications are not equal')
  }

  # block size
  k1 = ncol(d1)
  k2 = ncol(d2)

  # cat(paste0('\nParameters of First Design:\n',
  #              'treatments:',v1,',\n',
  #              'blocks: ',b1,',\n',
  #              'replication: ', r1,',\n',
  #              'block size: ', k1))
  # cat(paste0('\nParameters of Second Design:\n',
  #              'treatments: ',v2,',\n',
  #              'blocks: ',b2,',\n',
  #              'replication: ', r2,',\n',
  #              'block size: ', k2, '\n'))
  ## check the criteria
  if(v1==r2 &
     v2==r1 &
     b1==b2 &
     k1== k2){

    ### calculation N(v1*b1) matrix from D1 design
    n_mat = n_matrix(d1, v1, b1)
    # cat("\n")
    # print("Incidence Matrix (n_matrix)")
    # print(n_mat)

    ### Calculation A(v1*b1) matrix Calculation
    a_mat = a_matrix(d2, n_mat)
    # cat("\n")
    # print("a_matrix")
    # print(a_mat)

    ### rearrange the A matrix
    rearranged_a_matrix = rearrange_mat(a_mat, v2)
    # cat("\n")
    # print("Rearranged a_matrix")
    # print(rearranged_a_matrix)

    ### generate the final matrix
    if(rearranged_a_matrix [1]!='error'){
      final_matrix = matrix(, nrow = v1, ncol =v2)
      blk = ncol(rearranged_a_matrix)
      for (i in (1:v2)){
        for (j in (1:v1)){
          for (k in (1:blk)){
            if(rearranged_a_matrix[j,k]==i)
              break
          }
          final_matrix[j,i] = k
        }
      }
    #cat("\nCongratulations!! Final design has been successfully generated.\n")
      } else {
      return('Oops!! Design can'/'t be created. Some error occured.')
    }
  } else{
    return('Oops!! Design can'/'t be created. Criteria mis-matched error.')
  }
  #cat("\n")
  print("Final design(D)")
  print(final_matrix)

  return(final_matrix)
}

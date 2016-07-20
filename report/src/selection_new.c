__assume_aligned(aux_row_bitmap, MEM_LINE_SIZE);
for ( int at_row = 0; at_row < A_number_rows; ++at_row){
  aux_csc_values[at_row] =  0;
  non_zero = 0;
  field = (char*) g_quark_to_string ( at_row+1 );
  returned_strcmp = strcmp( field, comparation_key);
  if (
      ( (opp_code == LESS)  && (returned_strcmp < 0 ))
      ||
      ( (opp_code == LESS_EQ)  && (returned_strcmp <= 0 ))
      ||
      ( (opp_code == GREATER)  && (returned_strcmp > 0 ))
      ||
      ( (opp_code == GREATER_EQ)  && (returned_strcmp >= 0 ))
     ){
    returned_strcmp2 = strcmp( field, comparation_key2);
    if (
        ( (opp_code2 == LESS)  && (returned_strcmp2 < 0) )
        ||
        ( (opp_code2 == LESS_EQ)  && (returned_strcmp2 <= 0) )
        ||
        ( (opp_code2 == GREATER)  && (returned_strcmp2 > 0) )
        ||
        ( (opp_code2 == GREATER_EQ)  && (returned_strcmp2 >= 0) )
       ){
      aux_row_bitmap[at_row] =  1;
    }
  }
}

__assume_aligned(A_csc_values, MEM_LINE_SIZE);
__assume_aligned(A_row_ind, MEM_LINE_SIZE);
__assume_aligned(A_col_ptr, MEM_LINE_SIZE);

for ( int at_column = 0; at_column < A_number_columns; ++at_column){
  aux_csc_col_ptr[at_column] =  at_non_zero;
  const int iaa = A_row_ind[at_column] +1 ;
  non_zero = 0;
  if ( aux_row_bitmap[iaa] == 1 ){
    aux_csc_row_ind[at_non_zero] =  at_column;
    aux_csc_values[at_non_zero] =  A_csc_values[at_column];
    at_non_zero++;
  }
}


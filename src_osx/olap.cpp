/* ---------------------------------------------------------------------------
 **    Filename: olap_search.c
 **
 **     License: This file is part of OLAP PROJECT.
 **
 **              OLAP PROJECT is free software: you can redistribute it
 **              and/or modify it under the terms of the GNU General Public
 **              License as published by the Free Software Foundation,
 **              either version 3 of the License, or (at your option)
 **              any later version.
 **
 **              OLAP is distributed in the hope that it will be useful,
 **              but WITHOUT ANY WARRANTY; without even the implied warranty of
 **              MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 **              GNU General Public License for more details.
 **
 **              You should have received a copy of the GNU General Public
 **              License along with OLAP.
 **              If not, see <http://www.gnu.org/licenses/>.
 **
 ** Description: The proposed file focus on a typed linear algebra approach,
 **              encoding OLAP functionality solely in terms of Linear Algebra
 **              operations, represented in the CSR and BSR format, and
 **              recurring to Intel MKL and GLIB libraries.
 **
 **     Authors: Filipe Oliveira <a57816@alunos.uminho.pt>
 **          and Sérgio Caldas   <a57779@alunos.uminho.pt>
 **
 ** University of Minho, High Performance Computing Dpt. , April 2016
 ** -------------------------------------------------------------------------*/

#ifndef _olap_cpp
#define _olap_cpp 
#include "olap.hpp"

void col_read_csc (
    char* filename, int col_number
/*    ,
    int* n_nnz, int* n_rows, int* n_cols,
    vector<float>& A_csc_values,
    vector<int>& A_row_ind,
    vector<int>& A_col_ptr
 */   ){

  FILE* f = fopen(filename, "r" );
  fseek(f, 0, SEEK_END);
  long fileSize = ftell(f);
  vector<char> dev(fileSize);
  cout << " file has size " << fileSize << endl;
  fclose(f);
}

#endif

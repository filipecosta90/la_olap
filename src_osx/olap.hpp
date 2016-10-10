/* ---------------------------------------------------------------------------
 **    Filename: olap_search.h
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
 **
 ** University of Minho, High Performance Computing Dpt. , October 2016
 ** -------------------------------------------------------------------------*/

//Cache-Lines size is (typically) 32 bytes
#define MEM_LINE_SIZE 32
#define ARRAY_SIZE MEM_LINE_SIZE 
#define GROWTH_FACTOR 2
#define MAX_FIELD_SIZE 128
#define MAX_REG_SIZE 1024
#define LESS 1
#define LESS_EQ 2
#define GREATER 3
#define GREATER_EQ 4

#ifndef _olap_hpp
#define _olap_hpp
#include <string>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>

//GPU libs
//#include <thrust/device_vector.h>                                               
//#include <thrust/copy.h>                                                        
//#include <thrust/count.h> 

using namespace std;
map<string, map<string,vector<float> > > data_dict_mx_values;
map<string, map<string,vector<int> > > data_dict_mx_row_ind;
map<string, map<string,vector<int> > > data_dict_mx_col_ptr;
map<string, map<string,int > > data_dict_mx_n_rows;
map<string, map<string,int > > data_dict_mx_n_cols;
map<string, map<string,int > > data_dict_mx_n_nnz;


void col_read_csc (
    char* filename, int col_number
/*    ,
    int* n_nnz, int* n_rows, int* n_cols,
    vector<float>& A_csc_values,
    vector<int>& A_row_ind,
    vector<int>& A_col_ptr
 */   );




#endif

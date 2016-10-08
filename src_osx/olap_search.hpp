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
 **          and Sérgio Caldas   <a57779@alunos.uminho.pt>
 **
 ** University of Minho, High Performance Computing Dpt. , April 2016
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

using namespace std;
map<string, map<string,vector<float> > > data_dict_mx_values;
map<string, map<string,vector<int> > > data_dict_mx_row_ind;

////////////////////////////////////// AUX /////////////////////////////////////
////////////////////////////////////// AUX /////////////////////////////////////
////////////////////////////////////// AUX /////////////////////////////////////

//starts at position 1 
#endif

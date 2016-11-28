#include <cctype>
#include <fstream>
#include <cassert>
#include <fcntl.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>

//GPU libs
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/count.h>

#include "olap_driver.hh"
#include "olap_parser.hh"
#include "olap_scanner.hh"

struct is_newline_break{                                                                               
  __host__ __device__                                                           
    bool operator()(const char x)                                               
    {                                                                           
      return x == '\n';                                                           
    }                                                                           
};

struct parse_functor
{
  const char *source;
  int **dest_col_start;
  int **dest_col_size;
  const unsigned int *ind;
  const unsigned int *num_cols_to_parse;
  const char *field_separator;
  const int *src_newline_ind;
  const unsigned int *dest_len;

  parse_functor(
      const char* _source, int** _dest_col_start, int** _dest_col_size, 
      const unsigned int* _ind, const unsigned int* _cnt, 
      const char* _separator,
      const int* _src_ind, const unsigned int* _dest_len
      ):
    source(_source), dest_col_start(_dest_col_start), dest_col_size(_dest_col_size), ind(_ind), num_cols_to_parse(_cnt),  field_separator(_separator), src_newline_ind(_src_ind), dest_len(_dest_len) {}

  template <typename IndexType>
    __host__ __device__
    void operator()(const IndexType & i) {
      unsigned int current_column = 0, total_parsed_cols = 0, inline_pos = 0, in_col_nonzeros, line_start;
      line_start = src_newline_ind[i]+1;

      while(total_parsed_cols < *num_cols_to_parse){
        // if its the column we want to parse
        if(ind[total_parsed_cols] == current_column) { //process
          dest_col_start[total_parsed_cols][i]=line_start+inline_pos;
          in_col_nonzeros = 0;
          while(source[line_start+inline_pos] != *field_separator){
            if(source[line_start+inline_pos] != 0) {
              in_col_nonzeros++;
            };
            inline_pos++;
          };
          //save the size of the column
          dest_col_size[total_parsed_cols][i]=in_col_nonzeros;
          inline_pos++;
          total_parsed_cols++;
        }
        // ignore the current column
        else{
          while(source[line_start+inline_pos] != *field_separator) {
            inline_pos++;
          };
          inline_pos++;
        };
        current_column++;
      }
    }
};

namespace OLAP{
  OLAP_Driver::~OLAP_Driver(){
    delete(scanner);
    scanner = nullptr;
    delete(parser);
    parser = nullptr;
  }

  void OLAP_Driver::parse( const char * const filename ){
    assert( filename != nullptr );
    std::ifstream in_file( filename );
    if( ! in_file.good() ){
      exit( EXIT_FAILURE );
    }
    parse_helper( in_file );
    return;
  }

  void OLAP_Driver::parse( std::istream &stream ){
    if( ! stream.good()  && stream.eof() ){
      return;
    }
    //else
    parse_helper( stream );
    return;
  }

  void OLAP_Driver::parse_helper( std::istream &stream ){
    delete(scanner);
    try{
      scanner = new OLAP::OLAP_Scanner( &stream );
    }
    catch( std::bad_alloc &ba ){
      std::cerr << "Failed to allocate scanner: (" <<
        ba.what() << "), exiting!!\n";
      exit( EXIT_FAILURE );
    }

    delete(parser);
    try{
      parser = new OLAP::OLAP_Parser( (*scanner) /* scanner */,
          (*this) /* driver */ );
    }
    catch( std::bad_alloc &ba ){
      std::cerr << "Failed to allocate parser: (" <<
        ba.what() << "), exiting!!\n";
      exit( EXIT_FAILURE );
    }
    const int accept( 0 );
    if( parser->parse() != accept ){
      std::cerr << "Parse failed!!\n";
    }
    return;
  }

  void OLAP_Driver::load_matrix_csc ( std::string  filename, int col_number, int max_col_size ){
    std::cout << "loading column " << col_number << " from file " << filename << std::endl;
    std::clock_t start1 = std::clock();
    int fd;
    fd = open ( filename.c_str(), O_RDONLY );
    if ( fd == -1 ){
      perror("open");
    }

    off_t file_size;
    file_size = lseek(fd, 0, SEEK_END);
    thrust::device_vector<char> dev(file_size);
    char* mapped_file;

    mapped_file = (char*)mmap (0, file_size, PROT_READ, MAP_SHARED, fd, 0);                  

    if (mapped_file == MAP_FAILED){                                                        
      perror ("mmap");                                                            
    }                                                                             

    if (close (fd) == -1){                                                       
      perror ("close");                                                           
    }                                                                             

    thrust::copy(mapped_file, mapped_file+file_size, dev.begin());                                     
    int line_count = std::count(dev.begin(), dev.end(), '\n');                        
    std::cout << "There are " << line_count << " total lines in a file with size " << file_size << std::endl;    

    // find out the position of every newline
    thrust::device_vector<int> dev_newline_pos(line_count+1); 
    thrust::copy_if(
        thrust::make_counting_iterator((unsigned int) 0 ), 
        thrust::make_counting_iterator((unsigned int) file_size), 
        dev.begin(), dev_newline_pos.begin()+1, 
        is_newline_break()
        ); 

    // field position based on column number
    thrust::device_vector<unsigned int> field_index(1);
    field_index[0]=col_number;

    // field separator
    thrust::device_vector<char> field_separator(1);
    field_separator[0] = '|';

    thrust::device_vector<int> dev_col_start1(line_count); 
    thrust::device_vector<int> dev_col_size1(line_count); 
    thrust::fill(dev_col_start1.begin(), dev_col_start1.end(), 0); 
    thrust::fill(dev_col_size1.begin(), dev_col_size1.end(), 0); 

    thrust::device_vector<int*> dest_col_start(1);
    thrust::device_vector<int*> dest_col_size(1);
    dest_col_start[0] = thrust::raw_pointer_cast(dev_col_start1.data());
    dest_col_size[0] = thrust::raw_pointer_cast(dev_col_size1.data());

    //fields max lengths 
    thrust::device_vector<unsigned int> dest_len(1); 
    dest_len[0] = max_col_size;

    //fields count
    thrust::device_vector<unsigned int> ind_cnt(1); 
    ind_cnt[0] = 1;

    thrust::counting_iterator<unsigned int> begin(0);

    parse_functor ff(
        (const char*) thrust::raw_pointer_cast(dev.data()),
        (int**) thrust::raw_pointer_cast(dest_col_start.data()), 
        (int**) thrust::raw_pointer_cast(dest_col_size.data()), 
        thrust::raw_pointer_cast(field_index.data()),
        thrust::raw_pointer_cast(ind_cnt.data()), 
        thrust::raw_pointer_cast(field_separator.data()), 
        thrust::raw_pointer_cast(dev_newline_pos.data()), 
        thrust::raw_pointer_cast(dest_len.data())
        );
    thrust::for_each(begin, begin + line_count, ff);
/*
    for (int pos = 0; pos < line_count ; pos ++  ){
      std::string col ( &(mapped_file[dev_col_start1[pos]]), dev_col_size1[pos] );
      std::cout << dev_col_start1[pos] << " " << dev_col_size1[pos] << " : " << col<<  std::endl;
    }
*/
    std::cout << "leaving load matrix " << std::endl;
  }

  std::ostream& OLAP_Driver::print( std::ostream &stream ){
    stream << red  << "Debug info: " << norm << "\n";
    stream << blue << "OLAP: "  << norm << "\n";

    return(stream);
  }
}



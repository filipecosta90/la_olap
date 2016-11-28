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
      return x == 10;                                                           
    }                                                                           
};

struct parse_functor
{
  const char *source;
  char **dest;
  const unsigned int *ind;
  const unsigned int *cnt;
  const char *separator;
  const int *src_ind;
  const unsigned int *dest_len;

  parse_functor(
      const char* _source, char** _dest, 
      const unsigned int* _ind, const unsigned int* _cnt, 
      const char* _separator,
      const int* _src_ind, const unsigned int* _dest_len
      ):
    source(_source), dest(_dest), ind(_ind), cnt(_cnt),  separator(_separator), src_ind(_src_ind), dest_len(_dest_len) {}

  template <typename IndexType>
    __host__ __device__
    void operator()(const IndexType & i) {
      unsigned int curr_cnt = 0, dest_curr = 0, j = 0, t, pos;
      pos = src_ind[i]+1;

      while(dest_curr < *cnt) {
        if(ind[dest_curr] == curr_cnt) { //process
          t = 0;
          while(source[pos+j] != *separator) {
            if(source[pos+j] != 0) {
              dest[dest_curr][dest_len[dest_curr]*i+t] = source[pos+j];
              t++;
            };
            j++;
          };
          j++;
          dest_curr++;
        }
        else {
          while(source[pos+j] != *separator) {
            j++;
          };
          j++;
        };
        curr_cnt++;
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
    std::cout << "going to count" << std::endl;
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
    for ( int pos = 0; pos < dev_newline_pos.size(); pos++ ){
      std::cout << dev_newline_pos[pos] << std::endl;
    }
    // field position based on column number
    thrust::device_vector<unsigned int> field_index(1);
    field_index[0]=col_number;

    // field separator
    thrust::device_vector<char> field_separator(1);
    field_separator[0] = '|';

    thrust::device_vector<char> dev_res1(line_count*max_col_size); 
    thrust::fill(dev_res1.begin(), dev_res1.end(), 0); 

    thrust::device_vector<char*> dest(1);
    dest[0] = thrust::raw_pointer_cast(dev_res1.data());

    //fields max lengths 
    thrust::device_vector<unsigned int> dest_len(1); 
    dest_len[0] = max_col_size;

    //fields count
    thrust::device_vector<unsigned int> ind_cnt(1); 
    ind_cnt[0] = 1;

    thrust::counting_iterator<unsigned int> begin(0);

    parse_functor ff(
        (const char*) thrust::raw_pointer_cast(dev.data()),
        (char**) thrust::raw_pointer_cast(dest.data()), 
        thrust::raw_pointer_cast(field_index.data()),
        thrust::raw_pointer_cast(ind_cnt.data()), 
        thrust::raw_pointer_cast(field_separator.data()), 
        thrust::raw_pointer_cast(dev_newline_pos.data()), 
        thrust::raw_pointer_cast(dest_len.data())
        );
    thrust::for_each(begin, begin + line_count, ff);
    
    for (int pos = 0; pos < (line_count * max_col_size ); pos ++  ){
      std::cout << dev_res1[pos];
      if( (pos % (max_col_size+1) )== 0 ){
        std::cout << std::endl;
      }
    }
    std::cout << "leaving load matrix " << std::endl;
  }

  std::ostream& OLAP_Driver::print( std::ostream &stream ){
    stream << red  << "Debug info: " << norm << "\n";
    stream << blue << "OLAP: "  << norm << "\n";

    return(stream);
  }
}



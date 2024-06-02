package neat_matrix_library

import "core:fmt"
import "core:c"
import "core:c/libc"
import "core:os"
import "base:runtime"
import "core:strings"
import "core:math/rand"

// #define BLACK "\033[0;30m"
// #define RED "\033[0;31m"
// #define GREEN "\033[0;32m"
// #define BLUE "\033[0;34m"
// #define PURPLE "\033[0;35m"
// #define CYAN "\033[0;36m"
// #define WHITE "\033[0;37m"
// #define YELLOW "\033[0;33m"
// #define RESET "\033[0m"


BLACK  :: "\033[0;30m"
RED    :: "\033[0;31m"
GREEN  :: "\033[0;32m"
BLUE   :: "\033[0;34m"
PURPLE :: "\033[0;35m"
CYAN   :: "\033[0;36m"
WHITE  :: "\033[0;37m"
YELLOW :: "\033[0;33m"
RESET  :: "\033[0m"

// jnc function
get_csting_from_location :: proc ( loc : runtime.Source_Code_Location ) -> cstring {
    str_loc_tmp := fmt.aprintf( "%s  [ %d : %d ], %s ",
                         loc.file_path,
                                loc.line,
                                loc.column,
                                loc.procedure )
    cstr_loc := strings.clone_to_cstring( str_loc_tmp )
    delete( str_loc_tmp )
    return cstr_loc
}

// #define NP_CHECK(ptr) \
//         if (!(ptr)) { \
//             fprintf(stderr, "%s:%d NULL POINTER: %s n", \
//                 __FILE__, __LINE__, (#ptr)); \
//             exit(-1); \
//         } \

// loation is of type runtime.Source_Code_Location
m_np_check :: proc ( ptr: rawptr, location := #caller_location ) {
    if ptr == nil {
        cstr_loc := get_csting_from_location( location )
        fmt.fprintf( os.stderr, "NULL POINTER: \n %s  \n", cstr_loc )
        os.exit( -1 )
    }
}


// #define BUFFER_SIZE 4096
BUFFER_SIZE : uint : 4096


// #define NML_ERROR(fmt) \
//     if (DEBUG_TRUE) { \
//         nml_log(stderr, __FILE__, __LINE__, RED fmt RESET); \
//     } \


// void nml_log(
//   FILE* stream,
//   const char* file_name,
//   unsigned int line,
//   const char* format,
//   ...
// ) {
// #if DEBUG_TRUE
//   va_list argp;
//   va_start(argp, format);
//   nml_vlog(stream, file_name, line, format, argp);
//   va_end(argp);
// #endif
// }


// The following function is the previous one but in Odin, but it has errors because of the va_list
// In Odin one can't make a va_list out of the parameters of a function.

// nml_log :: proc ( stream    : ^libc.FILE,
//                   file_name : cstring,
//                   line      : uint,
//                   format    : cstring,
//                   args: ..any ) {

//     if DEBUG_TRUE {
//         argp : ^c.va_list = new( c.va_list )
//         libc.va_start( args, format )
//         nml_vlog( stream, file_name, line, format, argp )
//         libc.va_end( argp )
//     }
// }



// void nml_vlog(
//   FILE* stream,
//   const char *file_name,
//   unsigned int line,
//   const char *format,
//   va_list argp
// ){
// #if DEBUG_TRUE
//   char buffer[BUFFER_SIZE];
//   char* level;
//   int stop;

//   if (stderr == stream) {
//       level = "ERROR";
//   } else if (stdout == stream) {
//       level = "INFO";
//   }

//   // Formating string and
//   // Check if the the string has been completly written and
//   // no buffer overflow occured
//   stop = vsnprintf(buffer, BUFFER_SIZE, format, argp);
//   if (stop < BUFFER_SIZE && stop > 0) {
//     fprintf(stream, "[%s:%d] [%s] %s\n", file_name, line, level, buffer);
//   }
// #endif
// }

nml_vlog_jnc :: proc (
    stream    : ^libc.FILE,
    file_name : cstring,
    line      : uint,
    format    : cstring,
  ){
    if g_debug_true {
        // buffer : [ BUFFER_SIZE ]byte
        buffer : [ ^ ]byte = make( [ ^ ]byte, BUFFER_SIZE ) 
        level  : cstring
        stop   : i32
    
        if libc.stderr == stream {
            level = "ERROR";
        } else if libc.stdout == stream {
            level = "INFO";
        }
    
        // Formating string and
        // Check if the the string has been completly written and
        // no buffer overflow occured
        stop = libc.snprintf(buffer, BUFFER_SIZE, format )
        if stop < i32( BUFFER_SIZE ) && stop > 0 {
            
            libc.fprintf(stream, "[%s:%d] [%s] %s\n", file_name, line, level, buffer );

            // cstr_loc := get_csting_from_location( location )
            // libc.fprintf( stream, " [%s] [%s] %s\n", cstr_loc, level, buffer )
        }
    }
}
  

// double nml_rand_interval(double min, double max) {
//   double d;
//   d = (double) rand() / ((double) RAND_MAX + 1);
//   return (min + d * (max - min));
// }



nml_rand_interval :: proc ( min : f64, max : f64 ) -> f64 {
    d: f64 = rand.float64( & g_rand_generator )
    return ( min + d * ( max - min ) )
}


// void nml_log(FILE *stream, const char *file_name,
//     unsigned int line, const char *format, ...);
  
// void nml_vlog(FILE* stream, const char *file_name,
//     unsigned int line, const char *format, va_list argp);
  
// #define NML_FLOG(stream, fmt, ...) \
//     if (DEBUG_TRUE) { \
//         nml_log(stream, __FILE__, __LINE__, fmt, __VA_ARGS__); \
//     } \

// The user of this function has to format the string before calling this function. 
m_mnl_flog ::proc ( stream : ^libc.FILE, fmt : cstring, location := #caller_location) {
    if g_debug_true {
        loc := location
        cstr_file_path := strings.clone_to_cstring( loc.file_path )
        nml_vlog_jnc( stream, cstr_file_path, uint( loc.line ), fmt )
    }
}

// #define NML_FINFO(fmt, ...) \
//     if (DEBUG_TRUE) { \
//         nml_log(stdout, __FILE__, __LINE__, fmt, __VA_ARGS__); \
//     } \

// The user of this function has to format the string before calling this function. 
m_mnl_finfo ::proc ( fmt : cstring, location := #caller_location) {
    if g_debug_true {
        loc := location
        cstr_file_path := strings.clone_to_cstring( loc.file_path )
        nml_vlog_jnc( libc.stdout, cstr_file_path, uint( loc.line ), fmt )
    }
}

// #define NML_FERROR(fmt, ...) \
//     if (DEBUG_TRUE) { \
//         nml_log(stderr, __FILE__, __LINE__, RED fmt RESET, __VA_ARGS__); \
//     } \


// The user of this function has to format the string before calling this function. 
m_mnl_ferror ::proc ( fmt : cstring, location := #caller_location) {
    if g_debug_true {
        loc := location
        cstr_file_path := strings.clone_to_cstring( loc.file_path )
        nml_vlog_jnc( libc.stderr, cstr_file_path, uint( loc.line ), fmt )
    }
}

// #define NML_ERROR(fmt) \
//     if (DEBUG_TRUE) { \
//         nml_log(stderr, __FILE__, __LINE__, RED fmt RESET); \
//     } \

m_mnl_error ::proc ( fmt : cstring, location := #caller_location) {
    if g_debug_true {
        loc := location
        cstr_file_path := strings.clone_to_cstring( loc.file_path )
        nml_vlog_jnc( libc.stderr, cstr_file_path, uint( loc.line ), fmt )
    }
}




/**

Author of the
port to Odin   : João Carvalho
Date           : 2024.06.01
Description    : This is a port of the neat_matrix_library from C to
                 the Odin programming language.
License        : Apache License, Version 2.0
                 ( The same as the original library )   



Original code from:

Copyright 20201 Andrei N. Ciobanu

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Original code can be found at:

    GitHub - nomemory - neat-matrix-library 
    https://github.com/nomemory/neat-matrix-library

Great Original code blog article:

    Writing your own linear algebra matrix library in C 
    https://www.andreinc.net/2021/01/20/writing-your-own-linear-algebra-matrix-library-in-c


Have fun and enjoy the code! :)

*/

package neat_matrix_library

import "core:fmt"
import "core:math/rand"
import "core:c/libc"
import "core:strings"
import "core:math"

// Enable this to allow debugging messages
// #define DEBUG_TRUE 1

g_debug_true :: true

g_rand_generator := rand.create( 42 )

// #define NML_MIN_COEF 0.000000000000001
NML_MIN_COEF :: 0.000000000000001


//// Library structures

// typedef struct nml_mat_s {
//     unsigned int num_rows;
//     unsigned int num_cols;
//     double **data;
//     int is_square;
//   } nml_mat;
  
//   typedef struct nml_mat_lup_s {
//     nml_mat *L;
//     nml_mat *U;
//     nml_mat *P;
//     unsigned int num_permutations;
//   } nml_mat_lup;
  
//   typedef struct nml_mat_qr_s {
//     nml_mat *Q;
//     nml_mat *R;
//   } nml_mat_qr;


// Library structures

Nml_mat :: struct {
    num_rows  : uint,
    num_cols  : uint,
    data      : [^][^]f64,
    is_square : bool
}

Nml_mat_lup :: struct {
    L                : ^Nml_mat,
    U                : ^Nml_mat,
    P                : ^Nml_mat,
    num_permutations : uint
}

Nml_mat_qr :: struct {
    Q : ^Nml_mat,
    R : ^Nml_mat
}


// *****************************************************************************

// #define DEFAULT_VALUE 0.0
DEFAULT_VALUE :: 0.0

// #define CANNOT_ADD "Cannot add two matrices with different dimensions.\n"
CANNOT_ADD :: "Cannot add two matrices with different dimensions.\n"

// #define CANNOT_SUBTRACT "Cannot subctract two matrices with different dimensions.\n"
CANNOT_SUBTRACT :: "Cannot subctract two matrices with different dimensions.\n"

// #define CANNOT_MULITPLY \
//   "Cannot multiply two matrices where \
//   the number of columns of the first one \
//   is different than the number of rows of the second one.\n" \
CANNOT_MULITPLY :: "Cannot multiply two matrices where the number of columns of the first one is different than the number of rows of the second one.\n"

// #define CANNOT_REMOVE_COLUMN "Cannot remove matrix column %d. The value should be less than %d.\n" 
CANNOT_REMOVE_COLUMN :: "Cannot remove matrix column %d. The value should be less than %d.\n"

// #define CANNOT_REMOVE_ROW "Cannot remove matrix row %d. The value should be less than %d.\n" 
CANNOT_REMOVE_ROW :: "Cannot remove matrix row %d. The value should be less than %d.\n"

// #define INVALID_ROWS \
//   "Cannot create matrix with 0 number of rows. Aborting.\n" \
INVALID_ROWS :: "Cannot create matrix with 0 number of rows. Aborting.\n"

// #define INVALID_COLS \
//     "Cannot create matrix with 0 number of cols. Aborting.\n" \
INVALID_COLS :: "Cannot create matrix with 0 number of cols. Aborting.\n"

// #define CANNOT_TRACE \
//     "Cannot calculate trace. Matrix needs to be square.\n" \
CANNOT_TRACE :: "Cannot calculate trace. Matrix needs to be square.\n"

// #define CANNOT_CROUT \
//     "Cannot apply crout algorithm. Matrix needs to be square.\n" \
CANNOT_CROUT :: "Cannot apply crout algorithm. Matrix needs to be square.\n"

// #define CANNOT_SWAP_ROWS \
//      "Cannot swap rows (%d, %d) because the matrix number of rows is %d.\n" \
CANNOT_SWAP_ROWS :: "Cannot swap rows (%d, %d) because the matrix number of rows is %d.\n"

//#define CANNOT_SWAP_COLUMNS \
//      "Cannot swap columns (%d, %d) because the matrix number or columns is %d.\n" \
CANNOT_SWAP_COLUMNS :: "Cannot swap columns (%d, %d) because the matrix number or columns is %d.\n"

// #define CANNOT_ROW_MULTIPLY \
//       "Cannot multiply row (%d), maximum number of rows is %d.\n" \
CANNOT_ROW_MULTIPLY :: "Cannot multiply row (%d), maximum number of rows is %d.\n"

// #define CANNOT_COL_MULTIPLY "Cannot multiply col (%d), maximum number of columns is %d.\n" 
CANNOT_COL_MULTIPLY :: "Cannot multiply col (%d), maximum number of columns is %d.\n" 

// #define CANNOT_ADD_TO_ROW \
//       "Cannot add %2.2f x (row=%d) to row=%d. Total number of rows is: %d.\n" \
CANNOT_ADD_TO_ROW :: "Cannot add %2.2v x (row=%d) to row=%d. Total number of rows is: %d.\n"

// #define CANNOT_LU_MATRIX_SQUARE \
//      "Canot LU. Matrix (%d, %d) needs to be square.\n" \
CANNOT_LU_MATRIX_SQUARE :: "Canot LU. Matrix (%d, %d) needs to be square.\n"

// #define CANNOT_LU_MATRIX_DEGENERATE \
//       "Cannot LU. Matrix is degenerate or almost degenerate.\n" \
CANNOT_LU_MATRIX_DEGENERATE :: "Cannot LU. Matrix is degenerate or almost degenerate.\n"

// #define CANNOT_SOLVE_LIN_SYS_INVALID_B \
//       "Cannot solve system. b[%d][%d] should have size b[%d][%d].\n" \
CANNOT_SOLVE_LIN_SYS_INVALID_B :: "Cannot solve system. b[%d][%d] should have size b[%d][%d].\n"

// #define CANNOT_SET_DIAG \
//       "Cannot set diag with value(=%2.2f). Matrix is not square.\n" \
CANNOT_SET_DIAG :: "Cannot set diag with value(=%2.2v). Matrix is not square.\n"

// #define CANNOT_CONCATENATE_H \
//       "Cannot concatenate. Matrices have a different number of rows. Expected %d, found: %d.\n" \
CANNOT_CONCATENATE_H :: "Cannot concatenate. Matrices have a different number of rows. Expected %d, found: %d.\n"

// #define CANNOT_CONCATENATE_V \
//       "Cannot concatenate. Matrices have a different number of cols. Expected %d, found: %d.\n" \
CANNOT_CONCATENATE_V :: "Cannot concatenate. Matrices have a different number of cols. Expected %d, found: %d.\n"

// #define CANNOT_GET_COLUMN \
//      "Cannot get column (%d). The matrix has %d number of columns.\n" \
CANNOT_GET_COLUMN :: "Cannot get column (%d). The matrix has %d number of columns.\n"

// #define CANNOT_GET_ROW \
//      "Cannot get row (%d). The matrix has %d number of rows.\n" \
CANNOT_GET_ROW :: "Cannot get row (%d). The matrix has %d number of rows.\n"

// #define INCONSITENT_ARRAY \
//      "Cannot found element %d in the array (NULL). Expected a total of : %d elements.\n"  \
INCONSITENT_ARRAY :: "Cannot found element %d in the array (NULL). Expected a total of : %d elements.\n"

// #define INCONSITENT_VARGS \
//      "Cannot find element %d in the varargs. Expecteda total of : %d varargs.\n" \
INCONSITENT_VARGS :: "Cannot find element %d in the varargs. Expecteda total of : %d varargs.\n"

// #define CANNOT_REF_MATRIX_DEGENERATE \
//       "Cannot compute REF. Matrix is degenerate or near degenerate.\n" \
CANNOT_REF_MATRIX_DEGENERATE :: "Cannot compute REF. Matrix is degenerate or near degenerate.\n"

// #define CANNOT_OPEN_FILE "Cannot open file '%s'. Please check the path is correct and you have reading rights.\n"
CANNOT_OPEN_FILE ::  "Cannot open file '%v'. Please check the path is correct and you have reading rights.\n"

// #define INVALID_MATRIX_FILE \
//       "Invalid matrix file: %s. Cannot read data.\n" \
INVALID_MATRIX_FILE :: "Invalid matrix file: %v. Cannot read data.\n"

// #define VECTOR_J_DEGENERATE \
//       "Vector on colum %d is generate or near degenerate. Cannot proceed further.\n" \
VECTOR_J_DEGENERATE :: "Vector on colum %d is generate or near degenerate. Cannot proceed further.\n"

// #define CANNOT_QR_NON_SQUARE \
//       "We cannot QA non-square matrix[%d, %d].\n" \
CANNOT_QR_NON_SQUARE :: "We cannot QA non-square matrix[%d, %d].\n"

// #define CANNOT_COLUMN_L2NORM \
//       "Cannot get column (%d). The matrix has %d numbers of columns.\n" \
CANNOT_COLUMN_L2NORM :: "Cannot get column (%d). The matrix has %d numbers of columns.\n"

// #define CANNOT_VECT_DOT_DIMENSIONS \
//       "The two vectors have different dimensions: %d and %d.\n" \
CANNOT_VECT_DOT_DIMENSIONS :: "The two vectors have different dimensions: %d and %d.\n"



// *****************************************************************************
//
// Constructing and destroying a matrix struct
//
// *****************************************************************************

// // Dynamically allocates a new matrix struct
// nml_mat *nml_mat_new(unsigned int num_rows, unsigned int num_cols) {
//   if (num_rows == 0) {
//     NML_ERROR(INVALID_ROWS);
//     return NULL;
//   }
//   if (num_cols == 0) {
//     NML_ERROR(INVALID_COLS);
//     return NULL;
//   }
//   nml_mat *m = calloc(1, sizeof(*m));
//   NP_CHECK(m);
//   m->num_rows = num_rows;
//   m->num_cols = num_cols;
//   m->is_square = (num_rows == num_cols) ? 1 : 0;
//   m->data = calloc(m->num_rows, sizeof(*m->data));
//   NP_CHECK(m->data);
//   int i;
//   for(i = 0; i < m->num_rows; ++i) {
//     m->data[i] = calloc(m->num_cols, sizeof(**m->data));
//     NP_CHECK(m->data[i]);
//   }
//   return m;
// }

// Dynamically allocates a new matrix struct
nml_mat_new :: proc ( num_rows : uint, num_cols : uint ) -> ^Nml_mat {
    if num_rows == 0 {
        m_mnl_error( INVALID_ROWS )
        return nil
    }
    if num_cols == 0 {
        m_mnl_error( INVALID_COLS )
        return nil
    }
    // nml_mat *m = calloc(1, sizeof(*m));
    m : ^Nml_mat = new( Nml_mat )
    m_np_check( m )
    m^.num_rows = num_rows
    m^.num_cols = num_cols
    m^.is_square = (num_rows == num_cols) ? true : false
    // m->data = calloc(m->num_rows, sizeof(*m->data));
    m^.data = make( [ ^ ][ ^ ]f64, m^.num_rows )
    m_np_check( m^.data )
    for i : uint = 0; i < m^.num_rows; i += 1 {
        // m->data[i] = calloc(m->num_cols, sizeof(**m->data));
        m^.data[ i ] = make( [ ^ ]f64, m^.num_cols )
        m_np_check( m^.data[ i ] )
    }
    return m
}

// nml_mat *nml_mat_rnd(unsigned int num_rows, unsigned int num_cols, double min, double max) {
//   nml_mat *r = nml_mat_new(num_rows, num_cols);
//   int i, j;
//   for(i = 0; i < num_rows; i++) {
//     for(j = 0; j < num_cols; j++) {
//       r->data[i][j] = nml_rand_interval(min, max);
//     }
//   }
//   return r;
// }

nml_mat_rnd :: proc ( num_rows : uint, num_cols : uint, min : f64, max : f64 ) -> ^Nml_mat {
    r : ^Nml_mat = nml_mat_new(num_rows, num_cols);
    for i : uint = 0; i < num_rows; i += 1 {
        for j : uint = 0; j < num_cols; j += 1 {
            r^.data[ i ][ j ] = nml_rand_interval( min, max )
        }
    }
    return r
}

// // Dynamically alloctes a new matrix struct
// // The matrix is square (number of rows is equl with the number of cols)
// nml_mat *nml_mat_sqr(unsigned int size) {
//   return nml_mat_new(size, size);
// }

// Dynamically alloctes a new matrix struct
// The matrix is square (number of rows is equl with the number of cols)
nml_mat_sqr :: proc ( size : uint ) -> ^Nml_mat {
    return nml_mat_new( size, size )
}

// nml_mat *nml_mat_sqr_rnd(unsigned int size, double min, double max) {
//     return nml_mat_rnd(size, size, min, max);
//   }

nml_mat_sqr_rnd :: proc ( size : uint, min : f64, max : f64 ) -> ^Nml_mat {
    return nml_mat_rnd( size, size, min, max )
}

// // Dynamically allocates a a new matrix struct
// // The resulting matrix is an ID matrix (I)
// // I has 1.0 on the first diagonal
// // I is square
// nml_mat *nml_mat_eye(unsigned int size) {
//   nml_mat *r = nml_mat_new(size, size);
//   int i;
//   for(i = 0; i < r->num_rows; i++) {
//     r->data[i][i] = 1.0;
//   }
//   return r;
// }

// Dynamically allocates a a new matrix struct
// The resulting matrix is an ID matrix (I)
// I has 1.0 on the first diagonal
// I is square
nml_mat_eye :: proc ( size : uint ) -> ^Nml_mat {
    r : ^Nml_mat = nml_mat_new( size, size )
    for i : uint = 0; i < r^.num_rows; i += 1 {
        r^.data[ i ][ i ] = 1.0
    }
    return r;
}

// // Dynamically allocates a new matrix struct
// // Initialise the matrix by reading values from a vector
// nml_mat *nml_mat_from(unsigned int num_rows, unsigned int num_cols, unsigned int n_vals, double *vals) {
//   nml_mat *m = nml_mat_new(num_rows, num_cols);
//   int i, j, v_idx;
//   for(i = 0; i < m->num_rows; i++) {
//     for(j = 0; j < m->num_cols; j++) {
//       v_idx = i * m->num_cols + j;
//       m->data[i][j] = (v_idx < n_vals) ? vals[v_idx] : 0.0;
//     }
//   }
//   return m;
// }

// Dynamically allocates a new matrix struct
// Initialise the matrix by reading values from a vector
nml_mat_from ::  proc ( num_rows : uint, num_cols : uint, n_vals : uint, vals : [ ^ ]f64 ) -> ^Nml_mat {
    m : ^Nml_mat = nml_mat_new( num_rows, num_cols )
    for i : uint = 0; i < m^.num_rows; i += 1 {
        for j : uint = 0; j < m^.num_cols; j += 1 {
            v_idx : uint = i * m^.num_cols + j
            m^.data[ i ][ j ] = ( v_idx < n_vals ) ? vals[ v_idx ] : 0.0
        }
    }
    return m
}

// // Dynamically allocates a new Matrix
// // Initialise the matrix by copying another one
// nml_mat *nml_mat_cp(nml_mat *m) {
//   nml_mat *r  = nml_mat_new(m->num_rows, m->num_cols);
//   int i,j;
//   for(i = 0; i < r->num_rows; i++) {
//     for(j = 0; j < r->num_cols; j++) {
//       r->data[i][j] = m->data[i][j];
//     }
//   }
//   return r;
// }

// Dynamically allocates a new Matrix
// Initialise the matrix by copying another one
nml_mat_cp :: proc ( m : ^Nml_mat ) -> ^Nml_mat {
    r : ^Nml_mat = nml_mat_new( m^.num_rows, m^.num_cols )
    for i : uint = 0; i < r^.num_rows; i += 1 {
          for j : uint = 0; j < r^.num_cols; j += 1 {
            r^.data[ i ][ j ] = m^.data[ i ][ j ]
        }
    }
    return r
}
  
// nml_mat *nml_mat_fromfile(const char *file) {
//     FILE *m_file = fopen(file, "r");
//     if (NULL == m_file) {
//       NML_FERROR(CANNOT_OPEN_FILE, file);
//       return NULL;
//     }
//     nml_mat *r = nml_mat_fromfilef(m_file);
//     fclose(m_file);
//     return r;
//   }

nml_mat_fromfile :: proc ( file : cstring ) -> ^Nml_mat {
    m_file : ^libc.FILE = libc.fopen( file, "r" )
    if nil == m_file {
        // NML_FERROR( CANNOT_OPEN_FILE, file );
    
        str_error := fmt.aprintf( CANNOT_OPEN_FILE, file )
        cstr_error := strings.clone_to_cstring( str_error )
        delete( str_error )
    
        m_mnl_ferror( cstr_error )
        return nil
    }
    r : ^Nml_mat = nml_mat_fromfilef( m_file )
    libc.fclose( m_file )
    return r
}

// nml_mat *nml_mat_fromfilef(FILE *f) {
//     int i, j;
//     unsigned int num_rows = 0, num_cols = 0;
//     fscanf(f, "%d", &num_rows);
//     fscanf(f, "%d", &num_cols);
//     nml_mat *r = nml_mat_new(num_rows, num_cols);
//     for(i = 0; i < r->num_rows; i++) {
//       for(j = 0; j < num_cols; j++) {
//         fscanf(f, "%lf\t", &r->data[i][j]);
//       }
//     }
//     return r;
//   }

nml_mat_fromfilef :: proc ( f : ^libc.FILE ) -> ^Nml_mat {
    num_rows : int = 0
    num_cols : int = 0
    libc.fscanf( f, "%d", & num_rows )
    libc.fscanf( f, "%d", & num_cols )
    r : ^Nml_mat = nml_mat_new( uint( num_rows ), uint( num_cols ) )
    for i : uint = 0; i < r^.num_rows; i += 1 {
        for j : uint = 0; j < uint( num_cols ); j += 1 {
            libc.fscanf( f, "%lf\t", & ( r^.data[ i ][ j ] ) )
        }
    }
    return r
}

// // Frees a matrix structure
// void nml_mat_free(nml_mat *matrix) {
//   int i;
//   for(i = 0; i < matrix->num_rows; ++i) {
//     free(matrix->data[i]);
//   }
//   free(matrix->data);
//   free(matrix);
// }

// Frees a matrix structure
nml_mat_free :: proc ( matrix_m : ^Nml_mat )  {
    for i : uint = 0; i < matrix_m^.num_rows; i +=1 {
        free( matrix_m^.data[ i ] )
    }
    free( matrix_m^.data )
    free( matrix_m )
}

// *****************************************************************************
//
// Matrix Equality
//
// *****************************************************************************


// // Checks if two matrices have the same dimesions
// int nml_mat_eqdim(nml_mat *m1, nml_mat *m2) {
//     return (m1->num_cols == m2->num_cols) &&
//             (m1->num_rows == m2->num_rows);
// }

// Checks if two matrices have the same dimesions
nml_mat_eqdim :: proc ( m1 : ^Nml_mat, m2 : ^Nml_mat ) -> bool {
    res := ( m1^.num_cols == m2^.num_cols ) &&
           ( m1^.num_rows == m2^.num_rows )
    return res
}

// // Checks if two matrices have the same dimensions, and the elements
// // are all equal to each other with a given tolerance;
// // For exact equality use tolerance = 0.0
// int nml_mat_eq(nml_mat *m1, nml_mat *m2, double tolerance) {
//   if (!nml_mat_eqdim(m1, m2)) {
//     return 0;
//   }
//   int i, j;
//   for(i = 0; i < m1->num_rows; i++) {
//     for(j = 0; j < m1->num_cols; j++) {
//       if (fabs(m1->data[i][j] - m2->data[i][j]) > tolerance) {
//         return 0;
//       }
//     }
//   }
//   return 1;
// }


// Checks if two matrices have the same dimensions, and the elements
// are all equal to each other with a given tolerance;
// For exact equality use tolerance = 0.0
nml_mat_eq :: proc ( m1 : ^Nml_mat, m2 : ^Nml_mat, tolerance : f64 ) -> bool {
    if ! nml_mat_eqdim( m1, m2 ) {
        return false
    }
    for i : uint = 0; i < m1^.num_rows; i += 1 {
        for j : uint = 0; j < m1^.num_cols; j += 1 {
            if math.abs( m1^.data[ i ][ j ] - m2^.data[ i ][ j ] ) > tolerance {
            return false
            }
        }
    }
    return true
}

// *****************************************************************************
//
// Matrix printing
//
// *****************************************************************************

// // Prints the matrix on the stdout
// void nml_mat_print(nml_mat *matrix) {
//   nml_mat_printf(matrix, "%lf\t\t");
// }

// Prints the matrix on the stdout
nml_mat_print :: proc ( matrix_m  : ^Nml_mat ) {
    nml_mat_printf( matrix_m, "%lf\t\t" )
}

// // Prints the matrix on the stdout (with a custom formatting for elements)
// void nml_mat_printf(nml_mat *matrix, const char *d_fmt) {
//   int i, j;
//   fprintf(stdout, "\n");
//   for(i = 0; i < matrix->num_rows; ++i) {
//     for(j = 0; j < matrix->num_cols; ++j) {
//       fprintf(stdout, d_fmt, matrix->data[i][j]);
//     }
//     fprintf(stdout, "\n");
//   }
//   fprintf(stdout, "\n");
// }

// Prints the matrix on the stdout (with a custom formatting for elements)
nml_mat_printf :: proc ( matrix_m : ^Nml_mat, d_fmt : cstring ) {
    libc.fprintf( libc.stdout, "\n" )
    for i : uint = 0; i < matrix_m^.num_rows; i += 1 {
        for j : uint = 0; j < matrix_m^.num_cols; j += 1 {
            libc.fprintf( libc.stdout, d_fmt, matrix_m^.data[ i ][ j ] )
        }
        libc.fprintf( libc.stdout, "\n" )
    }
    libc.fprintf( libc.stdout, "\n" )
}


// *****************************************************************************
//
// Accessing and modifying matrix elements
//
// *****************************************************************************

// double nml_mat_get(nml_mat *matrix, unsigned int i, unsigned int j) {
//   return matrix->data[i][j];
// }

nml_mat_get :: #force_inline proc ( matrix_m : ^Nml_mat, i : uint, j : uint ) -> f64 {
    return matrix_m^.data[ i ][ j ]
}

// nml_mat *nml_mat_col_get(nml_mat *m, unsigned int col) {
//     if (col >= m->num_cols) {
//       NML_FERROR(CANNOT_GET_COLUMN, col, m->num_cols);
//       return NULL;
//     }
//     nml_mat *r = nml_mat_new(m->num_rows, 1);
//     int j;
//     for(j = 0; j < r->num_rows; j++) {
//       r->data[j][0] = m->data[j][col];
//     }
//     return r;
// }
  


nml_mat_col_get :: proc ( m : ^Nml_mat, col : uint ) -> ^Nml_mat {
    if col >= m^.num_cols {
        // NML_FERROR(CANNOT_GET_COLUMN, col, m->num_cols);
        
        str_error := fmt.aprintf( CANNOT_GET_COLUMN, col, m^.num_cols )
        cstr_error := strings.clone_to_cstring( str_error )
        delete( str_error )

        m_mnl_ferror( cstr_error )
        return nil
    }
    r : ^Nml_mat = nml_mat_new( m^.num_rows, 1 )
    for j : uint = 0; j < r^.num_rows; j += 1 {
        r^.data[ j ][ 0 ] = m^.data[ j ][ col ]
    }
    return r
}

// nml_mat *nml_mat_row_get(nml_mat *m, unsigned int row) {
//   if (row >= m->num_rows) {
//     NML_FERROR(CANNOT_GET_ROW, row, m->num_rows);
//     return NULL;
//   }
//   nml_mat *r = nml_mat_new(1, m->num_cols);
//   memcpy(r->data[0], m->data[row], m->num_cols * sizeof(*r->data[0]));
//   return r;
// }

nml_mat_row_get :: proc ( m : ^Nml_mat, row : uint ) -> ^Nml_mat {
  if row >= m^.num_rows {
    // NML_FERROR(CANNOT_GET_ROW, row, m->num_rows);

    str_error := fmt.aprintf( CANNOT_GET_ROW, row, m^.num_rows )
    cstr_error := strings.clone_to_cstring( str_error )
    delete( str_error )

    m_mnl_ferror( cstr_error )
    
    return nil
  }
  r : ^Nml_mat = nml_mat_new( 1, m^.num_cols )
  libc.memcpy( r^.data[ 0 ], m^.data[ row ], m^.num_cols * size_of( r^.data[ 0 ] ) )
  return r
}

// void nml_mat_set(nml_mat *matrix, unsigned int i, unsigned int j, double value) {
//     matrix->data[i][j] = value;
// }
  
nml_mat_set :: #force_inline proc ( matrix_m : ^Nml_mat, i : uint, j : uint, value : f64 ) {
  matrix_m^.data[ i ][ j ] = value
}

// // Sets all elements of a matrix to a given value
// void nml_mat_all_set(nml_mat *matrix, double value) {
//   int i, j;
//   for(i = 0; i < matrix->num_rows; i++) {
//     for(j = 0; j < matrix->num_cols; j++) {
//       matrix->data[i][j] = value;
//     }
//   }
// }

// Sets all elements of a matrix to a given value
nml_mat_all_set :: proc ( matrix_m : ^Nml_mat, value : f64 ) {
    for i : uint = 0; i < matrix_m^.num_rows; i +=1 {
        for j : uint = 0; j < matrix_m^.num_cols; j += 1 {
            matrix_m^.data[ i ] [ j ] = value
        }
    }
}

// // Sets all elements of the matrix to given value
// int nml_mat_diag_set(nml_mat *m, double value) {
//   if (!m->is_square) {
//     NML_FERROR(CANNOT_SET_DIAG, value);
//     return 0;
//   }
//   int i;
//   for(i = 0; i < m->num_rows; i++) {
//     m->data[i][i] = value;
//   }
//   return 1;
// }

// Sets all elements of the matrix to given value
nml_mat_diag_set :: proc ( m : ^Nml_mat, value : f64 ) -> bool {
    if m^.is_square  {
        // NML_FERROR(CANNOT_SET_DIAG, value);
    
        str_error := fmt.aprintf( CANNOT_SET_DIAG, value )
        cstr_error := strings.clone_to_cstring( str_error )
        delete( str_error )

        m_mnl_ferror( cstr_error )
        return false
    }
  
    for i : uint = 0; i < m^.num_rows; i += 1 {
        m^.data[ i ][ i ] = value
    }
    return true
}

// nml_mat *nml_mat_row_mult(nml_mat *m, unsigned int row, double num) {
//   nml_mat *r = nml_mat_cp(m);
//   if (!nml_mat_row_mult_r(r, row, num)) {
//     nml_mat_free(r);
//     return NULL;
//   }
//   return r;
// }

nml_mat_row_mult :: proc ( m : ^Nml_mat, row : uint, num : f64 ) -> ^Nml_mat {
    r : ^Nml_mat = nml_mat_cp( m )
    if !nml_mat_row_mult_r( r, row, num ) {
        nml_mat_free( r )
        return nil
    }
    return r
}

// int nml_mat_row_mult_r(nml_mat *m, unsigned int row, double num) {
//   if (row>= m->num_rows) {
//     NML_FERROR(CANNOT_ROW_MULTIPLY, row, m->num_rows);
//     return 0;
//   }
//   int i;
//   for(i=0; i < m->num_cols; i++) {
//     m->data[row][i] *= num;
//   }
//   return 1;
// }

nml_mat_row_mult_r :: proc ( m : ^Nml_mat, row : uint, num :  f64 ) -> bool {
    if row >= m^.num_rows {
        // NML_FERROR(CANNOT_ROW_MULTIPLY, row, m->num_rows);
    
        str_error := fmt.aprintf( CANNOT_ROW_MULTIPLY, row, m^.num_rows )
        cstr_error := strings.clone_to_cstring( str_error )
        delete( str_error )

        m_mnl_ferror( cstr_error )
    
        return false
    }
    for i : uint = 0; i < m^.num_cols; i += 1 {
        m^.data[ row ][ i ] *= num
    }
    return true
}

// nml_mat *nml_mat_col_mult(nml_mat *m, unsigned int col, double num) {
//   nml_mat *r = nml_mat_cp(m);
//   if (!nml_mat_col_mult_r(r, col, num)) {
//     nml_mat_free(r);
//     return NULL;
//   }
//   return r;
// }

nml_mat_col_mult :: proc ( m : ^Nml_mat, col : uint, num : f64 ) -> ^Nml_mat {
    r : ^Nml_mat = nml_mat_cp( m )
    if !nml_mat_col_mult_r( r, col, num ) {
        nml_mat_free( r )
        return nil
  }
  return r
}

// int nml_mat_col_mult_r(nml_mat *m, unsigned int col, double num) {
//     if (col>=m->num_cols) {
//       NML_FERROR(CANNOT_COL_MULTIPLY, col, m->num_cols);
//       return 0;
//     }
//     int i;
//     for(i = 0; i < m->num_rows; i++) {
//       m->data[i][col] *= num;
//     }
//     return 1;
//   }

nml_mat_col_mult_r :: proc ( m :  ^Nml_mat, col : uint, num: f64 ) -> bool {
    if col >= m^.num_cols {
        // NML_FERROR(CANNOT_COL_MULTIPLY, col, m->num_cols);
    
        str_error := fmt.aprintf( CANNOT_COL_MULTIPLY, col, m ^.num_cols )
        cstr_error := strings.clone_to_cstring( str_error )
        delete( str_error )

        m_mnl_ferror( cstr_error )
    
        return false
    }
    for i : uint = 0; i < m^.num_rows; i += 1 {
        m^.data[ i ][ col ] *= num
    }
    return true
}

// nml_mat *nml_mat_row_addrow(nml_mat *m, unsigned int where, unsigned int row, double multiplier) {
//   nml_mat *r = nml_mat_cp(m);
//   if (!nml_mat_row_addrow_r(m, where, row, multiplier)) {
//     nml_mat_free(r);
//     return NULL;
//   }
//   return r;
// }

nml_mat_row_addrow :: proc ( m : ^Nml_mat, where_m : uint, row : uint, multiplier : f64 ) -> ^Nml_mat {
    r : ^Nml_mat = nml_mat_cp( m )
    if !nml_mat_row_addrow_r( m, where_m, row, multiplier ) {
        nml_mat_free( r )
        return nil
    }
    return r;
}

// int nml_mat_row_addrow_r(nml_mat *m, unsigned int where, unsigned int row, double multiplier) {
//     if (where >= m->num_rows || row >= m->num_rows) {
//       NML_FERROR(CANNOT_ADD_TO_ROW, multiplier, row, where, m->num_rows);
//       return 0;
//     }
//     int i = 0;
//     for(i = 0; i < m->num_cols; i++) {
//       m->data[where][i] += multiplier * m->data[row][i];
//     }
//     return 1;
//   }

nml_mat_row_addrow_r :: proc ( m : ^Nml_mat, where_m : uint, row : uint, multiplier : f64 ) -> bool {
    if where_m >= m^.num_rows || row >= m^.num_rows {
        // NML_FERROR(CANNOT_ADD_TO_ROW, multiplier, row, where_m, m->num_rows);

        str_error := fmt.aprintf( CANNOT_ADD_TO_ROW, multiplier, row, where_m, m^.num_rows )
        cstr_error := strings.clone_to_cstring( str_error )
        delete( str_error )

        m_mnl_ferror( cstr_error )
    
        return false
    }
    for i : uint = 0; i < m^.num_cols; i += 1 {
        m^.data[ where_m ][ i ] += multiplier * m^.data[ row ][ i ]
    }
    return true
}

// nml_mat *nml_mat_smult(nml_mat *m, double num) {
//   nml_mat *r = nml_mat_cp(m);
//   nml_mat_smult_r(r, num);
//   return r;
// }

nml_mat_smult :: proc ( m : ^Nml_mat, num : f64 ) -> ^Nml_mat {
    r : ^Nml_mat = nml_mat_cp( m )
    nml_mat_smult_r( r, num )
    return r
}

// int nml_mat_smult_r(nml_mat *m, double num) {
//   int i, j;
//   for(i = 0; i < m->num_rows; i++) {
//     for(j = 0; j < m->num_cols; j++) {
//       m->data[i][j] *= num;
//     }
//   }
//   return 1;
// }

nml_mat_smult_r :: proc ( m : ^Nml_mat, num : f64 ) -> bool {
    for i : uint = 0; i < m^.num_rows; i += 1 {
        for j : uint = 0; j < m^.num_cols; j += 1 {
            m^.data[ i ][ j ] *= num
        }
    }
    return true
}

// *****************************************************************************
//
// Modifying the matrix structure
//
// *****************************************************************************


// nml_mat *nml_mat_col_rem(nml_mat *m, unsigned int column) {
//   if(column >= m->num_cols) {
//     NML_FERROR(CANNOT_REMOVE_COLUMN, column, m->num_cols);
//     return NULL;
//   }
//   nml_mat *r = nml_mat_new(m->num_rows, m->num_cols-1);
//   int i, j, k;
//   for(i = 0; i < m->num_rows; i++) {
//     for(j = 0, k=0; j < m->num_cols; j++) {
//       if (column!=j) {
//         r->data[i][k++] = m->data[i][j];
//       }
//     }
//   }
//   return r;
// }

nml_mat_col_rem :: proc ( m : ^Nml_mat, column : uint ) -> ^Nml_mat {
    if column >= m^.num_cols {
        // NML_FERROR(CANNOT_REMOVE_COLUMN, column, m->num_cols);

        str_error := fmt.aprintf( CANNOT_REMOVE_COLUMN, column, m^.num_cols )
        cstr_error := strings.clone_to_cstring( str_error )
        delete( str_error )

        m_mnl_ferror( cstr_error )
        return nil
    }
    r : ^Nml_mat = nml_mat_new( m^.num_rows, m^.num_cols - 1 )
    for i : uint = 0; i < m^.num_rows; i += 1 {
        k : uint = 0
        for j : uint = 0; j < m^.num_cols; j +=1  {
            if column != j {
                r^.data[ i ][ k ] = m^.data[ i ][ j ]
                k *= 1
            }
        }
    }
    return r;
}

// nml_mat *nml_mat_row_rem(nml_mat *m, unsigned int row) {
//   if (row >= m->num_rows) {
//     NML_FERROR(CANNOT_REMOVE_ROW, row, m->num_rows);
//     return NULL;
//   }
//   nml_mat *r = nml_mat_new(m->num_rows-1, m->num_cols);
//   int i, j, k;
//   for(i = 0, k = 0; i < m->num_rows; i++) {
//     if (row!=i) {
//       for(j = 0; j < m->num_cols; j++) {
//         r->data[k][j] = m->data[i][j];
//       }
//       k++;
//     }
//   }
//   return r;
// }

nml_mat_row_rem :: proc ( m : ^Nml_mat, row : uint ) -> ^Nml_mat {
    if row >= m^.num_rows {
        // NML_FERROR(CANNOT_REMOVE_ROW, row, m->num_rows);
 
        str_error := fmt.aprintf( CANNOT_REMOVE_ROW, row, m^.num_rows )
        cstr_error := strings.clone_to_cstring( str_error )
        delete( str_error )

        m_mnl_ferror( cstr_error )

        return nil
    }
    r : ^Nml_mat = nml_mat_new( m^.num_rows - 1, m^.num_cols )
    k : uint = 0;
    for i : uint = 0; i < m^.num_rows; i += 1 {
        if (row!=i) {
            for j : uint = 0; j < m^.num_cols; j += 1 {
                r^.data[ k ][ j ] = m^.data[ i ][ j ]
            }
            k += 1
        }
    }
    return r
}

// nml_mat *nml_mat_row_swap(nml_mat *m, unsigned int row1, unsigned int row2) {
//   nml_mat *r = nml_mat_cp(m);
//   if (!nml_mat_row_swap_r(r, row1, row2)) {
//     nml_mat_free(r);
//     return NULL;
//   }
//   return r;
// }

nml_mat_row_swap :: proc ( m : ^Nml_mat, row1 : uint, row2 : uint ) -> ^Nml_mat {
    r : ^Nml_mat = nml_mat_cp( m )
    if !nml_mat_row_swap_r( r, row1, row2 ) {
        nml_mat_free( r )
        return nil
    }
    return r
}

// int nml_mat_row_swap_r(nml_mat *m, unsigned int row1, unsigned int row2) {
//   if (row1 >= m->num_rows || row2 >= m->num_rows) {
//     NML_FERROR(CANNOT_SWAP_ROWS, row1, row2, m->num_rows);
//     return 0;
//   }
//   double *tmp = m->data[row2];
//   m->data[row2] = m->data[row1];
//   m->data[row1] = tmp;
//   return 1;
// }

nml_mat_row_swap_r :: proc ( m : ^Nml_mat, row1 : uint, row2 : uint ) -> bool {
    if row1 >= m^.num_rows || row2 >= m^.num_rows  {
        // NML_FERROR(CANNOT_SWAP_ROWS, row1, row2, m->num_rows);
        
        str_error := fmt.aprintf( CANNOT_SWAP_ROWS, row1, row2, m^.num_rows )
        cstr_error := strings.clone_to_cstring( str_error )
        delete( str_error )

        m_mnl_ferror( cstr_error )
        
        return false
    }
    tmp : ^f64 = m^.data[ row2 ]
    m^.data[ row2 ] = m^.data[ row1 ]
    m^.data[ row1 ] = tmp
    return true
}

// nml_mat *nml_mat_col_swap(nml_mat *m, unsigned int col1, unsigned int col2) {
//   nml_mat *r = nml_mat_cp(m);
//   if (!nml_mat_col_swap_r(r, col1, col2)) {
//     nml_mat_free(r);
//     return NULL;
//   }
//   return r;
// }


nml_mat_col_swap :: proc ( m : ^Nml_mat, col1 : uint, col2 : uint ) -> ^Nml_mat{
  r : ^Nml_mat = nml_mat_cp( m )
  if !nml_mat_col_swap_r( r, col1, col2 ) {
    nml_mat_free( r )
    return nil
  }
  return r
}

// int nml_mat_col_swap_r(nml_mat *m, unsigned int col1, unsigned int col2) {
//   if (col1 >= m->num_cols || col2 >= m->num_cols) {
//     NML_FERROR(CANNOT_SWAP_ROWS, col1, col2, m->num_cols);
//     return 0;
//   }
//   double tmp;
//   int j;
//   for(j = 0; j < m->num_rows; j++) {
//     tmp = m->data[j][col1];
//     m->data[j][col1] = m->data[j][col2];
//     m->data[j][col2] = tmp;
//   }
//   return 1;
// }

nml_mat_col_swap_r :: proc ( m : ^Nml_mat, col1 : uint, col2 : uint ) -> bool {
    if col1 >= m^.num_cols || col2 >= m^.num_cols {
        // NML_FERROR(CANNOT_SWAP_ROWS, col1, col2, m->num_cols);
    
        str_error := fmt.aprintf( CANNOT_SWAP_ROWS, col1, col2, m^.num_cols )
        cstr_error := strings.clone_to_cstring( str_error )
        delete( str_error )

        m_mnl_ferror( cstr_error )
   
        return false
    }
    for j : uint = 0; j < m^.num_rows; j += 1 {
        tmp : f64 = m^.data[ j ][ col1 ]
        m^.data[ j ][ col1 ] = m^.data[ j ][ col2 ]
        m^.data[ j ][ col2 ] = tmp
    }
    return true
}

// nml_mat *nml_mat_cath(unsigned int mnum, nml_mat **marr) {
//   if (0==mnum) {
//     return NULL;
//   }
//   if (1==mnum) {
//     // We just return the one matrix supplied as the first param
//     // no need for additional logic
//     return nml_mat_cp(marr[0]);
//   }
//   // We calculate the total number of columns to know how to allocate memory
//   // for the resulting matrix
//   int i,j,k,offset;
//   unsigned int lrow, ncols;
//   lrow = marr[0]->num_rows;
//   ncols = marr[0]->num_cols;
//   for(k = 1; k < mnum; k++) {
//     if (NULL == marr[k]) {
//       NML_FERROR(INCONSITENT_ARRAY, k, mnum);
//       return NULL;
//     }
//     if (lrow != marr[k]->num_rows) {
//       NML_FERROR(CANNOT_CONCATENATE_H, lrow, marr[k]->num_rows);
//       return NULL;
//     }
//     ncols+=marr[k]->num_cols;
//   }
//   // At this point we know how the resulting matrix looks like,
//   // we allocate memory for it accordingly
//   nml_mat *r = nml_mat_new(lrow, ncols);
//   for(i = 0; i < r->num_rows; i++) {
//     k = 0;
//     offset = 0;
//     for(j = 0; j < r->num_cols; j++) {
//       // If the column index of marr[k] overflows
//       if (j-offset == marr[k]->num_cols) {
//         offset += marr[k]->num_cols;
//         // We jump to the next matrix in the array
//         k++;
//       }
//       r->data[i][j] = marr[k]->data[i][j - offset];
//     }
//   }
//   return r;
// }

nml_mat_cath :: proc ( mnum : uint, marr : [ ^ ]^Nml_mat ) -> ^Nml_mat {
    if 0 == mnum {
        return nil
    }
    if 1 == mnum {
        // We just return the one matrix supplied as the first param
        // no need for additional logic
        return nml_mat_cp (marr[ 0 ] )
    }
    // We calculate the total number of columns to know how to allocate memory
    // for the resulting matrix
    i,j,k : uint
    offset : uint
    lrow, ncols : uint
    lrow  = marr[ 0 ]^.num_rows
    ncols = marr[ 0 ]^.num_cols
    for k = 1; k < mnum; k += 1 {
        if nil == marr[ k ] {
            // NML_FERROR(INCONSITENT_ARRAY, k, mnum);

            str_error := fmt.aprintf( INCONSITENT_ARRAY, k, mnum )
            cstr_error := strings.clone_to_cstring( str_error )
            delete( str_error )
    
            m_mnl_ferror( cstr_error )
       
            return nil
        }
        if lrow != marr[ k ]^.num_rows {
            // NML_FERROR(CANNOT_CONCATENATE_H, lrow, marr[k]->num_rows);
       
            str_error := fmt.aprintf( CANNOT_CONCATENATE_H, lrow, marr[ k ]^.num_rows )
            cstr_error := strings.clone_to_cstring( str_error )
            delete( str_error )
    
            m_mnl_ferror( cstr_error )

            return nil
        }
        ncols += marr[ k ]^.num_cols
    }
    // At this point we know how the resulting matrix looks like,
    // we allocate memory for it accordingly
    r : ^Nml_mat = nml_mat_new( lrow, ncols )
    for i = 0; i < r^.num_rows; i += 1 {
        k = 0
        offset = 0
        for j = 0; j < r^.num_cols; j += 1 {
            // If the column index of marr[k] overflows
            if j - offset == marr[ k ]^.num_cols {
                offset += marr[ k ]^.num_cols
                // We jump to the next matrix in the array
                k += 1
            }
            r^.data[ i ][ j ] = marr[ k ]^.data[ i ][ j - offset ]
        }
    }
    return r
}

// // Concatenates a variable number of matrices into one.
// // The concentation is done vertically this means the matrices need to have
// // the same number of columns, while the number of rows is allowed to
// // be variable
// nml_mat *nml_mat_catv(unsigned int mnum, nml_mat **marr) {
//   if (0 == mnum) {
//     return NULL;
//   }
//   if (1 == mnum) {
//     return nml_mat_cp(marr[0]);
//   }
//   // We check to see if the matrices have the same number of columns
//   int lcol, i, j, k, offset;
//   unsigned int numrows;
//   nml_mat *r;
//   lcol = marr[0]->num_cols;
//   numrows = 0;
//   for(i = 0; i < mnum; i++) {
//     if (NULL==marr[i]) {
//       NML_FERROR(INCONSITENT_ARRAY, i, mnum);
//       return NULL;
//     }
//     if (lcol != marr[i]->num_cols) {
//       NML_FERROR(CANNOT_CONCATENATE_V,lcol,marr[i]->num_cols);
//       return NULL;
//     }
//     // In the same time we calculate the resulting matrix number of rows
//     numrows+=marr[i]->num_rows;
//   }
//   // At this point we know the dimensions of the resulting Matrix
//   r = nml_mat_new(numrows, lcol);
//   // We start copying the values one by one
//   for(j = 0; j < r->num_cols; j++) {
//     offset = 0;
//     k = 0;
//     for(i = 0; i < r->num_rows; i++) {
//       if (i - offset == marr[k]->num_rows) {
//         offset += marr[k]->num_rows;
//         k++;
//       }
//       r->data[i][j] = marr[k]->data[i-offset][j];
//     }
//   }
//   nml_mat_print(r);
//   return r;
// }


// Concatenates a variable number of matrices into one.
// The concentation is done vertically this means the matrices need to have
// the same number of columns, while the number of rows is allowed to
// be variable
nml_mat_catv :: proc ( mnum : uint, marr : [ ^ ]^Nml_mat ) -> ^Nml_mat {
    if  0 == mnum {
        return nil
    }
    if 1 == mnum {
        return nml_mat_cp( marr[0] )
    }
    // We check to see if the matrices have the same number of columns
    lcol, i, j, k, offset : uint
    numrows : uint
    r : ^Nml_mat
    lcol = marr[ 0 ]^.num_cols
    numrows = 0
    for i = 0; i < mnum; i +=1 {
        if nil == marr[ i ] {
            // NML_FERROR(INCONSITENT_ARRAY, i, mnum);
            
            str_error := fmt.aprintf( INCONSITENT_ARRAY, i, mnum )
            cstr_error := strings.clone_to_cstring( str_error )
            delete( str_error )
    
            m_mnl_ferror( cstr_error )
       
            return nil
        }
        if lcol != marr[ i ]^.num_cols {
            // NML_FERROR(CANNOT_CONCATENATE_V,lcol,marr[i]->num_cols);
            
            str_error := fmt.aprintf( CANNOT_CONCATENATE_V,lcol,marr[ i ]^.num_cols )
            cstr_error := strings.clone_to_cstring( str_error )
            delete( str_error )
    
            m_mnl_ferror( cstr_error )

            return nil
        }
        // In the same time we calculate the resulting matrix number of rows
        numrows += marr[ i ]^.num_rows
    }
    // At this point we know the dimensions of the resulting Matrix
    r = nml_mat_new( numrows, lcol )
    // We start copying the values one by one
    for j = 0; j < r^.num_cols; j += 1 {
        offset = 0
        k = 0
        for i = 0; i < r^.num_rows; i +=1  {
            if  i - offset == marr[ k ]^.num_rows {
                offset += marr[ k ]^.num_rows
                k += 1
            }
            r^.data[ i ][ j ] = marr[ k ]^.data[ i - offset ][ j ]
        }
    }
    nml_mat_print( r )
    return r
}

// *****************************************************************************
//
// Matrix Operations
//
// *****************************************************************************
//
// Matrix Operations
//


// nml_mat *nml_mat_add(nml_mat *m1, nml_mat *m2) {
//   nml_mat *r = nml_mat_cp(m1);
//   if (!nml_mat_add_r(r, m2)) {
//     nml_mat_free(r);
//     return NULL;
//   }
//   return r;
// }


nml_mat_add :: proc ( m1 : ^Nml_mat, m2 : ^Nml_mat ) -> ^Nml_mat {
    r : ^Nml_mat = nml_mat_cp( m1 )
    if !nml_mat_add_r( r, m2 ) {
        nml_mat_free( r )
        return nil
    }
    return r
}

// int nml_mat_add_r(nml_mat *m1, nml_mat *m2) {
//   if (!nml_mat_eqdim(m1, m2)) {
//     NML_ERROR(CANNOT_ADD);
//     return 0;
//   }
//   int i, j;
//   for(i = 0; i < m1->num_rows; i++) {
//     for(j = 0; j < m1->num_cols; j++) {
//       m1->data[i][j] += m2->data[i][j];
//     }
//   }
//   return 1;
// }

nml_mat_add_r :: proc ( m1 : ^Nml_mat, m2 : ^Nml_mat ) -> bool {
    if  !nml_mat_eqdim( m1, m2 ) {
        // NML_ERROR(CANNOT_ADD);

        m_mnl_error( CANNOT_ADD )
        return false
    }
    for i : uint = 0; i < m1^.num_rows; i += 1 {
        for j : uint = 0; j < m1^.num_cols; j += 1 {
            m1^.data[ i ][ j ] += m2^.data[ i ][ j ]
        }
    }
    return true
}

// nml_mat *nml_mat_sub(nml_mat *m1, nml_mat *m2) {
//   nml_mat *r = nml_mat_cp(m2);
//   if (!nml_mat_sub_r(r, m2)) {
//     nml_mat_free(r);
//     return NULL;
//   }
//   return r;
// }



nml_mat_sub :: proc ( m1 : ^Nml_mat, m2 : ^Nml_mat ) -> ^Nml_mat {
    r : ^Nml_mat = nml_mat_cp( m2 )
    if !nml_mat_sub_r( r, m2 ) {
        nml_mat_free( r )
        return nil
    }
    return r
}

// int nml_mat_sub_r(nml_mat *m1, nml_mat *m2) {
//   if (!nml_mat_eqdim(m1, m2)) {
//     NML_ERROR(CANNOT_SUBTRACT);
//     return 0;
//   }
//   int i, j;
//   for(i = 0; i < m1->num_rows; i++) {
//     for(j = 0; j < m1->num_cols; j++) {
//       m1->data[i][j] -= m2->data[i][j];
//     }
//   }
//   return 1;
// }

nml_mat_sub_r :: proc ( m1 : ^Nml_mat, m2 : ^Nml_mat ) -> bool {
    if !nml_mat_eqdim( m1, m2 ) {
        // NML_ERROR(CANNOT_SUBTRACT);

        m_mnl_error( CANNOT_SUBTRACT )
        return false
    }
    for i : uint = 0; i < m1^.num_rows; i += 1 {
        for j : uint = 0; j < m1^.num_cols; j += 1 {
            m1^.data[ i ][ j ] -= m2^.data[ i ][ j ]
        }
    }
    return true
}

// nml_mat *nml_mat_dot(nml_mat *m1, nml_mat *m2) {
//   if (!(m1->num_cols == m2->num_rows)) {
//     NML_ERROR(CANNOT_MULITPLY);
//     return NULL;
//   }
//   int i, j, k;
//   nml_mat *r = nml_mat_new(m1->num_rows, m2->num_cols);
//   for(i = 0; i < r->num_rows; i++) {
//     for(j = 0; j < r->num_cols; j++) {
//       for(k = 0; k < m1->num_cols; k++) {
//         r->data[i][j] += m1->data[i][k] * m2->data[k][j];
//       }
//     }
//   }
//   return r;
// }

nml_mat_dot :: proc  ( m1 : ^Nml_mat, m2 : ^Nml_mat ) -> ^Nml_mat {
    if !( m1^.num_cols == m2^.num_rows ) {
        // NML_ERROR(CANNOT_MULITPLY);
        m_mnl_error( CANNOT_MULITPLY )
        return nil
    }
    r : ^Nml_mat = nml_mat_new( m1^.num_rows, m2^.num_cols )
    for i : uint = 0; i < r^.num_rows; i += 1 {
        for j : uint = 0; j < r^.num_cols; j += 1 {
            for k : uint = 0; k < m1^.num_cols; k += 1 {
                r^.data[ i ][ j ] += m1^.data[ i ][ k ] * m2^.data[ k ][ j ]
            }
        }
    }
    return r
}

// nml_mat *nml_mat_transp(nml_mat *m) {
//   int i, j;
//   nml_mat *r = nml_mat_new(m->num_cols, m->num_rows);
//   for(i = 0; i < r->num_rows; i++) {
//     for(j = 0; j < r->num_cols; j++) {
//       r->data[i][j] = m->data[j][i];
//     }
//   }
//   return r;
// }


nml_mat_transp :: proc ( m : ^Nml_mat ) -> ^Nml_mat {
    r : ^Nml_mat = nml_mat_new( m^.num_cols, m^.num_rows )
    for i : uint = 0; i < r^.num_rows; i += 1 {
        for j : uint = 0; j < r^.num_cols; j += 1 {
            r^.data[ i ][ j ] = m^.data[ j ][ i ]
        }
    }
    return r
}

// double nml_mat_trace(nml_mat* m) {
//   if (!m->is_square) {
//     NML_ERROR(CANNOT_TRACE);
//   }
//   int i;
//   double trace = 0.0;
//   for(i = 0; i < m->num_rows; i++) {
//     trace += m->data[i][i];
//   }
//   return trace;
// }

nml_mat_trace :: proc ( m : ^Nml_mat ) -> f64 {
    if !m^.is_square {
        // NML_ERROR(CANNOT_TRACE);
        m_mnl_error( CANNOT_TRACE )
    }
    trace : f64 = 0.0
    for i : uint = 0; i < m^.num_rows; i += 1 {
        trace += m^.data[ i ][ i ]
    }
    return trace
}

// *****************************************************************************
//
// Row Echelon
//
// *****************************************************************************


// // Finds the first non-zero element on the col column, under the row row.
// // This is used to determine the pivot in gauss Elimination
// // If not pivot is found, returns -1
// int _nml_mat_pivotidx(nml_mat *m, unsigned int col, unsigned int row) {
//   // No validations are made, this is an API Method
//   int i;
//   for(i = row; i < m->num_rows; i++) {
//     if (fabs(m->data[i][col]) > NML_MIN_COEF) {
//       return i;
//     }
//   }
//   return -1;
// }

// Finds the first non-zero element on the col column, under the row row.
// This is used to determine the pivot in gauss Elimination
// If not pivot is found, returns -1
_nml_mat_pivotidx :: proc ( m : ^Nml_mat, col : uint, row : uint ) -> int {
    // No validations are made, this is an API Method
    for i : uint = row; i < m ^.num_rows; i += 1 {
        if math.abs( m^.data[ i ][ col ] ) > NML_MIN_COEF {
          return int( i )
        }
    }
    return -1
}





// // Find the max element from the column "col" under the row "row"
// // This is needed to pivot in Gauss-Jordan elimination
// // If pivot is not found, return -1
// int _nml_mat_pivotmaxidx(nml_mat *m, unsigned int col, unsigned int row) {
//   int i, maxi;
//   double micol;
//   double max = fabs(m->data[row][col]);
//   maxi = row;
//   for(i = row; i < m->num_rows; i++) {
//     micol = fabs(m->data[i][col]);
//     if (micol>max) {
//       max = micol;
//       maxi = i;
//     }
//   }
//   return (max < NML_MIN_COEF) ? -1 : maxi;
// }


// Find the max element from the column "col" under the row "row"
// This is needed to pivot in Gauss-Jordan elimination
// If pivot is not found, return -1
_nml_mat_pivotmaxidx :: proc ( m : ^Nml_mat, col : uint, row : uint ) -> int {
    i    : uint
    maxi : uint
    micol : f64
    max : f64 = math.abs( m^.data[ row ][ col ] )
    maxi = row 
    for i : uint = row; i < m^.num_rows; i += 1 {
        micol = math.abs( m^.data[ i ][ col ] )
        if  micol > max {
            max = micol
            maxi = i
        }
    }
    return ( max < NML_MIN_COEF ) ? -1 : int( maxi )
}

// // Retrieves the matrix in Row Echelon form using Gauss Elimination
// nml_mat *nml_mat_ref(nml_mat *m) {
//   nml_mat *r = nml_mat_cp(m);
//   int i, j, k, pivot;
//   j = 0, i = 0;
//   while(j < r->num_cols && i < r->num_cols) {
//     // Find the pivot - the first non-zero entry in the first column of the matrix
//     pivot = _nml_mat_pivotidx(r, j, i);
//     if (pivot<0) {
//       // All elements on the column are zeros
//       // We move to the next column without doing anything
//       j++;
//       continue;
//     }
//     // We interchange rows moving the pivot to the first row that doesn't have
//     // already a pivot in place
//     if (pivot!=i) {
//       nml_mat_row_swap_r(r, i, pivot);
//     }
//     // Multiply each element in the pivot row by the inverse of the pivot
//     nml_mat_row_mult_r(r, i, 1/r->data[i][j]);
//     // We add multiplies of the pivot so every element on the column equals 0
//     for(k = i+1; k < r->num_rows; k++) {
//       if (fabs(r->data[k][j]) > NML_MIN_COEF) {
//         nml_mat_row_addrow_r(r, k, i, -(r->data[k][j]));
//       } 
//     }
//     i++;
//     j++;
//   }
//   return r;
// }

// Retrieves the matrix in Row Echelon form using Gauss Elimination
nml_mat_ref :: proc ( m : ^Nml_mat ) -> ^Nml_mat {
    r : ^Nml_mat = nml_mat_cp( m )
    i, j, k : uint
    pivot : int
    j = 0
    i = 0
    for j < r^.num_cols && i < r^.num_cols {
        // Find the pivot - the first non-zero entry in the first column of the matrix
        pivot = _nml_mat_pivotidx( r, j, i )
        if pivot < 0 {
            // All elements on the column are zeros
            // We move to the next column without doing anything
            j += 1
            continue
        }
        // We interchange rows moving the pivot to the first row that doesn't have
        // already a pivot in place
        if uint( pivot ) != i {
            nml_mat_row_swap_r( r, i, uint( pivot ) )
        }
        // Multiply each element in the pivot row by the inverse of the pivot
        nml_mat_row_mult_r( r, i, 1 / r^.data[ i ][ j ] )
        // We add multiplies of the pivot so every element on the column equals 0
        for k = i+1; k < r^.num_rows; k += 1 {
            if math.abs( r^.data[ k ][ j ] ) > NML_MIN_COEF {
                nml_mat_row_addrow_r( r, k, i, -(r^.data[ k ][ j ] ) )
            }
        }
        i += 1
        j += 1
    }
    return r
}

// // Retrieves the matrix in Reduced Row Echelon using Guass-Jordan Elimination
// nml_mat *nml_mat_rref(nml_mat *m) {
//   nml_mat* r = nml_mat_cp(m);
//   int i,j,k,pivot;
//   i = 0;
//   j = 0;
//   while(j < r->num_cols && i < r->num_rows) {
//     // We find the pivot, the maximum row id (fabs) in the column
//     pivot = _nml_mat_pivotmaxidx(r, j, i);
//     if (pivot<0) {
//       // No pivot, we change columns
//       j++;
//       continue;
//     }
//     // We interchange rows to out the pivot row into the 
//     // desired position
//     if (pivot!=i) {
//       nml_mat_row_swap_r(r, i, pivot);
//     }
//     // We create 1 in the pivot position
//     nml_mat_row_mult_r(r, i, 1/r->data[i][j]);
//      // We put zeros on the colum with the pivot
//     for(k = 0; k < r->num_rows; k++) {
//       if (!(k==i)) {
//         nml_mat_row_addrow_r(r, k, i, -(r->data[k][j]));
//       }
//     }
//     i++;
//     j++;
//   }
//   return r;
// }

// Retrieves the matrix in Reduced Row Echelon using Guass-Jordan Elimination
nml_mat_rref :: proc ( m : ^Nml_mat ) -> ^Nml_mat {
    r : ^Nml_mat = nml_mat_cp( m )
    i,j,k : uint
    pivot : int 
    i = 0
    j = 0
    for j < r^.num_cols && i < r^.num_rows {
        // We find the pivot, the maximum row id (fabs) in the column
        pivot = _nml_mat_pivotmaxidx( r, j, i )
        if pivot < 0 {
            // No pivot, we change columns
            j += 1
            continue
        }
        // We interchange rows to out the pivot row into the 
        // desired position
        if uint( pivot ) != i {
            nml_mat_row_swap_r( r, i, uint( pivot ) )
        }
        // We create 1 in the pivot position
        nml_mat_row_mult_r( r, i, 1 / r^.data[ i ][ j ] )
        // We put zeros on the colum with the pivot
        for k = 0; k < r^.num_rows; k += 1 {
            if !( k == i ) {
                nml_mat_row_addrow_r( r, k, i, -(r^.data[ k ][ j ] ) )
            }
        }
        i += 1
        j += 1
    }
    return r
}

// *****************************************************************************
//
// LUP Decomposition
//
// *****************************************************************************

// // Finds the maxid on the column (starting from k -> num_rows)
// // This method is used for pivoting in LUP decomposition
// int _nml_mat_absmaxr(nml_mat *m, unsigned int k) {
//   // Find max id on the column;
//   int i;
//   double max = m->data[k][k];
//   int maxIdx = k;
//   for(i = k+1; i < m->num_rows; i++) {
//     if (fabs(m->data[i][k]) > max) {
//       max = fabs(m->data[i][k]);
//       maxIdx = i;
//     }
//   }
//   return maxIdx;
// }

// Finds the maxid on the column (starting from k -> num_rows)
// This method is used for pivoting in LUP decomposition
_nml_mat_absmaxr :: proc ( m : ^Nml_mat, k :uint ) -> int {
    // Find max id on the column;
    max : f64 = m^.data[ k ][ k ]
    maxIdx : uint = k
    for i : uint = k+1; i < m^.num_rows; i += 1 {
        if  math.abs( m^.data[ i ][ k ] ) > max {
            max = math.abs( m^.data[ i ][ k ] )
            maxIdx = i
        }
    }
    return int( maxIdx )
}

// // Allocates memory for a new nml_mat_lup structure
// nml_mat_lup *nml_mat_lup_new(nml_mat *L, nml_mat *U, nml_mat *P, unsigned int num_permutations) {
//   nml_mat_lup *r = malloc(sizeof(*r));
//   NP_CHECK(r);
//   r->L = L;
//   r->U = U;
//   r->P = P;
//   r->num_permutations = num_permutations;
//   return r;
// }

// Allocates memory for a new nml_mat_lup structure
nml_mat_lup_new :: proc ( L : ^Nml_mat, U : ^Nml_mat, P : ^Nml_mat, num_permutations : uint ) -> ^Nml_mat_lup {
    // nml_mat_lup *r = malloc(sizeof(*r));  
    r : ^Nml_mat_lup = new( Nml_mat_lup )
    // NP_CHECK(r);
    m_np_check( r )
    r^.L = L
    r^.U = U
    r^.P = P
    r^.num_permutations = num_permutations
    return r
}

// void nml_mat_lup_free(nml_mat_lup* lu) {
//   nml_mat_free(lu->P);
//   nml_mat_free(lu->L);
//   nml_mat_free(lu->U);
//   free(lu);
// }

nml_mat_lup_free :: proc ( lu : ^Nml_mat_lup ) {
    nml_mat_free( lu^.P )
    nml_mat_free( lu^.L )
    nml_mat_free( lu^.U )
    free( lu )
}

// void nml_mat_lup_print(nml_mat_lup *lu) {
//   nml_mat_print(lu->L);
//   nml_mat_print(lu->U);
//   nml_mat_print(lu->P);
// }

nml_mat_lup_print :: proc ( lu : ^Nml_mat_lup) {
    nml_mat_print( lu^.L )
    nml_mat_print( lu^.U )
    nml_mat_print( lu^.P )
}

// void nml_mat_lup_printf(nml_mat_lup *lu, const char *fmt) {
//   nml_mat_printf(lu->L, fmt);
//   nml_mat_printf(lu->U, fmt);
//   nml_mat_printf(lu->P, fmt);
// }

nml_mat_lup_printf :: proc ( lu : ^Nml_mat_lup, fmt : cstring ) {
    nml_mat_printf( lu^.L, fmt )
    nml_mat_printf( lu^.U, fmt )
    nml_mat_printf( lu^.P, fmt )
}

// nml_mat_lup *nml_mat_lup_solve(nml_mat *m) {
//   if (!m->is_square) {
//     NML_FERROR(CANNOT_LU_MATRIX_SQUARE, m->num_rows, m-> num_cols);
//     return NULL;
//   }
//   nml_mat *L = nml_mat_new(m->num_rows, m->num_rows);
//   nml_mat *U = nml_mat_cp(m);
//   nml_mat *P = nml_mat_eye(m->num_rows);

//   int j,i, pivot;
//   unsigned int num_permutations = 0;
//   double mult;

//   for(j = 0; j < U->num_cols; j++) {
//     // Retrieves the row with the biggest element for column (j)
//     pivot = _nml_mat_absmaxr(U, j);
//     if (fabs(U->data[pivot][j]) < NML_MIN_COEF) {
//       NML_ERROR(CANNOT_LU_MATRIX_DEGENERATE);
//       return NULL;
//     }
//     if (pivot!=j) {
//       // Pivots LU and P accordingly to the rule
//       nml_mat_row_swap_r(U, j, pivot);
//       nml_mat_row_swap_r(L, j, pivot);
//       nml_mat_row_swap_r(P, j, pivot);
//       // Keep the number of permutations to easily calculate the
//       // determinant sign afterwards
//       num_permutations++;
//     }
//     for(i = j+1; i < U->num_rows; i++) {
//       mult = U->data[i][j] / U->data[j][j];
//       // Building the U upper rows
//       nml_mat_row_addrow_r(U, i, j, -mult);
//       // Store the multiplier in L
//       L->data[i][j] = mult;
//     }
//   }
//   nml_mat_diag_set(L, 1.0);
//
//   return nml_mat_lup_new(L, U, P, num_permutations);
// }


nml_mat_lup_solve :: proc ( m : ^Nml_mat ) -> ^Nml_mat_lup {
    if !m^.is_square {
        // NML_FERROR(CANNOT_LU_MATRIX_SQUARE, m->num_rows, m->num_cols);
        
        str_error := fmt.aprintf( CANNOT_LU_MATRIX_SQUARE, m^.num_rows, m^.num_cols )
        cstr_error := strings.clone_to_cstring( str_error )
        delete( str_error )
    
        m_mnl_ferror( cstr_error )

        return nil
    }
    L : ^Nml_mat = nml_mat_new( m^.num_rows, m^.num_rows )
    U : ^Nml_mat = nml_mat_cp( m )
    P : ^Nml_mat = nml_mat_eye( m^.num_rows )

    // j, i  : uint 
    pivot : int
    num_permutations : uint = 0
    mult  : f64

    for j : uint = 0; j < U^.num_cols; j += 1 {
        // Retrieves the row with the biggest element for column (j)
        pivot = _nml_mat_absmaxr( U, j )
        if math.abs( U^.data[ pivot ][ j ]) < NML_MIN_COEF {
            // NML_ERROR(CANNOT_LU_MATRIX_DEGENERATE);            
            m_mnl_error( CANNOT_LU_MATRIX_DEGENERATE )
            return nil
        }
        if uint( pivot ) != j {
            // Pivots LU and P accordingly to the rule
            nml_mat_row_swap_r( U, j, uint( pivot ) )
            nml_mat_row_swap_r( L, j, uint( pivot ) )
            nml_mat_row_swap_r( P, j, uint( pivot ) )
            // Keep the number of permutations to easily calculate the
            // determinant sign afterwards
            num_permutations += 1
        }
        for i : uint = j+1; i < U^.num_rows; i += 1 {
            mult = U^.data[ i ][ j ] / U^.data[ j ][ j ]
            // Building the U upper rows
            nml_mat_row_addrow_r( U, i, j, -mult )
            // Store the multiplier in L
            L^.data[ i ][ j ] = mult
        }
    }
    nml_mat_diag_set( L, 1.0 )

    return nml_mat_lup_new( L, U, P, num_permutations )
}

// // After the LU(P) factorisation the determinant can be easily calculated
// // by multiplying the main diagonal of matrix U with the sign.
// // the sign is -1 if the number of permutations is odd
// // the sign is +1 if the number of permutations is even
// double nml_mat_det(nml_mat_lup* lup) {
//   int k;
//   int sign = (lup->num_permutations%2==0) ? 1 : -1;
//   nml_mat *U = lup->U;
//   double product = 1.0;
//   for(k = 0; k < U->num_rows; k++) {
//     product *= U->data[k][k];
//   }
//   return product * sign;
// }


// After the LU(P) factorisation the determinant can be easily calculated
// by multiplying the main diagonal of matrix U with the sign.
// the sign is -1 if the number of permutations is odd
// the sign is +1 if the number of permutations is even
nml_mat_det :: proc ( lup :^Nml_mat_lup ) -> f64 {
    sign : int = (lup^.num_permutations % 2 == 0 ) ? 1 : -1
    U : ^Nml_mat = lup^.U
    product : f64 = 1.0
    for k : uint  = 0; k < U^.num_rows; k += 1 {
        product *= U^.data[ k ][ k ]
    }
    return product * f64( sign )
}


// // Returns LU matrix from a LUP structure
// nml_mat *nml_mat_lu_get(nml_mat_lup* lup) {
//   nml_mat *r = nml_mat_cp(lup->U);
//   // Copy L (without first diagonal in result)
//   int i, j;
//   for(i = 1; i < lup->L->num_rows; i++) {
//     for(j = 0; j < i; j++) {
//       r->data[i][j] = lup->L->data[i][j];
//     }
//   }
//   return r;
// }

// Returns LU matrix from a LUP structure
nml_mat_lu_get :: proc ( lup : ^Nml_mat_lup ) -> ^Nml_mat  {
    r : ^Nml_mat = nml_mat_cp( lup^.U )
    // Copy L (without first diagonal in result)
    for i : uint = 1; i < lup^.L^.num_rows; i += 1 {
        for j : uint = 0; j < i; j += 1 {
            r^.data[ i ][ j ] = lup^.L^.data[ i ][ j ]
        }
    }
    return r
}

// *****************************************************************************
//
// Solving linear systems of equations
//
// *****************************************************************************

// // Forward substitution algorithm
// // Solves the linear system L * x = b
// //
// // L is lower triangular matrix of size NxN
// // B is column matrix of size Nx1
// // x is the solution column matrix of size Nx1
// //
// // Note: In case L is not a lower triangular matrix, the algorithm will try to
// // select only the lower triangular part of the matrix L and solve the system
// // with it.
// //
// // Note: In case any of the diagonal elements (L[i][i]) are 0 the system cannot
// // be solved
// //
// // Note: This function is usually used with an L matrix from a LU decomposition
// nml_mat *nml_ls_solvefwd(nml_mat *L, nml_mat *b) {
//   nml_mat* x = nml_mat_new(L->num_cols, 1);
//   int i,j;
//   double tmp;
//   for(i = 0; i < L->num_cols; i++) {
//     tmp = b->data[i][0];
//     for(j = 0; j < i ; j++) {
//       tmp -= L->data[i][j] * x->data[j][0];
//     }
//     x->data[i][0] = tmp / L->data[i][i];
//   }
//   return x;
// }


// Forward substitution algorithm
// Solves the linear system L * x = b
//
// L is lower triangular matrix of size NxN
// B is column matrix of size Nx1
// x is the solution column matrix of size Nx1
//
// Note: In case L is not a lower triangular matrix, the algorithm will try to
// select only the lower triangular part of the matrix L and solve the system
// with it.
//
// Note: In case any of the diagonal elements (L[i][i]) are 0 the system cannot
// be solved
//
// Note: This function is usually used with an L matrix from a LU decomposition
nml_ls_solvefwd :: proc ( L : ^Nml_mat, b : ^Nml_mat ) -> ^Nml_mat {
    x : ^Nml_mat = nml_mat_new( L^.num_cols, 1 )
    tmp : f64
    for i : uint = 0; i < L^.num_cols; i += 1 {
        tmp = b^.data[ i ][ 0 ]
        for j : uint = 0; j < i ; j += 1 {
            tmp -= L^.data[ i ][ j ] * x^.data[ j ][ 0 ]
        }
        x^.data[ i ][ 0 ] = tmp / L^.data[ i ][ i ]
    }
    return x;
}

// // Back substition algorithm
// // Solves the linear system U *x = b
// //
// // U is an upper triangular matrix of size NxN
// // B is a column matrix of size Nx1
// // x is the solution column matrix of size Nx1
// //
// // Note in case U is not an upper triangular matrix, the algorithm will try to
// // select only the upper triangular part of the matrix U and solve the system
// // with it
// //
// // Note: In case any of the diagonal elements (U[i][i]) are 0 the system cannot
// // be solved
// nml_mat *nml_ls_solvebck(nml_mat *U, nml_mat *b) {
//   nml_mat *x = nml_mat_new(U->num_cols, 1);
//   int i = U->num_cols, j;
//   double tmp;
//   while(i-->0) {
//     tmp = b->data[i][0];
//     for(j = i; j < U->num_cols; j++) {
//       tmp -= U->data[i][j] * x->data[j][0];
//     }
//     x->data[i][0] = tmp / U->data[i][i];
//   }
//   return x;
// }


// Back substition algorithm
// Solves the linear system U *x = b
//
// U is an upper triangular matrix of size NxN
// B is a column matrix of size Nx1
// x is the solution column matrix of size Nx1
//
// Note in case U is not an upper triangular matrix, the algorithm will try to
// select only the upper triangular part of the matrix U and solve the system
// with it
//
// Note: In case any of the diagonal elements (U[i][i]) are 0 the system cannot
// be solved
nml_ls_solvebck :: proc ( U : ^Nml_mat, b : ^Nml_mat ) -> ^Nml_mat {
    x : ^Nml_mat = nml_mat_new( U^.num_cols, 1 )
    i : uint = U^.num_cols
    i -= 1
    for  i > 0 {
        tmp : f64 = b^.data[ i ][ 0 ]
        for j : uint = i; j < U^.num_cols; j += 1 {
            tmp -= U^.data[ i ][ j ] * x^.data[ j ][ 0 ]
        }
        x^.data[ i ][ 0 ] = tmp / U^.data[ i ][ i ]
        i -= 1
    }
    return x
}

// // A[n][n] is a square matrix
// // m contains matrices L, U, P for A[n][n] so that P*A = L*U
// //
// // The linear system is:
// // A*x=b  =>  P*A*x = P*b  =>  L*U*x = P*b  =>
// // (where b is a matrix[n][1], and x is a matrix[n][1])
// //
// // if y = U*x , we solve two systems:
// //    L * y = P b (forward substition)
// //    U * x = y (backward substition)
// //
// // We obtain and return x
// nml_mat *nml_ls_solve(nml_mat_lup *lu, nml_mat* b) {
//   if (lu->U->num_rows != b->num_rows || b->num_cols != 1) {
//     NML_FERROR(CANNOT_SOLVE_LIN_SYS_INVALID_B,
//       b->num_rows,
//       b->num_cols,
//       lu->U->num_rows,
//       1);
//       return NULL;
//   }
//   nml_mat *Pb = nml_mat_dot(lu->P, b);
//
//   // We solve L*y = P*b using forward substition
//   nml_mat *y = nml_ls_solvefwd(lu->L, Pb);
//
//   // We solve U*x=y
//   nml_mat *x = nml_ls_solvebck(lu->U, y);
//
//   nml_mat_free(y);
//   nml_mat_free(Pb);
//   return x;
// }


// A[n][n] is a square matrix
// m contains matrices L, U, P for A[n][n] so that P*A = L*U
//
// The linear system is:
// A*x=b  =>  P*A*x = P*b  =>  L*U*x = P*b  =>
// (where b is a matrix[n][1], and x is a matrix[n][1])
//
// if y = U*x , we solve two systems:
//    L * y = P b (forward substition)
//    U * x = y (backward substition)
//
// We obtain and return x
nml_ls_solve :: proc ( lu : ^Nml_mat_lup, b : ^Nml_mat ) -> ^Nml_mat {
    if lu^.U^.num_rows != b^.num_rows || b^.num_cols != 1 {
        // NML_FERROR(CANNOT_SOLVE_LIN_SYS_INVALID_B,
        //     b->num_rows,
        //     b->num_cols,
        //     lu->U->num_rows,
        //     1);

        str_error := fmt.aprintf( CANNOT_SOLVE_LIN_SYS_INVALID_B,
                           b^.num_rows,
                                  b^.num_cols,
                                  lu^.U^.num_rows,
                                  1 )
        cstr_error := strings.clone_to_cstring( str_error )
        delete( str_error )
    
        m_mnl_ferror( cstr_error )

        return nil
    }
    Pb : ^Nml_mat = nml_mat_dot( lu^.P, b )

    // We solve L*y = P*b using forward substition
    y : ^Nml_mat = nml_ls_solvefwd( lu^.L, Pb )

    // We solve U*x=y
    x : ^Nml_mat = nml_ls_solvebck( lu^.U, y )

    nml_mat_free( y )
    nml_mat_free( Pb )
    return x
}

// // Calculates the inverse of a matrix
// nml_mat *nml_mat_inv(nml_mat_lup *lup) {
//   unsigned n = lup->L->num_cols;
//   nml_mat *r = nml_mat_sqr(n);
//   nml_mat *I = nml_mat_eye(lup->U->num_rows);
//   nml_mat *invx;
//   nml_mat *Ix;
//   int i,j;
//   for(j =0; j < n; j++) {
//     Ix = nml_mat_col_get(I, j);
//     invx = nml_ls_solve(lup, Ix);
//     for(i = 0; i < invx->num_rows; i++) {
//       r->data[i][j] = invx->data[i][0];
//     }
//     nml_mat_free(invx);
//     nml_mat_free(Ix);
//   }
//   nml_mat_free(I);
//   return r;
// }

// Calculates the inverse of a matrix
nml_mat_inv :: proc ( lup : ^Nml_mat_lup ) -> ^Nml_mat {
    n : uint = lup^.L^.num_cols
    r : ^Nml_mat = nml_mat_sqr( n )
    I : ^Nml_mat = nml_mat_eye( lup^.U^.num_rows )
    invx : ^Nml_mat
    Ix   : ^Nml_mat
    for j : uint = 0; j < n; j += 1 {
        Ix = nml_mat_col_get( I, j )
        invx = nml_ls_solve(lup, Ix )
        for i : uint = 0; i < invx^.num_rows; i += 1 {
            r^.data[ i ][ j ] = invx^.data[ i ][ 0 ]
        }
        nml_mat_free( invx )
        nml_mat_free( Ix )
    }
    nml_mat_free( I )
    return r
}

// *****************************************************************************
//
// QR Decomposition
//
// *****************************************************************************




// // Useful for QR decomposition
// // Represents the (dot) product of two vectors:
// // vector1 = m1col column from m1
// // vector2 = m2col column from m2
// double nml_vect_dot(nml_mat *m1, unsigned int m1col, nml_mat *m2, unsigned m2col) {
//   if (m1->num_rows!=m2->num_rows) {
//     NML_FERROR(CANNOT_VECT_DOT_DIMENSIONS, m1->num_rows, m2->num_rows);
//   }
//   if (m1col >= m1->num_cols) {
//     NML_FERROR(CANNOT_GET_COLUMN, m1col, m1->num_cols);
//   }
//   if (m2col >= m2->num_cols) {
//     NML_FERROR(CANNOT_GET_COLUMN, m2col, m2->num_cols);
//   }
//   int i;
//   double dot = 0.0;
//   for(i = 0; i < m1->num_rows; i++) {
//     dot += m1->data[i][m1col] * m2->data[i][m2col];
//   }
//   return dot;
// }


// Useful for QR decomposition
// Represents the (dot) product of two vectors:
// vector1 = m1col column from m1
// vector2 = m2col column from m2
nml_vect_dot :: proc ( m1 : ^Nml_mat, m1col : uint, m2 : ^Nml_mat, m2col : uint ) -> f64 {
    if m1^.num_rows != m2^.num_rows {
        // NML_FERROR(CANNOT_VECT_DOT_DIMENSIONS, m1->num_rows, m2->num_rows);

        str_error := fmt.aprintf( CANNOT_VECT_DOT_DIMENSIONS, m1^.num_rows, m2^.num_rows )
        cstr_error := strings.clone_to_cstring( str_error )
        delete( str_error )

        m_mnl_ferror( cstr_error )
    }
    if m1col >= m1^.num_cols {
        // NML_FERROR(CANNOT_GET_COLUMN, m1col, m1->num_cols);

        str_error := fmt.aprintf( CANNOT_GET_COLUMN, m1col, m1^.num_cols )
        cstr_error := strings.clone_to_cstring( str_error )
        delete( str_error )

        m_mnl_ferror( cstr_error )
    }
    if m2col >= m2^.num_cols {
        // NML_FERROR(CANNOT_GET_COLUMN, m2col, m2->num_cols);

        str_error := fmt.aprintf( CANNOT_GET_COLUMN, m2col, m2^.num_cols )
        cstr_error := strings.clone_to_cstring( str_error )
        delete( str_error )

        m_mnl_ferror( cstr_error )
    }
    dot : f64 = 0.0;
    for i : uint = 0; i < m1^.num_rows; i += 1 {
        dot += m1^.data[ i ][ m1col ] * m2^.data[ i ][ m2col ]
    }
    return dot
}

// // Calculates the l2 norm for a colum in the matrix
// double nml_mat_col_l2norm(nml_mat *m, unsigned int col) {
//   if (col >= m->num_cols) {
//     NML_FERROR(CANNOT_COLUMN_L2NORM, col, m->num_cols);
//   }
//   double doublesum = 0.0;
//   int i;
//   for(i = 0; i < m->num_rows; i++) {
//     doublesum += (m->data[i][col]*m->data[i][col]);
//   }
//   return sqrt(doublesum);
// }


// Calculates the l2 norm for a colum in the matrix
nml_mat_col_l2norm :: proc ( m : ^Nml_mat, col : uint ) -> f64 {
    if col >= m^.num_cols {
        // NML_FERROR(CANNOT_COLUMN_L2NORM, col, m->num_cols);

        str_error := fmt.aprintf( CANNOT_COLUMN_L2NORM, col, m^.num_cols )
        cstr_error := strings.clone_to_cstring( str_error )
        delete( str_error )

        m_mnl_ferror( cstr_error )
    }
    doublesum : f64 = 0.0
    for i : uint = 0; i < m^.num_rows; i += 1 {
        doublesum += ( m^.data[ i ][ col ] * m^.data[ i ][ col ] )
    }
    return math.sqrt( doublesum )
}

// // Calculates the l2norm for each column
// // Keeps results into 1 row matrix
// nml_mat *nml_mat_l2norm(nml_mat *m) {
//   int i, j;
//   nml_mat *r = nml_mat_new(1, m->num_cols);
//   double square_sum;
//   for(j = 0; j < m->num_cols; j++) {
//     square_sum = 0.0;
//     for(i = 0; i < m->num_rows; i++) {
//       square_sum+=m->data[i][j]*m->data[i][j];
//     }
//     r->data[0][j] = sqrt(square_sum);
//   }
//   return r;
// }

// Calculates the l2norm for each column
// Keeps results into 1 row matrix
nml_mat_l2norm :: proc ( m : ^Nml_mat ) -> ^Nml_mat {
    r : ^Nml_mat = nml_mat_new( 1, m^.num_cols )
    for j : uint = 0; j < m^.num_cols; j += 1 {
        square_sum : f64 = 0.0
        for i : uint = 0; i < m^.num_rows; i += 1 {
            square_sum += m^.data[ i ][ j ] * m^.data[ i ][ j ]
        }
        r^.data[ 0 ][ j ] = math.sqrt( square_sum )
    }
    return r;
}

// nml_mat *nml_mat_normalize(nml_mat *m) {
//   nml_mat *r = nml_mat_cp(m);
//   if (!nml_mat_normalize_r(r)) {
//     nml_mat_free(r);
//     return NULL;
//   }
//   return r;
// }

nml_mat_normalize :: proc ( m : ^Nml_mat ) -> ^Nml_mat {
    r : ^Nml_mat = nml_mat_cp( m )
    if ! nml_mat_normalize_r( r ) {
        nml_mat_free( r )
        return nil
    }
    return r
}

// int nml_mat_normalize_r(nml_mat *m) {
//   nml_mat *l2norms = nml_mat_l2norm(m);
//   int j;
//   for(j = 0; j < m->num_cols; j++) {
//     if (l2norms->data[0][j] < NML_MIN_COEF) {
//       NML_FERROR(VECTOR_J_DEGENERATE, j);
//       nml_mat_free(l2norms);
//       return 0;
//     }
//     nml_mat_col_mult_r(m, j, 1/l2norms->data[0][j]);
//   }
//   nml_mat_free(l2norms);
//   return 1;
// }

nml_mat_normalize_r :: proc ( m : ^Nml_mat ) -> bool {
    l2norms : ^Nml_mat = nml_mat_l2norm( m )
    for j : uint  = 0; j < m^.num_cols; j += 1 {
        if l2norms^.data[ 0 ][ j ] < NML_MIN_COEF {
            // NML_FERROR(VECTOR_J_DEGENERATE, j);
            
            str_error := fmt.aprintf( VECTOR_J_DEGENERATE, j )
            cstr_error := strings.clone_to_cstring( str_error )
            delete( str_error )
    
            m_mnl_ferror( cstr_error )            
            
            nml_mat_free( l2norms )
            return false
        }
        nml_mat_col_mult_r( m, j, 1 / l2norms^.data[ 0 ][ j ] )
    }
    nml_mat_free( l2norms )
    return true
}

// nml_mat_qr *nml_mat_qr_new() {
//   nml_mat_qr *qr = malloc(sizeof(*qr));
//   NP_CHECK(qr);
//   return qr;
// }

nml_mat_qr_new :: proc () -> ^Nml_mat_qr {
    // nml_mat_qr *qr = malloc(sizeof(*qr));
    qr : ^Nml_mat_qr = new( Nml_mat_qr )
    // NP_CHECK(qr);
    m_np_check( qr )
    return qr
}

// void nml_mat_qr_free(nml_mat_qr *qr) {
//   nml_mat_free(qr->Q);
//   nml_mat_free(qr->R);
//   free(qr);
// }

nml_mat_qr_free :: proc ( qr : ^Nml_mat_qr ) {
    nml_mat_free( qr^.Q )
    nml_mat_free( qr^.R )
    free( qr )
}

// // M = QR
// nml_mat_qr *nml_mat_qr_solve(nml_mat *m) {

//   nml_mat_qr *qr = nml_mat_qr_new();
//   nml_mat *Q = nml_mat_cp(m);
//   nml_mat *R = nml_mat_new(m->num_rows, m->num_cols);

//   int j, k;
//   double l2norm;
//   double rkj;
//   nml_mat *aj;
//   nml_mat *qk;
//   for(j=0; j < m->num_cols; j++) {    
//     rkj = 0.0;
//     aj = nml_mat_col_get(m, j);
//     for(k = 0; k < j; k++) {
//        rkj = nml_vect_dot(m, j, Q, k);
//        R->data[k][j] = rkj;
//        qk = nml_mat_col_get(Q, k);
//        nml_mat_col_mult_r(qk, 0, rkj);
//        nml_mat_sub_r(aj, qk);
//        nml_mat_free(qk);
//     }
//     for(k = 0; k < Q->num_rows; k++) {
//       Q->data[k][j] = aj->data[k][0];
//     }
//     l2norm = nml_mat_col_l2norm(Q, j);
//     nml_mat_col_mult_r(Q, j, 1/l2norm);
//     R->data[j][j] = l2norm;
//     nml_mat_free(aj);
//   }
//   qr->Q = Q;
//   qr->R = R;
//   return qr;
// }


// M = QR
nml_mat_qr_solve :: proc ( m : ^Nml_mat ) -> ^Nml_mat_qr {

    qr : ^Nml_mat_qr = nml_mat_qr_new( )
    Q  : ^Nml_mat = nml_mat_cp( m )
    R  : ^Nml_mat = nml_mat_new( m^.num_rows, m^.num_cols )

    l2norm : f64
    rkj    : f64
    aj : ^Nml_mat
    qk : ^Nml_mat
    for j : uint = 0; j < m^.num_cols; j += 1 {    
        rkj = 0.0
        aj = nml_mat_col_get( m, j )
        for k : uint = 0; k < j; k += 1 {
            rkj = nml_vect_dot( m, j, Q, k )
            R^.data[ k ][ j ] = rkj
            qk = nml_mat_col_get( Q, k )
            nml_mat_col_mult_r( qk, 0, rkj )
            nml_mat_sub_r( aj, qk )
            nml_mat_free( qk )
        }
        for k : uint = 0; k < Q^.num_rows; k += 1 {
            Q^.data[ k ][ j ] = aj^.data[ k ][ 0 ]
        }
        l2norm = nml_mat_col_l2norm( Q, j )
        nml_mat_col_mult_r( Q, j, 1 / l2norm )
        R^.data[ j ][ j ] = l2norm
        nml_mat_free( aj )
    }
    qr^.Q = Q
    qr^.R = R
    return qr
}




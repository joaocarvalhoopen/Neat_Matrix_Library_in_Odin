package ex_backward_substitution

import m "./../../nml"

import "core:fmt"
import "core:c/libc"

import "core:os"

main :: proc ( ) {
    // path_to_matrix_s_files = "./../../examples_odin/data/matrix10_upper_triangular.data"
    path_to_matrix_s_files : cstring = "./examples_odin/data/matrix10_upper_triangular.data"
    input : ^libc.FILE = libc.fopen( path_to_matrix_s_files, "r")
    if input == nil {
        fmt.printf( "Error: could not open file\n" )
        os.exit( 1 )
    }

    A : ^m.Nml_mat = m.nml_mat_fromfilef( input )
    B : ^m.Nml_mat = m.nml_mat_fromfilef( input )
  
    x : ^m.Nml_mat = m.nml_ls_solvebck( A, B )

    m.nml_mat_print( A )
    m.nml_mat_print( B )
    m.nml_mat_print( x )

    m.nml_mat_free( A )
    m.nml_mat_free( B )
    m.nml_mat_free( x )

    libc.fclose( input )
}

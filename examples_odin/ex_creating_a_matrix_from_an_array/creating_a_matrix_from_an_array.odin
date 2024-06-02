package ex_creating_a_matrix_from_an_array

import m "./../../nml"

import "core:fmt"
import "core:c/libc"

import "core:os"
import "core:mem"

main :: proc ( ) {
    array := [?]f64 { 1.0, 0.2, 3.0, 4.0, 5.0, 3.1 }
    my : ^m.Nml_mat
    
    // 3 rows, 2 columns
    // read exactly 6 numbers from array[6]
    my = m.nml_mat_from( 3, 2, 6, & array[ 0 ] )
    m.nml_mat_print( my )
    m.nml_mat_free( my )

    // 4 rows, 2 columns
    // read exactly 3 numbers from array[6]
    my = m.nml_mat_from( 4, 2, 3, & array[ 0 ] )
    m.nml_mat_print( my )
    m.nml_mat_free( my )
}



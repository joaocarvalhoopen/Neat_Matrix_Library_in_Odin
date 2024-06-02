package ex_creating_a_matrix_from_user_input

import m "./../../nml"

import "core:fmt"
import "core:c/libc"

import "core:os"

main :: proc ( ) {
    from_file2 : ^m.Nml_mat = m.nml_mat_fromfilef( libc.stdin )
    m.nml_mat_print( from_file2 )
    m.nml_mat_free( from_file2 )
}
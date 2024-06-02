package ex_creating_a_matrix_from_file

import m "./../../nml"

import "core:fmt"
import "core:c/libc"

import "core:os"
import "core:mem"

main :: proc ( ) {
    f : cstring = "./examples_odin/data/matrix1.data";
    from_file : ^m.Nml_mat = m.nml_mat_fromfile( f )
    m.nml_mat_print( from_file )
    m.nml_mat_free( from_file )

    // Or if the file is already opened

    m_file : ^libc.FILE = libc.fopen( "./examples_odin/data/matrix2.data", "r" )
    from_file2 : ^m.Nml_mat = m.nml_mat_fromfilef( m_file )
    m.nml_mat_print( from_file2 )
    m.nml_mat_free( from_file2 )
    libc.fclose( m_file )
}
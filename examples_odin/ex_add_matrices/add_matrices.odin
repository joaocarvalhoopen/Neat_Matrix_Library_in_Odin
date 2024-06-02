package examples_odin

import m "./../../nml"

import "core:fmt"

main :: proc ( ) {
    m1 : ^m.Nml_mat = m.nml_mat_sqr_rnd( 4, 0.0, 10.0 )
    m2 : ^m.Nml_mat = m.nml_mat_sqr_rnd( 4, 0.0, 10.0 )
    fmt.printf( "m1=\n" )
    m.nml_mat_print( m1 )
    
    fmt.printf( "m2=\n" )
    m.nml_mat_print( m2 )

    // Add the matrices to, result is kept in m3
    // m1 and m2 remain unchanged
    m3 : ^m.Nml_mat = m.nml_mat_add( m1, m2 )
    fmt.printf( "m3=\n" )
    m.nml_mat_print( m3 )

    // Add the matrices, the result is kept in m1
    // m1 is modified, m2 remains unchanged
    m.nml_mat_add_r( m1, m2 )
    fmt.printf( "m1=\n" )
    m.nml_mat_print( m1 )

    m.nml_mat_free( m1 )
    m.nml_mat_free( m2 )
    m.nml_mat_free( m3 )
}
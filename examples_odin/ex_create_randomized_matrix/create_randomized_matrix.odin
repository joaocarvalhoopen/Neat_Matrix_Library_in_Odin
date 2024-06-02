package ex_create_randomized_matrix

import m "./../../nml"

import "core:fmt"
import "core:c/libc"

import "core:os"
import "core:mem"

import "core:math/rand"

main :: proc ( ) {
    // srand(time(NULL)); // Should be called once per program
    m.g_rand_generator = rand.create( 0 ) // Should be called once per program
    m1 : ^m.Nml_mat = m.nml_mat_rnd( 5, 5, -10.0, 10.0 )
    m.nml_mat_print( m1 )
    m.nml_mat_free( m1 )
}
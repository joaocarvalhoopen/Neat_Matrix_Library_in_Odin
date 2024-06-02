package ex_concatenation_matrices

import m "./../../nml"

import "core:fmt"
import "core:c/libc"

import "core:os"
import "core:mem"

main :: proc ( ) {
    I    : ^m.Nml_mat = m.nml_mat_eye( 3 )
    Ix2  : ^m.Nml_mat = m.nml_mat_smult( I, 2.0 )
    rndm : ^m.Nml_mat = m.nml_mat_rnd( 3, 4, 1.0, 5.0 )
  
    // nml_mat **ms = malloc(sizeof(*ms) * 2);
    ms : [^]^m.Nml_mat = make( [^]^m.Nml_mat, 2 )
    ms[ 0 ] = I
    ms[ 1 ] = Ix2
    
    concats1 : ^m.Nml_mat = m.nml_mat_cath( 2, ms )
  
    ms[ 0 ] = concats1
    ms[ 1 ] = rndm
  
    concats2 : ^m.Nml_mat = m.nml_mat_cath( 2, ms )
  
    fmt.printf( "\nConcatenate horizontally\n" )
    fmt.printf( "I=\n" )
    m.nml_mat_print( I )
    fmt.printf( "Ix2=\n" )
    m.nml_mat_print( Ix2 )
    fmt.printf( "rndm=\n" )
    m.nml_mat_print( rndm )
    fmt.printf( "concats1=\n" )
    m.nml_mat_print( concats1 )
    fmt.printf( "concats2=\n" )
    m.nml_mat_print( concats2 )
  
    free( ms )
    m.nml_mat_free( I )
    m.nml_mat_free( Ix2 )
    m.nml_mat_free( concats1 )
    m.nml_mat_free( concats2 )
    m.nml_mat_free( rndm )
  
    // -------------------------------------
    // Vertical concatenation
    // -------------------------------------
  
    A : ^m.Nml_mat = m.nml_mat_rnd( 3, 4, 1.0, 4.0 )
    B : ^m.Nml_mat = m.nml_mat_rnd( 5, 4, 10.0, 20.0 )
    C : ^m.Nml_mat = m.nml_mat_eye( 4 )
  
    // nml_mat **ABarr = malloc( sizeof( *ABarr ) * 2 );
    ABarr : [^]^m.Nml_mat = make( [^]^m.Nml_mat, 2 )

    ABarr[ 0 ] = A
    ABarr[ 1 ] = B
    ABCat : ^m.Nml_mat = m.nml_mat_catv( 2, ABarr )
  
    fmt.printf( "\nA=\n" )
    m.nml_mat_print( A )
    fmt.printf( "\nB=\n" )
    m.nml_mat_print( B )
    fmt.printf( "\nC=\n" )
    m.nml_mat_print( C )
    fmt.printf( "\nA concat B =\n" )
    m.nml_mat_print( ABCat )
  
    free( ABarr )
    m.nml_mat_free( A )
    m.nml_mat_free( B )
    m.nml_mat_free( C ) 
}
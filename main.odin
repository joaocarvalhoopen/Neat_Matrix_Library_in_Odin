package main

import "./nml"
import "core:fmt"
import "core:os"


/**********************************************************************
*
* IMPORTNAT NOTE:
*
* This file is only a scratchpad for testing the NML library!
* See the NML diretory for the actual implementation and 
* the NML documentation in the examples dir and in the Blog Article
* from the Original NML author.
**********************************************************************/



T_Matrix_Struct :: struct {
    rows : int,
    cols : int,
    data : [^][^]f64,
}

test_matrix_create :: proc ( rows: int, cols : int ) -> ^T_Matrix_Struct {   
    m : ^T_Matrix_Struct = new( T_Matrix_Struct )
    m^.rows = rows
    m^.cols = cols
    m^.data =  make( [ ^ ][ ^ ]f64, rows )
    if m^.data == nil {
        fmt.println( "Error allocating memory in teste_create() 1." )
        os.exit(1)
    }
    for i : int = 0; i < m.rows; i += 1 {
        m^.data[ i ] = make( [ ^ ]f64, cols )
        if m^.data[ i ] == nil {
            fmt.println( "Error allocating memory in teste_create() 2." )
            os.exit(1)
        }   
    }
    fmt.println( m )
    return m
}

test_matrix_set :: proc ( m : ^T_Matrix_Struct ) {
    val : f64 = 1
    for i : int = 0; i < m.rows; i += 1 {
        for j : int = 0; j < m.cols; j += 1 {
            m^.data[i][j] = val
            val += 1 
        }
    }
}

test_matrix_print :: proc ( m : ^T_Matrix_Struct ) {
    for i : int = 0; i < m.rows; i += 1 {
        for j : int = 0; j < m.cols; j += 1 {
            fmt.print( m^.data[i][j] )
            fmt.print( " " )
        }
        fmt.println( "" )
    }
}

test_run :: proc ( ) {
    matrix_m : ^T_Matrix_Struct = test_matrix_create( 2,  3)
    test_matrix_set( matrix_m )
    test_matrix_print( matrix_m )
}


main :: proc ( ) {
    fmt.println("Begin NML library ...")

    test_run()

    bla_ptr :^int
    nml.m_np_check( bla_ptr )


}
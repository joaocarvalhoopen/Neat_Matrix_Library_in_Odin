all:
	odin build . -out:nml.exe --debug
opti:
	odin build . -out:nml.exe -o:speed

clean:
	rm nml.exe

examples: ./examples_odin/ex_add_matrices/add_matrices.odin ./examples_odin/ex_backward_substitution/backward_substitution.odin ./examples_odin/ex_concatenate_matrices/concatenate_matrices.odin ./examples_odin/ex_create_randomized_matrix/create_randomized_matrix.odin ./examples_odin/ex_creating_a_matrix_from_an_array/creating_a_matrix_from_an_array.odin ./examples_odin/ex_creating_a_matrix_from_file/creating_a_matrix_from_file.odin ./examples_odin/ex_creating_a_matrix_from_user_input/creating_a_matrix_from_user_input.odin
	odin build ./examples_odin/ex_add_matrices/add_matrices.odin -file -out:./examples_odin/ex_add_matrices/add_matrices.exe --debug 
	odin build ./examples_odin/ex_backward_substitution/backward_substitution.odin -file -out:./examples_odin/ex_backward_substitution/backward_substitution.exe --debug 
	odin build ./examples_odin/ex_concatenate_matrices/concatenate_matrices.odin -file -out:./examples_odin/ex_concatenate_matrices/concatenate_matrices.exe --debug
	odin build ./examples_odin/ex_create_randomized_matrix/create_randomized_matrix.odin -file -out:./examples_odin/ex_create_randomized_matrix/create_randomized_matrix.exe --debug
	odin build ./examples_odin/ex_creating_a_matrix_from_an_array/creating_a_matrix_from_an_array.odin -file -out:./examples_odin/ex_creating_a_matrix_from_an_array/creating_a_matrix_from_an_array.exe --debug
	odin build ./examples_odin/ex_creating_a_matrix_from_file/creating_a_matrix_from_file.odin -file -out:./examples_odin/ex_creating_a_matrix_from_file/creating_a_matrix_from_file.exe --debug	
	odin build ./examples_odin/ex_creating_a_matrix_from_user_input/creating_a_matrix_from_user_input.odin -file -out:./examples_odin/ex_creating_a_matrix_from_user_input/creating_a_matrix_from_user_input.exe --debug



clean_examples:
	rm ./examples_odin/ex_add_matrices/add_matrices.exe
	rm ./examples_odin/ex_backward_substitution/backward_substitution.exe
	rm ./examples_odin/ex_concatenate_matrices/concatenate_matrices.exe
	rm ./examples_odin/ex_create_randomized_matrix/create_randomized_matrix.exe
	rm ./examples_odin/ex_creating_a_matrix_from_an_array/creating_a_matrix_from_an_array.exe
	rm ./examples_odin/ex_creating_a_matrix_from_file/creating_a_matrix_from_file.exe
	rm ./examples_odin/ex_creating_a_matrix_from_user_input/creating_a_matrix_from_user_input.exe

run_examples:
	./examples_odin/ex_add_matrices/add_matrices.exe
#	./examples_odin/ex_backward_substitution/backward_substitution.exe
#	./examples_odin/ex_concatenate_matrices/concatenate_matrices.exe
#	./examples_odin/ex_create_randomized_matrix/create_randomized_matrix.exe
#	./examples_odin/ex_creating_a_matrix_from_an_array/creating_a_matrix_from_an_array.exe
#	./examples_odin/ex_creating_a_matrix_from_file/creating_a_matrix_from_file.exe
#	./examples_odin/ex_creating_a_matrix_from_user_input/creating_a_matrix_from_user_input.exe
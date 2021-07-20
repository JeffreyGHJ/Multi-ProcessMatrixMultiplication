# Multi-ProcessMatrixMultiplication
Purpose:
	
	matrix.c performs matrix multiplication on two square matrices using multiple processes.
	The command line takes 3 filenames as arguments. The first 2 specify the input file names 
	of binary files that must both contain square matrices of the same size. The 3rd
	argument specifies the name of the output file for the resulting square matrix. This program
	requires that the first two lines of the matrix files be the dimensions of the matrix.
	(in this case both numbers will be the same because they must be square matrices)

To Compile:

	mpicc matrix.c -o matrix -lm

To Run:
	
	mpirun -np [num processes] ./matrix [input file 1] [input file 2] [output file]
	
		-[num processes] can be any integer 1 or geater
		-[input file 1] must be a properly formatted binary file containing a matrix
		-[input file 2] must be a properly formatted binary file containing a matrix
		-[output file] specifies the file name to store resulting product matrix

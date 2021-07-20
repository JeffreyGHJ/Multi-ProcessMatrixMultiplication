//matrix.c


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/time.h>
#include <mpi.h>
#include <math.h>

#define MASTER 0
#define NO_ERR 0
#define DATA_MSG 0
#define PROMPT_MSG 1
#define RESPONSE_MSG 2
#define MALLOC_ERROR -2

#define MIN(a,b)           ((a)<(b)?(a):(b))
#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)

typedef double datatype;

int genmat(int, char**);
int prtmat(char*);
int validate(int, char**);
double *my_malloc(int, int);
double read_row_striped_matrix(char*, double***, double**, int*, MPI_Comm);
void print_row_striped_matrix (double**, int, MPI_Comm);
void print_submatrix (double**, int, int);
void pass_blocks(MPI_Comm, double**, double**, double**, double***, double***, double***, int);
void create_row_striped_matrix(double***, double**, int, MPI_Comm);
void matrix_mult(double***, double***, double***, int, int, MPI_Comm);
void write_row_striped_matrix(char*, double**, int, MPI_Comm);

int main (int argc, char * argv[]) 
{
 	int numtasks;
	int taskId;
	int *n;
	double temp = 4.20;
	char *infileA;
	char *infileB;
	char *outfileC;
	
	//THE MATRICES ARE REALLY JUST 1-D ARRAYS
	double *matrixA;
	double *matrixB;
	double *matrixC;
	double **subsA;
	double **subsB;
	double **subsC;
	
	double       **lptr;
	double        *rptr; 
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskId);
	
	if (taskId == MASTER)
	{
		validate(argc, argv);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	infileA = argv[1];
	infileB = argv[2];
	outfileC = argv[3];

	MPI_Barrier(MPI_COMM_WORLD);
	
	//DISTRIBUTE MATRICES AMONGST PROCESSES
	read_row_striped_matrix(infileA, &subsA, &matrixA, n, MPI_COMM_WORLD);
	read_row_striped_matrix(infileB, &subsB, &matrixB, n, MPI_COMM_WORLD);
	create_row_striped_matrix(&subsC, &matrixC, *n, MPI_COMM_WORLD);
	
	MPI_Barrier(MPI_COMM_WORLD);
    
    if (!taskId)
    	printf("Input Matrix A:\n");
    
    print_row_striped_matrix(&matrixA, *n, MPI_COMM_WORLD);
    	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if (!taskId)
    	printf("Input Matrix B:\n");
	
	print_row_striped_matrix(&matrixB, *n, MPI_COMM_WORLD);
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	//ALL PROCESSES NOW PERFORM MATRIX MULT...
	pass_blocks(MPI_COMM_WORLD, &matrixA, &matrixB, &matrixC, &subsA, &subsB, &subsC, *n);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (!taskId)
    	printf("Result Matrix C:\n");
    	
    write_row_striped_matrix(outfileC, &matrixC, *n, MPI_COMM_WORLD);
    
    //print_row_striped_matrix(&matrixC, *n, MPI_COMM_WORLD);
    
    if(!taskId)
    	prtmat(outfileC);

	MPI_Finalize();
	return NO_ERR;
}


void write_row_striped_matrix(char* outfile, double** matrixC, int n, MPI_Comm comm)
{
	int p, id, i;
	FILE *foutptr;
	datatype *write_buffer;
	int max_block_size;
	int datum_size;
	int local_rows;

	double *recv_storage;
	
	MPI_Status status;
	
	MPI_Comm_size (comm, &p);
	MPI_Comm_rank (comm, &id);
	
	local_rows = BLOCK_SIZE(id, p, n);
	max_block_size = BLOCK_SIZE(p-1,p,n);
	datum_size = sizeof(double);
	
	if ( id == MASTER )
	{
		foutptr = fopen (outfile, "w");

		fwrite (&n, sizeof(int), 1, foutptr);
		fwrite (&n, sizeof(int), 1, foutptr);
		
		write_buffer = my_malloc (id, datum_size * n * max_block_size);

		fwrite ( *matrixC, datum_size, n * local_rows, foutptr);
		
		for ( i = 1; i < p; i++)
		{	
			local_rows = BLOCK_SIZE(i, p, n);
			
			MPI_Recv (write_buffer, local_rows * n, MPI_DOUBLE, 
				i, DATA_MSG, MPI_COMM_WORLD, &status);
				
			fwrite ( write_buffer, datum_size, n * local_rows, foutptr);
		}
			
		free(write_buffer);
		fclose (foutptr);
	}
	else
	{
		MPI_Send (*matrixC, local_rows * n, MPI_DOUBLE, MASTER, DATA_MSG, MPI_COMM_WORLD);
	}
}

void pass_blocks(MPI_Comm comm, double **matrixA, double **matrixB, double **matrixC,
	double*** subsA, double*** subsB, double*** subsC, int n)
{
	int p, id, i;
	int rotated_id, local_rows, recv_rows; 
	int rotation, dest, source; 
	int max_block, next_rows, next_id;
	double *recv_storage;
	double **subsR;
	int    max_block_size;
	int datum_size;
	int         prompt;          /* Dummy variable */
	MPI_Status status;
	
	MPI_Comm_size (comm, &p);
	MPI_Comm_rank (comm, &id);
	
	//DETERMINE SOURCE AND DESTIONATION PROCS
	if( id == 0 )
	{
		source = p - 1;
	}
	else
	{
		source = id - 1;
	}
	
	dest = (id + 1) % p;
	local_rows = BLOCK_SIZE(id,p,n);
	datum_size = sizeof(double);   
	max_block_size = BLOCK_SIZE(p-1,p,n);
		     
	recv_storage = calloc (1, max_block_size * n * datum_size);
		     
	subsR =(double **) my_malloc (id, max_block_size * datum_size);
		     
	subsR[0] = recv_storage;
		     
		     
     for (i = 1; i < max_block_size; i++) 
     {
        subsR[i] = subsR[i-1] + n * datum_size;
     }
     
	for (rotation = 0; rotation < p; rotation++)
	{
		rotated_id = (id + (rotation * ( p - 1 ))) % p;
		next_id = (id + ((rotation + 1) * ( p - 1))) % p;
		local_rows = BLOCK_SIZE(rotated_id, p, n);
		next_rows = BLOCK_SIZE(next_id, p, n);
	
		if ( id == 0 )
		{
			if( rotation == 0 )
			{
				matrix_mult(subsA, subsB, subsC, n, rotated_id, comm);
				
				MPI_Send (*matrixB, local_rows * n, MPI_DOUBLE, dest, DATA_MSG, MPI_COMM_WORLD);
					
				MPI_Recv (recv_storage, next_rows * n, MPI_DOUBLE, 
		        	source, DATA_MSG, MPI_COMM_WORLD, &status);
			
			}
			else
			{
				matrix_mult(subsA, &subsR, subsC, n, rotated_id, comm);
				
				MPI_Send(recv_storage, local_rows * n, MPI_DOUBLE, dest, DATA_MSG, MPI_COMM_WORLD);
				
				MPI_Recv (recv_storage, next_rows * n, MPI_DOUBLE, source, 
   					DATA_MSG, MPI_COMM_WORLD, &status);
			}
		}
		else
		{
			if( rotation == 0)
			{
				matrix_mult(subsA, subsB, subsC, n, rotated_id, comm);
			
				MPI_Send (*matrixB, local_rows * n, MPI_DOUBLE, dest, DATA_MSG, MPI_COMM_WORLD);
				
				MPI_Recv (recv_storage, next_rows * n, MPI_DOUBLE, 
		        	source, DATA_MSG, MPI_COMM_WORLD, &status);
			}
			else
			{
				matrix_mult(subsA, &subsR, subsC, n, rotated_id, comm);
			
				MPI_Send(recv_storage, local_rows * n, MPI_DOUBLE, dest, DATA_MSG, MPI_COMM_WORLD);
				
				MPI_Recv (recv_storage, next_rows * n, MPI_DOUBLE, source, 
   					DATA_MSG, MPI_COMM_WORLD, &status);
			}
		}
	}
}

void matrix_mult(double*** subsA, double*** subsB, double*** subsC, 
	int n, int rotated_id, MPI_Comm comm)
{
	int p, id, local_rows, rot_local_rows;
	int row, column, k;
	int row_offset, col_offset;
	
	double **subptrA = &(*subsA[0]);
	double **subptrB = &(*subsB[0]);
	double **subptrC = &(*subsC[0]);
	
	MPI_Comm_size (comm, &p);
	MPI_Comm_rank (comm, &id);
	
	local_rows = BLOCK_SIZE(id,p,n);
	col_offset = BLOCK_LOW(rotated_id, p, n);
	rot_local_rows = BLOCK_SIZE(rotated_id,p,n);
	
	for ( row = 0; row < local_rows ; row++ )
	{
		row_offset = n * row;
		
		for ( column = 0; column < n; column++ )
		{
			for ( k = 0; k < rot_local_rows; k++)
			{
				
				*( *(subptrC) + (row_offset) + column) += 
					*( *(subptrA) + (row_offset) + col_offset + k) *
					*( *(subptrB) + (k * n) + column);
					
			}
		}
	}
}


//PROBABLY COMPLETE
double read_row_striped_matrix (
	char        *infile,   /* IN - File name */
	double      ***subs,     /* OUT - 2D submatrix indices */
	double       **storage,  /* OUT - Submatrix stored here */
	int         *n,        /* OUT - Matrix cols */
	MPI_Comm     comm)     /* IN - Communicator */
{
	int datum_size = sizeof(double);
	int          i;
	int          id;           /* Process rank */
	FILE        *infileptr;    /* Input file pointer */
	int          local_rows;   /* Rows on this proc */
	double       **lptr;         /* Pointer into 'subs' */
	int          p;            /* Number of processes */
	double        *rptr;         /* Pointer into 'storage' */
	MPI_Status   status;       /* Result of receive */
	int          x;            /* Result of read */

	MPI_Comm_size (comm, &p);
	MPI_Comm_rank (comm, &id);

	/* Process p-1 opens file, reads size of matrix,
	and broadcasts matrix dimensions to other procs */
	if (id == (p-1)) 
	{
		infileptr = fopen (infile, "r");
		
	 	fread (n, sizeof(int), 1, infileptr);
	 	fread (n, sizeof(int), 1, infileptr);
	}

	MPI_Bcast (n, 1, MPI_INT, p-1, comm);

	local_rows = BLOCK_SIZE(id,p,*n);

	/* Dynamically allocate matrix. Allow double subscripting
	through 'a'. */
	*storage = my_malloc(id, local_rows * (*n) * datum_size);
	
	*subs = (double**) my_malloc (id, local_rows * sizeof(double*));
	
	lptr = &(*subs[0]);
	rptr = *storage;
	
	for (i = 0; i < local_rows; i++) 
	{
		*(lptr++)= rptr;
		rptr += *n * datum_size;
	}
	
	/* Process p-1 reads blocks of rows from file and
	sends each block to the correct destination process.
	The last block it keeps. */
	if (id == (p-1)) 
	{
		for (i = 0; i < p-1; i++)
		{
	 		x = fread (*storage, datum_size, BLOCK_SIZE(i,p,*n) * *n, infileptr);

	 		MPI_Send (*storage, BLOCK_SIZE(i,p,*n) * *n, MPI_DOUBLE, i, DATA_MSG, comm);
		}
	
		x = fread (*storage, datum_size, local_rows * *n, infileptr);
	
		fclose (infileptr);
	
	} 
	else
	{
		MPI_Recv (*storage, local_rows * *n, MPI_DOUBLE, p-1, DATA_MSG, comm, &status);
	}
	
}

void create_row_striped_matrix(double*** subs, double** storage, int n, MPI_Comm comm)
{
	int datum_size = sizeof(double);
	int          i;
	int          id;           /* Process rank */
	int          local_rows;   /* Rows on this proc */
	double       **lptr;         /* Pointer into 'subs' */
	int          p;            /* Number of processes */
	double       *rptr;         /* Pointer into 'storage' */
	MPI_Status   status;       /* Result of receive */

	MPI_Comm_size (comm, &p);
	MPI_Comm_rank (comm, &id);

	local_rows = BLOCK_SIZE(id,p,n);

	*storage = calloc( 1 , local_rows * n * datum_size);
	
	*subs = (double**) my_malloc (id, local_rows * sizeof(double*));
	
	lptr = &(*subs[0]);
	rptr = *storage;
	
	for (i = 0; i < local_rows; i++) 
	{
		*(lptr++)= rptr;
		rptr += n * datum_size;
	}
}

//PROBABLY COMPLETE
void print_row_striped_matrix (double **a, int n, MPI_Comm comm)
{
   MPI_Status  status;          /* Result of receive */
   double      *bstorage;        /* Elements received from
                                   another process */
   double      **b;               /* 2D array indexing into
                                   'bstorage' */
   int         datum_size;      /* Bytes per element */
   int         i;
   int         id;              /* Process rank */
   int         local_rows;      /* This proc's rows */
   int         max_block_size;  /* Most matrix rows held by
                                   any process */
   int         prompt;          /* Dummy variable */
   int         p;               /* Number of processes */

   MPI_Comm_rank (comm, &id);
   MPI_Comm_size (comm, &p);
   
   local_rows = BLOCK_SIZE(id,p,n);
   
	if (!id) //MASTER TASK
	{	
     	print_submatrix (a, local_rows, n);
      
      
		  if (p > 1) 
		  {
		     datum_size = sizeof(double);
		     
		     max_block_size = BLOCK_SIZE(p-1,p,n);
		     
		     bstorage = my_malloc (id, max_block_size * n * datum_size);
		     
		     b =(double **) my_malloc (id, max_block_size * datum_size);
		     
		     b[0] = bstorage;
		     
		     for (i = 1; i < max_block_size; i++) 
		     {
		        b[i] = b[i-1] + n * datum_size;
		     }
		     
		     for (i = 1; i < p; i++) 
		     {
		        MPI_Send (&prompt, 1, MPI_INT, i, PROMPT_MSG, MPI_COMM_WORLD);
		        
		        MPI_Recv (bstorage, BLOCK_SIZE(i,p,n)*n, MPI_DOUBLE, 
		        	i, RESPONSE_MSG, MPI_COMM_WORLD, &status);
		        	
		        print_submatrix (b, BLOCK_SIZE(i,p,n), n);
		     }
		     
		     free (b);
		     free (bstorage);
		  }
      
      putchar ('\n');
      
   } 
   else 
   {
   	
      MPI_Recv (&prompt, 1, MPI_INT, 0, PROMPT_MSG, MPI_COMM_WORLD, &status);
      
      MPI_Send (*a, local_rows * n, MPI_DOUBLE, 0, RESPONSE_MSG, MPI_COMM_WORLD);
      
   }
}

//PROBABLY COMPLETE
void print_submatrix (
   double       **a,     /* OUT - Doubly-subscripted array */
   int          rows,    /* OUT - Matrix rows */
   int          cols)    /* OUT - Matrix cols */
{
   int i, j;

   for (i = 0; i < rows; i++) 
   {
      for (j = 0; j < cols; j++) 
      {
         	printf ("[%6.2f] ", *( *(a) + i * cols + j));    
      }
      putchar ('\n');
   }
}

//PROBABLY COMPLETE
double *my_malloc(int id, int bytes)
{
	double *buffer;
	
	if ((buffer = malloc ((size_t) bytes)) == NULL) 
	{
		printf ("Error: Malloc failed for process %d\n", id);
		fflush (stdout);
		
		MPI_Abort (MPI_COMM_WORLD, MALLOC_ERROR);
	}
	
	return buffer;
}

//ADD CHECK THAT OUTFILE C CAN BE USED
int validate(int argc, char** argv)
{
	FILE *finptrA;
	FILE *finptrB;
	
	int m, n, l, k;
	
	if (argc != 4)
	{
		fprintf(stderr, "Usage: %s [Matrix A] [Matrix B] [Matrix C]\n", argv[0]);
		MPI_Finalize();
		exit(1);
	}
	
	finptrA = fopen (argv[1], "r");
	finptrB = fopen (argv[2], "r");
	
	if (!finptrA) 
	{
		perror("ERROR: can't open matrix A file\n");
		MPI_Finalize();
		exit(2);
	}
	
	if (!finptrB)
	{
		perror("ERROR: can't open matrix B file\n");
		MPI_Finalize();
		exit(3);
	}
	
	if(	fread(&m, sizeof(int), 1, finptrA) != 1 ||
	 	fread(&n, sizeof(int), 1, finptrA) != 1) 
	{
		fprintf(stderr, "Error reading matrix file A\n");
		MPI_Finalize();
		exit(4);
	}
	
	if(	fread(&l, sizeof(int), 1, finptrB) != 1 ||
	 	fread(&k, sizeof(int), 1, finptrB) != 1) 
	{
		fprintf(stderr, "Error reading matrix file B\n");
		MPI_Finalize();
		exit(5);
	}
	
	if ( m != n )
	{
		fprintf(stderr, "ERROR: Matrix A must be a square matrix\n");
		MPI_Finalize();
		exit(6);
	}
	
	if ( l != k )
	{
		fprintf(stderr, "ERROR: Matrix B must be a square matrix\n");
		MPI_Finalize();
		exit(7);
	}
	
	if ( m != l )
	{
		fprintf(stderr, "ERROR: Matrix A and B must have the same dimensions\n");
		MPI_Finalize();
		exit(8);
	}
	
	fclose (finptrA);
	fclose (finptrB);
	
	return NO_ERR;	
}

//GIVEN FUNCTIONS
int prtmat(char* infile)
{
	int i, j;
	int m, n;
	FILE *finptr;
	double *a;
	int blocksize;

	finptr = fopen (infile, "r");

	if(	fread(&m, sizeof(int), 1, finptr) != 1 ||
	 	fread(&n, sizeof(int), 1, finptr) != 1) 
	{
		perror("Error reading matrix file");
		return 3;
	}
	
	//printf("Matrix[%d][%d]:\n", m, n);

	a = (double*)malloc(n*sizeof(double));
	
	for (i = 0; i < m; i++) 
	{
		if(fread(a, sizeof(datatype), n, finptr) != n) 
		{
	  		perror("Error reading matrix file");
	  		return 3;
		}

		for(j=0; j<n; j++)
			printf(" %6.2lf  ", (double)a[j]);
		
		printf("\n");
	}

	free(a);
	fclose (finptr);
	return 0;
}


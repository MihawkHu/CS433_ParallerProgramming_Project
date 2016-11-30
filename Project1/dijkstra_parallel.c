/*
Purpose:The implementation of the paraller dijkstra algorithm using mpi for
        solving theshortest path problem: find the length of the shortest path
        between a specified vertex and all other vertices in a directed graph.

Input:  n, the number of vertices in the digraph
        mat, the adjacency matrix of the digraph

Output: The submatrix assigned to each process and the complete matrix printed
        from process 0.  Both print "i" instead of 1000000 for infinity.

Compile:mpicc -g -Wall -o dijkstra_parallel dijkstra_parallel.c
Run:    mpiexec -n <p> ./dijkstra_parallel (on lab machines)
        csmpiexec -n <p> ./dijkstra_parallel (on the penguin cluster)

Date:   2016-10
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <time.h>

// #define DEBUG

#define INFINITY 1000000
#define MAX_STRING 10000

// input functions, they are copied from ./mpi_io.c
int Read_n(int my_rank, MPI_Comm comm, FILE *fp);
MPI_Datatype Build_blk_col_type(int n, int loc_n);
void Read_matrix(int loc_mat[], int n, int loc_n,
        MPI_Datatype blk_col_mpi_t, int my_rank, MPI_Comm comm, FILE *fp);
void Print_local_matrix(int loc_mat[], int n, int loc_n, int my_rank);
void Print_matrix(int loc_mat[], int n, int loc_n,
        MPI_Datatype blk_col_mpi_t, int my_rank, MPI_Comm comm);

void Find_min_dist(int* my_min,int loc_n, int *loc_dist, int *loc_known,
        int my_rank, MPI_Comm comm);

// output functions
void Print_dists(int dist[], int n, int option, FILE *fp);
void Print_paths(int pred[], int n, int option, FILE *fp);

int main(int argc, char const *argv[])
{
    int *loc_mat, *loc_dist, *loc_pred, *loc_known, *glo_known, *glo_dist, *glo_pred;
    int n, loc_n, p, my_rank, u, loc_new_dist, i, v;
    int my_min[2], glo_min[2];
    
    // input file
    FILE *fpin;
    char *ch = argv[1];
    if ((fpin = fopen(ch, "r")) == 0) {
        printf("Open input file failed. %d\n", my_rank);        
    }
    else {
        // printf("Read file done\n");
    }
    
    MPI_Comm comm;
    MPI_Datatype blk_col_mpi_t;

    MPI_Init(&argc, &argv);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &my_rank);

    n = Read_n(my_rank, comm, fpin);
    loc_n = n / p;

    //global variable
    glo_known = malloc(n * sizeof(int));
    glo_dist = malloc(n * sizeof(int));
    glo_pred = malloc(n * sizeof(int));

    loc_mat = malloc(n * loc_n * sizeof(int)); // part of the matrix
    loc_dist = malloc(loc_n * sizeof(int)); //store part of the distance
    loc_pred = malloc(loc_n * sizeof(int)); //store part of the pred info

    blk_col_mpi_t = Build_blk_col_type(n, loc_n);
    Read_matrix(loc_mat, n, loc_n, blk_col_mpi_t, my_rank, comm, fpin);
    
    clock_t t1 = clock();
    
    //local initialization
    loc_known = malloc(loc_n * sizeof(int)); //store part of the known info
    for(v = 0; v < loc_n; v++) {
        loc_pred[v] = 0;
        loc_known[v] = 0;
        loc_dist[v] = loc_mat[v];
    }
    //local initialization end

    //initialize the known[0] to be true(1)
    if(my_rank == 0) {
        loc_known[0] = 1;
    }

    //update the globale known/pred/distance(gather and broadcast)
    MPI_Allgather(loc_known, loc_n, MPI_INT, glo_known, loc_n, MPI_INT, comm);
    MPI_Allgather(loc_pred, loc_n, MPI_INT, glo_pred, loc_n, MPI_INT, comm);
    MPI_Allgather(loc_dist, loc_n, MPI_INT, glo_dist, loc_n, MPI_INT, comm);
    //allgather end

    for(i = 1; i < n; i++){
        Find_min_dist(my_min, loc_n, loc_dist, loc_known, my_rank, comm);

        //calculate the cheapest edge and the vertex
        MPI_Allreduce(my_min, glo_min, 1, MPI_2INT, MPI_MINLOC, comm);
        //get the vertex to be added
        u = glo_min[1];
        //update glo_known
        glo_known[u] = 1;

        //for process where u locate, update local known info
        if(my_rank == u / loc_n) {
            loc_known[u%loc_n] = 1;
        }//no need to allgather, because we have update glo_known

        //update local distance and pred info
        for(v = 0; v < loc_n; v++) {
            if(loc_known[v] == 0) {
                loc_new_dist = glo_dist[u] + loc_mat[u*loc_n + v];
                if(loc_new_dist < loc_dist[v]) {
                    loc_dist[v] = loc_new_dist;
                    loc_pred[v] = u;
                }
            }
        }

        //update global distance and pred info
        MPI_Allgather(loc_pred, loc_n, MPI_INT, glo_pred, loc_n, MPI_INT, comm);
        MPI_Allgather(loc_dist, loc_n, MPI_INT, glo_dist, loc_n, MPI_INT, comm);
    }
    
    //output info
    if (my_rank == 0) {
        clock_t t2 = clock();
        printf("Paraller time used: %f s\n", (double)(t2-t1)/CLOCKS_PER_SEC);
        
        if (argc == 2) {
            Print_dists(glo_dist, n, 0, NULL);
            Print_paths(glo_pred, n, 0, NULL);
        }
        
        else if (argc == 3) {
            char *ch = argv[2];
            FILE *fpout;
            if (fpout = fopen(ch, "w")) {
                Print_dists(glo_dist, n, 1, fpout);
                Print_paths(glo_pred, n, 1, fpout);
                fclose(fpout);
            }
            else {
                printf("Open output file failed\n");
            } 
        }
    }
    //free space
    free(loc_mat);
    free(loc_dist);
    free(loc_pred);
    free(loc_known);
    free(glo_known);
    free(glo_dist);
    free(glo_pred);
    MPI_Type_free(&blk_col_mpi_t);
    fclose(fpin);

    MPI_Finalize();
    return 0;
}

// Find the minimam number of loc_dist[] whose distance is unknown
// each process will be distributed loc_n = n / p vertices
void Find_min_dist(int *my_min, int loc_n, int *loc_dist, int *loc_known,
        int my_rank, MPI_Comm comm) {
    int v;
    int loc_u = -1, loc_min_dist = INFINITY;

    for (v = 0; v < loc_n; ++v) {
        if (loc_known[v] == 0) {
            if (loc_dist[v] < loc_min_dist) {
                loc_u = v;
                loc_min_dist = loc_dist[v];
            }
        }
    }

    my_min[0] = loc_min_dist;
    my_min[1] = my_rank*loc_n + loc_u;
}

// print dist to console or file
void Print_dists(int dist[], int n, int option, FILE *fp) {
    int v;
    
    // print to console
    if (option == 0) {
        printf("  v    dist 0->v\n");
        printf("----   ---------\n");
        
        for (v = 1; v < n; v++)
        printf("%3d       %4d\n", v, dist[v]);
        printf("\n");
    }
    
    // print to file
    if (option == 1) {
        fprintf(fp, "  v    dist 0->v\n");
        fprintf(fp, "----   ---------\n");

        for (v = 1; v < n; v++)
            fprintf(fp, "%3d       %4d\n", v, dist[v]);
        fprintf(fp, "\n");
    }
}

// print path to console or file
void Print_paths(int pred[], int n, int option, FILE *fp) {
    int v, w, *path, count, i;

    path =  malloc(n*sizeof(int));
    
    // print to console
    if (option == 0) {
        printf("  v     Path 0->v\n");
        printf("----    ---------\n");
        for (v = 1; v < n; v++) {
            printf("%3d:    ", v);
            count = 0;
            w = v;
            while (w != 0) {
                path[count] = w;
                count++;
                w = pred[w];
            }
            printf("0 ");
            for (i = count-1; i >= 0; i--)
            printf("%d ", path[i]);
            printf("\n");
        }
    }
   
   // print to file
   if (option == 1) {
       fprintf(fp, "  v     Path 0->v\n");
       fprintf(fp, "----    ---------\n");
       
       for (v = 1; v < n; v++) {
           fprintf(fp, "%3d:    ", v);
           count = 0;
           w = v;
           while (w != 0) {
               path[count] = w;
               count++;
               w = pred[w];
           }
           fprintf(fp, "0 ");
           for (i = count-1; i >= 0; i--)
               fprintf(fp, "%d ", path[i]);
           fprintf(fp, "\n");
      }
   }
   
   free(path);
}

/*---------------------------------------------------------------------
 * Function:  Read_n
 * Purpose:   Read in the number of rows in the matrix on process 0
 *            and broadcast this value to the other processes
 * In args:   my_rank:  the calling process' rank
 *            comm:  Communicator containing all calling processes
 * Ret val:   n:  the number of rows in the matrix
 */
int Read_n(int my_rank, MPI_Comm comm, FILE *fpin) {
    int n;

    if (my_rank == 0)
        fscanf(fpin, "%d", &n);
    MPI_Bcast(&n, 1, MPI_INT, 0, comm);
    return n;
}  /* Read_n */


/*---------------------------------------------------------------------
 * Function:  Build_blk_col_type
 * Purpose:   Build an MPI_Datatype that represents a block column of
 *            a matrix
 * In args:   n:  number of rows in the matrix and the block column
 *            loc_n = n/p:  number cols in the block column
 * Ret val:   blk_col_mpi_t:  MPI_Datatype that represents a block
 *            column
 */
MPI_Datatype Build_blk_col_type(int n, int loc_n) {
    MPI_Aint lb, extent;
    MPI_Datatype block_mpi_t;
    MPI_Datatype first_bc_mpi_t;
    MPI_Datatype blk_col_mpi_t;

    MPI_Type_contiguous(loc_n, MPI_INT, &block_mpi_t);
    MPI_Type_get_extent(block_mpi_t, &lb, &extent);

    MPI_Type_vector(n, loc_n, n, MPI_INT, &first_bc_mpi_t);
    MPI_Type_create_resized(first_bc_mpi_t, lb, extent,
        &blk_col_mpi_t);
    MPI_Type_commit(&blk_col_mpi_t);

    MPI_Type_free(&block_mpi_t);
    MPI_Type_free(&first_bc_mpi_t);

    return blk_col_mpi_t;
}  /* Build_blk_col_type */

/*---------------------------------------------------------------------
 * Function:  Read_matrix
 * Purpose:   Read in an nxn matrix of ints on process 0, and
 *            distribute it among the processes so that each
 *            process gets a block column with n rows and n/p
 *            columns
 * In args:   n:  the number of rows in the matrix and the submatrices
 *            loc_n = n/p:  the number of columns in the submatrices
 *            blk_col_mpi_t:  the MPI_Datatype used on process 0
 *            my_rank:  the caller's rank in comm
 *            comm:  Communicator consisting of all the processes
 * Out arg:   loc_mat:  the calling process' submatrix (needs to be
 *               allocated by the caller)
 */
void Read_matrix(int loc_mat[], int n, int loc_n,
        MPI_Datatype blk_col_mpi_t, int my_rank, MPI_Comm comm, FILE *fp) {
    int* mat = NULL, i, j;

    if (my_rank == 0) {
        mat = malloc(n*n*sizeof(int));
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                fscanf(fp, "%d", &mat[i*n + j]);
    }

    MPI_Scatter(mat, 1, blk_col_mpi_t,
        loc_mat, n*loc_n, MPI_INT, 0, comm);

    if (my_rank == 0) free(mat);
}  /* Read_matrix */


/*---------------------------------------------------------------------
 * Function:  Print_local_matrix
 * Purpose:   Store a process' submatrix as a string and print the
 *            string.  Printing as a string reduces the chance
 *            that another process' output will interrupt the output.
 *            from the calling process.
 * In args:   loc_mat:  the calling process' submatrix
 *            n:  the number of rows in the submatrix
 *            loc_n:  the number of cols in the submatrix
 *            my_rank:  the calling process' rank
 */
void Print_local_matrix(int loc_mat[], int n, int loc_n, int my_rank) {
    char temp[MAX_STRING];
    char *cp = temp;
    int i, j;

    sprintf(cp, "Proc %d >\n", my_rank);
    cp = temp + strlen(temp);
    for (i = 0; i < n; i++) {
        for (j = 0; j < loc_n; j++) {
            if (loc_mat[i*loc_n + j] == INFINITY)
                sprintf(cp, " i ");
            else
                sprintf(cp, "%2d ", loc_mat[i*loc_n + j]);
            cp = temp + strlen(temp);
        }
        sprintf(cp, "\n");
        cp = temp + strlen(temp);
    }

    printf("%s\n", temp);
}  /* Print_local_matrix */


/*---------------------------------------------------------------------
 * Function:  Print_matrix
 * Purpose:   Print the matrix that's been distributed among the
 *            processes.
 * In args:   loc_mat:  the calling process' submatrix
 *            n:  number of rows in the matrix and the submatrices
 *            loc_n:  the number of cols in the submatrix
 *            blk_col_mpi_t:  MPI_Datatype used on process 0 to
 *               receive a process' submatrix
 *            my_rank:  the calling process' rank
 *            comm:  Communicator consisting of all the processes
 */
void Print_matrix(int loc_mat[], int n, int loc_n,
        MPI_Datatype blk_col_mpi_t, int my_rank, MPI_Comm comm) {
    int* mat = NULL, i, j;

    if (my_rank == 0) mat = malloc(n*n*sizeof(int));
    MPI_Gather(loc_mat, n*loc_n, MPI_INT,
            mat, 1, blk_col_mpi_t, 0, comm);
    if (my_rank == 0) {
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++)
                if (mat[i*n + j] == INFINITY)
                    printf(" i ");
                else
                    printf("%2d ", mat[i*n + j]);
            printf("\n");
        }
        free(mat);
    }
}  /* Print_matrix */
/* File:     dijkstra.c
 * Purpose:  Implement Dijkstra's algorithm for solving the single-source 
 *           shortest path problem:  find the length of the shortest path 
 *           between a specified vertex and all other vertices in a 
 *           directed graph.
 *
 * Input:    n, the number of vertices in the digraph
 *           mat, the adjacency matrix of the digraph
 * Output:   A list showing the cost of the shortest path
 *           from vertex 0 to every other vertex in the graph.
 *
 * Compile:  gcc -g -Wall -o dijkstra dijkstra.c
 * Run:      ./dijkstra
 *           For large matrices, put the matrix into a file with n as
 *           the first line and run with ./dijkstra < large_matrix
 *
 * Notes:
 * 1.  Edge lengths should be nonnegative.
 * 2.  The distance from v to w may not be the same as the distance from
 *     w to v.
 * 3.  If there is no edge between two vertices, the length is the constant
 *     INFINITY.  So input edge length should be substantially less than
 *     this constant.
 * 4.  The cost of travelling from a vertex to itself is 0.  So the adjacency
 *     matrix has zeroes on the main diagonal.
 * 5.  No error checking is done on the input.
 * 6.  The adjacency matrix is stored as a 1-dimensional array and subscripts
 *     are computed using the formula:  the entry in the ith row and jth
 *     column is mat[i*n + j]
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

const int INFINITY = 1000000;

void Read_matrix(int mat[], int n, FILE *fp);
void Print_matrix(int mat[], int n);
void Print_dists(int dist[], int n, FILE *fp);
void Print_paths(int pred[], int n, FILE *fp);
int  Find_min_dist(int dist[], int known[], int n);
void Dijkstra(int mat[], int dist[], int pred[], int n);

int main(int argc, char const *argv[]) {

   int  n;
   int *mat, *dist, *pred;
   
   // input file 
   FILE *fpin;
   char *ch = argv[1];
   if (fpin = fopen(ch, "r")){
    //    printf("Read file done\n");
   }
   else {
       printf("Open input file failed\n");
   }
   

   // printf("How many vertices?\n");
   fscanf(fpin, "%d", &n);
   mat = malloc(n*n*sizeof(int));
   dist = malloc(n*sizeof(int));
   pred = malloc(n*sizeof(int));

   // printf("Enter the matrix\n");
   Read_matrix(mat, n, fpin);
   
   clock_t t1 = clock();
   
   Dijkstra(mat, dist, pred, n);
   
   clock_t t2 = clock();
   printf("Serial time used: %f s\n", (double)(t2-t1)/CLOCKS_PER_SEC);

   // printf("The distance from 0 to each vertex is:\n");
   // printf("The shortest path from 0 to each vertex is:\n");
   
   FILE *fpout;
   ch = argv[2];
   if (fpout = fopen(ch, "w")) {
       Print_dists(dist, n, fpout);
       Print_paths(pred, n, fpout);
   }
   else {
       printf("Output file failed.\n");
   }
   
   free(mat);
   free(dist);
   free(pred);
   fclose(fpin);
   fclose(fpout);
   
   return 0;
}  /* main */

/*-------------------------------------------------------------------
 * Function:  Read_matrix
 * Purpose:   Read in the adjacency matrix
 * In arg:    n
 * Out arg:   mat
 */
void Read_matrix(int mat[], int n, FILE *fp) {
   int i, j;

   for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
         fscanf(fp, "%d", &mat[i*n+j]);
}  /* Read_matrix */

/*-------------------------------------------------------------------
 * Function:  Print_matrix
 * Purpose:   Print the contents of the matrix
 * In args:   mat, n
 */
void Print_matrix(int mat[], int n) {
   int i, j;

   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++)
         if (mat[i*n+j] == INFINITY)
            printf("i ");
         else
            printf("%d ", mat[i*n+j]);
      printf("\n");
   }
}  /* Print_matrix */

/*-------------------------------------------------------------------
 * Function:    Dijkstra
 * Purpose:     Apply Dijkstra's algorithm to the matrix mat
 * In args:     n:  the number of vertices
 *              mat:  adjacency matrix for the graph
 * Out args:    dist:  dist[v] = distance 0 to v.
 *              pred:  pred[v] = predecessor of v on a 
 *                  shortest path 0->v.
 */
void Dijkstra(int mat[], int dist[], int pred[], int n) {
   int i, u, v, *known, new_dist;

   /* known[v] = true, if the shortest path 0->v is known */
   /* known[v] = false, otherwise                         */
   known = malloc(n*sizeof(int));

   /* Initialize d and p */
   dist[0] = 0; pred[0] = 0; known[0] = 1; 
   for (v = 1; v < n; v++) {
      dist[v] = mat[0*n + v];
      pred[v] = 0;
      known[v] = 0;
   }

#     ifdef DEBUG
      printf("i = 0\n");
      Print_dists(dist, n);
#     endif

   /* On each pass find an additional vertex */
   /* whose distance to 0 is known           */
   for (i = 1; i < n; i++) {
      u = Find_min_dist(dist, known, n);

      known[u] = 1;

      for (v = 1; v < n; v++) 
         if (!known[v]) {
            new_dist = dist[u] + mat[u*n + v];
            if (new_dist < dist[v]) {
               dist[v] = new_dist;
               pred[v] = u;
            }
         }

#     ifdef DEBUG
      printf("i = %d\n", i);
      Print_dists(dist, n);
#     endif
   } /* for i */

   free(known);
}  /* Dijkstra */

/*-------------------------------------------------------------------
 * Function:    Find_min_dist
 * Purpose:     Find the vertex u with minimum distance to 0
 *              (dist[u]) among the vertices whose distance 
 *              to 0 is not known.
 * In args:     dist:  dist[v] = current estimate of distance
 *                 0->v
 *              known:  whether the minimum distance 0-> is
 *                 known
 *              n:  the total number of vertices
 * Ret val:     The vertex u whose distance to 0, dist[u]
 *              is a minimum among vertices whose distance
 *              to 0 is not known.
 */
int Find_min_dist(int dist[], int known[], int n) {
   int v, u, best_so_far = INFINITY;

   for (v = 1; v < n; v++)
      if (!known[v])
         if (dist[v] < best_so_far) {
            u = v;
            best_so_far = dist[v];
         }

   return u;
}  /* Find_min_dist */


/*-------------------------------------------------------------------
 * Function:    Print_dists
 * Purpose:     Print the length of the shortest path from 0 to each
 *              vertex
 * In args:     n:  the number of vertices
 *              dist:  distances from 0 to each vertex v:  dist[v]
 *                 is the length of the shortest path 0->v
 */
void Print_dists(int dist[], int n, FILE *fp) {
   int v;

   fprintf(fp, "  v    dist 0->v\n");
   fprintf(fp, "----   ---------\n");
                  
   for (v = 1; v < n; v++)
      fprintf(fp, "%3d       %4d\n", v, dist[v]);
   fprintf(fp, "\n");
} /* Print_dists */  


/*-------------------------------------------------------------------
 * Function:    Print_paths
 * Purpose:     Print the shortest path from 0 to each vertex
 * In args:     n:  the number of vertices
 *              pred:  list of predecessors:  pred[v] = u if
 *                 u precedes v on the shortest path 0->v
 */
void Print_paths(int pred[], int n, FILE *fp) {
   int v, w, *path, count, i;

   path =  malloc(n*sizeof(int));

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

   free(path);
}  /* Print_paths */
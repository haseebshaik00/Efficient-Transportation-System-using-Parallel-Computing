#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#define MAX_STRING 10000

#define INFINITY 1000000
 
 int Read_n(int my_rank, MPI_Comm comm);
 MPI_Datatype Build_blk_col_type(int n, int loc_n);
 void Read_matrix(int loc_mat[], int n, int loc_n,MPI_Datatype blk_col_mpi_t, int my_rank, MPI_Comm comm);
 void Print_local_matrix(int loc_mat[], int n, int loc_n, int my_rank);
 void Print_matrix(int loc_mat[], int n, int loc_n,MPI_Datatype blk_col_mpi_t, int my_rank, MPI_Comm comm);
 void Dijkstra(int mat[], int dist[], int pred[], int n, int loc_n, int my_rank,MPI_Comm comm);
 void Initialize_matrix(int mat[], int loc_dist[], int loc_pred[], int known[], int loc_n, int my_rank);

 int Find_min_dist(int loc_dist[], int known[], int loc_n, int my_rank, MPI_Comm comm);
 int Global_vertex(int loc_u, int loc_n, int my_rank);
 void Print_dists(int loc_dist[], int n, int loc_n, int my_rank, MPI_Comm comm);
 void Print_paths(int loc_pred[], int n, int loc_n, int my_rank, MPI_Comm comm);


   int main(int argc, char* argv[]) 
    {
     int *loc_mat, *loc_dist, *loc_pred;
     int n, loc_n, p, my_rank;
     MPI_Comm comm;
     MPI_Datatype blk_col_mpi_t;
     MPI_Init(&argc, &argv);
     comm = MPI_COMM_WORLD;
     MPI_Comm_size(comm, &p);
     MPI_Comm_rank(comm, &my_rank);
     n = Read_n(my_rank, comm);  
     loc_n = n/p;
     loc_mat  = malloc(n*loc_n*sizeof(int));
     loc_dist = malloc(n*loc_n*sizeof(int));
     loc_pred = malloc(n*loc_n*sizeof(int));
     blk_col_mpi_t = Build_blk_col_type(n, loc_n);
     Read_matrix(loc_mat, n, loc_n, blk_col_mpi_t, my_rank, comm); 
     Dijkstra(loc_mat, loc_dist, loc_pred, n, loc_n, my_rank, comm);
     Print_dists(loc_dist, n, loc_n, my_rank, comm);
     Print_paths(loc_pred, n, loc_n, my_rank, comm);
     free(loc_mat);
     free(loc_dist);
     free(loc_pred);
     MPI_Type_free(&blk_col_mpi_t);
     MPI_Finalize();
     return 0;
   }
   int Read_n(int my_rank, MPI_Comm comm) 
   {  int n;
      if (my_rank == 0)
         printf("Enter number of vertices in the matrix: \n"); scanf("%d", &n);
      MPI_Bcast(&n, 1, MPI_INT, 0, comm); return n;
   }
  MPI_Datatype Build_blk_col_type(int n, int loc_n) 
  {
    MPI_Aint lb, extent;
    MPI_Datatype block_mpi_t;
    MPI_Datatype first_bc_mpi_t;
    MPI_Datatype blk_col_mpi_t;
    MPI_Type_contiguous(loc_n, MPI_INT, &block_mpi_t);
    MPI_Type_get_extent(block_mpi_t, &lb, &extent);
    MPI_Type_vector(n, loc_n, n, MPI_INT, &first_bc_mpi_t);
    MPI_Type_create_resized(first_bc_mpi_t, lb, extent,&blk_col_mpi_t);
    MPI_Type_commit(&blk_col_mpi_t); MPI_Type_free(&block_mpi_t);
    MPI_Type_free(&first_bc_mpi_t);
    return blk_col_mpi_t;
  }
 void Read_matrix(int loc_mat[], int n, int loc_n,MPI_Datatype blk_col_mpi_t, int my_rank, MPI_Comm comm)
 { 
   int* mat = NULL, i, j;
   if (my_rank == 0) 
   {
     mat = malloc(n*n*sizeof(int)); 
     for (i = 0; i < n; i++)
     for (j = 0; j < n; j++)
       scanf("%d", &mat[i*n + j]);
   }
   MPI_Scatter(mat, 1, blk_col_mpi_t,loc_mat, n*loc_n, MPI_INT, 0, comm);

   if (my_rank == 0) 
      free(mat);
  }
  void Print_local_matrix(int loc_mat[], int n, int loc_n, int my_rank) 
 { 
   char temp[MAX_STRING];
   char *cp = temp; int i, j;
   sprintf(cp, "Proc %d >\n", my_rank); cp = temp + strlen(temp);
   for (i = 0; i < n; i++) 
   {
     for (j = 0; j < loc_n; j++) 
     {
       if (loc_mat[i*loc_n + j] == INFINITY) 
           sprintf(cp, " i ");
       else
           sprintf(cp, "%2d ", loc_mat[i*loc_n + j]); cp = temp + strlen(temp);
     }
     sprintf(cp, "\n");
     cp = temp + strlen(temp);
   }
   printf("%s\n", temp);
 }
  void Print_matrix(int loc_mat[], int n, int loc_n, MPI_Datatype blk_col_mpi_t, int my_rank, MPI_Comm comm) 
  { 
    int* mat = NULL, i, j;
    if (my_rank == 0) 
       mat = malloc(n*n*sizeof(int)); 
    MPI_Gather(loc_mat, n*loc_n, MPI_INT, mat, 1, blk_col_mpi_t, 0, comm); 
    if (my_rank == 0) 
    {
      for (i = 0; i < n; i++) 
     { 
         for (j = 0; j < n; j++)
            if (mat[i*n + j] == INFINITY) 
               printf(" i ");
            else
               printf("%2d ", mat[i*n + j]); printf("\n");
     }
     free(mat);
    }
  }
  void Dijkstra(int mat[], int loc_dist[], int loc_pred[], int n, int loc_n, int my_rank, MPI_Comm comm) 
 {
    int i, u, *known, new_dist; int loc_u, loc_v;
    known = malloc(loc_n*sizeof(int));
    Initialize_matrix(mat, loc_dist, loc_pred, known, loc_n, my_rank); 
    for (i = 1; i < n; i++) 
    {
      loc_u = Find_min_dist(loc_dist, known, loc_n, my_rank, comm); int my_min[2], glbl_min[2];
      int g_min_dist;
      if (loc_u < INFINITY) 
      { 
        my_min[0] = loc_dist[loc_u];
        my_min[1] = Global_vertex(loc_u, loc_n, my_rank);
      } 
      else 
      {
        my_min[0] = INFINITY; my_min[1] = INFINITY;
      }
     MPI_Allreduce(my_min, glbl_min, 1, MPI_2INT, MPI_MINLOC, comm); u = glbl_min[1];
     g_min_dist = glbl_min[0]; if (u/loc_n == my_rank) {
     loc_u = u % loc_n; known[loc_u] = 1;
    }
    for (loc_v = 0; loc_v < loc_n; loc_v++) 
       if (!known[loc_v]) 
       {
         new_dist = g_min_dist + mat[u*loc_n + loc_v]; 
         if (new_dist < loc_dist[loc_v]) 
         {
           loc_dist[loc_v] = new_dist; loc_pred[loc_v] = u;
         }
       }
   }
   free(known);
  }
  void Initialize_matrix(int mat[], int loc_dist[], int loc_pred[], int known[], int loc_n, int my_rank) 
 {
      for (int v = 0; v < loc_n; v++) 
      { 
         loc_dist[v] = mat[0*loc_n + v]; 
         loc_pred[v] = 0;
         known[v] = 0;
      }
      if (my_rank == 0) { known[0] = 1;
  }
 }
  int Find_min_dist(int loc_dist[], int loc_known[], int loc_n, int my_rank, MPI_Comm comm) 
 {
    int loc_v, loc_u;
    int loc_min_dist = INFINITY; 
    loc_u = INFINITY;
    for (loc_v = 0; loc_v < loc_n; loc_v++) 
       if (!loc_known[loc_v])
          if (loc_dist[loc_v] < loc_min_dist) 
          { 
             loc_u = loc_v;
             loc_min_dist = loc_dist[loc_v];
          }
       return loc_u;
  }
  int Global_vertex(int loc_u, int loc_n, int my_rank) 
 { 
    int global_u = loc_u + my_rank*loc_n;
    return global_u;
 }
  void Print_dists(int loc_dist[], int n, int loc_n, int my_rank, MPI_Comm comm) 
 { 
    int v;
    int* dist = NULL; 
    if (my_rank == 0) 
    {
      dist = malloc(n*sizeof(int));
    }
    MPI_Gather(loc_dist, loc_n, MPI_INT, dist, loc_n, MPI_INT, 0, comm); 
    if (my_rank == 0) 
    {
       printf("The distance from 0 to each vertex is:\n"); 
       printf(" v dist 0->v\n");
       printf("	      \n");
       for (v = 1; v < n; v++)
           printf("%3d	  %4d\n", v, dist[v]); printf("\n");
       free(dist);
    }
 }
  void Print_paths(int loc_pred[], int n, int loc_n, int my_rank, MPI_Comm comm) 
 { 
    int v, w, *path, count, i;
    int* pred = NULL; 
    if (my_rank == 0) 
    {
      pred = malloc(n*sizeof(int));
    }
    MPI_Gather(loc_pred, loc_n, MPI_INT, pred, loc_n, MPI_INT, 0, comm); 
    if (my_rank == 0) 
    {
       path = malloc(n*sizeof(int));
       printf("The shortest path from 0 to each vertex is:\n"); 
       printf("  v	  Path 0->v\n");
       printf("	      \n");
       for (v = 1; v < n; v++) 
       { 
          printf("%3d:	    ", v); 
          count = 0;
          w = v;
          while (w != 0) 
          {  
             path[count] = w; 
             count++;
             w = pred[w];
          }
          printf("0 ");
          for (i = count-1; i >= 0; i--)
             printf("%d ", path[i]); printf("\n");
      }
    free(path); free(pred);
   }
 }

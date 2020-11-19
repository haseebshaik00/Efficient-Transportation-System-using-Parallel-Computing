#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<limits.h>
#define V 4

int minDistance(int dist[], int sptSet[])
{
    int min = INT_MAX, min_index;

    for (int v = 0; v < V; v++)
        if (sptSet[v] == 0 && dist[v] <= min)
            min = dist[v], min_index = v;
    return min_index;
}

void printSolution(int dist[])
{
    printf("Vertex \t\t Distance from Source\n");
    for (int i = 0; i < V; i++)
        printf("%d \t\t %d\n", i, dist[i]);
}

void dijkstra(int graph[V][V], int src)
{
    int dist[V];
    int sptSet[V];
    for (int i = 0; i < V; i++)
        dist[i] = INT_MAX, sptSet[i] = 0;
    dist[src] = 0;
    for (int count = 0; count < V - 1; count++) {
        int u = minDistance(dist, sptSet);
        sptSet[u] =    1;
        for (int v = 0; v < V; v++)
            if (!sptSet[v] && graph[u][v] && dist[u] != INT_MAX
                && dist[u] + graph[u][v] < dist[v])
                dist[v] = dist[u] + graph[u][v];
    }
    printSolution(dist);
}

int main()
{
    double time_spent = 0.0;
	clock_t begin = clock();
    int graph[V][V] = {{0,5,0,10},{0,0,3,0},{0,0,0,1},{0,0,0,0}};
    printf("Adjency Matrix: \n0\t5\t-\t10 \n-\t0\t3\t-\t\n-\t-\t0\t1\n-\t-\t-\t0\n\n");
    dijkstra(graph,0);
    clock_t end = clock();
    time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
	printf("\nTime elapsed is %f seconds\n", time_spent);
    return 0;
}

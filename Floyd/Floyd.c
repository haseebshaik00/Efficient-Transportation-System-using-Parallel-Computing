#include<stdio.h>
#include <time.h>
#define V 4
#define INF 99999

void printSolution(int dist[][V])
{
	for (int i = 0; i < V; i++)
    {
        for (int j = 0; j < V; j++)
        {
            if (dist[i][j] == INF)
                printf("%s\t", "-");
            else
                printf ("%d\t", dist[i][j]);
        }
        printf("\n");
    }
}

void floydWarshall (int graph[][V])
{
    int dist[V][V], i, j, k;
    for (i = 0; i < V; i++)
        for (j = 0; j < V; j++)
            dist[i][j] = graph[i][j];
    for (k = 0; k < V; k++)
    {
        for (i = 0; i < V; i++)
        {
            for (j = 0; j < V; j++)
            {
                if (dist[i][k] + dist[k][j] < dist[i][j])
                    dist[i][j] = dist[i][k] + dist[k][j];
            }
        }
    }
    printf ("\nMatrix after Floyd Warshall Algorithm\n");
    printSolution(dist);
}

int main()
{
    double time_spent = 0.0;
	clock_t begin = clock();
    int graph[V][V] = {{0,5,INF,10},{INF,0,3,INF},{INF,INF,0,1},{INF,INF,INF,0}};
    printf("Initial Matrix\n");
    printSolution(graph);
    floydWarshall(graph);
    clock_t end = clock();
    time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Time elpased is %f seconds", time_spent);
	return 0;
}

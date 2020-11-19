#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <limits.h>

struct Edge
{
    int source, destination, weight;
};

struct Graph
{
    int V, E;
    struct Edge* edge;
};

struct Graph* createGraph(int V, int E)
{
    struct Graph* graph = (struct Graph*) malloc( sizeof(struct Graph));
    graph->V = V;
    graph->E = E;
    graph->edge = (struct Edge*) malloc( graph->E * sizeof( struct Edge ) );
    return graph;
}

void FinalSolution(int dist[], int n)
{
    printf("Vertex distance from Source vertex :\n");
    int i;
    for (i = 0; i < n; ++i){
		printf("%d \t\t %d\n", i, dist[i]);
	}
}
void BellmanFord(struct Graph* graph, int source)
{
    int V = graph->V;
    int E = graph->E;
    int StoreDistance[V];
    int i,j;
    for (i = 0; i < V; i++)
        StoreDistance[i] = INT_MAX;
    StoreDistance[source] = 0;
    for (i = 1; i <= V-1; i++){
        for (j = 0; j < E; j++){
            int u = graph->edge[j].source;
            int v = graph->edge[j].destination;
            int weight = graph->edge[j].weight;
            if (StoreDistance[u] + weight < StoreDistance[v])
                StoreDistance[v] = StoreDistance[u] + weight;
        }
    }
    for (i = 0; i < E; i++){
        int u = graph->edge[i].source;
        int v = graph->edge[i].destination;
        int weight = graph->edge[i].weight;
        if (StoreDistance[u] + weight < StoreDistance[v])
            printf("This graph contains negative edge cycle\n");
    }
    FinalSolution(StoreDistance, V);
    return;
}

int main()
{
    double time_spent = 0.0;
	clock_t begin = clock();
    int V=4,E=4,S=0;
    struct Graph* graph = createGraph(V, E);
    int i,k=0;
    int a[]={0,3,10,0,1,5,1,2,3,2,3,1};
    for(i=0;i<E;i++)
    {
        graph->edge[i].source=a[k++];
        graph->edge[i].destination=a[k++];
        graph->edge[i].weight=a[k++];
    }
    BellmanFord(graph, S);
    clock_t end = clock();
    time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
	printf("\nTime elapsed is %f seconds\n", time_spent);
    return 0;
}

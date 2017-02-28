// Author : Somay Jain
#include <stdio.h>
#include <stdlib.h>

#define IDX(x, y, t, X, Y) (((t) * X * Y) + ((y) * X) + (x))

void simulateRow(int y, int X, int numRows, int currentTimestep, int *buffer);
int isInside(int x, int y, int X, int Y);
int getNeighbours(int x, int y, int X, int Y, int currentTimestep, int *buffer);

int main(int argc, char** argv)
{
    if(argc < 5)
    {
        printf("Too few command line arguments\n");
        printf("Usage: life <input file name> <# of generations> <X_limit> <Y_limit>\n");
        exit(1);
    }

    int X, Y, nIter;
    int currentTimestep = 0;

    nIter = atoi(argv[2]);
    X = atoi(argv[3]);
    Y = atoi(argv[4]);
    
    int *buffer = (int*) calloc(Y * X * 2, sizeof(int));

    FILE *fp = fopen(argv[1], "r");
    int x, y;
    while(fscanf(fp, "%d %d", &x, &y) != EOF)
    {
    	buffer[IDX(x, y, currentTimestep, X, Y)] = 1;
    }

    for(int iter = 0; iter < nIter; iter++)
    {
        for(int y = 0; y < Y; y++)
        {
            simulateRow(y, X, Y, currentTimestep, buffer);
        }
        currentTimestep = 1 - currentTimestep;
    }

    // Print the ans.
    for(int x = 0; x < X; x++)
    {
        for(int y = 0; y < Y; y++)
        {
            if(buffer[IDX(x, y, currentTimestep, X, Y)] == 1)
                printf("%d %d\n", x, y);
        }
    }

	return 0;
}

int getNeighbours(int x, int y, int X, int Y, int currentTimestep, int *buffer)
{
    int num_neighbours = 0;
    if(x > 0)
    {
        // get left
        num_neighbours += buffer[IDX(x-1, y, currentTimestep, X, Y)];

        if(y > 0)
            // get left down
            num_neighbours += buffer[IDX(x-1, (y-1), currentTimestep, X, Y)];

        if(y < Y-1)
            // get left up
            num_neighbours += buffer[IDX(x-1, (y+1), currentTimestep, X, Y)];
    }
    if(x < X-1)
    {
        // get right
        num_neighbours += buffer[IDX(x+1, y, currentTimestep, X, Y)];

        if(y > 0)
            // get right down
            num_neighbours += buffer[IDX(x+1, (y-1), currentTimestep, X, Y)];

        if(y < Y-1)
            // get right up
            num_neighbours += buffer[IDX(x+1, (y+1), currentTimestep, X, Y)];
    }
    if(y > 0)
        // get bottom
        num_neighbours += buffer[IDX(x, (y-1), currentTimestep, X, Y)];
    if(y < Y-1)
        // get top
        num_neighbours += buffer[IDX(x, (y+1), currentTimestep, X, Y)];
    return num_neighbours;
}

int isInside(int x, int y, int X, int Y)
{
    if (x < 0 || x >= X)
        return 0;
    if(y < 0 || y >= Y)
        return 0;
    return 1;
}

void simulateRow(int y, int X, int numRows, int currentTimestep, int *buffer)
{
    for(int x=0; x<X; x++)
    {
        int num_neighbours = getNeighbours(x, y, X, numRows, currentTimestep, buffer);
  
        if(buffer[IDX(x, y, currentTimestep, X, numRows)] == 1)
        {
            // Cell was alive.
            if(num_neighbours == 2 || num_neighbours == 3)
                buffer[IDX(x, y, 1 - currentTimestep, X, numRows)] = 1;
            else
                buffer[IDX(x, y, 1 - currentTimestep, X, numRows)] = 0;
        }
        else
        {
            if(num_neighbours == 3)
                buffer[IDX(x, y, 1 - currentTimestep, X, numRows)] = 1;
            else
                buffer[IDX(x, y, 1 - currentTimestep, X, numRows)] = 0;
        }
    }
}
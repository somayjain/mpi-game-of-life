 // Author : Somay Jain
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <vector>

using namespace std;

// t: timestep (0 or 1)
// X: width of the array
// Y: height of the array
// x, y : the element index desired.
#define IDX(x, y, t, X, Y) (((t) * X * Y) + ((y) * X) + (x))

typedef struct SubGrid
{
    int leftX, leftY;       // Top left corner of the original grid
    int rightX, rightY;     // Bottom right corner of the original grid
    int numX, numY;         // Actual number of rows, cols (including padding)
    int pid;                // Process ID
    int Px, Py;             // Location of this process in the grid
} SubGrid;

void getWorkingDimHelper(int pidDir, int numProcessorsDir, int sizeDir, int *lo, int *high);
void getGridWorkingDims(SubGrid *grid, int coord[2], int dims[2], int X, int Y, int pid);
int simulatePartialRow(int y, int fromX, int tillX, int maxX, int maxY, int currentTimestep, int *buffer);
void simulateHelper(int x, int y, int maxX, int maxY, int currentTimestep, int *buffer);
bool isInsideGrid(int* coord, SubGrid grid);
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
    MPI_Status status;
    int tag = 12345;

    int dirX[8] = {0,  1,  0, -1,  1,  1, -1, -1};
    int dirY[8] = {1,  0, -1,  0,  1, -1, -1,  1};
    int neighbourPids[8];


    int X, Y, nIter, idx;
    int numProcessors, pid;
    int currentTimestep = 0;
    SubGrid grid;

    MPI_Init(&argc,&argv);
    
    nIter = atoi(argv[2]);
    X = atoi(argv[3]);
    Y = atoi(argv[4]);
    
    MPI_Comm_size(MPI_COMM_WORLD, &numProcessors);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    int dims[2] = {0, 0};
    int periods[2] = {1, 1};
    MPI_Dims_create(numProcessors, 2, dims);
    
    if(dims[0] > dims[1])
    {
        int temp = dims[0];
        dims[0] = dims[1];
        dims[1] = temp;
    }

    // if(pid == 0)
        // printf("Dims %d %d\n", dims[0], dims[1]);

    MPI_Comm new_comm;

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &new_comm);

    // Get process cartesian coordinates for the current process.
    int coord[2];
    MPI_Cart_coords(new_comm, pid, 2, coord);

    // Get the details of the SubGrid this process is working on.
    getGridWorkingDims(&grid, coord, dims, X, Y, pid);
    // if(pid == 1)
        // printf("ID %d x, y (%d, %d) working on (%d, %d)->(%d, %d) \n", pid, coord[0], coord[1],
        // grid.leftX, grid.leftY, grid.rightX, grid.rightY);
    // printf("Area %d\n", (grid.rightX - grid.leftX + 1) *  (grid.rightY - grid.leftY + 1));
    // printf("ID %d (%d, %d)\n", pid, grid.numX, grid.numY);

    int *buffer = (int*) calloc(grid.numX * grid.numY * 2, sizeof(int));

    if(pid == 0)
    {
        // Get the working region of all other processes so that we know whom to
        // send the data.
        SubGrid* allGrids = (SubGrid*) malloc(sizeof(SubGrid) * numProcessors);
        int ctr = 0;
        vector< vector< int > > input_data;
        vector<int> empty_row;
        for(int i=0; i<dims[0]; i++)
        {
            for(int j=0; j<dims[1]; j++)
            {
                coord[0] = i;
                coord[1] = j;
                int rank;
                MPI_Cart_rank(new_comm, coord, &rank);
                getGridWorkingDims(&allGrids[ctr++], coord, dims, X, Y, rank);

                input_data.push_back(empty_row);
            }
        }

        FILE *fp = fopen(argv[1], "r");
        while(fscanf(fp, "%d%d", &coord[0], &coord[1]) != EOF)
        {
            // Check if this coord belongs to the current process.
            if(isInsideGrid(coord, grid))
            {
                // printf("(%d, %d) belongs to current process, putting at (%d,%d)\n", coord[0], coord[1], coord[0] - grid.leftX + 1, coord[1] - grid.leftY + 1);
                buffer[IDX(coord[0] - grid.leftX + 1, coord[1] - grid.leftY + 1, 0, grid.numX, grid.numY)] = 1;
            }
            else
            {
                for(int i=1; i<numProcessors; i++)
                {
                    if(isInsideGrid(coord, allGrids[i]))
                    {
                        input_data[i].push_back(coord[0]);
                        input_data[i].push_back(coord[1]);
                        break;
                    }
                }
            }
        }

        for(int i=1; i<numProcessors; i++)
        {
            int size = input_data[i].size();
            MPI_Send(&size, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
            
            if(size > 0)
            {
                // Send the vector.
                MPI_Send(&input_data[i].front(), input_data[i].size(), MPI_INT, i, tag, MPI_COMM_WORLD);
                // printf("Sending %d elements to %d\n", size, i);
            }
        }
        fclose(fp);
    }
    else
    {
        int size;
        int *input_data;
        MPI_Recv(&size, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
        
        
        if(size > 0)
        {
            input_data = (int*) malloc(sizeof(int) * size);
            MPI_Recv(input_data, size, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
            // printf("Pid %d Received %d elements\n", pid, size);
        }

        for(int i=0; i<size; i+=2)
        {
            // printf("(%d, %d) belongs to process %d, putting at (%d,%d)\n", input_data[i], input_data[i+1], pid, input_data[i] - grid.leftX + 1, input_data[i+1] - grid.leftY + 1);
            idx = IDX(input_data[i] - grid.leftX + 1, input_data[i+1] - grid.leftY + 1, 0, grid.numX, grid.numY);
            buffer[idx] = 1;
        }
        if(size > 0)
            free(input_data);
    }

    // Reading complete, should synchronize here?
    MPI_Barrier(MPI_COMM_WORLD);

    // Populate pids of the neighbours
    for(int i=0; i<8; i++)
    {
        coord[0] = grid.Px + dirX[i];
        coord[1] = grid.Py + dirY[i];
        if(coord[0] < 0 || coord[0] >= dims[0])
            neighbourPids[i] = -1;
        else if(coord[1] < 0 || coord[1] >= dims[1])
            neighbourPids[i] = -1;
        else
            MPI_Cart_rank(new_comm, coord, &neighbourPids[i]);
    }

    MPI_Request request[8] = {0}, send_req[8];
    MPI_Status recv_status[8];
    int err_code;
    
    // Generate the data type for the column
    MPI_Datatype colSlice;
    MPI_Type_vector(grid.numY - 2, 1, grid.numX, MPI_INT, &colSlice);
    MPI_Type_commit(&colSlice);

    for(int iter = 0; iter < nIter; iter++)
    {
        for(int i=0; i<8; i++)
            request[i] = NULL;

        if(neighbourPids[0] != -1)
        {
            // Up
            idx = IDX(1, grid.numY - 1, currentTimestep, grid.numX, grid.numY);
            err_code = MPI_Irecv(&buffer[idx], grid.numX - 2, MPI_INT, 
                                 neighbourPids[0], tag, MPI_COMM_WORLD, &request[0]);

            idx = IDX(1, grid.numY - 2, currentTimestep, grid.numX, grid.numY);
            err_code = MPI_Isend(&buffer[idx], grid.numX - 2, MPI_INT, 
                                 neighbourPids[0], tag, MPI_COMM_WORLD, &send_req[0]);
        }

        if(neighbourPids[1] != -1)
        {
            // Right
            idx = IDX(grid.numX - 1, 1, currentTimestep, grid.numX, grid.numY);
            err_code = MPI_Irecv(&buffer[idx], 1, colSlice,
                                 neighbourPids[1], tag, MPI_COMM_WORLD, &request[1]);

            idx = IDX(grid.numX - 2, 1, currentTimestep, grid.numX, grid.numY);
            err_code = MPI_Isend(&buffer[idx], 1, colSlice,
                                 neighbourPids[1], tag, MPI_COMM_WORLD, &send_req[1]);
        }

        if(neighbourPids[2] != -1)
        {
            // Down
            idx = IDX(1, 0, currentTimestep, grid.numX, grid.numY);
            err_code = MPI_Irecv(&buffer[idx], grid.numX - 2, MPI_INT, 
                                 neighbourPids[2], tag, MPI_COMM_WORLD, &request[2]);

            idx = IDX(1, 1, currentTimestep, grid.numX, grid.numY);
            err_code = MPI_Isend(&buffer[idx], grid.numX - 2, MPI_INT, 
                                 neighbourPids[2], tag, MPI_COMM_WORLD, &send_req[2]);

        }

        if(neighbourPids[3] != -1)
        {
            // Left
            idx = IDX(0, 1, currentTimestep, grid.numX, grid.numY);
            err_code = MPI_Irecv(&buffer[idx], 1, colSlice,
                                 neighbourPids[3], tag, MPI_COMM_WORLD, &request[3]);

            idx = IDX(1, 1, currentTimestep, grid.numX, grid.numY);
            err_code = MPI_Isend(&buffer[idx], 1, colSlice,
                             neighbourPids[3], tag, MPI_COMM_WORLD, &send_req[3]);
        }
        
        if(neighbourPids[4] != -1)
        {
            // Top right
            idx = IDX(grid.numX - 1, grid.numY - 1, currentTimestep, grid.numX, grid.numY);
            err_code = MPI_Irecv(&buffer[idx], 1, MPI_INT, 
                                 neighbourPids[4], tag, MPI_COMM_WORLD, &request[4]);

            idx = IDX(grid.numX - 2, grid.numY - 2, currentTimestep, grid.numX, grid.numY);
            err_code = MPI_Isend(&buffer[idx], 1, MPI_INT, 
                                 neighbourPids[4], tag, MPI_COMM_WORLD, &send_req[4]);
        }

        if(neighbourPids[5] != -1)
        {
            // Bottom right
            idx = IDX(grid.numX - 1, 0, currentTimestep, grid.numX, grid.numY);
            err_code = MPI_Irecv(&buffer[idx], 1, MPI_INT, 
                                 neighbourPids[5], tag, MPI_COMM_WORLD, &request[5]);

            
            idx = IDX(grid.numX - 2, 1, currentTimestep, grid.numX, grid.numY);
            err_code = MPI_Isend(&buffer[idx], 1, MPI_INT, 
                                 neighbourPids[5], tag, MPI_COMM_WORLD, &send_req[5]);
        }

        if(neighbourPids[6] != -1)
        {
            // Bottom left
            idx = IDX(0, 0, currentTimestep, grid.numX, grid.numY);
            err_code = MPI_Irecv(&buffer[idx], 1, MPI_INT, 
                                 neighbourPids[6], tag, MPI_COMM_WORLD, &request[6]);

            idx = IDX(1, 1, currentTimestep, grid.numX, grid.numY);
            err_code = MPI_Isend(&buffer[idx], 1, MPI_INT, 
                                 neighbourPids[6], tag, MPI_COMM_WORLD, &send_req[6]);
        }

        if(neighbourPids[7] != -1)
        {
            // Top left
            idx = IDX(0, grid.numY - 1, currentTimestep, grid.numX, grid.numY);
            err_code = MPI_Irecv(&buffer[idx], 1, MPI_INT, 
                                 neighbourPids[7], tag, MPI_COMM_WORLD, &request[7]);

            idx = IDX(1, grid.numY - 2, currentTimestep, grid.numX, grid.numY);
            err_code = MPI_Isend(&buffer[idx], 1, MPI_INT, 
                                 neighbourPids[7], tag, MPI_COMM_WORLD, &send_req[7]);
            
        }

        /*************************************************************/

        int num_processed = 0;
        // Do computation for [2 to end - 3]
        for(int y = 2; y < grid.numY - 2; y++)
            num_processed += simulatePartialRow(y, 2, grid.numX - 2, grid.numX, grid.numY, currentTimestep, buffer);

        // MPI wait all
        for(int i=0; i<8; i++)
        {
            if(neighbourPids[i] != -1)
                MPI_Wait(&request[i], &recv_status[i]);
        }


        // Do computation for remaining elements.
        num_processed += simulatePartialRow(1, 1, grid.numX - 1, grid.numX, grid.numY, currentTimestep, buffer);
        num_processed += simulatePartialRow(grid.numY - 2, 1, grid.numX - 1, grid.numX, grid.numY, currentTimestep, buffer);
        
        // x = {1, grid.numX-2}, y = [1- grid.numY-2]
        for(int y=2; y < grid.numY - 2; y++)
        {
            simulateHelper(1, y, grid.numX, grid.numY, currentTimestep, buffer);
            simulateHelper(grid.numX - 2, y, grid.numX, grid.numY, currentTimestep, buffer);
            num_processed += 2;
        }

        // printf("PID %d processed %d/%d\n", pid, num_processed, (grid.numX-2) * (grid.numY-2));

        currentTimestep = 1 - currentTimestep;
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Print the ans.
    for(int x = 1; x < grid.numX - 1; x++)
    {
        for(int y=1; y < grid.numY - 1; y++)
        {
            if(buffer[IDX(x, y, currentTimestep, grid.numX, grid.numY)] == 1)
            {
                printf("%d %d\n", x + grid.leftX - 1, y + grid.leftY - 1);
                fflush(stdout);
            }
        }
    }

    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;

}

void getGridWorkingDims(SubGrid *grid, int coord[2], int dims[2], int X, int Y, int pid)
{
    getWorkingDimHelper(coord[0], dims[0], X, &((*grid).leftX), (&(*grid).rightX));
    getWorkingDimHelper(coord[1], dims[1], Y, &((*grid).leftY), (&(*grid).rightY));

    ((*grid).numX) = ((*grid).rightX) - ((*grid).leftX) + 3;
    ((*grid).numY) = ((*grid).rightY) - ((*grid).leftY) + 3;

    ((*grid).pid) = pid;

    ((*grid).Px) = coord[0];
    ((*grid).Py) = coord[1];
}

void simulateHelper(int x, int y, int maxX, int maxY, int currentTimestep, int *buffer)
{
    int num_neighbours = getNeighbours(x, y, maxX, maxY, currentTimestep, buffer);
    
    if(buffer[IDX(x, y, currentTimestep, maxX, maxY)] == 1)
    {
        // Cell was alive.
        if(num_neighbours == 2 || num_neighbours == 3)
            buffer[IDX(x, y, 1 - currentTimestep, maxX, maxY)] = 1;
        else
            buffer[IDX(x, y, 1 - currentTimestep, maxX, maxY)] = 0;
    }
    else
    {
        if(num_neighbours == 3)
            buffer[IDX(x, y, 1 - currentTimestep, maxX, maxY)] = 1;
        else
            buffer[IDX(x, y, 1 - currentTimestep, maxX, maxY)] = 0;
    }
}

int simulatePartialRow(int y, int fromX, int tillX, int maxX, int maxY, int currentTimestep, int *buffer)
{
    int num_processed = 0;
    for(int x = fromX; x < tillX; x++)
    {
        simulateHelper(x, y, maxX, maxY, currentTimestep, buffer);
        num_processed += 1;
    }
    return num_processed;
}

void getWorkingDimHelper(int pidDir, int numProcessorsDir, int sizeDir, int *lo, int *high)
{
    // Assumption: number of processes < sizeDir.
    int temp1 = sizeDir / numProcessorsDir;
    int temp2 = sizeDir % numProcessorsDir;
    if (pidDir < temp2)
    {    
        // return temp1 + 1;
        *lo = pidDir * (temp1 + 1);
        (*high) = (pidDir + 1) * (temp1 + 1) - 1;
        return;
    }
    else
    {
        *lo = temp2 * (temp1 + 1) + (pidDir - temp2) * temp1;
        (*high) = temp2 * (temp1 + 1) + ((pidDir + 1 - temp2) * temp1) - 1;
        return;
    }
}

bool isInsideGrid(int* coord, SubGrid grid)
{
    if(coord[0] < grid.leftX || coord[0] > grid.rightX)
        return false;
    if(coord[1] < grid.leftY || coord[1] > grid.rightY)
        return false;
    return true;
}

int isInside(int x, int y, int X, int Y)
{
    if (x < 0 || x >= X)
        return 0;
    if(y < 0 || y >= Y)
        return 0;
    return 1;
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

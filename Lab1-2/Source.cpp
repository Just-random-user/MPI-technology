#include <mpi.h>
#include <stdio.h>
#include <iostream>
using namespace std;

// Variant 0: ring
int main(int argc, char* argv[])
{
	int ProcNum, ProcRank, Value = 256;
	MPI_Status Status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	int num = 0;
	if (ProcRank == 0)
	{
		cout << "Enter the number of iteraions: ";
		while (num < 1)
			cin >> num;
	}
	MPI_Bcast(&num, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	for (int i = 0; i < num; i++)
	{
		if (ProcRank == 0)
			MPI_Send(&Value, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);

		MPI_Recv(&Value, 1, MPI_INT, (ProcRank + ProcNum - 1) % ProcNum,
			MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
		cout << "Process " << ProcRank << " recieved the data" << endl;
		if (ProcRank != 0)
			MPI_Send(&Value, 1, MPI_INT, (ProcRank + 1) % ProcNum, 0, MPI_COMM_WORLD);
	}
	MPI_Finalize();
	return 0;
}
#include <mpi.h>
#include <stdio.h>
#include <iostream>
using namespace std;

// Variant 0: ring
int main(int argc, char* argv[])
{
	int ProcNum, ProcRank, Value = 1;
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
		for (int j = 0; j < ProcNum; j++)
		{
			if (ProcRank == j)
				cout << "Process " << ProcRank << " is now sending a data" << endl;
			MPI_Bcast(&Value, 1, MPI_INT, j, MPI_COMM_WORLD);
			if (Value == ProcRank)
			{
				Value = (Value + 1) % ProcNum;
				cout << "Process " << ProcRank << " has gained a data!";
				cout << endl;
			}
		}
	}
	MPI_Finalize();
	return 0;
}
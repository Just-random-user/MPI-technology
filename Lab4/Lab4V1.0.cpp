#include <mpi.h>
#include <stdio.h>
#include <complex>
#include <iostream>
#include <time.h>
#include <random>
#include <vector>
using namespace std;

static bool debugging = true;

complex<double>** MatrixMul(complex<double>** A, complex<double>** B, complex<double>** C, int size);
complex<double>** MatrixFill(complex<double>** matrix, int size);
complex<double>** MatrixCopy(complex<double>** from, complex<double>** to, int size);
void PrintMatrix(complex<double>** matrix, int size);
void Synchronyze(complex<double>*** matrixes, int size, int amount, int ProcNum, MPI_Datatype MY_CVECTOR);

// Variant 0: complex matrix
int main(int argc, char* argv[])
{
	int ProcNum, ProcRank, Value = 1;
	//MPI_Status Status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	//main data initializing
	int size = 0, amount = 0;

	srand(ProcRank);
	if (ProcRank == 0)
	{
		cout << "Enter the size of the matrix and their amount: ";
		fflush(stdout);
		cin >> size >> amount;
	}
	MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&amount, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//memory initializing
	complex<double>*** matrArray;
	matrArray = new complex<double>**[amount];
	for (int i = 0; i < amount; i++)
	{
		matrArray[i] = new complex<double>*[size];
		for (int j = 0; j < size; j++)
		{
			matrArray[i][j] = new complex<double>[size];
		}
	}


	//matrix filling
	for (int i = 0; i < amount; i++)
	{
		if (ProcRank == i % ProcNum)
		{
			matrArray[i] = MatrixFill(matrArray[i], size);
			cout << "Process " << ProcRank << " has finished the matrix " << i << endl;
			fflush(stdout);
		}
	}

	//type commiting
	MPI_Datatype MY_CVECTOR;
	MPI_Type_contiguous(size, MPI_DOUBLE_COMPLEX, &MY_CVECTOR);
	MPI_Type_commit(&MY_CVECTOR);

	Synchronyze(matrArray, size, amount, ProcNum, MY_CVECTOR);

	//output
	if (ProcRank == 0)
	{
		for (int i = 0; i < amount; i++)
		{
			cout << endl << "Matrix " << i << endl;
			PrintMatrix(matrArray[i], size);
		}
	}
	fflush(stdout);

	int currAmount = amount;
	int oldAmount = amount;

	//recursive matrix multiplication (by pairs)
	while (currAmount != 1)
	{
		currAmount = ceil((double)currAmount / 2);
		complex<double>*** matrArray1;
		matrArray1 = new complex<double>**[currAmount];
		for (int i = 0; i < currAmount; i++)
		{
			matrArray1[i] = new complex<double>*[size];
			for (int j = 0; j < size; j++)
			{
				matrArray1[i][j] = new complex<double>[size];
			}
		}
		int counter = 0;
		for (int i = 0; i < currAmount; i++)
		{
			if (oldAmount % 2 != 0 && counter + 1 == oldAmount)
			{
				if (ProcRank == i % ProcNum)
				{
					MatrixCopy(matrArray[counter], matrArray1[i], size);
				}
				break;
			}
			if (ProcRank == i % ProcNum)
			{
				matrArray1[i] = MatrixMul(matrArray[counter], matrArray[counter + 1], matrArray1[i], size);
			}
			counter += 2;
		}

		for (int i = 0; i < oldAmount; i++)
		{
			for (int j = 0; j < size; j++)
			{
				delete[](matrArray[i][j]);
			}
			delete[](matrArray[i]);
		}

		oldAmount = currAmount;
		matrArray = matrArray1;

		Synchronyze(matrArray, size, currAmount, ProcNum, MY_CVECTOR);

		//every-iteration output for answer check
		if (ProcRank == 0 && currAmount > 1 && debugging == true)
		{
			for (int i = 0; i < currAmount; i++)
			{
				cout << endl << "middle matrix " << i << ": " << endl;
				PrintMatrix(matrArray[i], size);
			}
		}
	}

	//result output
	if (ProcRank == 0)
	{
		cout << endl << "result: " << endl;
		PrintMatrix(matrArray[0], size);
	}

	for (int i = 0; i < currAmount; i++)
	{
		for (int j = 0; j < size; j++)
		{
			delete[](matrArray[i][j]);
		}
		delete[](matrArray[i]);
	}
	delete[](matrArray);
	MPI_Type_free(&MY_CVECTOR);
	MPI_Finalize();
	return 0;
}

complex<double>** MatrixMul(complex<double>** A, complex<double>** B,
	complex<double>** C, int size)
{
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
		{
			C[i][j] = 0;
			for (int k = 0; k < size; k++)
				C[i][j] += A[i][k] * B[k][j];
		}
	return C;
}

complex<double>** MatrixCopy(complex<double>** from, complex<double>** to, int size)
{
	for (int j = 0; j < size; j++)
	{
		for (int k = 0; k < size; k++)
			to[j][k] = from[j][k];
	}
	return from;
}

complex<double>** MatrixFill(complex<double>** matrix, int size)
{
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
		{
			matrix[i][j] = complex<double>(int(rand() % 100), int(rand() % 100));
		}
	return matrix;
}

void PrintMatrix(complex<double>** matrix, int size)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			cout << matrix[i][j].real();
			if (matrix[i][j].imag() >= 0)
				cout << "+";
			cout << matrix[i][j].imag() << "i ";
		}
		cout << endl;
	}
	return;
}

void Synchronyze(complex<double>*** matrixes, int size, int amount, int ProcNum, MPI_Datatype MY_CVECTOR)
{
	for (int i = 0; i < amount; i++)
	{
		for (int j = 0; j < size; j++)
			MPI_Bcast(matrixes[i][j], 1, MY_CVECTOR, i % ProcNum, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	}
}
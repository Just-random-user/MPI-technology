#include "mpi.h"
#include "fft.h"
#include "Header.h"

#include <iostream>
#include <stdio.h>
#include <locale.h>
#include <time.h>
#include <iomanip>
#include <complex>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>

using namespace std;


vector<int> SchonhageStrassenMultiplication(std::vector<int> source_vec_1, std::vector<int> source_vec_2)
{
	reverse(source_vec_1.begin(), source_vec_1.end());
	reverse(source_vec_2.begin(), source_vec_2.end());
	while (source_vec_2.size() < source_vec_1.size())
	{
		source_vec_2.push_back(0);
	}
	int Vec_size = source_vec_1.size() + source_vec_2.size() + 1;
	for (int i = source_vec_1.size(); i < Vec_size; i++)
	{
		source_vec_1.push_back(0);
		source_vec_2.push_back(0);
	}
	std::vector<std::complex<double>> vec_1(source_vec_1.begin(), source_vec_1.end());
	std::vector<std::complex<double>> vec_2(source_vec_2.begin(), source_vec_2.end());
	std::vector<std::complex<double>> v;

	std::vector<std::complex<double>> fft_1 = fft::fft(vec_1);
	std::vector<std::complex<double>> fft_2 = fft::fft(vec_2);
	for (int i = 0; i < fft_2.size(); ++i) {
		v.push_back(fft_1[i] * fft_2[i]);
	}
	std::vector<std::complex<double>> linearConvolution = fft::ifft(v);
	std::vector<int> round_to_int = fft::round(linearConvolution);
	int res = 0;
	large_num result(0);
	large_num POW(1);

	for (int i = 0; i < round_to_int.size(); ++i) {
		if (i == 0)
		{
			result = result + (large_num(round_to_int[i]) * large_num((int)pow(10, i)));
		}
		else {
			POW = POW * 10;
			result = result + (large_num(round_to_int[i]) * POW);
		}
	}
	vector<int> result_vector = result.resVector();
	return result_vector;
}


int main(int argc, char* argv[]) {
	int ProcNum, ProcRank, RecvRank;
	MPI_Status Status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	int world_proc_size_copy = ProcNum;

	MPI_Group new_group;
	MPI_Group original_group;
	int ranks_for_new_group[10] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	MPI_Comm_group(MPI_COMM_WORLD, &original_group);
	MPI_Group_incl(original_group, 10, ranks_for_new_group, &new_group);

	MPI_Comm new_comm;
	MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_comm);
	int periods = 0;
	MPI_Cart_create(MPI_COMM_WORLD, 1, &world_proc_size_copy, &periods, 0, &new_comm);
	
	std::vector<int> source_vec_1{ 1, 2, 3, 4 };
	std::vector<int> source_vec_2{ 4, 3, 2, 1 };
	std::vector<int> vec_buff;
	std::vector<int> vec_buff_1;
	std::vector<int> vec_buff_2;
	std::vector<int> vec_buff_3;
	std::vector<int> vec_buff_4;
	std::vector<int> vec_buff_5;
	std::vector<int> vec_buff_6;
	std::vector<int> vec_buff_7;
	std::vector<int> vec_buff_8;
	std::vector<int> vec_buff_9;

	vec_buff_1 = SchonhageStrassenMultiplication(source_vec_1, source_vec_2);

	MPI_Datatype type_1; // создание нового типа
	MPI_Type_contiguous(vec_buff_1.size(), MPI_INT, &type_1);
	MPI_Type_commit(&type_1);

	int *buff_2 = new int[vec_buff_1.size()];
	int *buff_1 = new int[vec_buff_1.size()];
	for (int i = 0; i < vec_buff_1.size(); i++)
		buff_1[i] = vec_buff_1[i];

	for (int i = 0; i < ProcNum; i++)
	{
		if ((ProcRank == i) && (i % 2 != 0) && (ProcRank != 0))
			MPI_Send(buff_1, 1, type_1, ProcRank - 1, 1, new_comm);
		else if ((ProcRank == i) && (i % 2 == 0))
			MPI_Recv(buff_2, 1, type_1, ProcRank + 1, MPI_ANY_TAG, new_comm, &Status);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (ProcRank % 2 == 0)
	{
		for (int i = 0; i < vec_buff_1.size(); i++)
			vec_buff_2.push_back(buff_2[i]);
		vec_buff_3 = SchonhageStrassenMultiplication(vec_buff_1, vec_buff_2);// 0  2  4  6  8 
	}

	MPI_Barrier(MPI_COMM_WORLD);

	int *buff_3 = new int[vec_buff_3.size()];
	int *buff_4 = new int[vec_buff_3.size()];


	MPI_Barrier(MPI_COMM_WORLD);

	if (ProcRank % 2 == 0)
	{
		for (int i = 0; i < vec_buff_3.size(); i++)
		{
			buff_3[i] = vec_buff_3[i];
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Datatype type_2; // создание нового типа
	MPI_Type_contiguous(vec_buff_3.size(), MPI_INT, &type_2);
	MPI_Type_commit(&type_2);

	MPI_Barrier(MPI_COMM_WORLD);

	for (int i = 0; i < ProcNum; i++)
	{
		if (((ProcRank == i) && (ProcRank == 2)) || ((ProcRank == i) && (ProcRank == 6)))
			MPI_Send(buff_3, 1, type_2, ProcRank - 2, 1, new_comm);
		else if (((ProcRank == i) && (ProcRank == 0)) || ((ProcRank == i) && (ProcRank == 4)))
			MPI_Recv(buff_4, 1, type_2, ProcRank + 2, MPI_ANY_TAG, new_comm, &Status);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (ProcRank == 0 || ProcRank == 4)
	{
		for (int i = 0; i < vec_buff_3.size(); i++)
			vec_buff_4.push_back(buff_4[i]);
		vec_buff_5 = SchonhageStrassenMultiplication(vec_buff_3, vec_buff_4);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Datatype type_3; // создание нового типа
	MPI_Type_contiguous(vec_buff_5.size(), MPI_INT, &type_3);
	MPI_Type_commit(&type_3);

	int *buff_5 = new int[vec_buff_5.size()];
	int *buff_6 = new int[vec_buff_5.size()];

	MPI_Barrier(MPI_COMM_WORLD);

	if (ProcRank == 4 )
	{
		for (int i = 0; i < vec_buff_5.size(); i++)
			buff_5[i] = vec_buff_5[i];
		MPI_Send(buff_5, 1, type_3, 0, 1, new_comm);
	}
	if (ProcRank == 0)
	{
		MPI_Recv(buff_6, 1, type_3, 4, MPI_ANY_TAG, new_comm, &Status);
		for (int i = 0; i < vec_buff_5.size(); i++)
			vec_buff_6.push_back(buff_6[i]);
		vec_buff_7 = SchonhageStrassenMultiplication(vec_buff_5, vec_buff_6);// 0 4 8
	}

	MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Datatype type_4; // создание нового типа
	MPI_Type_contiguous(vec_buff_7.size(), MPI_INT, &type_4);
	MPI_Type_commit(&type_4);


	int *buff_7 = new int[vec_buff_3.size()];
	int *buff_8 = new int[vec_buff_3.size()];


	MPI_Barrier(MPI_COMM_WORLD);

	if (ProcRank == 8)
	{
		for (int i = 0; i < vec_buff_3.size(); i++)
			buff_7[i] = vec_buff_3[i];
		MPI_Send(buff_7, 1, type_2, 0, 1, new_comm);
	}

	if (ProcRank == 0)
	{
		MPI_Recv(buff_8, 1, type_2, 8, MPI_ANY_TAG, new_comm, &Status);
		for (int i = 0; i < vec_buff_3.size(); i++)
		{
			vec_buff_8.push_back(buff_8[i]);
		}
	
		vec_buff_9 = SchonhageStrassenMultiplication(vec_buff_7, vec_buff_8);// 0 4 8
		cout << "Result: ";
		for (int i = 0; i < vec_buff_9.size(); i++)
			cout << vec_buff_9[i];
	}

	delete[] buff_1;
	delete[] buff_2;
	delete[] buff_3;
	delete[] buff_4;
	delete[] buff_5;
	delete[] buff_6;
	delete[] buff_7;
	delete[] buff_8;

	MPI_Type_free(&type_1);
	MPI_Group_free(&new_group);
	MPI_Comm_free(&new_comm);

	MPI_Finalize();

	return 0;
}


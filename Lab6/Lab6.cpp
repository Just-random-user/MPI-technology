#include "mpi.h"
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

static double constexpr pi = 3.14159265358979323846;

int upper_log2(const int x)
{
	for (int i = 0; i < 30; ++i)
	{
		const int temp = 1 << i;
		if (temp >= x)
		{
			return i;
		}
	}
	return 30;
}


const std::vector<std::complex<double>> phase_vec(const int len)
{
	std::vector<std::complex<double>> res(len);
	const double radius = 1;
	for (int i = 0; i < len; ++i)
	{
		const double phase = -2 * pi * i / len;
		res[i] = std::polar(radius, phase);
	}
	return res;
}


void forward(std::vector<std::complex<double>>& prev, std::vector<std::complex<double>>& temp,
	const std::vector<std::complex<double>>& phases, const int turn, const int n_bits)
{
	if (turn == n_bits)
	{
		return;
	}

	const int group_size = 1 << (turn + 1); // size of butterfly group
	const int num_groups = prev.size() / group_size;
	const int phase_angular_freq = num_groups;
	for (int i_group = 0; i_group < num_groups; ++i_group)
	{
		const int base_index = i_group * group_size;
		// iterate through within the butterfly group
		for (int j = 0; j < group_size / 2; ++j)
		{
			const int x0_index = base_index + j;
			const int x1_index = base_index + group_size / 2 + j;
			prev[x1_index] *= phases[j * phase_angular_freq];
			temp[x0_index] = prev[x0_index] + prev[x1_index];
			temp[x1_index] = prev[x0_index] - prev[x1_index];
		}
	}

	forward(temp, prev, phases, turn + 1, n_bits);
}


void bit_reversal_permutation(std::vector<std::complex<double>>& vec, const int n_bits)
{
	if (vec.size() <= 2)
	{
		return;
	}
	if (vec.size() == 4)
	{
		std::swap(vec[1], vec[3]);
		return;
	}
	std::vector<int> bit_rerversal(vec.size());
	// initialize the first 4 elements
	bit_rerversal[0] = 0;
	bit_rerversal[1] = 1 << (n_bits - 1);
	bit_rerversal[2] = 1 << (n_bits - 2);
	bit_rerversal[3] = bit_rerversal[1] + bit_rerversal[2];


	for (int k = 3; k <= n_bits; ++k)
	{
		const int nk = (1 << k) - 1;          // n_k = 2^k - 1
		const int n_km1 = (1 << (k - 1)) - 1; // n_(k-1) = 2^(k-1) - 1
		bit_rerversal[nk] = bit_rerversal[n_km1] + (1 << (n_bits - k));
		for (int i = 1; i <= n_km1; ++i)
		{
			bit_rerversal[nk - i] = bit_rerversal[nk] - bit_rerversal[i];
		}
	}

	// permute the input vector according to bit_rerversal[]
	for (int i = 0; i < vec.size(); ++i)
	{
		if (bit_rerversal[i] > i)
		{
			std::swap(vec[i], vec[bit_rerversal[i]]);
		}
	}
}


std::vector<std::complex<double>> fft(const std::vector<std::complex<double>>& inputs)
{
	if (inputs.empty())
	{
		return {};
	}
	const int n_bits = upper_log2(inputs.size());
	const int len = 1 << n_bits;
	const std::vector<std::complex<double>> phases = phase_vec(len);
	std::vector<std::complex<double>> prev(len);
	std::vector<std::complex<double>> temp(len);
	std::copy(inputs.begin(), inputs.end(), prev.begin());
	bit_reversal_permutation(prev, n_bits);
	// butterfly forwarding from input to output, Cooley Tukey algorithm
	forward(prev, temp, phases, 0, n_bits);
	return (n_bits % 2 == 1) ? temp : prev;
}
std::vector<std::complex<double>> ifft(const std::vector<std::complex<double>>& inputs)
{
	std::vector<std::complex<double>> reverse_freq_spectrum(inputs);
	std::reverse(std::next(reverse_freq_spectrum.begin()), reverse_freq_spectrum.end());

	const double len = reverse_freq_spectrum.size();
	std::transform(reverse_freq_spectrum.begin(), reverse_freq_spectrum.end(), reverse_freq_spectrum.begin(),
		[len](const std::complex<double>& num) { return num / len; });
	// fft
	return fft(reverse_freq_spectrum);
}
std::vector<int> round(const std::vector<std::complex<double>>& vec)
{
	std::vector<int> res(vec.size());
	std::transform(vec.begin(), vec.end(), res.begin(), [](const std::complex<double>& num) -> int { return std::round(num.real()); });
	return res;
}
std::vector<double> real(const std::vector<std::complex<double>>& vec)
{
	std::vector<double> res(vec.size());
	std::transform(vec.begin(), vec.end(), res.begin(), [](const std::complex<double>& num) -> double { return num.real(); });
	return res;
}

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

	std::vector<std::complex<double>> fft_1 = fft(vec_1);
	std::vector<std::complex<double>> fft_2 = fft(vec_2);
	for (int i = 0; i < fft_2.size(); ++i) {
		v.push_back(fft_1[i] * fft_2[i]);
	}
	std::vector<std::complex<double>> linearConvolution = ifft(v);
	std::vector<int> round_to_int = round(linearConvolution);
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


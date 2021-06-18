#include <mpi.h>
#include <stdio.h>
#include <complex>
#include <iostream>
#include <time.h>
#include <random>
#include <vector>
#include <string>
using namespace std;

class my_long
{
public:
	string digit;
	my_long()
	{
		digit = "";
	}
	my_long(string pattern)
	{
		digit = pattern;
	}
	~my_long() {}
	my_long operator+ (const my_long& summand)
	{
		my_long ans;
		int curr_digit = 0;
		int carry_flag = 0;
		string temp1 = digit, temp2 = summand.digit;

		string* str_ptr_grater, *str_ptr_lesser;
		if (temp1.size() > temp2.size())
		{
			str_ptr_grater = &temp1;
			str_ptr_lesser = &temp2;
		}
		else
		{
			str_ptr_grater = &temp2;
			str_ptr_lesser = &temp1;
		}
		//to the same size
		reverse(str_ptr_lesser->begin(), str_ptr_lesser->end());
		while (str_ptr_lesser->size() < str_ptr_grater->size())
			str_ptr_lesser->push_back('0');
		reverse(str_ptr_lesser->begin(), str_ptr_lesser->end());

		for (int i = 1; i <= temp1.size(); i++)
		{
			curr_digit = temp1.at(temp1.size() - i) + temp2.at(temp2.size() - i) - '0';
			curr_digit += carry_flag;
			carry_flag = (curr_digit - '0') / 10;
			curr_digit -= carry_flag * 10;

			ans.digit.push_back(curr_digit);
		}
		if (carry_flag != 0)
			ans.digit.push_back(carry_flag + '0');
		reverse(ans.digit.begin(), ans.digit.end());
		return ans;
	}
	my_long operator- (const my_long& summand)
	{
		my_long ans;
		int curr_digit = 0;
		int carry_flag = 0;
		string temp1 = digit, temp2 = summand.digit;

		string* str_ptr_grater, * str_ptr_lesser;
		if (temp1.size() > temp2.size())
		{
			str_ptr_grater = &temp1;
			str_ptr_lesser = &temp2;
		}
		else
		{
			str_ptr_grater = &temp2;
			str_ptr_lesser = &temp1;
		}
		//to the same size
		reverse(str_ptr_lesser->begin(), str_ptr_lesser->end());
		while (str_ptr_lesser->size() < str_ptr_grater->size())
			str_ptr_lesser->push_back('0');
		reverse(str_ptr_lesser->begin(), str_ptr_lesser->end());

		for (int i = 1; i <= temp1.size(); i++)
		{
			curr_digit = temp1.at(temp1.size() - i) - temp2.at(temp2.size() - i) + '0';
			curr_digit -= carry_flag;
			if (curr_digit < '0')
			{
				curr_digit += 10;
				carry_flag = 1;
			}
			else
				carry_flag = 0;

			ans.digit.push_back(curr_digit);
		}
		if (carry_flag != 0)
			ans.digit.push_back(carry_flag + '0');
		reverse(ans.digit.begin(), ans.digit.end());
		return ans;
	}
	void to_size(size_t size)
	{
		reverse(digit.begin(), digit.end());
		while (digit.size() < size)
			digit.push_back('0');
		reverse(digit.begin(), digit.end());
	}
};

my_long Karatsuba(my_long LA, my_long LB);
my_long multiply_my_long(my_long LA, my_long LB);
string long_multiplication(string a, string b);
my_long clear_zeroes_my_long(my_long digit);
void synchronyze(my_long* arr, size_t size, int ProcRank, int ProcNum);
void synchronyzefrom(my_long* arr, size_t size, int ProcRank, int root); 
void string_to_array(string st, char* arr);
string array_to_string(char* arr, int amount);

// Variant 0: Karatsuba
int main(int argc, char* argv[])
{
	int ProcNum, ProcRank, Value = 1;
	//MPI_Status Status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Datatype MPI_BIGINT;

	my_long *a = nullptr, *temp, c;
	size_t amount = 0, temp_amount;

	if (ProcRank == 0)
	{
		std::cout << "Enter the amount of digits: ";
		cin >> amount;
	}
	MPI_Bcast(&amount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	a = new my_long[amount];
	if (ProcRank == 0)
	{
		for (int i = 0; i < amount; i++)
		{
			std::cout << "Enter the " << i << " digit: ";
			cin >> a[i].digit;
		}
	}

	synchronyzefrom(a, amount, ProcRank, 0);


	//if (ProcRank == 2)
	//	cout << a[0].digit << endl;

	while (amount > 1)
	{
		if (ProcRank == 0)
			cout << "=============================" << endl;
		temp_amount = ceil((double)amount / 2);
		temp = new my_long[temp_amount];
		
		size_t counter = 0;
		for (int i = 0; i < temp_amount; i++)
		{
			if (counter + 1 == amount && ProcRank == i % ProcNum)
			{
				temp[i].digit = a[counter].digit;
				cout << "process " << i << " copied digit N" << i << endl;
				break;
			}
			if (ProcRank == i % ProcNum)
			{
				temp[i] = multiply_my_long(a[counter], a[counter + 1]);
				cout << "process " << i << " finished digit N" << i << endl;
			}
			counter += 2; 
		}

		delete[] a;
		a = temp;
		amount = temp_amount;
		synchronyze(a, amount, ProcRank, ProcNum);
	}

	if (ProcRank == 0)
		cout << "=============================" << endl;

	c = a[0];
	delete[] a;
	c = clear_zeroes_my_long(c);
	if (ProcRank == 0)
		std::cout << "Answer: " << c.digit << endl;

	MPI_Finalize();
	return 0;
}

my_long Karatsuba(my_long LA, my_long LB)
{
	size_t max = LA.digit.size() > LB.digit.size() ? LA.digit.size() : LB.digit.size();
	LA.to_size(max);
	LB.to_size(max);
	size_t m = LA.digit.size() / 2;
	my_long A0, A1, B0, B1;

	A0.digit.assign(LA.digit.begin(), LA.digit.end() - m);
	A1.digit.assign(LA.digit.end() - m, LA.digit.end());
	B0.digit.assign(LB.digit.begin(), LB.digit.end() - m);
	B1.digit.assign(LB.digit.end() - m, LB.digit.end());

	my_long A0_B0 = multiply_my_long(A0, B0);
	my_long A1_B1 = multiply_my_long(A1, B1);
	my_long mid = multiply_my_long((A0 + A1), (B0 + B1)) - A0_B0 - A1_B1;

	for (int i = 0; i < m; i++)
		mid.digit.push_back('0');

	my_long end;
	end.digit = A0_B0.digit;

	for (int i = 0; i < 2 * m; i++)
		end.digit.push_back('0');

	my_long ans;
	ans = A1_B1 + mid + end;

	return ans;
}

string long_multiplication(string a, string b)
{
	if (a.size() != 1 || b.size() != 1)
		return "";
	return to_string(stoi(a) * stoi(b));
}

my_long multiply_my_long(my_long LA, my_long LB)
{
	my_long ans;
	if (LA.digit.size() == 1 && LB.digit.size() == 1)
		ans.digit = long_multiplication(LA.digit, LB.digit);
	else
		ans = Karatsuba(LA, LB);

	return ans;
}

my_long clear_zeroes_my_long(my_long digit)
{
	reverse(digit.digit.begin(), digit.digit.end());
	while (digit.digit.at(digit.digit.size() - 1) == '0')
		digit.digit.pop_back();
	reverse(digit.digit.begin(), digit.digit.end());
	return digit;
}

void synchronyze(my_long* arr, size_t size, int ProcRank, int ProcNum)
{
	MPI_Datatype MPI_BIGINT;
	for (int i = 0; i < size; i++)
	{
		int SZ;
		if (ProcRank == i % ProcNum)
		{
			SZ = arr[i].digit.size();
		}
		MPI_Bcast(&SZ, 1, MPI_INT, i % ProcNum, MPI_COMM_WORLD);

		MPI_Type_contiguous(SZ, MPI_CHAR, &MPI_BIGINT);
		MPI_Type_commit(&MPI_BIGINT);
		char* buff = new char[SZ];
		if (ProcRank == i % ProcNum)
			string_to_array(arr[i].digit, buff);

		MPI_Bcast(buff, 1, MPI_BIGINT, i % ProcNum, MPI_COMM_WORLD);

		arr[i].digit = array_to_string(buff, SZ);
		delete[] buff;
		MPI_Type_free(&MPI_BIGINT);
	}
}

void synchronyzefrom(my_long* arr, size_t size, int ProcRank, int root)
{
	MPI_Datatype MPI_BIGINT;
	for (int i = 0; i < size; i++)
	{
		int SZ;
		if (ProcRank == root)
		{
			SZ = arr[i].digit.size();
		}
		MPI_Bcast(&SZ, 1, MPI_INT, root, MPI_COMM_WORLD);

		MPI_Type_contiguous(SZ, MPI_CHAR, &MPI_BIGINT);
		MPI_Type_commit(&MPI_BIGINT);
		char* buff = new char[SZ];
		if (ProcRank == root)
			string_to_array(arr[i].digit, buff);

		MPI_Bcast(buff, 1, MPI_BIGINT, root, MPI_COMM_WORLD);

		arr[i].digit = array_to_string(buff, SZ);
		delete[] buff;
		MPI_Type_free(&MPI_BIGINT);
	}
}

void string_to_array(string st, char* arr)
{
	for (int i = 0; i < st.size(); i++)
		arr[i] = st.at(i);
}

string array_to_string(char* arr, int amount)
{
	string st;
	for (int i = 0; i < amount; i++)
		st.push_back(arr[i]);
	return st;
}

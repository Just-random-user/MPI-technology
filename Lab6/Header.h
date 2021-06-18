#pragma once
#pragma once
#include <string>
#include <vector>
#include <stdexcept>
#include <exception>

using namespace std;

class large_num
{
public:
	large_num();
	large_num(const string&);
	large_num(const vector<int>&);
	large_num(const int&);
	large_num(const long long int&);

	vector<int> resVector() const noexcept;

	large_num(const large_num&);
	large_num& operator = (const large_num&);

	const large_num operator +() const;
	const large_num operator -() const;

	bool operator ==(const large_num&) const noexcept;
	bool operator != (const large_num&) const noexcept;
	bool operator <(const large_num&) const noexcept;
	bool operator <= (const large_num&) const noexcept;
	bool operator >(const large_num&) const noexcept;
	bool operator >= (const large_num&) const noexcept;

	large_num operator +(const large_num&)const noexcept;
	large_num operator -(const large_num&) const noexcept;
	large_num operator *(const large_num&)const noexcept;
	large_num operator /(const large_num&) const;

	large_num operator%(const large_num&)const noexcept;

	large_num gcd(const large_num&);
	large_num random();
	large_num inverse_number() const noexcept;

	void print()const noexcept;

private:
	class DivideByZeroException : public exception {};
	vector<int> digits;
	bool is_positive;

	void rm_zeros();
	void shift_right();
	void shift_left();

	bool odd()const noexcept;///checks if the current number is even
	bool even()const noexcept;///checks if the current number is odd

	void extend_vector(vector<int>& a, int n);
	void fft(vector<int> & a, bool invert);

	large_num pow(large_num)const;
	large_num sqrt_bi();
	large_num to_binary() const noexcept;
	large_num log2_ceil() const;

	vector<int> mul(const vector<int>& a, const vector<int>& b);
	void vect_rm_zeros(vector<int>& a);
	vector<int> to_vec(const string& a);
	vector<int> pow(const vector<int>& a, long long n);

	large_num schonhage_strassen(const large_num&) noexcept;
};
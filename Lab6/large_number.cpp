#pragma once
#include "Header.h"
#include <algorithm>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <math.h>

#define PI 3.14

const int MIN_LEN_ = 216;

large_num::large_num() {
	this->is_positive = true;
}

large_num::large_num(const string& number) {
	string str = number;
	if (str.length() == 0) {
		this->is_positive = true;
	}
	else {
		if (str[0] == '-') {
			str = str.substr(1);
			this->is_positive = false;
		}
		else this->is_positive = true;

		for (int i = str.length() - 1; i >= 0; --i) {
			this->digits.push_back((int)str[i] - '0');
		}

		this->rm_zeros();
	}
}

large_num::large_num(const vector<int>& vec) {  //not reversed
	this->digits = vec;
	this->is_positive = true;
}

large_num::large_num(const large_num& other) {
	this->digits = other.digits;
	this->is_positive = other.is_positive;
}

large_num::large_num(const int& n) {
	if (n < 0) this->is_positive = false;
	else this->is_positive = true;

	int d = n;
	while (d != 0) {
		this->digits.push_back(d % 10);
		d /= 10;
	}
}

large_num::large_num(const long long& n) {
	if (n < 0) this->is_positive = false;
	else this->is_positive = true;

	long long d = n;
	while (d != 0) {
		this->digits.push_back(d % 10);
		d /= 10;
	}
}
std::vector<int> large_num::resVector() const noexcept {  //not reversed
	std::vector<int> vec;
	for (int i = (int)this->digits.size() - 1; i >= 0; i--) {
		vec.push_back(this->digits[i]);
	}
	return vec;
}

large_num& large_num::operator = (const large_num& other) {
	this->digits = other.digits;
	this->is_positive = other.is_positive;
	return (*this);
}

void large_num::rm_zeros() {
	while (this->digits.size() > 1 && this->digits.back() == 0) {
		this->digits.pop_back();
	}

	if (this->digits.size() == 1 && this->digits[0] == 0) this->is_positive = true;
}

void large_num::print() const noexcept {
	if (!this->is_positive) { cout << "-"; }

	for (int i = (int)this->digits.size() - 1; i >= 0; i--) {
		cout << this->digits[i];
	}
}

bool large_num::operator ==(const large_num& other) const noexcept {
	if (this->is_positive != other.is_positive) return false;

	if (this->digits.empty()) {
		if (other.digits.empty() || (other.digits.size() == 1 && other.digits[0] == 0))
		{
			return true;
		}
		else  return false;
	}

	if (other.digits.empty()) {
		if (this->digits.size() == 1 && this->digits[0] == 0)
		{
			return true;
		}
		else  return false;
	}

	if (this->digits.size() != other.digits.size()) return false;

	for (unsigned i = 0; i < this->digits.size(); i++) {
		if (this->digits[i] != other.digits[i]) { return false; }
	}

	return true;
}

bool large_num::operator !=(const large_num& other) const noexcept {
	return !(*this == other);
}

//return copy number
const large_num large_num::operator +() const {
	return large_num(*this);
}

//return reversed number (with reversed is_positive)
const large_num large_num::operator -() const {
	large_num copy(*this);
	copy.is_positive = !copy.is_positive;
	return copy;
}

large_num large_num::operator + (const large_num& other) const noexcept {
	if (!this->is_positive) {
		if (!other.is_positive) return -(-(*this) + (-other));
		else return (other - (-(*this)));
	}
	else if (!other.is_positive) return ((*this) - (-other));

	vector<int> res;
	int carry = 0;

	for (unsigned i = 0; i < max(this->digits.size(), other.digits.size()) || carry != 0; ++i) {
		int total = carry + (i < this->digits.size() ? this->digits[i] : 0) + (i < other.digits.size() ? other.digits[i] : 0);

		if (total > 9) {
			carry = 1;
			res.push_back(total - 10);
		}
		else {
			carry = 0;
			res.push_back(total);
		}
	}

	large_num result(res);
	result.rm_zeros();
	return result;
}

bool large_num::operator < (const large_num& other) const noexcept {
	if (*this == other) return false;

	if (!this->is_positive) {
		if (!other.is_positive) return ((-other) < (-(*this)));
		else return true;
	}

	else  if (!other.is_positive) { return false; }
	else {
		if (this->digits.size() != other.digits.size()) 
		{
			return this->digits.size() < other.digits.size();
		}
		else 
		{
			for (unsigned i = this->digits.size() - 1; i >= 0; --i) {
				if (this->digits[i] != other.digits[i]) {
					return this->digits[i] < other.digits[i];
				}
			}
			return false;
		}
	}
}

bool large_num::operator <= (const large_num& other) const noexcept {
	return (*this < other || *this == other);
}

bool large_num::operator >= (const large_num& other) const noexcept {
	return !(*this < other);
}

bool large_num::operator > (const large_num& other) const noexcept {
	return !(*this <= other);
}

large_num large_num::operator - (const large_num& other) const noexcept {
	if (!other.is_positive) return (*this + (-other));
	else if (!this->is_positive) return -((-(*this)) + other);
	else if (*this < other) return -(other - (*this));

	vector<int> res;
	int carry = 0;

	for (unsigned i = 0; i < max(other.digits.size(), this->digits.size()) || carry != 0; i++) {
		int total = this->digits[i] - (carry + (i < other.digits.size() ? other.digits[i] : 0));
		if (total < 0) {
			carry = 1;
			res.push_back(total + 10);
		}
		else {
			carry = 0;
			res.push_back(total);
		}
	}

	large_num r(res);
	r.rm_zeros();
	return r;
}

large_num large_num::operator * (const large_num& other) const noexcept {
	large_num result;
	result.digits.resize(this->digits.size() + other.digits.size() + 1, 0);

	for (unsigned i = 0; i < this->digits.size(); i++) {
		int carry = 0;

		for (unsigned j = 0; j < other.digits.size() || carry != 0; ++j) {
			int total = result.digits[i + j] + this->digits[i] * (j < other.digits.size() ? other.digits[j] : 0) + carry;
			result.digits[i + j] = static_cast<int>(total % 10);
			carry = static_cast<int>(total / 10);
		}
	}

	result.is_positive = this->is_positive == other.is_positive;
	result.rm_zeros();
	return result;
}

//shifts all digits by 1 to the right (multiplies by 10)
void large_num::shift_right() {
	if (this->digits.size() == 0) {
		this->digits.push_back(0);
		return;
	}

	this->digits.push_back(this->digits[this->digits.size() - 1]);

	for (unsigned i = this->digits.size() - 2; i > 0; --i) {
		this->digits[i] = this->digits[i - 1];
	}

	this->digits[0] = 0;
}

void large_num::shift_left() {
	if (this->digits.size() == 0) {
		this->digits.push_back(0);
		return;
	}

	this->digits.erase(digits.begin());
	// this->digits.push_back(0);
}


large_num large_num::operator / (const large_num& other) const {
	if (other == large_num("0")) throw large_num::DivideByZeroException();

	large_num b = other;
	b.is_positive = true;

	large_num result, current;
	result.digits.resize(this->digits.size());

	for (long long i = static_cast<long long>(this->digits.size()) - 1; i >= 0; --i) {
		current.shift_right();
		current.digits[0] = this->digits[i];
		current.rm_zeros();
		int x(0), l(0), r = 10;

		while (l <= r) {
			int m = (l + r) / 2;
			large_num mm(to_string(m));
			large_num t = b * mm;
			if (t <= current) {
				x = m;
				l = m + 1;
			}
			else r = m - 1;
		}
		result.digits[i] = x;
		large_num z(to_string(x));
		current = current - b * z;
	}

	result.is_positive = this->is_positive == other.is_positive;
	result.rm_zeros();
	return result;
}

bool large_num::odd() const noexcept {
	if (this->digits.size() == 0) return false;
	return this->digits[0] & 1;
}

bool large_num::even() const noexcept {
	return !this->odd();
}

large_num large_num::pow(large_num n) const {
	large_num a(*this), result(1);
	while (n != 0) {
		if (n.odd()) result = result * a;
		a = a * a;
		n = n / large_num(2);
	}
	return result;
}

large_num large_num::operator % (const large_num& mod) const noexcept {
	if (*this < mod) { return *this; }
	if (*this == mod) {
		large_num null("0");
		return null;
	}

	large_num res;
	large_num q = *this / mod;
	res = *this - q * mod;
	return res;
}


large_num large_num::gcd(const large_num& other) {
	large_num a = *this, b = other;
	while (a != b) {
		if (a > b) {
			large_num t = a;
			a = b;
			b = t;
		}
		b = b - a;
	}
	return a;
}

large_num large_num::random() {
	large_num one("1");
	large_num subb = *this - one;
	vector<int> p = subb.digits;

	for (int i = 0; i < p.size(); i++) {
		if (p[i] == 0) continue;
		int a = 1 + rand() % p[i];
		p[i] = a;
	}

	large_num res(p);
	res.rm_zeros();
	return res;
}


void large_num::extend_vector(vector<int>& a, int n) {
	while (n & (n - 1)) { ++n; }
	a.resize(n);
}


large_num large_num::sqrt_bi() {
	large_num one("1");
	large_num two("2");
	large_num l("0");
	large_num r = (*this) / two;
	while (l <= r) {
		large_num m = (r + l) / two;
		large_num t = m * m;
		if (t <= (*this)) {
			l = m + one;
		}
		else r = m - one;
	}
	return l - one;
}

large_num large_num::log2_ceil() const {
	int log_2 = (static_cast<int> (this->to_binary().digits.size())) - 1;
	if (log_2 < 0) { throw DivideByZeroException(); } // log of 0 is +oo
	cout << log_2 << endl;
	int log_e = ceil((double)log_2 * 0.69314718056);
	large_num result(log_e);

	return result;
}

large_num large_num::to_binary() const noexcept {
	large_num zero("0");
	vector<int> binary_representation;
	large_num convert_num = *this;

	while (convert_num > zero) {
		binary_representation.push_back(convert_num.digits[0] % 2);
		convert_num = convert_num / 2;
	}

	if (!this->is_positive) {
		binary_representation.at(0) = 0;
	}
	return large_num(binary_representation);
}
large_num large_num::inverse_number() const noexcept {
	if (!this->is_positive) return (-*this).inverse_number();
	large_num zero("0");
	large_num one("1");
	large_num two("2");
	large_num four("4");

	if (*this == zero) { return zero; }

	large_num b((-1) + (int)this->to_binary().digits.size());
	large_num r = two.pow(b);
	large_num s = r;
iteration_lable:
	r = two * r - (*this) * ((r*r) / two.pow(b)) / two.pow(b);
	if (r <= s) goto sharp_result;
	s = r;
	goto iteration_lable;
sharp_result:
	large_num y(four.pow(b) - (*this)*r);
	while (y < zero) {
		r = r - one;
		y = y + *this;
	}
	return r;
}

large_num large_num::schonhage_strassen(const large_num &other) noexcept {
	vector<int> a = this->digits, b = other.digits;
	size_t n = max(a.size(), b.size());
	a.resize(n), b.resize(n);
	fft(a, false), fft(b, false);

	for (size_t i = 0; i < n; ++i) {
		a[i] += b[i];
	}

	fft(a, true);

	large_num res(a);
	return res;
}

void large_num::fft(std::vector<int>& a, bool invert) {
	int n = (int)a.size();
	if (n == 1)  return;

	vector<int> a0(n / 2), a1(n / 2);

	for (int i = 0, j = 0; i < n; i += 2, ++j) {
		a0[j] = a[i];
		a1[j] = a[i + 1];
	}

	fft(a0, invert);
	fft(a1, invert);

	double ang = 2 * PI / n * (invert ? -1 : 1);

	int w(1);
	for (int i = 0; i < n / 2; ++i) {
		a[i] = a0[i] + w * a1[i];
		a[i + n / 2] = a0[i] - w * a1[i];
		if (invert)
			a[i] /= 2, a[i + n / 2] /= 2;
		w *= w;
	}
}


vector<int> large_num::mul(const vector<int>& a, const vector<int>& b) {
	vector<int> res(a.size() + b.size() + 1, 0);

	for (unsigned i = 0; i < a.size(); i++) {
		int carry = 0;
		for (unsigned j = 0; j < b.size() || carry != 0; ++j) {
			int cur = res[i + j] + a[i] * (j < b.size() ? b[j] : 0) + carry;
			res[i + j] = static_cast<int>(cur % 10);
			carry = static_cast<int>(cur / 10);
		}
	}
	vect_rm_zeros(res);
	return res;
}


void large_num::vect_rm_zeros(vector<int>& a) {
	while (a.size() > 1 && a.back() == 0) {
		a.pop_back();
	}
}

vector<int> large_num::to_vec(const string& a) {
	vector<int> res;

	for (int i = a.size() - 1; i >= 0; i--) {
		res.push_back(a[i] - 48);
	}
	return res;
}

vector<int> large_num::pow(const vector<int>& a, long long n) {
	vector<int> b = a, res = { 1 };

	while (n != 0) {
		if (n % 2 == 1) res = mul(res, b);
		b = mul(b, b);
		n /= 2;
	}
	return res;
}

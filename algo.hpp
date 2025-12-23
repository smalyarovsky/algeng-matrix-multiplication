#pragma once
#include <bits/stdc++.h>
using namespace std;

struct Algo {
    virtual ~Algo() = default;
    virtual string name() const = 0;
    virtual bool multiply(const double* A, const double* B, double* C, int n) = 0;
    virtual pair<unsigned long long, unsigned long long> ops(int n) = 0;
};

unique_ptr<Algo> make_naive();
unique_ptr<Algo> make_strassen();
unique_ptr<Algo> make_alphaevolve();


#include "algo.hpp"
using namespace std;

static unsigned long long G_MULS = 0;
static unsigned long long G_ADDS = 0;

template<bool Count>
static inline double op_add(double a, double b) {
    if constexpr (Count) ++G_ADDS;
    return a + b;
}

template<bool Count>
static inline double op_mul(double a, double b) {
    if constexpr (Count) ++G_MULS;
    return a * b;
}

template<bool Count>
static void naive_mul(const double* A, const double* B, double* C, int n) {
    for (int i = 0; i < n; ++i) {
        double* crow = C + (size_t)i * n;
        for (int k = 0; k < n; ++k) {
            double aik = A[(size_t)i * n + k];
            const double* brow = B + (size_t)k * n;
            for (int j = 0; j < n; ++j) {
                crow[j] = op_add<Count>(crow[j], op_mul<Count>(aik, brow[j]));
            }
        }
    }
}

struct Naive : Algo {
    string name() const override { return "naive"; }

    bool multiply(const double* A, const double* B, double* C, int n) override {
        naive_mul<false>(A, B, C, n);
        return true;
    }

    pair<unsigned long long, unsigned long long> ops(int n) override {
        G_MULS = 0;
        G_ADDS = 0;

        vector<double> A((size_t)n * n, 1.0);
        vector<double> B((size_t)n * n, 1.0);
        vector<double> C((size_t)n * n, 0.0);

        naive_mul<true>(A.data(), B.data(), C.data(), n);

        return {G_MULS, G_ADDS};
    }
};

unique_ptr<Algo> make_naive() {
    return make_unique<Naive>();
}

#include "algo.hpp"
using namespace std;

static unsigned long long G_MULS = 0, G_ADDS = 0;

template<bool Count> static inline double op_add(double a, double b){ if constexpr(Count) ++G_ADDS; return a+b; }
template<bool Count> static inline double op_sub(double a, double b){ if constexpr(Count) ++G_ADDS; return a-b; }
template<bool Count> static inline double op_mul(double a, double b){ if constexpr(Count) ++G_MULS; return a*b; }

static int next_pow2(int n){ int p=1; while(p<n) p<<=1; return p; }

template<bool Count>
static void add_mat(const double* A, const double* B, double* C, int n, int sA, int sB, int sC){
    for(int i=0;i<n;++i){
        const double* a=A+(size_t)i*sA;
        const double* b=B+(size_t)i*sB;
        double* c=C+(size_t)i*sC;
        for(int j=0;j<n;++j) c[j]=op_add<Count>(a[j],b[j]);
    }
}

template<bool Count>
static void sub_mat(const double* A, const double* B, double* C, int n, int sA, int sB, int sC){
    for(int i=0;i<n;++i){
        const double* a=A+(size_t)i*sA;
        const double* b=B+(size_t)i*sB;
        double* c=C+(size_t)i*sC;
        for(int j=0;j<n;++j) c[j]=op_sub<Count>(a[j],b[j]);
    }
}

template<bool Count>
static void naive_mul_add(const double* A, const double* B, double* C, int n, int sA, int sB, int sC){
    for(int i=0;i<n;++i){
        double* crow=C+(size_t)i*sC;
        for(int k=0;k<n;++k){
            double aik=A[(size_t)i*sA+(size_t)k];
            const double* brow=B+(size_t)k*sB;
            for(int j=0;j<n;++j){
                crow[j]=op_add<Count>(crow[j], op_mul<Count>(aik, brow[j]));
            }
        }
    }
}

static size_t ws_need_doubles(int n){
    if(n<=1) return 0;
    int m=n/2;
    size_t m2=(size_t)m*(size_t)m;
    return 9ull*m2 + ws_need_doubles(m);
}

template<bool Count>
static bool strassen_rec(const double* A, const double* B, double* C,
                         int n, int sA, int sB, int sC,
                         double* ws){
    if(n==1){
        C[0]=op_add<Count>(C[0], op_mul<Count>(A[0], B[0]));
        return true;
    }
    if(n&1) return false;

    int m=n/2;
    size_t m2=(size_t)m*(size_t)m;

    const double* A11=A;
    const double* A12=A+m;
    const double* A21=A+(size_t)m*sA;
    const double* A22=A+(size_t)m*sA+m;

    const double* B11=B;
    const double* B12=B+m;
    const double* B21=B+(size_t)m*sB;
    const double* B22=B+(size_t)m*sB+m;

    double* C11=C;
    double* C12=C+m;
    double* C21=C+(size_t)m*sC;
    double* C22=C+(size_t)m*sC+m;

    double* M1=ws; ws+=m2;
    double* M2=ws; ws+=m2;
    double* M3=ws; ws+=m2;
    double* M4=ws; ws+=m2;
    double* M5=ws; ws+=m2;
    double* M6=ws; ws+=m2;
    double* M7=ws; ws+=m2;
    double* S =ws; ws+=m2;
    double* T =ws; ws+=m2;

    memset(M1,0,m2*sizeof(double));
    memset(M2,0,m2*sizeof(double));
    memset(M3,0,m2*sizeof(double));
    memset(M4,0,m2*sizeof(double));
    memset(M5,0,m2*sizeof(double));
    memset(M6,0,m2*sizeof(double));
    memset(M7,0,m2*sizeof(double));

    add_mat<Count>(A11,A22,S,m,sA,sA,m);
    add_mat<Count>(B11,B22,T,m,sB,sB,m);
    if(!strassen_rec<Count>(S,T,M1,m,m,m,m,ws)) return false;

    add_mat<Count>(A21,A22,S,m,sA,sA,m);
    if(!strassen_rec<Count>(S,B11,M2,m,m,sB,m,ws)) return false;

    sub_mat<Count>(B12,B22,T,m,sB,sB,m);
    if(!strassen_rec<Count>(A11,T,M3,m,sA,m,m,ws)) return false;

    sub_mat<Count>(B21,B11,T,m,sB,sB,m);
    if(!strassen_rec<Count>(A22,T,M4,m,sA,m,m,ws)) return false;

    add_mat<Count>(A11,A12,S,m,sA,sA,m);
    if(!strassen_rec<Count>(S,B22,M5,m,m,sB,m,ws)) return false;

    sub_mat<Count>(A21,A11,S,m,sA,sA,m);
    add_mat<Count>(B11,B12,T,m,sB,sB,m);
    if(!strassen_rec<Count>(S,T,M6,m,m,m,m,ws)) return false;

    sub_mat<Count>(A12,A22,S,m,sA,sA,m);
    add_mat<Count>(B21,B22,T,m,sB,sB,m);
    if(!strassen_rec<Count>(S,T,M7,m,m,m,m,ws)) return false;

    for(int i=0;i<m;++i){
        double* c11=C11+(size_t)i*sC;
        double* c12=C12+(size_t)i*sC;
        double* c21=C21+(size_t)i*sC;
        double* c22=C22+(size_t)i*sC;

        const double* p1=M1+(size_t)i*m;
        const double* p2=M2+(size_t)i*m;
        const double* p3=M3+(size_t)i*m;
        const double* p4=M4+(size_t)i*m;
        const double* p5=M5+(size_t)i*m;
        const double* p6=M6+(size_t)i*m;
        const double* p7=M7+(size_t)i*m;

        for(int j=0;j<m;++j){
            c11[j]=op_add<Count>(op_sub<Count>(op_add<Count>(p1[j],p4[j]),p5[j]),p7[j]);
            c12[j]=op_add<Count>(p3[j],p5[j]);
            c21[j]=op_add<Count>(p2[j],p4[j]);
            c22[j]=op_add<Count>(op_add<Count>(op_sub<Count>(p1[j],p2[j]),p3[j]),p6[j]);
        }
    }

    return true;
}

struct Strassen : Algo {
    string name() const override { return "strassen"; }

    bool multiply(const double* A, const double* B, double* C, int n) override {
        int p=next_pow2(n);
        size_t pp=(size_t)p*(size_t)p;

        vector<double> Ap(pp,0.0), Bp(pp,0.0), Cp(pp,0.0);
        for(int i=0;i<n;++i){
            memcpy(Ap.data()+(size_t)i*p, A+(size_t)i*n, (size_t)n*sizeof(double));
            memcpy(Bp.data()+(size_t)i*p, B+(size_t)i*n, (size_t)n*sizeof(double));
        }

        size_t need=ws_need_doubles(p);
        vector<double> ws(need);
        bool ok=strassen_rec<false>(Ap.data(),Bp.data(),Cp.data(),p,p,p,p,ws.data());
        if(!ok) return false;

        for(int i=0;i<n;++i) memcpy(C+(size_t)i*n, Cp.data()+(size_t)i*p, (size_t)n*sizeof(double));
        return true;
    }

    pair<unsigned long long, unsigned long long> ops(int n) override {
        int p=next_pow2(n);
        size_t pp=(size_t)p*(size_t)p;

        vector<double> Ap(pp), Bp(pp), Cp(pp,0.0);
        mt19937 rng(777);
        uniform_real_distribution<double> dist(-0.5,0.5);
        for(size_t i=0;i<pp;++i){ Ap[i]=dist(rng); Bp[i]=dist(rng); }

        G_MULS=0; G_ADDS=0;

        size_t need=ws_need_doubles(p);
        vector<double> ws(need);
        bool ok=strassen_rec<true>(Ap.data(),Bp.data(),Cp.data(),p,p,p,p,ws.data());
        if(!ok) return {0,0};

        return {G_MULS,G_ADDS};
    }
};

unique_ptr<Algo> make_strassen(){ return make_unique<Strassen>(); }

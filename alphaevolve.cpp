#include "algo.hpp"
#include <cmath>
using namespace std;

static unsigned long long G_MULS = 0, G_ADDS = 0;

template<bool Count> static inline double add(double a, double b){ if constexpr(Count) ++G_ADDS; return a + b; }
template<bool Count> static inline double sub(double a, double b){ if constexpr(Count) ++G_ADDS; return a - b; }
template<bool Count> static inline double neg(double a){ if constexpr(Count) ++G_ADDS; return -a; }

template<bool Count> static inline double mul_rank(double a, double b){ if constexpr(Count) ++G_MULS; return a * b; }

static inline double div2(double x){ return ldexp(x, -1); }
static inline double div4(double x){ return ldexp(x, -2); }
static inline double div8(double x){ return ldexp(x, -3); }

template<bool Count>
static inline void alpha48_mul4x4_add(const double* A, const double* B, double* C)
{
    const double a11 = A[0];
    const double a12 = A[1];
    const double a13 = A[2];
    const double a14 = A[3];
    const double a21 = A[4];
    const double a22 = A[5];
    const double a23 = A[6];
    const double a24 = A[7];
    const double a31 = A[8];
    const double a32 = A[9];
    const double a33 = A[10];
    const double a34 = A[11];
    const double a41 = A[12];
    const double a42 = A[13];
    const double a43 = A[14];
    const double a44 = A[15];

    const double b11 = B[0];
    const double b12 = B[1];
    const double b13 = B[2];
    const double b14 = B[3];
    const double b21 = B[4];
    const double b22 = B[5];
    const double b23 = B[6];
    const double b24 = B[7];
    const double b31 = B[8];
    const double b32 = B[9];
    const double b33 = B[10];
    const double b34 = B[11];
    const double b41 = B[12];
    const double b42 = B[13];
    const double b43 = B[14];
    const double b44 = B[15];

    double s[48];
    double t[48];
    double p[48];

    {
        const double x16 = add<Count>(a13,a24);
        const double x17 = sub<Count>(a11,a22);
        const double x18 = add<Count>(a34,a43);
        const double x19 = sub<Count>(a31,a42);
        const double x20 = sub<Count>(a32,a41);
        const double x21 = add<Count>(a14,a23);
        const double x22 = add<Count>(a33,a44);
        const double x23 = sub<Count>(a12,a21);
        const double x24 = sub<Count>(neg<Count>(a12),a22);
        const double x25 = sub<Count>(a13,a23);
        const double x26 = sub<Count>(a32,a42);
        const double x27 = sub<Count>(neg<Count>(a34),a44);
        const double x28 = add<Count>(a33,a43);
        const double x29 = sub<Count>(a14,a24);
        const double x30 = sub<Count>(a31,a41);
        const double x31 = add<Count>(a11,a21);
        const double x32 = sub<Count>(neg<Count>(x16),x21);
        const double x33 = add<Count>(x17,x23);
        const double x34 = sub<Count>(x18,x22);
        const double x35 = sub<Count>(x19,x20);
        const double x36 = sub<Count>(x17,x23);
        const double x37 = sub<Count>(x16,x21);
        const double x38 = add<Count>(x19,x20);
        const double x39 = add<Count>(x18,x22);
        const double x40 = sub<Count>(x26,x30);
        const double x41 = sub<Count>(x31,x24);
        const double x42 = add<Count>(x27,x28);
        const double x43 = add<Count>(x25,x29);
        const double x44 = sub<Count>(neg<Count>(x40),x32);

        const double s24 = add<Count>(x33,x42);
        const double x46 = add<Count>(a32,a42);
        const double s8  = add<Count>(x35,x43);
        const double s34 = sub<Count>(x41,x34);
        const double x49 = sub<Count>(a34,a44);
        const double x50 = add<Count>(a31,a41);
        const double x51 = sub<Count>(a33,a43);
        const double x52 = add<Count>(a14,a24);
        const double x53 = add<Count>(a13,a23);
        const double x54 = sub<Count>(a12,a22);
        const double x55 = sub<Count>(a11,a21);
        const double s9  = add<Count>(x30,x31);
        const double s12 = add<Count>(x25,x28);
        const double x58 = add<Count>(x19,x22);
        const double x59 = sub<Count>(x35,x32);
        const double x62 = sub<Count>(neg<Count>(x38),x37);
        const double x63 = add<Count>(add<Count>(x51,x36),x49);
        const double x64 = add<Count>(x16,x17);
        const double s42 = sub<Count>(x30,x31);
        const double x66 = add<Count>(x21,x23);
        const double s36 = sub<Count>(x24,x26);
        const double s22 = sub<Count>(x27,x29);
        const double x71 = add<Count>(x18,x20);
        const double s38 = add<Count>(x26,x24);
        const double x73 = sub<Count>(x38,x37);
        const double x74 = add<Count>(add<Count>(x50,x37),x46);
        const double x75 = sub<Count>(x39,x36);
        const double x76 = sub<Count>(x33,x34);
        const double x77 = add<Count>(x32,x35);
        const double x78 = add<Count>(x33,x34);
        const double x79 = sub<Count>(add<Count>(x53,x38),x52);
        const double x80 = add<Count>(x36,x39);
        const double x81 = sub<Count>(add<Count>(x39,x55),x54);

        const double s2  = sub<Count>(x25,x28);
        const double s6  = add<Count>(x29,x27);
        const double s0  = sub<Count>(x44,x63);
        const double s1  = sub<Count>(x30,x55);
        const double s3  = add<Count>(x32,x42);
        const double s4  = sub<Count>(neg<Count>(x74),s24);
        const double s5  = add<Count>(x76,x62);
        const double s7  = sub<Count>(x62,x76);
        const double s10 = sub<Count>(x73,x78);
        const double s11 = add<Count>(x66,x71);
        const double s13 = sub<Count>(x53,x28);
        const double s14 = sub<Count>(s42,s2);
        const double s15 = sub<Count>(neg<Count>(x73),x78);
        const double s16 = sub<Count>(x41,x35);
        const double s17 = sub<Count>(x59,x80);
        const double s18 = sub<Count>(x66,x71);
        const double s19 = add<Count>(x75,x77);
        const double s20 = sub<Count>(x25,x51);
        const double s21 = add<Count>(x58,x64);
        const double s23 = add<Count>(s8,x81);
        const double s25 = add<Count>(x63,x44);
        const double s26 = sub<Count>(x46,x24);
        const double s27 = neg<Count>(x44);
        const double s28 = sub<Count>(x81,s8);
        const double s29 = sub<Count>(x52,x27);
        const double s30 = add<Count>(x26,x54);
        const double s31 = sub<Count>(x58,x64);
        const double s32 = sub<Count>(x74,s24);
        const double s33 = add<Count>(x34,x43);
        const double s35 = add<Count>(x59,x80);
        const double s37 = sub<Count>(neg<Count>(s22),s36);
        const double s39 = add<Count>(x33,x40);
        const double s40 = sub<Count>(x77,x75);
        const double s41 = add<Count>(x29,x49);
        const double s43 = add<Count>(s9,s12);
        const double s44 = sub<Count>(x79,s34);
        const double s45 = add<Count>(x79,s34);
        const double s46 = sub<Count>(x50,x31);
        const double s47 = sub<Count>(s6,s38);

        s[0]=s0; s[1]=s1; s[2]=s2; s[3]=s3; s[4]=s4; s[5]=s5; s[6]=s6; s[7]=s7;
        s[8]=s8; s[9]=s9; s[10]=s10; s[11]=s11; s[12]=s12; s[13]=s13; s[14]=s14; s[15]=s15;
        s[16]=s16; s[17]=s17; s[18]=s18; s[19]=s19; s[20]=s20; s[21]=s21; s[22]=s22; s[23]=s23;
        s[24]=s24; s[25]=s25; s[26]=s26; s[27]=s27; s[28]=s28; s[29]=s29; s[30]=s30; s[31]=s31;
        s[32]=s32; s[33]=s33; s[34]=s34; s[35]=s35; s[36]=s36; s[37]=s37; s[38]=s38; s[39]=s39;
        s[40]=s40; s[41]=s41; s[42]=s42; s[43]=s43; s[44]=s44; s[45]=s45; s[46]=s46; s[47]=s47;
    }

    {
        const double y16 = sub<Count>(b11,b14);
        const double y17 = sub<Count>(b21,b23);
        const double y18 = sub<Count>(b31,b34);
        const double y19 = sub<Count>(b41,b43);
        const double y20 = sub<Count>(neg<Count>(b12),b13);
        const double t17 = add<Count>(y17,y19);
        const double y22 = add<Count>(b32,b33);
        const double t2  = add<Count>(b12,b32);
        const double t6  = sub<Count>(b42,b22);
        const double t7  = add<Count>(y18,y16);
        const double t14 = sub<Count>(b12,b31);
        const double t15 = sub<Count>(y16,y18);
        const double t16 = add<Count>(y17,y20);
        const double t19 = sub<Count>(y19,y17);
        const double t22 = sub<Count>(sub<Count>(t17,b24),b44);
        const double t33 = sub<Count>(sub<Count>(y18,b42),b44);
        const double t35 = sub<Count>(y20,y22);
        const double t38 = sub<Count>(b21,b41);
        const double t39 = sub<Count>(sub<Count>(neg<Count>(y16),b22),b24);
        const double t42 = add<Count>(b11,b31);
        const double t43 = sub<Count>(sub<Count>(sub<Count>(b13,b34),y16),y22);
        const double t47 = add<Count>(b22,b41);
        const double v48 = sub<Count>(t22,t6);
        const double v49 = sub<Count>(t19,t35);
        const double v50 = sub<Count>(t33,t39);
        const double v51 = add<Count>(t16, div2(v49));
        const double v52 = sub<Count>(div2(v48), t47);
        const double v53 = add<Count>(t14,t43);
        const double v54 = sub<Count>(neg<Count>(t33),t39);
        const double v55 = sub<Count>(t14,t43);
        const double g56 = div2(v54);
        const double v56 = add<Count>(v51,g56);
        const double g57 = div2(v53);
        const double v57 = add<Count>(g57,v52);
        const double v58 = sub<Count>(v51,g56);
        const double v59 = sub<Count>(v52,g57);
        const double t24 = add<Count>(v48,v55);
        const double t3  = add<Count>(t16,v49);
        const double t8  = sub<Count>(v48,v55);
        const double t21 = add<Count>(v50,v49);
        const double v65 = sub<Count>(t7, div2(t21));
        const double v66 = sub<Count>(t2,t14);
        const double v67 = add<Count>(t15,v58);
        const double v68 = add<Count>(t47,t6);
        const double v69 = sub<Count>(neg<Count>(t38),t47);
        const double v71 = add<Count>(t42,t14);
        const double v72 = sub<Count>(v56,t17);
        const double g73 = div2(t24);
        const double g74 = div2(t8);
        const double g75 = div2(add<Count>(add<Count>(t35,v50),t19));
        const double t13 = sub<Count>(v66,t3);
        const double t11 = v58*2;
        const double t30 = sub<Count>(v69,t39);
        const double t31 = sub<Count>(v49,v50);
        const double t1  = add<Count>(t39,v71);
        const double t18 = v56*2;
        const double t25 = add<Count>(g75,v57);
        const double t20 = sub<Count>(v66,t33);
        const double t4  = add<Count>(v65,g73);
        const double t46 = add<Count>(t16,v71);
        const double t29 = sub<Count>(v68,t3);
        const double t5  = sub<Count>(v54,t15);
        const double t10 = sub<Count>(v50,t7);
        const double t12 = sub<Count>(v55,t2);
        const double t34 = v59*2;
        const double t28 = add<Count>(v72,g74);
        const double t23 = sub<Count>(g75,g74);
        const double t9  = add<Count>(t42,v53);
        const double t32 = sub<Count>(v67,g73);
        const double t37 = sub<Count>(v48,t47);
        const double t27 = v57*2;
        const double t45 = sub<Count>(v65,v59);
        const double t44 = sub<Count>(v67,v59);
        const double t26 = add<Count>(t16,v69);
        const double t36 = sub<Count>(v52*2,t38);
        const double t0  = add<Count>(v57,v72);
        const double t41 = add<Count>(t33,v68);
        const double t40 = sub<Count>(t17, v51*2);

        t[0]=t0; t[1]=t1; t[2]=t2; t[3]=t3; t[4]=t4; t[5]=t5; t[6]=t6; t[7]=t7;
        t[8]=t8; t[9]=t9; t[10]=t10; t[11]=t11; t[12]=t12; t[13]=t13; t[14]=t14; t[15]=t15;
        t[16]=t16; t[17]=t17; t[18]=t18; t[19]=t19; t[20]=t20; t[21]=t21; t[22]=t22; t[23]=t23;
        t[24]=t24; t[25]=t25; t[26]=t26; t[27]=t27; t[28]=t28; t[29]=t29; t[30]=t30; t[31]=t31;
        t[32]=t32; t[33]=t33; t[34]=t34; t[35]=t35; t[36]=t36; t[37]=t37; t[38]=t38; t[39]=t39;
        t[40]=t40; t[41]=t41; t[42]=t42; t[43]=t43; t[44]=t44; t[45]=t45; t[46]=t46; t[47]=t47;
    }

    for(int i=0;i<48;++i) p[i] = mul_rank<Count>(s[i], t[i]);

    {
        const double p0=p[0],p1=p[1],p2=p[2],p3=p[3],p4=p[4],p5=p[5],p6=p[6],p7=p[7],p8=p[8],p9=p[9],p10=p[10],p11=p[11],p12=p[12],p13=p[13],p14=p[14],p15=p[15],p16=p[16],p17=p[17],p18=p[18],p19=p[19],p20=p[20],p21=p[21],p22=p[22],p23=p[23],p24=p[24],p25=p[25],p26=p[26],p27=p[27],p28=p[28],p29=p[29],p30=p[30],p31=p[31],p32=p[32],p33=p[33],p34=p[34],p35=p[35],p36=p[36],p37=p[37],p38=p[38],p39=p[39],p40=p[40],p41=p[41],p42=p[42],p43=p[43],p44=p[44],p45=p[45],p46=p[46],p47=p[47];

        const double e5  = add<Count>(sub<Count>(p20,p41),p33);
        const double v61 = sub<Count>(p34,p44);
        const double e7  = add<Count>(sub<Count>(p9*2,p40),p11);
        const double e8  = sub<Count>(sub<Count>(p26,p46),p16);
        const double e12 = add<Count>(add<Count>(p30,p39),p1);
        const double e15 = sub<Count>(p24,p32);
        const double v62 = add<Count>(p0,p27);
        const double v71 = (p46-p13)*2;
        const double v70 = div4(p0+p44);
        const double v69 = div4(p32-p28);
        const double v68 = (p41-p30)*2;
        const double v65 = (p20-p1)*2;
        const double v63 = (p29-p26)*2;
        const double e28 = sub<Count>(v61,p18);
        const double e29 = sub<Count>(sub<Count>(v61,p15*2),e15);
        const double v54 = add<Count>(e8,e12);
        const double v51 = v69 + div2(e12-p47-p14-e8);
        const double e35 = p36*2+p5-p15+v62+e28+e7;
        const double e38 = (p43-p9+p12)*2+v54;
        const double e39 = p43+p37+v54;
        const double e40 = p36-p9+e39;
        const double e43 = div2(p14-p47-e38)-v70-e5;
        const double z21 = p21+p45+p11-p4+v68+v63+e29;
        const double z8  = p8+p28+v62+e29+(e7-e39)*2;
        const double z6  = p6+v69-e43;
        const double z3  = p3+e5-p29-p13+e38;
        const double z2  = p2+v69+e43;
        const double w43 = div2(p42+v70-v51);
        const double w42 = div2(p38+v70+v51);
        const double w41 = div8(p31+p45+p4+p5*2+v71+e15+v65+e28);
        const double w40 = div4(p25+p45-(p9+p36)*2);
        const double w39 = div4(p23+e40*2-p4);
        const double w38 = div8(z21);
        const double w37 = div8(z8);
        const double w36 = div2(z6);
        const double w35 = div4(z3);
        const double w34 = -div2(z2);
        const double e51 = div4(v68-p7-p15)-w43-w38;
        const double e52 = div4(p17+e35)-w37;
        const double c31 = w42+w37;
        const double c11 = w37-w42;
        const double c32 = -w37-w36;
        const double e54 = w35-w41;
        const double w24 = w38+w35;
        const double w29 = div8(z8-z21);
        const double c12 = div2(z6-z3)-w37;
        const double w23 = -e52;
        const double c21 = e52-w43;
        const double e55 = div4(p40+v63-p35+z3)+w29-w34;
        const double e56 = e54+div4(p5-p10+v65)-w36;
        const double e57 = div4(p19+v71+e35)-w41-c31;
        const double e60 = w23+div2(p22+p12+e40);
        const double c41 = w23-w43;
        const double w20 = e60;
        const double w17 = e54-e60;
        const double c42 = w34-w20;
        const double c22 = w20-div2(z3+z2);
        const double c34 = w29-e56;
        const double c14 = w29+e56;
        const double c33 = w24-w40-e57;
        const double c13 = w24-w39+e57;
        const double c44 = e51-w17;
        const double c24 = e51+w17;
        const double c43 = w40-w41+e55;
        const double c23 = w41+w39+e55;

        C[0]  = add<Count>(C[0],  c11);
        C[1]  = add<Count>(C[1],  c12);
        C[2]  = add<Count>(C[2],  c13);
        C[3]  = add<Count>(C[3],  c14);
        C[4]  = add<Count>(C[4],  c21);
        C[5]  = add<Count>(C[5],  c22);
        C[6]  = add<Count>(C[6],  c23);
        C[7]  = add<Count>(C[7],  c24);
        C[8]  = add<Count>(C[8],  c31);
        C[9]  = add<Count>(C[9],  c32);
        C[10] = add<Count>(C[10], c33);
        C[11] = add<Count>(C[11], c34);
        C[12] = add<Count>(C[12], c41);
        C[13] = add<Count>(C[13], c42);
        C[14] = add<Count>(C[14], c43);
        C[15] = add<Count>(C[15], c44);
    }
}

struct AlphaEvolve : Algo {
    string name() const override { return "alphaevolve"; }

    bool multiply(const double* A, const double* B, double* C, int n) override {
        if (n != 4) return false;
        alpha48_mul4x4_add<false>(A, B, C);
        return true;
    }

    pair<unsigned long long, unsigned long long> ops(int n) override {
        if (n != 4) return {0,0};
        G_MULS = 0; G_ADDS = 0;
        double A[16], B[16], C[16];
        mt19937 rng(777);
        uniform_real_distribution<double> dist(-0.5,0.5);
        for(int i=0;i<16;++i){ A[i]=dist(rng); B[i]=dist(rng); C[i]=0.0; }
        alpha48_mul4x4_add<true>(A,B,C);
        return {G_MULS, G_ADDS};
    }
};

unique_ptr<Algo> make_alphaevolve(){ return make_unique<AlphaEvolve>(); }

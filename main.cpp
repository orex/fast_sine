#include <iostream>
#include <iomanip>
#include <cmath>
#include <array>
#include <bitset>
#include <vector>
#include <boost/timer/timer.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/predef.h>
#include <random>
#include <fstream>

namespace bm = boost::math;
namespace bmp = boost::multiprecision;

typedef bmp::number<bmp::mpfr_float_backend<30> > bmp_prec_t;

using namespace std;

typedef union { uint32_t a[2]; uint64_t i; double x; } mynumber;

vector<double> fc;

// Very straightforward sin
double sin_e1(double x) {
  double result = 0;
  int sign = 1;
  for(int i = 1; i < 25; i += 2) {
    result += sign * pow(x, i) / bm::unchecked_factorial<double>(i);
    sign = -sign;
  }
  return result;
}

// Factorial and power optimization
double sin_e2(double x) {
  double result = 0;
  int sign = 1;
  double xx = x * x;
  double pw = x;
  double fti = 1.0;
  for(int i = 1; i < 25; i += 2) {
    fti /= i;
    result += sign * pw * fti;
    fti /= ( i + 1 );
    sign = -sign;
    pw  *= xx;
  }
  return result;
}

// Reversed calculation of sin_e1
double sin_e3(double x) {
  double result = 0;
  for(int i = 25; i >= 1; i -= 2) {
    result += (((i - 1) % 4 == 0) ? 1 : -1 ) * pow(x, i) / bm::unchecked_factorial<double>(i);
  }
  return result;
}

// Horner's method direct
double sin_e4(double x) {
  double xx = x * x;
  double res = fc[25];
  for(int i = 23; i >= 1; i -= 2) {
    res = fc[i] + xx * res;
  }
  return x * res;
}


#pragma GCC push_options
#pragma GCC optimize ("O2")
// Horner's method improved accuracy
double sin_e5(double x) {
  double xx = x * x;
  double res = fc[27];
  for(int i = 25; i >= 3; i -= 2) {
    res = fc[i] + xx * res;
  }
  return x + xx * (x * res);
}

#pragma GCC pop_options


#if BOOST_ARCH_X86_64
// Obsolete x87 FPU code.
inline
double fsin(double x) {
  double result;
  asm ("fsin" :"=t" (result) : "0" (x));
  return result;
}

#endif

// Fixed point calculations block

#pragma GCC push_options
#pragma GCC optimize ("Ofast")

#if BOOST_ARCH_X86_64
/*
inline uint64_t mul2(const uint64_t a, const uint64_t b) {
  uint64_t res;
  asm("mulq %2" : "=d"(res) : "a"(a), "rm"(b) : "cc");
  return res;
}*/

inline uint64_t mul2(const uint64_t a, const uint64_t b) {
  uint64_t res;
  asm("mulx %2, %0, %0" : "=r"(res) : "d"(a), "rm"(b));
  return res;
}

#elif BOOST_ARCH_ARM
inline uint64_t mul2(const uint64_t a, const uint64_t b) {
  uint64_t res;
  asm("UMULH %0, %1, %2" : "=r"(res) : "r"(a), "r"(b));
  return res;
}
#else

inline uint64_t mul2(const uint64_t a, const uint64_t b) {
  union {
    __uint128_t o;
    uint64_t u[2];
  } res;
  __builtin_mul_overflow(a, b, &res.o);
  return res.u[1];
}

#endif

constexpr array<uint64_t, 18> tsx = { // 2^64/i!
    0x0000000000000000LL,
    0x0000000000000000LL,
    0x8000000000000000LL,
    0x2aaaaaaaaaaaaaaaLL,
    0x0aaaaaaaaaaaaaaaLL,
    0x0222222222222222LL,
    0x005b05b05b05b05bLL,
    0x000d00d00d00d00dLL,
    0x0001a01a01a01a01LL,
    0x00002e3bc74aad8eLL,
    0x0000049f93edde27LL,
    0x0000006b99159fd5LL,
    0x00000008f76c77fcLL,
    0x00000000b092309dLL,
    0x000000000c9cba54LL,
    0x0000000000d73f9fLL,
    0x00000000000d73f9LL,
    0x000000000000ca96LL
};

constexpr mynumber toint    = {{0x00000000, 0x43F00000}};  /*  18446744073709551616 = 2^64     */
constexpr mynumber todouble = {{0x00000000, 0x3BF00000}};  /*  ~5.42101086242752217003726400434E-20 = 2^-64     */


// Fixed point calculations without factorial table (slow)
inline
double sin_e6(double xd) {
  uint64_t x = xd * toint.x;
  uint64_t xx = mul2(x, x);
  constexpr uint64_t half = uint64_t(1) << 63;
  uint64_t res =  half - xx / (2 * 17 * 16);
  for(int i = 15; i >= 3; i -= 2) {
    res = half - mul2(res, xx) / (i * (i - 1));
  }
  res = mul2(x, 2 * res);
  return res * todouble.x;
}

// Fixed points calculations with table (fast)
inline
double sin_e7(double xd) {
  uint64_t x = xd * toint.x;
  uint64_t xx = mul2(x, x);
  uint64_t res = tsx[19];
  for(int i = 17; i >= 3; i -= 2) {
    res = tsx[i] - mul2(res, xx);
  }
  res = mul2(res, x);
  res = x - mul2(xx, res);
  return res * todouble.x;
}
#pragma GCC pop_options

#define TEST_LOOP 1

#define SIN(a) sin_e7(a)
// ^^ Define function for the test here. ^^

#pragma GCC push_options
#pragma GCC optimize ("Ofast")
inline
void sin_ev(const vector<double> &x, vector<double> &y) {
  cout << "sin_e.." << endl;
  boost::timer::auto_cpu_timer at;
  for(int j = 0; j < TEST_LOOP; j++) {
    for (int i = 0; i < x.size(); i++) {
      y[i] = SIN(x[i]);
    }
  }
}
#pragma GCC pop_options

#pragma GCC push_options
#pragma GCC optimize ("O2")
void sin_iv(const vector<double> &x, vector<double> &y) {
  cout << "sin" << endl;
  boost::timer::auto_cpu_timer at;
  for(int j = 0; j < TEST_LOOP; j++) {
    for (int i = 0; i < x.size(); i++) {
      y[i] = sin(x[i]);
    }
  }
}
#pragma GCC pop_options

int main() {
  fc.resize(28);
  bmp_prec_t ft = 1;
  fc[1] = 1.0; //3 * 5;
  for(int i = 2; i < fc.size(); i++) {
    ft *= i;
    // factorial with sign for Taylor series
    fc[i] = double(1 / ft) * (( (i - 2) % 4 < 2) ? -1 : 1);
  }
  vector<double> xv, ye, yi;
  xv.resize(8 * 2000000);
  //xv.resize(50000);
  ye.resize(xv.size());
  yi.resize(xv.size());
  // Linear filling of input values
  for (int i = 0; i < xv.size(); i++) {
    xv[i] = 0.126 + (0.856  - 0.126) * i / double(xv.size());
  }
  //shuffle (xv.begin(), xv.end(), std::default_random_engine(200));
  //reverse (xv.begin(), xv.end());

  sin_ev(xv, ye);
  sin_iv(xv, yi);

  int co = 0, cn = 0;
  // Use mpfr library as "true" value
  bmp_prec_t avg = 0.0, div = 0.0;
  double co_max = 0, cn_max = 0;
  //fstream fs("out.txt", fstream::out);
  for(int i = 0; i < xv.size(); i++) {
    mynumber dqs, dxv, dold, dnew;
    dxv.x = xv[i];
    dold.x = yi[i];
    dnew.x = ye[i];
    bmp_prec_t q = bmp::sin(bmp_prec_t(xv[i])); // <= True value of sin
    bmp_prec_t dd = bmp_prec_t(dnew.x) - q;
    // Average and std deviation
    div += dd * dd;
    avg += dd;
    dqs.x = double(q);
    double ulp = 0;
    {
      mynumber t;
      t.x = xv[i];
      t.i ^= 1;
      ulp = abs(t.x - xv[i]);
    }

    //fs << std::scientific << xv[i] << " " << double(ye[i] - q) << endl;

    // Bitwise compare of internal sin(double) function and rounded to double "True" value of sin
    if( dold.i != dqs.i )
      co++;

    // Bitwise compare of tested sin(double) function and rounded to double "True" value of sin
    if( dqs.i != dnew.i )
      cn++;
    // Maximum ulp
    co_max = max<double>(co_max, double(dold.x - q) / ulp);
    cn_max = max<double>(cn_max, double(dnew.x - q) / ulp);
  }
  avg /= xv.size();
  div /= xv.size();

  // The number of bitwise wrong results for libm sin(double)
  cout << "libm bitwise error: " <<  co << " / " << xv.size() << "(" << 100.0 * co / xv.size() << "%)" << endl;
  cout << "libm max ULP error: " <<  co_max << endl;

  // The number of bitwise wrong results for tested function sin(double)
  cout << "New ULP error: " <<  cn << " / " << xv.size() << "(" << 100.0 * cn / xv.size() << "%)" << endl;
  cout << "New sin max ULP error: " << cn_max << endl;

  // Average value of deviation and std of deviation for tested function
  cout << "  Avg / std new: " << double(avg) << " / " << double(sqrt( div - avg * avg )) << endl;
  return 0;
}
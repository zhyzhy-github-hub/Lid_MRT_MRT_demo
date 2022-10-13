// --------------------------------------------------------
// -----*****  A program about Lid-Driven, Writen by Zhang Yu
// ---------- Harbin institude of techbology
// --------------------------------------------------------
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using std::cout;
using std::endl;
using std::vector;

constexpr unsigned int Q{9};
constexpr unsigned int Q9{9};
constexpr unsigned int NX{101};
constexpr unsigned int NY{101};
constexpr unsigned int NSAVE{10000};
constexpr unsigned int NMSG{100};
// constexpr unsigned int Max_Iteration{30000};
constexpr unsigned int Max_Iteration{int(30000 / 0.577)};

const size_t Mesh_Size = sizeof(double) * NX * NY;
const size_t Mesh_Size_Scalar = sizeof(double) * NX * NY;
const size_t Mesh_Size_Population = sizeof(double) * NX * NY * Q;

const double U_top{0.0577};
const double U_bottom{0.0};
const double U_right{0.0};
const double U_left{0.0};
const double V_top{0.0};
const double V_bottom{0.0};
const double V_right{0.0};
const double V_left{0.0};
const double Relative_Error{1e-6};

const double dx{1.0};
const double dy{1.0};
const double dt{1.0 * 1};
const double c = dx / dt;
const double Cs{c * 1.0 / sqrt(3.0)};
const double Cs2 = Cs * Cs;

const double p_1 = 1.0 * c;

const double C_x[Q]{0, p_1, 0, -p_1, 0, p_1, -p_1, -p_1, p_1};
const double C_y[Q]{0, 0, p_1, 0, -p_1, p_1, p_1, -p_1, -p_1};

const int l_Cx[Q]{0, 1, 0, -1, 0, 1, -1, -1, 1};
const int l_Cy[Q]{0, 0, 1, 0, -1, 1, 1, -1, -1};

const double w9c{4.0 / 9.0};
const double w9s{1.0 / 9.0};
const double w9d{1.0 / 36.0};
const double w9[Q]{4.0 / 9.0,
                   1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
                   1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};
const int Opp[Q]{0, 3, 4, 1, 2, 7, 8, 5, 6};

const double rho0 = 1.0;

const double Re{1000};

// if TRT -----------------------
const double Magic{1.0 / 4.0};
// TRT -----------------------

double uTop = 0.1;
const double nu = U_top * NX / Re;
const double tau = nu / (Cs2 * dt) + 0.5;

const double s_nu = 1.0 / tau;
const double s_e = s_nu;

const double s_q = 8 * (2 - s_nu) / (8 - s_nu);
const double s_eps = 2.0 / (6 * nu + 1);

// ----------------- s2 and s3 can be 1???
// const double s_Q9[Q9] = {0, s_e, s_eps,
//  0, s_q, 0,
//  s_q, s_nu, s_nu};

// const double s_Q9[Q9] = {1, s_e, s_eps,
//                          1, s_q, 1,
//                          s_q, s_nu, s_nu};
const double s_Q9[Q9] = {1, 0.2, 0.2,
                         1, 1.2, 1,
                         1.2, s_nu, s_nu};

inline size_t scalar_indexQ9(unsigned int x, unsigned int y)
{
    return (NX * y + x);
}
inline size_t fQ9_index(unsigned int x, unsigned int y, unsigned int d)
{
    return (Q * (NX * y + x) + d);
}

inline double Compute_all(double *f, size_t max)
{
    double f_all{0};
    for (size_t i = 0; i < max; ++i)
    {
        f_all += f[i];
    }
    return f_all;
}

inline void meqD2Q9(double *meq, const double &rho, const double &ux, const double &uy)
{
    const double ux2 = ux * ux;
    const double uy2 = uy * uy;

    meq[0] = 1.0 * rho;
    meq[1] = rho * (3 * ux2 + 3 * uy2 - 2 * c * c) / (c * c);
    meq[2] = rho * (-3 * ux2 - 3 * uy2 + 1 * c * c) * 1. / (c * c);
    meq[3] = 1.0 * (rho * ux) * 1. / (c);
    meq[4] = 1.0 * (-rho * ux) * 1. / (c);
    meq[5] = 1.0 * (rho * uy) * 1. / (c);
    meq[6] = 1.0 * (-rho * uy) * 1. / (c);
    meq[7] = 1.0 * rho * (ux2 - uy2) * 1. / (c * c);
    meq[8] = 1.0 * rho * ux * uy * 1. / (c * c);

    // meq[0] = rho;
    // meq[1] = rho * (3 * ux2 + 3 * uy2 - 2);
    // meq[2] = rho * (-3 * ux2 - 3 * uy2 + 1);
    // meq[3] = rho * ux;
    // meq[4] = -rho * ux;
    // meq[5] = rho * uy;
    // meq[6] = -rho * uy;
    // meq[7] = rho * (ux2 - uy2);
    // meq[8] = rho * ux * uy;
}

inline void m2f(double *f, const double *m)
{
    const double m0 = m[0];
    const double m1 = m[1];
    const double m2 = m[2];
    const double m3 = m[3];
    const double m4 = m[4];
    const double m5 = m[5];
    const double m6 = m[6];
    const double m7 = m[7];
    const double m8 = m[8];
    f[0] = m0 / 9 - m1 / 9 + m2 / 9;
    f[1] = m0 / 9 - m1 / 36 - m2 / 18 + m3 / 6 - m4 / 6 + m7 / 4;
    f[2] = m0 / 9 - m1 / 36 - m2 / 18 + m5 / 6 - m6 / 6 - m7 / 4;
    f[3] = m0 / 9 - m1 / 36 - m2 / 18 - m3 / 6 + m4 / 6 + m7 / 4;
    f[4] = m0 / 9 - m1 / 36 - m2 / 18 - m5 / 6 + m6 / 6 - m7 / 4;
    f[5] = m0 / 9 + m1 / 18 + m2 / 36 + m3 / 6 + m4 / 12 + m5 / 6 + m6 / 12 + m8 / 4;
    f[6] = m0 / 9 + m1 / 18 + m2 / 36 - m3 / 6 - m4 / 12 + m5 / 6 + m6 / 12 - m8 / 4;
    f[7] = m0 / 9 + m1 / 18 + m2 / 36 - m3 / 6 - m4 / 12 - m5 / 6 - m6 / 12 + m8 / 4;
    f[8] = m0 / 9 + m1 / 18 + m2 / 36 + m3 / 6 + m4 / 12 - m5 / 6 - m6 / 12 - m8 / 4;
}

inline void f2m(double *m, const double *f)
{
    const double f0 = f[0];
    const double f1 = f[1];
    const double f2 = f[2];
    const double f3 = f[3];
    const double f4 = f[4];
    const double f5 = f[5];
    const double f6 = f[6];
    const double f7 = f[7];
    const double f8 = f[8];
    m[0] = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
    m[1] = -4 * f0 - f1 - f2 - f3 - f4 + 2 * f5 + 2 * f6 + 2 * f7 + 2 * f8;
    m[2] = 4 * f0 - 2 * f1 - 2 * f2 - 2 * f3 - 2 * f4 + f5 + f6 + f7 + f8;
    m[3] = f1 - f3 + f5 - f6 - f7 + f8;
    m[4] = -2 * f1 + 2 * f3 + f5 - f6 - f7 + f8;
    m[5] = f2 - f4 + f5 + f6 - f7 - f8;
    m[6] = -2 * f2 + 2 * f4 + f5 + f6 - f7 - f8;
    m[7] = f1 - f2 + f3 - f4;
    m[8] = f5 - f6 + f7 - f8;
}

void BottomLeftBC_HWBB(double *f0, const double *f1, const double ux, const double uy)
{
    f0[2] = f1[4];
    f0[5] = f1[7];
    f0[1] = f1[3];
    f0[8] = f1[8];
    f0[6] = f1[6];
}
void BottomRightBC_HWBB(double *f0, const double *f1, const double ux, const double uy)
{
    f0[2] = f1[4];
    f0[3] = f1[1];
    f0[6] = f1[8];
    f0[5] = f1[5];
    f0[7] = f1[7];
}
void TopLeftBC_HWBB(double *f0, const double *f1, const double ux, const double uy)
{
    f0[8] = f1[6];
    f0[1] = f1[3];
    f0[4] = f1[2];
    f0[5] = f1[5];
    f0[7] = f1[7];
}
void TopLeftBC_HWBB(double *f0, const double *f1, const double ux, const double uy, const double rho)
{
    f0[1] = f1[3] - 2 * w9[3] * rho * (C_x[3] * ux + C_y[3] * uy) / Cs2;
    f0[4] = f1[2] - 2 * w9[2] * rho * (C_x[2] * ux + C_y[2] * uy) / Cs2;
    f0[8] = f1[6] - 2 * w9[6] * rho * (C_x[6] * ux + C_y[6] * uy) / Cs2;
    f0[5] = f1[5];
    f0[7] = f1[7];
}
void TopRightBC_HWBB(double *f0, const double *f1, const double ux, const double uy)
{
    f0[7] = f1[5];
    f0[3] = f1[1];
    f0[4] = f1[2];
    f0[8] = f1[6];
    f0[6] = f1[8];
}
void TopRightBC_HWBB(double *f0, const double *f1, const double ux, const double uy, const double rho)
{
    f0[7] = f1[5] - 2 * w9[5] * rho * (C_x[5] * ux + C_y[5] * uy) / Cs2;
    f0[3] = f1[1] - 2 * w9[1] * rho * (C_x[1] * ux + C_y[1] * uy) / Cs2;
    f0[4] = f1[2] - 2 * w9[2] * rho * (C_x[2] * ux + C_y[2] * uy) / Cs2;
    f0[8] = f1[8];
    f0[6] = f1[6];
}

void BottomBC_HWBB(double *f0, const double *f1, const double ux, const double uy)
{
    const double sixth = 1.0 / 6.0;
    f0[2] = f1[4];
    f0[5] = f1[7]; // + sixth * ux;
    f0[6] = f1[8]; // - sixth * ux;
}
void TopBC_HWBB(double *f0, const double *f1, const double ux, const double uy, const double rho)
{
    f0[4] = f1[2] - 2 * w9[2] * rho * (C_x[2] * ux + C_y[2] * uy) / Cs2;
    f0[7] = f1[5] - 2 * w9[5] * rho * (C_x[5] * ux + C_y[5] * uy) / Cs2;
    f0[8] = f1[6] - 2 * w9[6] * rho * (C_x[6] * ux + C_y[6] * uy) / Cs2;
}

void TopBC_HWBB(double *f0, const double *f1, const double ux, const double uy)
{
    const double sixth = 1.0 / 6.0;
    f0[4] = f1[2];
    f0[7] = f1[5] - sixth * ux;
    f0[8] = f1[6] + sixth * ux;
}
void LeftBC_HWBB(double *f0, const double *f1, const double ux, const double uy)
{
    const double sixth = 1.0 / 6.0;
    f0[1] = f1[3];
    f0[5] = f1[7]; // - sixth * ux;
    f0[8] = f1[6]; // + sixth * ux;
}
void RightBC_HWBB(double *f0, const double *f1, const double ux, const double uy)
{
    const double sixth = 1.0 / 6.0;
    f0[3] = f1[1];
    f0[7] = f1[5]; // - sixth * ux;
    f0[6] = f1[8]; // + sixth * ux;
}

void Initialzation(double *rho, double *ux, double *uy)
{
    for (size_t i = 0; i < NX; ++i)
    {
        for (size_t j = 0; j < NY; ++j)
        {
            const size_t index = scalar_indexQ9(i, j);
            ux[index] = 0.0;
            uy[index] = 0.0;
            rho[index] = rho0;
            // if(i < 60 && i > 40 && j < 60 && j > 40){
            //     ux[index]  = 0.1;
            //     uy[index]  = 0.1;
            //     //rho[index] = 0.01;
            // }
        }
    }
}

void Init_f(double *f0, double *f1, double *ux, double *uy, double *rho)
{
    for (size_t i = 0; i < NX; ++i)
    {
        for (size_t j = 0; j < NY; ++j)
        {

            const size_t index = scalar_indexQ9(i, j);
            const size_t fQ9 = fQ9_index(i, j, 0);

            double feqQ9[Q9], meqQ9[Q9];
            double rho_ = rho[index];
            double ux_ = ux[index];
            double uy_ = uy[index];

            meqD2Q9(meqQ9, rho_, ux_, uy_);
            m2f(feqQ9, meqQ9);
            for (size_t k = 0; k < Q9; ++k)
            {
                const size_t fQ9_k = fQ9 + k;
                f1[fQ9_k] = meqQ9[k];
                f0[fQ9_k] = feqQ9[k];
            }
        }
    }
}

inline void Comp_Macro(const double *f0, double &rho, double &ux, double &uy)
{
    double rho_ = 0;
    double uy_ = 0;
    double ux_ = 0;
    for (size_t k = 0; k < Q9; ++k)
    {
        rho_ += f0[k];
        ux_ += f0[k] * C_x[k];
        uy_ += f0[k] * C_y[k];
    }
    rho = rho_;
    ux = ux_ / rho_;
    uy = uy_ / rho_;
}

inline void StreamingQ9(size_t i, size_t j, double *f0, double *f1)
{
    for (size_t k = 0; k < Q9; ++k)
    {
        int ip = (NX + i - l_Cx[k]) % NX;
        int jp = (NY + j - l_Cy[k]) % NY;
        size_t f_p_index = fQ9_index(ip, jp, k);
        f0[k] = f1[f_p_index];
    }
}

inline void CollisionQ9_MRT(double *m1, const double *m0, const double *meq0, const double *s, const double dt_)
{
    m1[0] = m0[0] - s[0] * (m0[0] - meq0[0]);
    m1[1] = m0[1] - s[1] * (m0[1] - meq0[1]);
    m1[2] = m0[2] - s[2] * (m0[2] - meq0[2]);
    m1[3] = m0[3] - s[3] * (m0[3] - meq0[3]);
    m1[4] = m0[4] - s[4] * (m0[4] - meq0[4]);
    m1[5] = m0[5] - s[5] * (m0[5] - meq0[5]);
    m1[6] = m0[6] - s[6] * (m0[6] - meq0[6]);
    m1[7] = m0[7] - s[7] * (m0[7] - meq0[7]);
    m1[8] = m0[8] - s[8] * (m0[8] - meq0[8]);
    // m1[0] = m0[0];
    // m1[1] = m0[1] - s_e * (m0[1] - meq0[1]);
    // m1[2] = m0[2] - s_eps * (m0[2] - meq0[2]);
    // m1[3] = m0[3];
    // m1[4] = m0[4] - s_q * (m0[4] - meq0[4]);
    // m1[5] = m0[5];
    // m1[6] = m0[6] - s_q * (m0[6] - meq0[6]);
    // m1[7] = m0[7] - s_nu * (m0[7] - meq0[7]);
    // m1[8] = m0[8] - s_nu * (m0[8] - meq0[8]);
}

void NSE(double *f0, double *f1, double *rho, double *ux, double *uy)
{

    double *f_post = (double *)malloc(Mesh_Size_Population);
    for (size_t j = 0; j < NY; ++j)
    {
        for (size_t i = 0; i < NX; ++i)
        {

            const size_t index = scalar_indexQ9(i, j);
            const size_t fQ9 = fQ9_index(i, j, 0);

            // ------  Compute the conserved quantities delta_rho and u

            double rho_ = rho[index];
            double ux_ = ux[index];
            double uy_ = uy[index];

            Comp_Macro(f0 + fQ9, rho_, ux_, uy_);
            rho[index] = rho_;
            ux[index] = ux_;
            uy[index] = uy_;

            double m0[Q9];
            double meq0[Q9];
            f2m(m0, f0 + fQ9);
            meqD2Q9(meq0, rho_, ux_, uy_);
            CollisionQ9_MRT(f1 + fQ9, m0, meq0, s_Q9, dt);
            m2f(f_post + fQ9, f1 + fQ9);
        }
    }

    for (size_t j = 0; j < NY; ++j)
    {
        for (size_t i = 0; i < NX; ++i)
        {
            const size_t index = scalar_indexQ9(i, j);
            const size_t fQ9 = fQ9_index(i, j, 0);
            StreamingQ9(i, j, f0 + fQ9, f_post);
        }
    }
    for (size_t j = 0; j < NY; ++j)
    {
        for (size_t i = 0; i < NX; ++i)
        {
            const size_t index = scalar_indexQ9(i, j);
            const size_t fQ9 = fQ9_index(i, j, 0);
            double ux0, uy0, ux1, uy1;
            // if(j == 0){
            if (j == 0 && i != 0 && i != NX - 1)
            {
                const size_t fQ91 = fQ9_index(i, j + 1, 0);
                const size_t id1 = scalar_indexQ9(i, j + 1);
                double feq0[Q9], feq1[Q9];
                ux0 = 0;
                uy0 = 0;
                ux1 = ux[id1];
                uy1 = uy[id1];
                BottomBC_HWBB(f0 + fQ9, f_post + fQ9, 0, 0);
            }
            else if (i == 0 && j == 0)
            {
                BottomLeftBC_HWBB(f0 + fQ9, f_post + fQ9, 0, 0);
                // BottomLeftBC_HWBB(f0 + fQ9, f0 + fQ9, 0, 0);
            }
            else if (i == NX - 1 && j == 0)
            {
                BottomRightBC_HWBB(f0 + fQ9, f_post + fQ9, 0, 0);
                // BottomRightBC_HWBB(f0 + fQ9, f0 + fQ9, 0, 0);
            }
            else if (i == 0 && j != 0 && j != NY - 1)
            {
                const size_t fQ91 = fQ9_index(i + 1, j, 0);
                const size_t id1 = scalar_indexQ9(i + 1, j);
                double feq0[Q9], feq1[Q9];
                LeftBC_HWBB(f0 + fQ9, f_post + fQ9, 0, 0);
            }
            else if (i == NX - 1 && j != 0 && j != NY - 1)
            {
                const size_t fQ91 = fQ9_index(i - 1, j, 0);
                const size_t id1 = scalar_indexQ9(i - 1, j);
                double feq0[Q9], feq1[Q9];
                RightBC_HWBB(f0 + fQ9, f_post + fQ9, 0, 0);
            }
            else if (i == 0 && j == NY - 1)
            {
                // TopLeftBC_HWBB(f0 + fQ9, f_post + fQ9, 1 * U_top, 0, rho[index]);
                TopLeftBC_HWBB(f0 + fQ9, f_post + fQ9, 0 * U_top, 0);
            }
            else if (i == NX - 1 && j == NY - 1)
            {
                // TopRightBC_HWBB(f0 + fQ9, f_post + fQ9, 1 * U_top, 0, rho[index]);
                TopRightBC_HWBB(f0 + fQ9, f_post + fQ9, 1 * U_top, 0);
            }
            // else if(j == NY - 1){
            else if (j == NY - 1 && i != 0 && i != NX - 1)
            {
                const size_t fQ91 = fQ9_index(i, j - 1, 0);
                const size_t id1 = scalar_indexQ9(i, j - 1);
                double feq0[Q9], feq1[Q9];
                double rho_b = rho[index];
                TopBC_HWBB(f0 + fQ9, f_post + fQ9, 1 * U_top, 0, rho_b);
                // TopBC_HWBB(f0 + fQ9, f_post + fQ9, 1 * U_top, 0, rho0);
            }
        }
    }
    free(f_post);
}

void outputTec(int m, double *f0, double *f1, double *rho, double *ux, double *uy)
{
    std::ostringstream name;
    name << "Lid_Re_" << Re << "_" << c << "_" << Cs2 << "_" << m << ".dat";
    std::ofstream out(name.str().c_str());
    out << "Title= \"Lid Driven Flow\"\n";
    out << "VARIABLES = \"X\", \"Y\", \"p\",\"rho\", \"U\", \"V\", \"UV\"\n";
    out << "ZONE T= \"BOX\",I=" << NX << ",J=" << NY << ",F=	POINT" << endl;

    for (size_t j = 0; j < NY; ++j)
    {
        for (size_t i = 0; i < NX; ++i)
        {
            const size_t index = scalar_indexQ9(i, j);
            const size_t fQ9 = fQ9_index(i, j, 0);
            out << std::fixed << std::setprecision(10)
                << " " << i << " " << j
                << " " << (rho[index] - rho0) * Cs2
                << " " << rho[index]
                << " " << ux[index] << " " << uy[index]
                << " " << sqrt(ux[index] * ux[index] + uy[index] * uy[index])
                << endl;
        }
    }
}
void Print_lattice_para(std::ostream &out)
{
    using std::endl;

    out << " Q      = " << Q << endl;
    out << " Q9     = " << Q9 << endl;
    out << " NX     = " << NX << endl;
    out << " NY     = " << NY << endl;
    out << " NSAVE  = " << NSAVE << endl;
    out << " NMSG   = " << NMSG << endl;
    out << " Max_Iteration  = " << Max_Iteration << endl;
    out << " dx  = " << dx << endl;
    out << " dy  = " << dy << endl;
    out << " dt  = " << dt << endl;
    out << " c   = " << c << endl;
    out << " Cs  = " << Cs << endl;
    out << " Cs2 = " << Cs2 << endl;
    out << " nu    = " << nu << endl;
    out << " tau   = " << tau << endl;
    out << " s_nu  = " << s_nu << endl;
    out << " s_e   = " << s_e << endl;
    out << " s_q   = " << s_q << endl;
    out << " s_eps = " << s_eps << endl;
    out << " s_0 = " << s_Q9[0] << endl;
    out << " s_1 = " << s_Q9[1] << endl;
    out << " s_2 = " << s_Q9[2] << endl;
    out << " s_3 = " << s_Q9[3] << endl;
    out << " s_4 = " << s_Q9[4] << endl;
    out << " s_5 = " << s_Q9[5] << endl;
    out << " s_6 = " << s_Q9[6] << endl;
    out << " s_7 = " << s_Q9[7] << endl;
    out << " s_8 = " << s_Q9[8] << endl;
}

int main()
{
    double *f0 = (double *)malloc(Mesh_Size_Population);
    double *f1 = (double *)malloc(Mesh_Size_Population);
    double *ux = (double *)malloc(Mesh_Size_Scalar);
    double *uy = (double *)malloc(Mesh_Size_Scalar);
    double *rho = (double *)malloc(Mesh_Size_Scalar);

    double dx1 = 0.5;
    double dx2 = 0.126;

    size_t k = dx1 / dx2;
    cout << "k = " << k << endl;

    Initialzation(rho, ux, uy);
    Init_f(f0, f1, ux, uy, rho);

    double rho_all = Compute_all(rho, NX * NY);
    double uy_all = Compute_all(uy, NX * NY);
    double ux_all = Compute_all(ux, NX * NY);

    std::ostringstream name;
    name << "Re_" << Re << "_" << c << "_"
         << ".txt";
    std::ofstream out(name.str().c_str());
    Print_lattice_para(out);

    cout << "nu = " << nu << ", tau = " << tau << endl;
    for (size_t i = 0; i < Q9; ++i)
    {
        cout << "S_" << i << " = " << s_Q9[i] << "\n";
    }
    cout << "----------------------------------------------\n";

    cout << "n = " << 0 << " , "
         << "ux_all = " << ux_all << " , "
         << "uy_all = " << uy_all << " , "
         << "rho_all = " << rho_all << " , "
         << endl;

    outputTec(0, f0, f1, rho, ux, uy);

    size_t n = 0;
    // for (n = 0; n < 1000000; ++n)
    for (n = 0; n < Max_Iteration; ++n)
    {
        NSE(f0, f1, rho, ux, uy);

        double rho_all = Compute_all(rho, NX * NY);
        double uy_all = Compute_all(uy, NX * NY);
        double ux_all = Compute_all(ux, NX * NY);

        size_t ndt = floor(n * dt);
        size_t nmsg = NMSG / dt;
        size_t nsave = NSAVE / dt;

        if ((n) % nmsg == 0 && n >= nmsg)
        {
            cout << "----------------------------------------------\n";
            cout
                << "n*dt = " << ndt << " , "
                << "n = " << n << " , "
                << "ux_all = " << ux_all << " , "
                << "uy_all = " << uy_all << " , "
                << "rho_all = " << rho_all << " , "
                << endl;
        }
        if ((n) % nsave == 0 && n >= nsave)
        {
            outputTec(ndt, f0, f1, rho, ux, uy);
        }
        if (ndt >= 100000)
        {
            break;
        }
    }
    outputTec(n, f0, f1, rho, ux, uy);

    out.close();

    return 0;
}

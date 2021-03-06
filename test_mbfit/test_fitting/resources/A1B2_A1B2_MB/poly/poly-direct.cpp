#include "poly-model.h"
    
namespace mb_system {

double poly_model::eval_direct(const double a[24], const double x[15])
{
    double p[24];
    p[0] = x[13] + x[11] + x[10] + x[12];
    p[1] = x[6] + x[4] + x[5] + x[3];
    p[2] = x[0];
    p[3] = x[1]*x[13] + x[7]*x[11] + x[2]*x[11] + x[2]*x[10] + x[8]*x[10] + x[7]*x[13] + x[8]*x[12] + x[1]*x[12];
    p[4] = x[1]*x[11] + x[2]*x[13] + x[7]*x[12] + x[1]*x[10] + x[8]*x[11] + x[2]*x[12] + x[7]*x[10] + x[8]*x[13];
    p[5] = x[2]*x[6] + x[4]*x[8] + x[1]*x[5] + x[3]*x[7];
    p[6] = x[13]*x[14] + x[11]*x[14] + x[9]*x[12] + x[9]*x[11] + x[12]*x[14] + x[9]*x[10] + x[10]*x[14] + x[9]*x[13];
    p[7] = x[6]*x[9] + x[3]*x[14] + x[4]*x[14] + x[5]*x[9];
    p[8] = x[0]*x[4] + x[0]*x[5] + x[0]*x[3] + x[0]*x[6];
    p[9] = x[3]*x[9] + x[6]*x[14] + x[4]*x[9] + x[5]*x[14];
    p[10] = x[12]*x[13] + x[10]*x[11] + x[11]*x[13] + x[10]*x[12];
    p[11] = x[0]*x[0];
    p[12] = x[5]*x[13] + x[4]*x[12] + x[3]*x[11] + x[6]*x[11] + x[5]*x[12] + x[3]*x[13] + x[4]*x[10] + x[6]*x[10];
    p[13] = x[11]*x[11] + x[13]*x[13] + x[12]*x[12] + x[10]*x[10];
    p[14] = x[0]*x[7] + x[0]*x[2] + x[0]*x[1] + x[0]*x[8];
    p[15] = x[5]*x[5] + x[6]*x[6] + x[3]*x[3] + x[4]*x[4];
    p[16] = x[3]*x[8] + x[2]*x[5] + x[1]*x[6] + x[4]*x[7];
    p[17] = x[11]*x[12] + x[10]*x[13];
    p[18] = x[1]*x[3] + x[5]*x[8] + x[5]*x[7] + x[2]*x[3] + x[6]*x[7] + x[6]*x[8] + x[2]*x[4] + x[1]*x[4];
    p[19] = x[3]*x[4] + x[5]*x[6];
    p[20] = x[0]*x[12] + x[0]*x[11] + x[0]*x[10] + x[0]*x[13];
    p[21] = x[0]*x[14] + x[0]*x[9];
    p[22] = x[4]*x[5] + x[3]*x[5] + x[3]*x[6] + x[4]*x[6];
    p[23] = x[3]*x[12] + x[4]*x[11] + x[5]*x[11] + x[6]*x[12] + x[6]*x[13] + x[5]*x[10] + x[4]*x[13] + x[3]*x[10];

    double energy(0);
    for(int i = 0; i < 24; ++i)
        energy += p[i]*a[i];

    return energy;

}
} // namespace mb_system
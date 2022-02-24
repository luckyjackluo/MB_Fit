#include "mbnrg_1b_A1B2X2_A1B2X2_deg4_v1.h"

////////////////////////////////////////////////////////////////////////////////

namespace mbnrg_A1B2X2_A1B2X2_deg4 {

mbnrg_A1B2X2_A1B2X2_deg4_v1::mbnrg_A1B2X2_A1B2X2_deg4_v1(const std::string mon1) {

    // =====>> BEGIN SECTION CONSTRUCTOR <<=====
    // =>> PASTE RIGHT BELOW THIS LINE <==


    // =====>> END SECTION CONSTRUCTOR <<=====
}

//----------------------------------------------------------------------------//

double mbnrg_A1B2X2_A1B2X2_deg4_v1::f_switch(const double r, double& g)
{
    if (r > m_ro) {
        g = 0.0;
        return 0.0;
    } else if (r > m_ri) {
        const double t1 = M_PI/(m_ro - m_ri);
        const double x = (r - m_ri)*t1;
        g = - std::sin(x)*t1/2.0;
        return (1.0 + std::cos(x))/2.0;
    } else {
        g = 0.0;
        return 1.0;
    }
}

//----------------------------------------------------------------------------//

 std::vector<double> mbnrg_A1B2X2_A1B2X2_deg4_v1::eval(const double *xyz1, const size_t n) {
    std::vector<double> energies(n,0.0);
    std::vector<double> energies_sw(n,0.0);

    std::vector<double> xyz(18);
    double sw = 0.0;
    polynomial my_poly;

    for (size_t j = 0; j < n; j++) {
        const double *mon1 = xyz1 + 18*j;


        if (false ) {
             continue;
        }

        std::copy(mon1, mon1 + 18, xyz.begin() + 0);


        const double* coords_A_1_a = xyz.data() + 0;

        const double* coords_B_1_a = xyz.data() + 3;

        const double* coords_B_2_a = xyz.data() + 6;

        const double* coords_A_2_b = xyz.data() + 9;

        const double* coords_B_3_b = xyz.data() + 12;

        const double* coords_B_4_b = xyz.data() + 15;

        double coords_X_1_a[3];
        std::fill(coords_X_1_a,coords_X_1_a + 3, 0.0);

        double coords_X_2_a[3];
        std::fill(coords_X_2_a,coords_X_2_a + 3, 0.0);

        double coords_X_3_b[3];
        std::fill(coords_X_3_b,coords_X_3_b + 3, 0.0);

        double coords_X_4_b[3];
        std::fill(coords_X_4_b,coords_X_4_b + 3, 0.0);


        double w12 =     -9.721486914088159e-02;  //from MBpol
        double w13 =     -9.721486914088159e-02;
        double wcross =   9.859272078406150e-02;

    
        monomer m1;
        m1.setup(coords_A_1_a, w12, wcross, coords_X_1_a, coords_X_2_a);

        variable vs[31];

        double xs[31];

        xs[0] = vs[0].v_coul(m_k_x_inter_A_A_0, coords_A_1_a, coords_A_2_b);
        xs[1] = vs[1].v_exp(m_k_x_intra_A_B_1, coords_A_1_a, coords_B_1_a);
        xs[2] = vs[2].v_exp(m_k_x_intra_A_B_1, coords_A_1_a, coords_B_2_a);
        xs[3] = vs[3].v_coul(m_k_x_inter_A_B_0, coords_A_1_a, coords_B_3_b);
        xs[4] = vs[4].v_coul(m_k_x_inter_A_B_0, coords_A_1_a, coords_B_4_b);
        xs[5] = vs[5].v_exp(m_k_x_inter_A_X_0, coords_A_1_a, coords_X_3_b);
        xs[6] = vs[6].v_exp(m_k_x_inter_A_X_0, coords_A_1_a, coords_X_4_b);
        xs[7] = vs[7].v_coul(m_k_x_inter_A_B_0, coords_A_2_b, coords_B_1_a);
        xs[8] = vs[8].v_coul(m_k_x_inter_A_B_0, coords_A_2_b, coords_B_2_a);
        xs[9] = vs[9].v_exp(m_k_x_intra_A_B_1, coords_A_2_b, coords_B_3_b);
        xs[10] = vs[10].v_exp(m_k_x_intra_A_B_1, coords_A_2_b, coords_B_4_b);
        xs[11] = vs[11].v_exp(m_k_x_inter_A_X_0, coords_A_2_b, coords_X_1_a);
        xs[12] = vs[12].v_exp(m_k_x_inter_A_X_0, coords_A_2_b, coords_X_2_a);
        xs[13] = vs[13].v_exp(m_k_x_intra_B_B_1, coords_B_1_a, coords_B_2_a);
        xs[14] = vs[14].v_coul(m_k_x_inter_B_B_0, coords_B_1_a, coords_B_3_b);
        xs[15] = vs[15].v_coul(m_k_x_inter_B_B_0, coords_B_1_a, coords_B_4_b);
        xs[16] = vs[16].v_exp(m_k_x_inter_B_X_0, coords_B_1_a, coords_X_3_b);
        xs[17] = vs[17].v_exp(m_k_x_inter_B_X_0, coords_B_1_a, coords_X_4_b);
        xs[18] = vs[18].v_coul(m_k_x_inter_B_B_0, coords_B_2_a, coords_B_3_b);
        xs[19] = vs[19].v_coul(m_k_x_inter_B_B_0, coords_B_2_a, coords_B_4_b);
        xs[20] = vs[20].v_exp(m_k_x_inter_B_X_0, coords_B_2_a, coords_X_3_b);
        xs[21] = vs[21].v_exp(m_k_x_inter_B_X_0, coords_B_2_a, coords_X_4_b);
        xs[22] = vs[22].v_exp(m_k_x_intra_B_B_1, coords_B_3_b, coords_B_4_b);
        xs[23] = vs[23].v_exp(m_k_x_inter_B_X_0, coords_B_3_b, coords_X_1_a);
        xs[24] = vs[24].v_exp(m_k_x_inter_B_X_0, coords_B_3_b, coords_X_2_a);
        xs[25] = vs[25].v_exp(m_k_x_inter_B_X_0, coords_B_4_b, coords_X_1_a);
        xs[26] = vs[26].v_exp(m_k_x_inter_B_X_0, coords_B_4_b, coords_X_2_a);
        xs[27] = vs[27].v_exp(m_k_x_inter_X_X_0, coords_X_1_a, coords_X_3_b);
        xs[28] = vs[28].v_exp(m_k_x_inter_X_X_0, coords_X_1_a, coords_X_4_b);
        xs[29] = vs[29].v_exp(m_k_x_inter_X_X_0, coords_X_2_a, coords_X_3_b);
        xs[30] = vs[30].v_exp(m_k_x_inter_X_X_0, coords_X_2_a, coords_X_4_b);


        sw = 1.0;

        energies[j] = my_poly.eval(xs,coefficients.data());
        energies_sw[j] = energies[j]*sw;

    }
    return energies_sw;

}

//----------------------------------------------------------------------------//

std::vector<double> mbnrg_A1B2X2_A1B2X2_deg4_v1::eval(const double *xyz1, double *grad1 , const size_t n, std::vector<double> *virial) {
    std::vector<double> energies(n,0.0);
    std::vector<double> energies_sw(n,0.0);

    std::vector<double> xyz(18);
    double sw = 0.0;
    polynomial my_poly;

    for (size_t j = 0; j < n; j++) {
        const double *mon1 = xyz1 + 18*j;


        if (false ) {
             continue;
        }

        std::vector<double> gradients(18,0.0);

        std::copy(mon1, mon1 + 18, xyz.begin() + 0);
        const double* coords_A_1_a = xyz.data() + 0;

        const double* coords_B_1_a = xyz.data() + 3;

        const double* coords_B_2_a = xyz.data() + 6;

        const double* coords_A_2_b = xyz.data() + 9;

        const double* coords_B_3_b = xyz.data() + 12;

        const double* coords_B_4_b = xyz.data() + 15;

        double coords_X_1_a[3];
        std::fill(coords_X_1_a,coords_X_1_a + 3, 0.0);

        double coords_X_2_a[3];
        std::fill(coords_X_2_a,coords_X_2_a + 3, 0.0);

        double coords_X_3_b[3];
        std::fill(coords_X_3_b,coords_X_3_b + 3, 0.0);

        double coords_X_4_b[3];
        std::fill(coords_X_4_b,coords_X_4_b + 3, 0.0);


        double* coords_A_1_a_g = gradients.data() + 0;

        double* coords_B_1_a_g = gradients.data() + 3;

        double* coords_B_2_a_g = gradients.data() + 6;

        double* coords_A_2_b_g = gradients.data() + 9;

        double* coords_B_3_b_g = gradients.data() + 12;

        double* coords_B_4_b_g = gradients.data() + 15;

        double coords_X_1_a_g[3];
        std::fill(coords_X_1_a_g,coords_X_1_a_g + 3, 0.0);

        double coords_X_2_a_g[3];
        std::fill(coords_X_2_a_g,coords_X_2_a_g + 3, 0.0);

        double coords_X_3_b_g[3];
        std::fill(coords_X_3_b_g,coords_X_3_b_g + 3, 0.0);

        double coords_X_4_b_g[3];
        std::fill(coords_X_4_b_g,coords_X_4_b_g + 3, 0.0);



        double w12 =     -9.721486914088159e-02;  //from MBpol
        double w13 =     -9.721486914088159e-02;
        double wcross =   9.859272078406150e-02;

    
        monomer m1;
        m1.setup(coords_A_1_a, w12, wcross, coords_X_1_a, coords_X_2_a);

        variable vs[31];

        double xs[31];


        double gxs[31];

        xs[0] = vs[0].v_coul(m_k_x_inter_A_A_0, coords_A_1_a, coords_A_2_b);
        xs[1] = vs[1].v_exp(m_k_x_intra_A_B_1, coords_A_1_a, coords_B_1_a);
        xs[2] = vs[2].v_exp(m_k_x_intra_A_B_1, coords_A_1_a, coords_B_2_a);
        xs[3] = vs[3].v_coul(m_k_x_inter_A_B_0, coords_A_1_a, coords_B_3_b);
        xs[4] = vs[4].v_coul(m_k_x_inter_A_B_0, coords_A_1_a, coords_B_4_b);
        xs[5] = vs[5].v_exp(m_k_x_inter_A_X_0, coords_A_1_a, coords_X_3_b);
        xs[6] = vs[6].v_exp(m_k_x_inter_A_X_0, coords_A_1_a, coords_X_4_b);
        xs[7] = vs[7].v_coul(m_k_x_inter_A_B_0, coords_A_2_b, coords_B_1_a);
        xs[8] = vs[8].v_coul(m_k_x_inter_A_B_0, coords_A_2_b, coords_B_2_a);
        xs[9] = vs[9].v_exp(m_k_x_intra_A_B_1, coords_A_2_b, coords_B_3_b);
        xs[10] = vs[10].v_exp(m_k_x_intra_A_B_1, coords_A_2_b, coords_B_4_b);
        xs[11] = vs[11].v_exp(m_k_x_inter_A_X_0, coords_A_2_b, coords_X_1_a);
        xs[12] = vs[12].v_exp(m_k_x_inter_A_X_0, coords_A_2_b, coords_X_2_a);
        xs[13] = vs[13].v_exp(m_k_x_intra_B_B_1, coords_B_1_a, coords_B_2_a);
        xs[14] = vs[14].v_coul(m_k_x_inter_B_B_0, coords_B_1_a, coords_B_3_b);
        xs[15] = vs[15].v_coul(m_k_x_inter_B_B_0, coords_B_1_a, coords_B_4_b);
        xs[16] = vs[16].v_exp(m_k_x_inter_B_X_0, coords_B_1_a, coords_X_3_b);
        xs[17] = vs[17].v_exp(m_k_x_inter_B_X_0, coords_B_1_a, coords_X_4_b);
        xs[18] = vs[18].v_coul(m_k_x_inter_B_B_0, coords_B_2_a, coords_B_3_b);
        xs[19] = vs[19].v_coul(m_k_x_inter_B_B_0, coords_B_2_a, coords_B_4_b);
        xs[20] = vs[20].v_exp(m_k_x_inter_B_X_0, coords_B_2_a, coords_X_3_b);
        xs[21] = vs[21].v_exp(m_k_x_inter_B_X_0, coords_B_2_a, coords_X_4_b);
        xs[22] = vs[22].v_exp(m_k_x_intra_B_B_1, coords_B_3_b, coords_B_4_b);
        xs[23] = vs[23].v_exp(m_k_x_inter_B_X_0, coords_B_3_b, coords_X_1_a);
        xs[24] = vs[24].v_exp(m_k_x_inter_B_X_0, coords_B_3_b, coords_X_2_a);
        xs[25] = vs[25].v_exp(m_k_x_inter_B_X_0, coords_B_4_b, coords_X_1_a);
        xs[26] = vs[26].v_exp(m_k_x_inter_B_X_0, coords_B_4_b, coords_X_2_a);
        xs[27] = vs[27].v_exp(m_k_x_inter_X_X_0, coords_X_1_a, coords_X_3_b);
        xs[28] = vs[28].v_exp(m_k_x_inter_X_X_0, coords_X_1_a, coords_X_4_b);
        xs[29] = vs[29].v_exp(m_k_x_inter_X_X_0, coords_X_2_a, coords_X_3_b);
        xs[30] = vs[30].v_exp(m_k_x_inter_X_X_0, coords_X_2_a, coords_X_4_b);


        sw = 1.0;

        energies[j] = my_poly.eval(xs,coefficients.data(),gxs);
        energies_sw[j] = energies[j]*sw;

        for (size_t i = 0; i < 31; i++) {
            gxs[i] *= sw;
        }

        vs[0].grads(gxs[0], coords_A_1_a_g, coords_A_2_b_g, coords_A_1_a, coords_A_2_b);
        vs[1].grads(gxs[1], coords_A_1_a_g, coords_B_1_a_g, coords_A_1_a, coords_B_1_a);
        vs[2].grads(gxs[2], coords_A_1_a_g, coords_B_2_a_g, coords_A_1_a, coords_B_2_a);
        vs[3].grads(gxs[3], coords_A_1_a_g, coords_B_3_b_g, coords_A_1_a, coords_B_3_b);
        vs[4].grads(gxs[4], coords_A_1_a_g, coords_B_4_b_g, coords_A_1_a, coords_B_4_b);
        vs[5].grads(gxs[5], coords_A_1_a_g, coords_X_3_b_g, coords_A_1_a, coords_X_3_b);
        vs[6].grads(gxs[6], coords_A_1_a_g, coords_X_4_b_g, coords_A_1_a, coords_X_4_b);
        vs[7].grads(gxs[7], coords_A_2_b_g, coords_B_1_a_g, coords_A_2_b, coords_B_1_a);
        vs[8].grads(gxs[8], coords_A_2_b_g, coords_B_2_a_g, coords_A_2_b, coords_B_2_a);
        vs[9].grads(gxs[9], coords_A_2_b_g, coords_B_3_b_g, coords_A_2_b, coords_B_3_b);
        vs[10].grads(gxs[10], coords_A_2_b_g, coords_B_4_b_g, coords_A_2_b, coords_B_4_b);
        vs[11].grads(gxs[11], coords_A_2_b_g, coords_X_1_a_g, coords_A_2_b, coords_X_1_a);
        vs[12].grads(gxs[12], coords_A_2_b_g, coords_X_2_a_g, coords_A_2_b, coords_X_2_a);
        vs[13].grads(gxs[13], coords_B_1_a_g, coords_B_2_a_g, coords_B_1_a, coords_B_2_a);
        vs[14].grads(gxs[14], coords_B_1_a_g, coords_B_3_b_g, coords_B_1_a, coords_B_3_b);
        vs[15].grads(gxs[15], coords_B_1_a_g, coords_B_4_b_g, coords_B_1_a, coords_B_4_b);
        vs[16].grads(gxs[16], coords_B_1_a_g, coords_X_3_b_g, coords_B_1_a, coords_X_3_b);
        vs[17].grads(gxs[17], coords_B_1_a_g, coords_X_4_b_g, coords_B_1_a, coords_X_4_b);
        vs[18].grads(gxs[18], coords_B_2_a_g, coords_B_3_b_g, coords_B_2_a, coords_B_3_b);
        vs[19].grads(gxs[19], coords_B_2_a_g, coords_B_4_b_g, coords_B_2_a, coords_B_4_b);
        vs[20].grads(gxs[20], coords_B_2_a_g, coords_X_3_b_g, coords_B_2_a, coords_X_3_b);
        vs[21].grads(gxs[21], coords_B_2_a_g, coords_X_4_b_g, coords_B_2_a, coords_X_4_b);
        vs[22].grads(gxs[22], coords_B_3_b_g, coords_B_4_b_g, coords_B_3_b, coords_B_4_b);
        vs[23].grads(gxs[23], coords_B_3_b_g, coords_X_1_a_g, coords_B_3_b, coords_X_1_a);
        vs[24].grads(gxs[24], coords_B_3_b_g, coords_X_2_a_g, coords_B_3_b, coords_X_2_a);
        vs[25].grads(gxs[25], coords_B_4_b_g, coords_X_1_a_g, coords_B_4_b, coords_X_1_a);
        vs[26].grads(gxs[26], coords_B_4_b_g, coords_X_2_a_g, coords_B_4_b, coords_X_2_a);
        vs[27].grads(gxs[27], coords_X_1_a_g, coords_X_3_b_g, coords_X_1_a, coords_X_3_b);
        vs[28].grads(gxs[28], coords_X_1_a_g, coords_X_4_b_g, coords_X_1_a, coords_X_4_b);
        vs[29].grads(gxs[29], coords_X_2_a_g, coords_X_3_b_g, coords_X_2_a, coords_X_3_b);
        vs[30].grads(gxs[30], coords_X_2_a_g, coords_X_4_b_g, coords_X_2_a, coords_X_4_b);

        m1.grads(coords_X_1_a_g, coords_X_2_a_g, w12, wcross, coords_A_1_a_g);


        for (size_t i = 0; i < 3; i++) {
            gradients[0 + i] += 0.0 ;
        }


        for (size_t i = 0; i < 18; i++) {
            grad1[i + j*18] += gradients[0 + i];
        }

        
        if (virial != 0) {
        
            (*virial)[0] += -coords_A_1_a[0]*coords_A_1_a_g[0]
                        -coords_B_1_a[0]*coords_B_1_a_g[0]
                        -coords_B_2_a[0]*coords_B_2_a_g[0]
                        -coords_A_2_b[0]*coords_A_2_b_g[0]
                        -coords_B_3_b[0]*coords_B_3_b_g[0]
                        -coords_B_4_b[0]*coords_B_4_b_g[0];

        
            (*virial)[1] += -coords_A_1_a[0]*coords_A_1_a_g[1]
                        -coords_B_1_a[0]*coords_B_1_a_g[1]
                        -coords_B_2_a[0]*coords_B_2_a_g[1]
                        -coords_A_2_b[0]*coords_A_2_b_g[1]
                        -coords_B_3_b[0]*coords_B_3_b_g[1]
                        -coords_B_4_b[0]*coords_B_4_b_g[1];

        
            (*virial)[2] += -coords_A_1_a[0]*coords_A_1_a_g[2]
                        -coords_B_1_a[0]*coords_B_1_a_g[2]
                        -coords_B_2_a[0]*coords_B_2_a_g[2]
                        -coords_A_2_b[0]*coords_A_2_b_g[2]
                        -coords_B_3_b[0]*coords_B_3_b_g[2]
                        -coords_B_4_b[0]*coords_B_4_b_g[2];

        
            (*virial)[4] += -coords_A_1_a[1]*coords_A_1_a_g[1]
                        -coords_B_1_a[1]*coords_B_1_a_g[1]
                        -coords_B_2_a[1]*coords_B_2_a_g[1]
                        -coords_A_2_b[1]*coords_A_2_b_g[1]
                        -coords_B_3_b[1]*coords_B_3_b_g[1]
                        -coords_B_4_b[1]*coords_B_4_b_g[1];

        
            (*virial)[5] += -coords_A_1_a[1]*coords_A_1_a_g[2]
                        -coords_B_1_a[1]*coords_B_1_a_g[2]
                        -coords_B_2_a[1]*coords_B_2_a_g[2]
                        -coords_A_2_b[1]*coords_A_2_b_g[2]
                        -coords_B_3_b[1]*coords_B_3_b_g[2]
                        -coords_B_4_b[1]*coords_B_4_b_g[2];

        
            (*virial)[8] += -coords_A_1_a[2]*coords_A_1_a_g[2]
                        -coords_B_1_a[2]*coords_B_1_a_g[2]
                        -coords_B_2_a[2]*coords_B_2_a_g[2]
                        -coords_A_2_b[2]*coords_A_2_b_g[2]
                        -coords_B_3_b[2]*coords_B_3_b_g[2]
                        -coords_B_4_b[2]*coords_B_4_b_g[2];

            (*virial)[3] = (*virial)[1];
            (*virial)[6] = (*virial)[2];
            (*virial)[7] = (*virial)[5];


        }


    }
    return energies_sw;

}

//----------------------------------------------------------------------------//
} // namespace mbnrg_A1B2X2_A1B2X2_deg4

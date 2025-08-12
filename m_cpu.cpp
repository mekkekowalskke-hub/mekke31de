#include <vector>
#include <complex>
#include <algorithm>
#include <cmath>

void GRAD(int SD, int NEN, const double* X,
          const std::vector<std::vector<double>>& DPHI,
          std::vector<std::vector<double>>& FX);
void INV(int SD, const std::vector<std::vector<double>>& F,
         std::vector<std::vector<double>>& INVF, double& DETF);
void ADJ(int SD, const std::vector<std::vector<double>>& F,
         std::vector<std::vector<double>>& ADJF);
double DET(int SD, const std::vector<std::vector<double>>& F);
void MAT_TRANS(const std::vector<std::vector<double>>& A, int N, int M,
               std::vector<std::vector<double>>& AT);
void MATVEC_PRO(const std::vector<std::vector<double>>& A, int N, int M,
                const std::vector<double>& B, std::vector<double>& C);
void BASIS_IN(const int SD, const int OD, const int* DIM_CP, const int* DIM_KP,
              const int* DIM_OP, const int NEN, const int NUMEL,
              const int IND, const int IB, const int D123, const int* SPAN,
              const double* CWEI, double* DPHI_G, double* DERIV,
              double* TEMP, double* WP, const double* BSS_BASIS);

void INTEGRAL(
    const int SD, const int OD, const int D1, const int D2, const int D3, 
    const int NUMDR, const int NUMKP, const int NUMEL, const int HFKS, const int NEN,
    const int IND, const int LDGG, const int* NUMEIG, const double* IW, const double* JAC,
    const double* OCC, const std::complex<double>* PSI, const double* KWEI, const int SUM,
    const double* X, const double* DPHI_G, double* RHO, double* NUM0, double* CHECK,
    double* THECK, std::complex<double>* PHECK) {
    
    std::vector<std::vector<double>> DPHI
        (NEN, std::vector<double>(SD));
    std::vector<double> PHI
        (NEN);
    std::vector<std::vector<double>> GRADPHI
        (NEN, std::vector<double>(SD));
    std::vector<std::vector<double>> INVFXT 
        (SD, std::vector<double>(SD));
    std::vector<std::vector<double>> INVFX
        (SD, std::vector<double>(SD));
    std::vector<std::vector<std::complex<double>>> PSI0
        (SUM, std::vector<std::complex<double>>(NUMKP));
    std::vector<std::vector<std::vector<std::complex<double>>>> DPSI0
        (SUM, std::vector<std::vector<std::complex<double>>>(NUMKP, std::vector<std::complex<double>>(SD)));

    std::vector<double> Q1
        (SD);
    std::vector<double> Q2
        (SD);
    std::vector<double> RHO0
        (NUMDR, 0.0);

    for (int l = 0; l < HFKS; l++) NUM0[l] = 0.0;
    
    double DETFX = 1.0;
    double DETJ = 0.0;
    std::vector<std::vector<double>> 
        FX(SD, std::vector<double>(SD));

    int I, J, K, L, P;
    int K0, K1, K2, K3;
    int IB;
    int POS;

    POS = IND * D1 * D2 * D3 * NEN;
    IB = 0;

    for (K1 = 0; K1 < D1; K1++){
        for (K2 = 0; K2 < D2; K2++){
            for (K3 = 0; K3 < D3; K3++){
                IB = IB + 1;
                for (int i = 0; i < NEN; i++){
                    PHI[i] = DPHI_G[(POS + (IB - 1) * NEN + i) * (SD + 1) + 0];
                    for (int d = 0; d < SD; d++) {
                        DPHI[i][d] = DPHI_G[(POS + (IB - 1) * NEN + i) * (SD + 1) + (d + 1)];
                    }
                }
            
        
            for (int j = 0; j < NEN; j++){
                for (int d = 0; d < SD; d++){
                    DPHI[j][d] *= JAC[d];
                }
            }

            GRAD(SD, NEN, X, DPHI, FX);

            if (LDGG == 2) {
                INV(SD, FX, INVFX, DETFX);
                MAT_TRANS(INVFX, SD, SD, INVFXT);
                for (int j = 0; j < NEN; j++) {
                    for (int d = 0; d < SD; d++) 
                        Q1[d] = DPHI[j][d];
                    
                    MATVEC_PRO(INVFXT, SD, SD, Q1, Q2);
                    for (int d = 0; d < SD; d++) 
                        GRADPHI[j][d] = Q2[d];
                    
                }
            }

            for (int n = 0; n < SUM; n++) 
                for (int p = 0; p < NUMKP; p++) 
                    PSI0[n][p] = std::complex<double>(0.0, 0.0);
                
            

            for (int i = 0; i < NEN; i++) 
                for (int n = 0; n < SUM; n++) 
                    for (int p = 0; p < NUMKP; p++) 
                        PSI0[n][p] += PHI[i] * PSI[i * SUM * NUMKP + n * NUMKP + p];
                    
                
            

            if (LDGG == 2) {
                for (int n = 0; n < SUM; n++) 
                    for (int p = 0; p < NUMKP; p++) 
                       for (int d = 0; d < SD; d++) 
                            DPSI0[n][p][d] = std::complex<double>(0.0, 0.0);
                        
                    
                
                for (int i = 0; i < NEN; ++i) 
                    for (int k = 0; k < SD; ++k) 
                        for (int n = 0; n < SUM; ++n) 
                            for (int p = 0; p < NUMKP; ++p)
                                DPSI0[n][p][k] += GRADPHI[i][k] * PSI[i * SUM * NUMKP + n * NUMKP + p];                   
                
            }
            std::fill(RHO0.begin(), RHO0.end(), 0.0);
            for (int l = 0; l < HFKS; ++l) {
                int K0 = 0;
                if (l == 1) K0 = NUMEIG[0];
                for (int k = K0; k < K0 + NUMEIG[l]; k++)
                    for (int p = 0; p < NUMKP; p++) 
                        RHO0[l] += KWEI[p] * OCC[k * NUMKP + p] * std::norm(PSI0[k][p]);
                    
                

                if (LDGG == 2) {
                    for (int k = K0; k < K0 + NUMEIG[l]; k++) 
                        for (int p = 0; p < NUMKP; p++) 
                            for (int s = 0; s < SD; s++) 
                                RHO0[HFKS + l * SD + s] += 
                                KWEI[p] * OCC[k * NUMKP + p] * 2.0 * (std::real(PSI0[k][p]) * std::real(DPSI0[k][p][s]) + std::imag(PSI0[k][p]) * std::imag(DPSI0[k][p][s]));
                            
                        
                    }
                }
            for (int r = 0; r < NUMDR; r++) {
                RHO[r * (D1 * D2 * D3) + (K1 * D2 * D3) + (K2 * D3) + K3] = RHO0[r];
            }
            DETJ = IW[K1] * DETFX;
            if (OD == 2) DETJ *= IW[K2];
            DETJ *= IW[K3];

            for (int n = 0; n < HFKS; n++) {
                NUM0[n] += DETJ * RHO0[n];
            }

            for (int s = 0; s < SUM; s++) {
                for(int p = 0; p < NUMKP; p++) {
                    for(int d = 0; d <= SD; d++) { 
                        PHECK[((IND * D1 * D2 * D3 * SUM) + (IB - 1) * SUM + s) * NUMKP * (SD + 1) + p * (SD + 1) + d] = 
                        (d == 0 ? PSI0[s][p] : DPSI0[s][p][d - 1]);
                    }
                }
            }
            for(int i = 0; i < NEN; i++){
               CHECK[(POS + (IB-1)*NEN + i) * (SD + 1) + 0] = PHI[i];
                for(int d = 0; d < SD; d++){
                    CHECK[(POS + (IB - 1)* NEN + i) * (SD + 1) + (d + 1)] = DPHI[i][d];
                }
            }
            for(int d1 = 0; d1 < SD; d1++){
                for(int d2 = 0; d2 < SD; d2++){
                    THECK[(IND * D1 * D2 * D3 * SD + (IB - 1) * SD + d1) * SD + d2] = INVFX[d1][d2];
                }
            }







            }  
        }
    }
}


void GRAD(int SD, int NEN, const double* X,
          const std::vector<std::vector<double>>& DPHI,
          std::vector<std::vector<double>>& FX) {
    for (int j1 = 0; j1 < SD; ++j1) {
        for (int j2 = 0; j2 < SD; ++j2) {
            double sum = 0.0;
            for (int k = 0; k < NEN; ++k) {
                sum += X[k * SD + j1] * DPHI[k][j2];
            }
            FX[j1][j2] = sum;
        }
    }
}

void INV(int SD, const std::vector<std::vector<double>>& F,
         std::vector<std::vector<double>>& INVF, double& DETF) {
    std::vector<std::vector<double>> ADJF(SD, std::vector<double>(SD));
    ADJ(SD, F, ADJF);
    DETF = DET(SD, F);
    for (int j1 = 0; j1 < SD; ++j1) {
        for (int j2 = 0; j2 < SD; ++j2) {
            INVF[j1][j2] = ADJF[j2][j1] / DETF;
        }
    }
}

void ADJ(int SD, const std::vector<std::vector<double>>& F,
         std::vector<std::vector<double>>& ADJF) {
    if (SD == 1) {
        ADJF[0][0] = 1.0;
    } else if (SD == 2) {
        ADJF[0][0] =  F[1][1];
        ADJF[0][1] = -F[1][0];
        ADJF[1][0] = -F[0][1];
        ADJF[1][1] =  F[0][0];
    } else if (SD == 3) {
        ADJF[0][0] = F[1][1]*F[2][2]-F[2][1]*F[1][2];
        ADJF[0][1] = F[2][0]*F[1][2]-F[1][0]*F[2][2];
        ADJF[0][2] = F[1][0]*F[2][1]-F[2][0]*F[1][1];
        ADJF[1][0] = F[0][2]*F[2][1]-F[0][1]*F[2][2];
        ADJF[1][1] = F[0][0]*F[2][2]-F[0][2]*F[2][0];
        ADJF[1][2] = F[0][1]*F[2][0]-F[0][0]*F[2][1];
        ADJF[2][0] = F[0][1]*F[1][2]-F[0][2]*F[1][1];
        ADJF[2][1] = F[0][2]*F[1][0]-F[1][2]*F[0][0];
        ADJF[2][2] = F[0][0]*F[1][1]-F[1][0]*F[0][1];
    }
}

double DET(int SD, const std::vector<std::vector<double>>& F) {
    if (SD == 1) {
        return F[0][0];
    } else if (SD == 2) {
        return F[0][0]*F[1][1] - F[0][1]*F[1][0];
    } else if (SD == 3) {
        return F[0][0]*( F[1][1]*F[2][2]-F[2][1]*F[1][2] ) +
               F[0][1]*( F[2][0]*F[1][2]-F[1][0]*F[2][2] ) +
               F[0][2]*( F[1][0]*F[2][1]-F[2][0]*F[1][1] );
    }
    return 0.0;
}

void MAT_TRANS(const std::vector<std::vector<double>>& A, int N, int M,
               std::vector<std::vector<double>>& AT) {
    for (int j1 = 0; j1 < N; ++j1) {
        for (int j2 = 0; j2 < M; ++j2) {
            AT[j2][j1] = A[j1][j2];
        }
    }
}

void MATVEC_PRO(const std::vector<std::vector<double>>& A, int N, int M,
                const std::vector<double>& B, std::vector<double>& C) {
    std::fill(C.begin(), C.end(), 0.0);
    for (int j = 0; j < N; ++j) {
        for (int k = 0; k < M; ++k) {
            C[j] += A[j][k] * B[k];
        }
    }
}



void BASIS(
    const int SD, const int OD, const int D1, const int D2, const int D3,
    const int* DIM_CP, const int* DIM_KP, const int* DIM_OP, const int NEN, 
    const int NUMEL, const int IND, const int D123, const int* SPAN,const double* CWEI,
    double* DPHI_G, double* DERIV, double* TEMP, double* WP, const double* BSS_BASIS
)
{
    int IB = 0;
    for (int K1 = 0; K1 < D1; K1++) {
        for (int K2 = 0; K2 < D2; K2++) {
            for (int K3 = 0; K3 < D3; K3++){
                IB++;
                BASIS_IN(
                    SD, OD, DIM_CP, DIM_KP, DIM_OP, NEN, NUMEL, IND, IB, D123, SPAN,
                    CWEI, DPHI_G, DERIV, TEMP, WP, BSS_BASIS
                );
            }
        }
    }
}



void BASIS_IN(
    const int SD, const int OD, const int* DIM_CP, const int* DIM_KP,
    const int* DIM_OP, const int NEN, const int NUMEL, const int IND,
    const int IB, const int D123, const int* SPAN, const double* CWEI,
    double* DPHI_G, double* DERIV, double* TEMP, double* WP,
    const double* BSS_BASIS)
{
    const int BSS_IND = IND * D123 + IB - 1;
    const int COLS = 2 * (DIM_OP[0] + 1) + 2 * (DIM_OP[1] + 1) + 2 * (DIM_OP[2] + 1);

    for (int no = 1; no <= SD; no++){
        DERIV[IND*SD + (no - 1)] = double((IND + 1) % SD) * double(no);
    }

    for (int j = 0; j <= SD; j++) {
        WP[IND * (SD + 1) + j] = 0.0;
    }

    for (int J1 = 0; J1 <= DIM_OP[0]; J1++) {
        const int K1 = SPAN[IND * SD] - DIM_OP[0] + J1;
        for (int J2 = 0; J2 <= DIM_OP[1]; J2++) {
            const int K2 = (OD == 2) ? (SPAN[IND * SD + 1] - DIM_OP[1] + J2) : 0;
            for (int J3 = 0; J3 <= DIM_OP[2]; J3++) {
                const int K3 = SPAN[IND * SD + 2] - DIM_OP[2] + J3; 

                for (int j = 0; j <= SD; ++j) TEMP[IND*(SD+1)+j] = 0.0;

                TEMP[IND * (SD + 1)] = BSS_BASIS[BSS_IND * COLS + J1];
                if (OD == 2){
                    TEMP[IND * (SD + 1)] *= BSS_BASIS[BSS_IND * COLS + (2 * (DIM_OP[0] +1 ) + J2)];
                }
                TEMP[IND * (SD + 1)] *= BSS_BASIS[BSS_IND * COLS + (2 * (DIM_OP[0]+1) + 2 * (DIM_OP[1]+1) + J3)];
                TEMP[IND * (SD + 1) + 1] = BSS_BASIS[BSS_IND * COLS + ((DIM_OP[0] + 1) + J1)];
                if (OD == 2) {
                    TEMP[IND * (SD + 1) + 1] *= BSS_BASIS[BSS_IND*COLS + (2 * (DIM_OP[0]+1) + J2)];
                    TEMP[IND * (SD + 1) + OD]  =
                    BSS_BASIS[BSS_IND*COLS + (J1)] * BSS_BASIS[BSS_IND * COLS + (2 * (DIM_OP[0]+1) + (DIM_OP[1]+1) + J2)];
                }
                for (int d = 1; d <= OD; d++) {
                    TEMP[IND * (SD+1) + d] *= BSS_BASIS[BSS_IND*COLS + (2 * (DIM_OP[0]+1) + 2 * (DIM_OP[1]+1) + J3)];
                }
                TEMP[IND*(SD+1) + SD] =
                BSS_BASIS[BSS_IND * COLS +(J1)] * BSS_BASIS[BSS_IND * COLS + (2 * (DIM_OP[0] + 1) + J2)] * BSS_BASIS[BSS_IND*COLS + (2*(DIM_OP[0]+1) + 2*(DIM_OP[1]+1) + (DIM_OP[2]+1) + J3)];

                for (int j = 0; j <= SD; j++) {
                    TEMP[IND * (SD + 1) + j] *=
                        CWEI[K1 * (DIM_CP[1] + 1) * (DIM_CP[2] + 1) +
                             K2 * (DIM_CP[2] + 1) + K3];
                    WP[IND * (SD + 1) + j] += TEMP[IND * (SD + 1) + j];
                }

                int no_index = (J1 * (DIM_OP[1] + 1) + J2) * (DIM_OP[2] + 1) + J3;
                DPHI_G[(((IND * D123) + (IB - 1)) * NEN + no_index) * (SD + 1) + 0] =
                    TEMP[IND * (SD + 1) + 0] / WP[IND * (SD + 1) + 0];
                for (int d = 1; d <= SD; d++) {
                    DPHI_G[(((IND * D123) + (IB - 1)) * NEN + no_index) * (SD + 1) + d] =
                        (TEMP[IND * (SD + 1) + d] -
                         WP[IND * (SD + 1) + d] *
                             DPHI_G[(((IND * D123) + (IB - 1)) * NEN + no_index) * (SD + 1) + 0]) /
                        WP[IND * (SD + 1) + 0];
                }
            }
        }
    }
    for (int no = 1; no <= SD; no++) {
        for (int t = 0; t < NEN; t++) {
            DPHI_G[((IND * D123 * NEN) + ((IB - 1) * NEN) + t) * (SD + 1) + no] *=
                DERIV[IND * SD + (no - 1)];
        }
    }
}

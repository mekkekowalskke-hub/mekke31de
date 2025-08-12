#include <algorithm>
#include <complex>
#include <iostream>
#include <numeric>
#include <vector>
#include "cpu.hpp"


const int SD = 3;
const int OD = 2;
const int NUMEL = 1000;
const int NUMIP = 4;

const int NUMDR = 4;
const int HFKS = 1;
const int NUMKP = 1;
const int NEN = int(pow((3 + 1), 3));

const int LDGG = 2;
const int LAGRANGE = 0;
const int N = 1;

int main(int argc, char* argv[]) {
    std::vector<int> NUMEIG(HFKS);
    int D1, D2, D3;
    int D123;
    int IB;
    int MOP;

    std::vector<int> DIM_CP(SD);
    std::vector<int> DIM_KP(SD);
    std::vector<int> DIM_OP(SD);
    std::vector<std::vector<int>> SPAN(NUMEL, std::vector<int>(SD));

    std::vector<double> NUM(HFKS);
    std::vector<double> NUM0(HFKS);
    std::vector<double> KWEI(NUMKP);
    std::vector<double> JAC(SD);

    std::vector<double> IW;

    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> RHO_S;

    std::vector<std::vector<std::vector<std::vector<double>>>> RHO;
    std::vector<std::vector<double>> OCC;
    std::vector<std::vector<std::vector<std::complex<double>>>> PSI;
    std::vector<std::vector<double>> X;
    std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>> DPHI_G;

    std::vector<std::vector<double>> DERIV;
    std::vector<std::vector<std::vector<double>>> CWEI;

    std::vector<std::vector<double>> TEMP;
    std::vector<std::vector<double>> WP;

    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> BSS_BASIS;
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> BAK_KNOT0;
    std::vector<double> BAK_KNOT;

    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> BSS_A;
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> BSS_NDU;
    std::vector<std::vector<std::vector<double>>> LEFT, RIGHT;

    std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>> CHECK;
    std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>> THECK;
    std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>>> PHECK;

    int IND;
    int I, J, K;
    int K1, K2, K3;

    D1 = NUMIP;
    D2 = std::max(1, NUMIP * (OD - 1));
    D3 = std::max(1, NUMIP * (SD - OD));

    D123 = D1 * D2 * D3;

    DIM_CP = {10, 10, 10};
    DIM_KP = {14, 14, 14};
    DIM_OP = {3, 3, 3};

    MOP = *std::max_element(DIM_OP.begin(), DIM_OP.end());

    NUMEIG[0] = 1 * 256 - 1;
    KWEI[0] = 2.0;
    int SUM = std::accumulate(NUMEIG.begin(), NUMEIG.end(), 0);

    std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>> CHECK
        (NUMEL, std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>(D1, std::vector<std::vector<std::vector<std::vector<double>>>>(D2, std::vector<std::vector<std::vector<double>>>(D3, std::vector<std::vector<double>>(NEN, std::vector<double>(SD + 1))))));
    std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>> THECK
        (NUMEL, std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>(D1, std::vector<std::vector<std::vector<std::vector<double>>>>(D2, std::vector<std::vector<std::vector<double>>>(D3, std::vector<std::vector<double>>(SD, std::vector<double>(SD))))));
    std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>>> PHECK
        (NUMEL, std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>>>(D1, std::vector<std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>>(D2, std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>(D3, std::vector<std::vector<std::vector<std::complex<double>>>>(SUM, std::vector<std::vector<std::complex<double>>>(NUMKP, std::vector<std::complex<double>>(SD + 1)))))));

    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>RHO_S
        (NUMEL, std::vector<std::vector<std::vector<std::vector<double>>>>(D1,std::vector<std::vector<std::vector<double>>>(D2,std::vector<std::vector<double>>(3, std::vector<double>(NUMDR)))));
    std::vector<std::vector<std::vector<std::vector<double>>>> RHO
        (NUMDR, std::vector<std::vector<std::vector<double>>>(D1,std::vector<std::vector<double>>(D2, std::vector<double>(D3))));
    std::vector<std::vector<double>> OCC
        (SUM, std::vector<double>(NUMKP));
    std::vector<std::vector<std::vector<std::complex<double>>>> PSI
        (NEN, std::vector<std::vector<std::complex<double>>>(SUM,std::vector<std::complex<double>>(NUMKP)));
    std::vector<std::vector<double>> X
        (NEN, std::vector<double>(SD));
    std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>> DPHI_G
        (NUMEL,std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>(D1,std::vector<std::vector<std::vector<std::vector<double>>>>(D2,std::vector<std::vector<std::vector<std::vector<double>>>>(D3,std::vector<std::vector<std::vector<double>>>(NEN,std::vector<std::vector<double>>(SD+1,std::vector<double>(SD+1)))))));

    std::vector<std::vector<double>> DERIV
        (NUMEL,std::vector<double>(SD));
    std::vector<std::vector<std::vector<double>>> CWEI
        (DIM_CP[0]+1,std::vector<std::vector<double>>(DIM_CP[1]+1,std::vector<double>(DIM_CP[2]+1)));

    std::vector<std::vector<double>> TEMP
        (NUMEL,std::vector<double>(SD+1));
    std::vector<std::vector<double>> WP
        (NUMEL,std::vector<double>(SD+1));

    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> BSS_BASIS
        (NUMEL, std::vector<std::vector<std::vector<std::vector<double>>>>(D1, std::vector<std::vector<std::vector<double>>>(D2, std::vector<std::vector<double>>(D3, std::vector<double>(2 * (DIM_OP[0] + 1) + 2 * (DIM_OP[1] + 1) + 2 * (DIM_OP[2] + 1))))));
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> BAK_KNOT0
        (NUMEL, std::vector<std::vector<std::vector<std::vector<double>>>>(D1, std::vector<std::vector<std::vector<double>>>(D2, std::vector<std::vector<double>>(D3, std::vector<double>(SD)))));
    std::vector<double> BAK_KNOT
        ((DIM_KP[0] + 1) + (DIM_KP[1] + 1) + (DIM_KP[2] + 1));
    
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> BSS_A
        (NUMEL, std::vector<std::vector<std::vector<std::vector<double>>>>(D1, std::vector<std::vector<std::vector<double>>>(D2, std::vector<std::vector<double>>(D3, std::vector<double>(2 * (MOP + 1))))));
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> BSS_NDU
        (NUMEL, std::vector<std::vector<std::vector<std::vector<double>>>>(D1, std::vector<std::vector<std::vector<double>>>(D2, std::vector<std::vector<double>>(D3, std::vector<double>((MOP + 1) * (MOP + 1))))));
    std::vector<std::vector<std::vector<double>>> LEFT
        (NUMEL, std::vector<std::vector<double>>(D123, std::vector<double>(MOP)));
    std::vector<std::vector<std::vector<double>>> RIGHT
        (NUMEL, std::vector<std::vector<double>>(D123, std::vector<double>(MOP)));



    for (int i = 0; i < NUMEL; i++) {
        if ((i + 1) % 4 == 1) {
            SPAN[i][0] = 10;
            SPAN[i][1] = 5;
            SPAN[i][2] = 10;
    }   else if ((i + 1) % 4 == 2) {
            SPAN[i][0]= 4;
            SPAN[i][1] = 5;
            SPAN[i][2] = 7;
    }   else if ((i + 1) % 4 == 3) {
            SPAN[i][0] = 1;
            SPAN[i][1] = 5;
            SPAN[i][2] = 3;
    }   else {
            SPAN[i][0] = 8;
            SPAN[i][1] = 5;
            SPAN[i][2] = 2;
    }
    }

    for (int i = 0; i <= DIM_CP[0]; i++)
        for (int j = 0; j <= DIM_CP[1]; j++)
            for (int k = 0; k <= DIM_CP[2]; k++)
                CWEI[i][j][k] = i * 0.001 + j * 0.01 + k * 0.001;
    
    for (int i = 0; i < NUMEL; i++) {
        for (int j = 0; j < D1; j++) {
            for (int k = 0; k < D2; k++) {
                for (int l = 0; l < D3; l++) {
                    for (int m = 0; m < SD; m++){
                        BAK_KNOT0[i][j][k][l][m] = double(((i + 1) % (D1 * D2 * D3))) / double(NUMEL);
                    }
                }
            }
        }
    }

    for (int i = 0; i < (DIM_KP[0] + 1) + (DIM_KP[1] + 1) + (DIM_KP[2] + 1); i++) {
        BAK_KNOT[i] = double(i + 1) / double(DIM_KP[0] + 1);
    }

    for (int i = 0; i < NEN; i++) {
        for (int j = 0; j < SUM; j++) {
            for (int k = 0; k < NUMKP; k++) {
                PSI[i][j][k] = std::complex<double>(double(j + 1)/double(SUM), double(k + 1)/double(NUMKP));
            }
        }
    }

    for (int i = 0; i < SUM; i++) {
        OCC[i][0] = double(i + 1) / double(SUM);
    }

    for (int i = 0; i < NEN; i++) {
        for (int j = 0; j < SD; j++) {
            X[i][j] = 0.001 * double(j + 1) / double(i + 1);
        }
    }

    IW.resize(4);
    IW = {1.0, 2.0, 3.0, 4.0};
    JAC.resize(3);
    JAC = {1.0, 2.0, 3.0};

    std::fill(NUM.begin(), NUM.end(), 0.0);

    for (int IND = 0; IND < NUMEL; IND++) {
    if (LAGRANGE == 0) {
        IB = 0;
        for (int K1 = 0; K1 < D1; K1++) {
            for (int K2 = 0; K2 < D2; K2++) {
                for (int K3 = 0; K3 < D3; K3++) {
                    IB++;
                    DERS_BASIS_FUNS(SD, OD, NUMEL, IND, IB, D123, N, MOP,
                        DIM_KP, SPAN, BAK_KNOT0, DIM_OP, BAK_KNOT, BSS_BASIS,
                        BSS_A, BSS_NDU, LEFT, RIGHT);
                    }
                }
            }
        }
    }

    for (int IND = 0; IND < NUMEL; IND++) {
        BASIS(SD, OD, D1, D2, D3, DIM_CP, DIM_KP, DIM_OP, NEN, NUMEL, IND, D123,
        SPAN, CWEI, DPHI_G, DERIV, TEMP, WP, BSS_BASIS);
    }

    for (int IND = 0; IND < NUMEL; IND++) {
        INTEGRAL(SD, OD, D1, D2, D3, NUMDR, NUMKP, NUMEL, HFKS, NEN, IND, LDGG,
        NUMEIG, IW, JAC, OCC, PSI, KWEI, X, DPHI_G, RHO, NUM0, CHECK, THECK, PHECK);

        for (int K1 = 0; K1 < D1; K1++) {
            for (int K2 = 0; K2 < D2; K2++) {
                for (int K3 = 0; K3 < D3; K3++) {
                    for (int m = 0; m < NUMDR; m++) {
                        RHO_S[IND][K1][K2][K3][m] = RHO[m][K1][K2][K3];
                    }
                }
            }
        }
    }
    for (int l = 0; l < HFKS; l++) {
        NUM[l] += NUM0[l];
    }

    for (int k = 0; k < D3; k++) {
        std::cout << "RHO: ";
        for (int m = 0; m < NUMDR; m++) {
            std::cout << RHO_S[99][D1 - 1][D2 - 1][k][m] << " ";
        }
    std::cout << "\n";
    }
    std::cout << "\n";
    std::cout << "NUM: ";
    for (int l = 0; l < HFKS; l++) {
        std::cout << NUM[l] << " ";
    }
    std::cout << "\n";

    return 0;
}


    




    
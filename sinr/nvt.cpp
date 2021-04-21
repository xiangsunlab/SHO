#include "nvt.h"
int main(){
    std::cout;
    cout.precision(16);
    normal_distribution<Real> normal_dist(0.0,1.0);
    INIT_HAMILTONIAN();gen.seed(153);
    for (int t = 0; t < 10000; t++){
        position[t] = 0;
        velocity[t] = 0;
    }
    for (int t = 0; t < sampling; t++){
        PrepOneTrajVV(t);
        PropOneTrajVV(t);
    }
    for (int t = 0; t < 10000; t++){
        cout<< position[t] << "  "<< velocity[t] << endl;
    }
    return 0;
}
void INIT_HAMILTONIAN(){

    R     .resize(DOFn,0);
    P     .resize(DOFn,0);
    V     .resize(DOFn,0);
    F     .resize(DOFn,0);
    mass  .resize(DOFn,0);
    ovmass.resize(DOFn,0);
    Req   .resize(DOFn,0);
    Omega .resize(DOFn,0);
    E     .resize(sampling * nsteps,0);
    U     .resize(sampling * nsteps,0);
    KE    .resize(sampling * nsteps,0);
    //position.resize(10000,0);  //0.01 as a grid and ranging from -50 to +50
    //velocity.resize(10000,0);  //0.01 as a grid and ranging from -50 to +50

    DT  = 0.002;
    DT2 = DT * 0.5;
    beta = 1.0;
    temp = 1.0 / beta;

    for (int j = 0; j < DOFn; j++) {
        mass  [j] = 2.0;
        ovmass[j] = 1.0 / mass[j];
    }
        for (int j = 0; j < DOFn; j++) {
        Req  [j] = 2.33;
        Omega[j] = 1.0;
    }

    return;
}
void UPDATE_F(){
    for (int j = 0; j < DOFn; j++){
        F[j] = - mass[j] * Omega[j] * Omega[j] * (R[j] - Req[j]);
        PEne = 0.5 * mass[j] * Omega[j] * Omega[j] * (R[j] - Req[j]) * (R[j] - Req[j]);
    }
    return;
}
void INIT_SAMPLING_Cl(){
    normal_distribution<Real> normal_dist(0.0,1.0);
    vector<Real> sigma_R;    //Standard derivation of each R[j] (STD = STandard Derivation)
    vector<Real> sigma_P;    //STD of each P[j];
    sigma_R.resize(DOFn,0);
    sigma_P.resize(DOFn,0);

    for (int j = 0; j < DOFn; j++) {
        sigma_R[j] = sqrt(1.0 / (beta * Omega[j] * Omega[j] * mass[j]) );
        sigma_P[j] = sqrt(temp * mass[j]);
    }
    for (int j = 0; j < DOFn; j++) {
        R[j] = normal_dist(gen) * sigma_R[j] + Req[j];
        P[j] = normal_dist(gen) * sigma_P[j];
        V[j] = P[j] * ovmass[j];
    }
    return;
}
void PropOneTrajVV(int t){
    UPDATE_F();
    for(int i = 0; i < nsteps; i++){
        MOVE_Nose_V(DT2);
        MOVE_Nose_R(DT2);
        UPDATE_F();
        MOVE_NoseOU(DT);
        MOVE_Nose_R(DT2);
        UPDATE_F();
        MOVE_Nose_V(DT2);
        Energy = KEne + PEne;
        if (abs(R[0]) >= 50 || abs(V[0]) >= 25){
        }
        else{
            int q = int((R[0] - Req[0] + 50 ) * 100);
            int p = int((V[0] * mass[0] + 50) * 100);
            //cout<<"R[0] = "<<R[0]<<endl;
            //cout<<"q = " <<q<<" p = "<<p<<endl;
            velocity[p] += 1;
            position[q] += 1;
            }
        }
    return;
}
void PrepOneTrajVV(int t){
    INIT_Nose();
    INIT_SAMPLING_Cl();
    UPDATE_F();
    for (int i = 0; i < nsteps; i++){
        MOVE_Nose(DT2);
        MOVE_Nose_V(DT2);
        MOVE_Nose_R(DT2);
        UPDATE_F();
        MOVE_NoseOU(DT);

        MOVE_Nose_R(DT2);
        UPDATE_F();
        MOVE_Nose_V(DT2);
        MOVE_Nose(DT2);

    }
    return;
}
void INIT_Nose() {
    int k;   // NHL counting parameters
    int j;   // DOFn counting parameters
    int indND;  //indND = k * DOFn + j;
    Q1.    resize(DOFn, 0);
    ovQ1.  resize(DOFn, 0);
    v1.    resize(NHL * DOFn, 0);
    Q2.    resize(DOFn, 0);
    ovQ2.  resize(DOFn, 0);
    v2.    resize(NHL * DOFn, 0);
    fric.  resize(NHL * DOFn, 0);
    wsy.  resize(NHL, 0);
    tsy.  resize(nsy, 0);    
    
    wsy[0] = 1.0 / (2.0 - pow (2.0, (1.0 / 3.0)));
    wsy[2] = 1.0 / (2.0 - pow (2.0, (1.0 / 3.0)));
    wsy[1] = 1.0 - wsy[0] - wsy[2];
    
    for (k = 0; k < NHL; k++){
        if (k == 0){
            for (j = 0; j < DOFn; j++){
                indND = k * DOFn + j;
                v1[indND] = 0.1;
                v2[indND] = 1;
                Q1[j] = temp * (DT * 100) *  (DT * 100);
                ovQ1[j] = 1.0 / Q1[j];
                Q2[j] = (DT * 100) * (DT * 100) * temp;
                ovQ2[j] = 1.0/ Q2[j];
                fric[j] = 0.1;
            }
        }
        else {
            for (j = 0; j < DOFn; j++){
                indND = k * DOFn + j;
                v1[indND] = 0.1;
                v2[indND] = 1;
                Q1[j] =  DT * 100 * temp *  (DT * 100);
                ovQ1[j] = 1.0 / Q1[j];
                Q2[j] = temp *  (DT * 100) *  (DT * 100);
                ovQ2[j] = 1.0/ Q2[j];
                fric[indND] = 0.1;
            }
        }
    }
    return;
}
void MOVE_Nose(Real dt){
    Real g(0);
    Real h(0);
    Real h_sum(0);
    Real dti(0);
    Real ev2(0);
    int indND;
    for (int j = 0; j < DOFn; j++){
        for (int ires = 0; ires < nres; ires++){
            for (int isy = 0; isy < nsy; isy++) {
                h_sum = 0;
                for (int k = 0; k < NHL; k++){
                    indND = k * DOFn + j;
                    g = (Q1[j] * v1[indND] * v1[indND] - temp ) * ovQ2[j]; 
                    v2[indND] += (0.5 * tsy[isy] * g);
                    ev2 = exp(- 2.0 * tsy[isy] * v2[indND]);
                    h_sum = h_sum + (Q1[j] * v1[indND] * v1[indND] * ev2);
                }
                h = sqrt((static_cast<Real>(NHL)) * temp / (V[j] * V[j] * mass[j] + h_sum) );
                for (int k = 0; k < NHL; k++){
                    indND = k * DOFn + j;
                    v1[indND] *= (h * exp(- tsy[isy] * v2[indND]));
                }
                for (int k = 0; k < NHL; k++){
                    indND = k * DOFn + j;
                    g = (Q1[j] * v1[indND] * v1[indND] - temp ) * ovQ2[j];
                    v2[indND] += (0.5 * tsy[isy] * g);
                }
            }
        }
    }
    return;
}
void MOVE_Nose_V(Real dt){
    Real arg(0);
    Real a(0);
    Real b(0);
    Real st(0);
    Real rb(0);
    Real sdot(0);
    Real g(0);
    int indND;
    KEne = 0;
    g = (static_cast<Real> (NHL)) * temp;
    for (int j = 0; j < DOFn; j++) {
        a = F[j] * V[j] / g;
        b = F[j] * F[j] * ovmass[j] / g;
        rb = sqrt(b);
        arg = (dt * 2.0) * rb * 0.5;
        if (arg > 0.00001){
            for (int k = 0; k < NHL; k++){
                indND =  k * DOFn + j;
                st = (1.0 / rb ) * sinh(arg) + ( a / b )* (cosh(arg) - 1.0);
                sdot = cosh(arg) + (a / rb) * sinh(arg);
            }
        }
        else {
            for (int k = 0; k < NHL; k++){
                indND =  k * DOFn + j;
                st = ((((b * a / 24.0) * DT2 + b / 6.0) * DT2 + 0.5 * a ) * DT2 + 1.0);
                sdot = (((b * a / 6.0) + b * 0.5) * DT2 + a) * DT2 + 1.0;
            }
        }
        V[j] = (V[j] + F[j] * ovmass[j] * st) / sdot;
        KEne += (0.5 * mass[j] * V[j] * V[j]);
        for (int k = 0; k < NHL; k++){
            indND =  k * DOFn + j;
            v1[indND] /= sdot;
        }
    }
    return;
}
void MOVE_Nose_R(Real dt){
    PEne = 0;
    for (int j = 0; j < DOFn; j++){
        R[j] += V[j] * DT2;
        PEne += 0.5 * mass[j] * Omega[j] * Omega[j] * (R[j] - Req[j]) * (R[j] - Req[j]);
    }
    return;
}
void MOVE_NoseOU(Real dt){
    normal_distribution<Real> normal_dist(0.0,1.0);
    Real rand_gauss;

    for (int j = 0; j < DOFn; j++){
        for (int k = 0; k < NHL; k++){
            int indND = k * DOFn + j;
            rand_gauss = normal_dist(gen);
            v2[indND] = exp( - fric[indND] * DT) * v2[indND] + sqrt((1.0 - exp( - 2.0 * fric[indND] * DT)) * temp * 2.0 * fric[indND] / ( Q2[j] )) * rand_gauss;
        }
    }
    return;
}


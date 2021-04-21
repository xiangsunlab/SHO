#include "nvt.h"
int main(){
    std::cout;
    cout.precision(16);
    normal_distribution<Real> normal_dist(0.0,1.0);
<<<<<<< HEAD
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
=======
    INIT_HAMILTONIAN();
    //For Check Function above 
    /* 
    for (int j = 0; j < DOFn; j++) {
        cout<<"mass["<<j<<"] = "<<mass[j]<<endl; 
    }
    for (int j = 0; j < DOFn; j++) {
        cout<<"ovmass["<<j<<"] = "<<ovmass[j]<<endl;
    }
    for (int j = 0; j < DOFn; j++) {
       cout<<"Omega["<<j<<"] = "<<Omega[j]<<endl;
    }
    for (int j = 0; j < DOFn; j++) {
       cout<<"Req["<<j<<"] = "<<Req[j]<<endl;
    }
    for (int j = 0; j < 2; j++){
        cout<< normal_dist(gen) <<endl;
    }
    */
    gen.seed(153);
    for (int t = 0; t < sampling; t++){
        /*for (int j = 0; j < 200; j++){
            cout<< normal_dist(gen) <<"  "<< normal_dist(gen) <<endl;
        }*/
        //INIT_SAMPLING_Cl();
    
        
        //for (int j = 0; j < DOFn; j++) {
        //    cout<<P[j]<<"   "<<R[j]<<endl;
        //}
        //for (int j = 0; j < DOFn; j++) {
        //    cout<<"P["<<j<<"] = "<<P[j]<<endl;
        PrepOneTrajVV(t); 
        PropOneTrajVV(t);
    }
    for (int t = 0; t < 10000; t++){
        //cout<< position[t] << "  "<< velocity[t] << endl; 
>>>>>>> origin/main
    }
    return 0;
}
void INIT_HAMILTONIAN(){
<<<<<<< HEAD
=======
    //normal_distribution<Real> normal_dist(0.0,1.0);
>>>>>>> origin/main

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
<<<<<<< HEAD
    //position.resize(10000,0);  //0.01 as a grid and ranging from -50 to +50
    //velocity.resize(10000,0);  //0.01 as a grid and ranging from -50 to +50

    DT  = 0.002;
=======
    position.resize(10000,0);  //0.01 as a grid and ranging from -50 to +50
    velocity.resize(10000,0);  //0.01 as a grid and ranging from -50 to +50
    
    DT  = 0.01;
>>>>>>> origin/main
    DT2 = DT * 0.5;
    beta = 1.0;
    temp = 1.0 / beta;

<<<<<<< HEAD
=======
    //Initial mass for each particle
>>>>>>> origin/main
    for (int j = 0; j < DOFn; j++) {
        mass  [j] = 2.0;
        ovmass[j] = 1.0 / mass[j];
    }
<<<<<<< HEAD
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
=======
    //Here we have a shifted SHO and Req is the equilibrated position
    //Omega here serves as a set of spring constants
    //For SB model we have a certain relationship between Req and Omega
    for (int j = 0; j < DOFn; j++) {
        Req  [j] = 2.33;
        Omega[j] = 1.0;
    }
    
    return;    
}
void MOVE_VV_V(Real dt){
    KEne = 0;
    for (int j = 0; j < DOFn; j++){
        P[j] += F[j] * dt;
        KEne += 0.5 * P[j] * P[j] * ovmass[j];
    }
    return;
}
void MOVE_LF_V(Real dt){
    KEne = 0;
    for (int j = 0; j < DOFn; j++){
        P[j] += F[j] * dt;
        KEne += 0.5 * P[j] * P[j] * ovmass[j];
    }
    return;
}
void UPDATE_F(){
    //PEne = 0;
    for (int j = 0; j < DOFn; j++){
        //PEne += 0.5 * mass[j] * Omega[j] * Omega[j] * (R[j] - Req[j]) * (R[j] - Req[j]);
        F[j] = - mass[j] * Omega[j] * Omega[j] * (R[j] - Req[j]);
    }
    return;
}

void MOVE_VV_R(Real dt){
    PEne = 0;
    for (int j = 0; j < DOFn; j++){
        V[j] =  P[j] * ovmass[j];
        R[j] += V[j] * dt;
        PEne += 0.5 * mass[j] * Omega[j] * Omega[j] * (R[j] - Req[j]) * (R[j] - Req[j]);
    }
    return;
}
void MOVE_LF_R(Real dt){
    PEne = 0;
    for (int j = 0; j < DOFn; j++){
        V[j] =  P[j] * ovmass[j];
        R[j] += V[j] * dt;
        PEne += 0.5 * mass[j] * Omega[j] * Omega[j] * (R[j] - Req[j]) * (R[j] - Req[j]);
    }
    return;
}
void INIT_SAMPLING_Wigner(){
    //C++ 11 Gaussian Distribution Number Generator
    normal_distribution<Real> normal_dist(0.0,1.0);
    //[1] Initialise sampling distribution width (sigma) for the Nuke DOF 
    vector<Real> sigma_R;    //Standard derivation of each R[j] (STD = STandard Derivation)
    vector<Real> sigma_P;    //STD of each P[j];

    sigma_R.resize(DOFn,0);
    sigma_P.resize(DOFn,0);
    
    for (int j = 0; j < DOFn; j++) {
        sigma_R[j] = sqrt(hbar/ (2.0 * Omega[j] * tanh(0.5 * beta * hbar * Omega[j])) );
        sigma_P[j] = sigma_R[j] * Omega[j];
    }
    //code for examine
    /*
    for (int j = 0; j < DOFn; j++) {
        cout<<"sigma_R["<<j<<"] = "<<sigma_R[j]<<endl;
    }
    for (int j = 0; j < DOFn; j++) {
        cout<<"sigma_P["<<j<<"] = "<<sigma_P[j]<<endl;
    }
    */
    //[2] generate gaussian random number distribution 
    for (int j = 0; j < DOFn; j++) {
        R[j] = normal_dist(gen) * sigma_R[j] + Req[j];
	P[j] = normal_dist(gen) * sigma_P[j];
    }
    //Code for test
    /*
    for (int j = 0; j < DOFn; j++) {
        cout<<"R["<<j<<"] = "<<R[j]<<endl;
    }
    for (int j = 0; j < DOFn; j++) {
        cout<<"P["<<j<<"] = "<<P[j]<<endl;
    }
    */
    return;
}
>>>>>>> origin/main
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
<<<<<<< HEAD
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
=======
    //Test Code
    /*
    for (int j = 0; j < DOFn; j++) {
        cout<<"sigma_R["<<j<<"] = "<<sigma_R[j]<<endl;
    }
    for (int j = 0; j < DOFn; j++) {
        cout<<"sigma_P["<<j<<"] = "<<sigma_P[j]<<endl;
    }*/
    
    for (int j = 0; j < DOFn; j++) {
        R[j] = normal_dist(gen) * sigma_R[j] + Req[j];
        P[j] = normal_dist(gen) * sigma_P[j];
    }
    //Test Code
    /*
    for (int j = 0; j < DOFn; j++) {
        cout<<"R["<<j<<"] = "<<R[j]<<endl;
    }
    for (int j = 0; j < DOFn; j++) {
        cout<<"P["<<j<<"] = "<<P[j]<<endl;
    }*/

    return;
}
void PropOneTrajVV(int t){
    //INIT_HAMILTONIAN();
    //INIT_Langevin();  //Langevin dynamics propagation
    UPDATE_F();
    
    //Check Code
    //cout<<"init force updated!"<<endl;
    /*for (int j = 0; j < DOFn; j++) {
        //cout<<"R["<<j<<"] = "
        cout<<R[j]<<"  "<<P[j]<<endl;
    }
    for (int j = 0; j < DOFn; j++) {
        cout<<"P["<<j<<"] = "<<P[j]<<endl;
    }
    for (int j = 0; j < DOFn; j++) {
        cout<<"F["<<j<<"] = "<<F[j]<<endl;
    }*/
    //cout<<sigma[0]<<"<-sigma  fric->"<<fric[0]<<endl;

    for (int i = 0; i < nsteps; i++){
        int ind = i + nsteps * t;
        //cout<<"updated here1"<<endl;
        MOVE_Nose(DT2);
        //cout<<"updated here1"<<endl;
        MOVE_Nose_V(DT2);
        //cout<<"updated here2"<<endl;
        /*for (int j = 0; j < DOFn; j++) {
            cout<<"R["<<j<<"] = "<<R[j]<<endl;
        }
        for (int j = 0; j < DOFn; j++) {
            cout<<"P["<<j<<"] = "<<P[j]<<endl;
        }*/
        
        MOVE_Nose_R(DT2);
        //cout<<"updated here3"<<endl;
        /*
        for (int j = 0; j < DOFn; j++) {
            cout<<"R["<<j<<"] = "<<R[j]<<endl;
        }*/
        UPDATE_F(); 
        MOVE_NoseOU(DT);
        
        //MOVE_OU(DT);
        /*
        for (int j = 0; j < DOFn; j++) {
            cout<<"F["<<j<<"] = "<<F[j]<<endl;
        }*/
        
        MOVE_Nose_R(DT2);
        UPDATE_F();
        MOVE_Nose_V(DT2);
        
        //for (int j = 0; j < DOFn; j++) {
        //    cout<<R[j]<<"  "<<P[j]<<endl;}
	/*
        for (int j = 0; j < DOFn; j++) {
            cout<<"P["<<j<<"] = "<<P[j]<<endl;
        }
        cout<<"---------C'est la "<<i<<"th step--------------------"<<endl;*/
        MOVE_Nose(DT2);
        Energy = KEne + PEne;
	//cout<<ind<<endl;
        E[ind] = Energy;
        U[ind] = PEne;
        KE[ind]= KEne;
        
        int q = int((R[0] - Req[0] + 50 ) * 100 );
        int p = int((P[0] + 50) * 100);
        velocity[p] += 1;
        position[q] += 1;
        //cout<<velocity[p]<<endl;
        //cout<<Energy<<"  "<<PEne<<"  "<<KEne<<endl;
        
    }
>>>>>>> origin/main
    return;
}
void PrepOneTrajVV(int t){
    INIT_Nose();
    INIT_SAMPLING_Cl();
    UPDATE_F();
<<<<<<< HEAD
    for (int i = 0; i < nsteps; i++){
        MOVE_Nose(DT2);
        MOVE_Nose_V(DT2);
        MOVE_Nose_R(DT2);
        UPDATE_F();
        MOVE_NoseOU(DT);

=======

    for (int i = 0; i < nsteps; i++){
        int ind = i + nsteps * t;
        MOVE_Nose(DT2);
        MOVE_Nose_V(DT2);
        MOVE_Nose_R(DT2);


        UPDATE_F();
        MOVE_NoseOU(DT);
>>>>>>> origin/main
        MOVE_Nose_R(DT2);
        UPDATE_F();
        MOVE_Nose_V(DT2);
        MOVE_Nose(DT2);

    }
    return;
}
<<<<<<< HEAD
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
    
=======
void MOVE_OU(Real dt){
    normal_distribution<Real> normal_dist(0.0,1.0);
    for (int j = 0; j < DOFn; j++){
        F[j] = - fric[j] * P[j] + sigma[j] * normal_dist(gen);
        P[j] += F[j] * dt; 
    }
    return;
}
void INIT_Langevin(){
    fric    .resize(DOFn, 0);
    Q       .resize(DOFn, 0);
    sigma   .resize(DOFn, 0);
    for (int j = 0; j < DOFn; j++ ){
        fric [j] = 0.01;
        Q[j]     = 0;
        sigma[j] = sqrt(2.0 * mass[j] * fric[j] * temp / DT);
    }
    return;
}

void INIT_NHC(){
    int k = 0;
    Q.    resize(NHL, 0);
    ovQ.  resize(NHL, 0);
    v.    resize(NHL, 0);
    G.    resize(NHL, 0);
    wsy.  resize(NHL, 0);
    tsy.  resize(nsy, 0);
    
    k = 0;
    Q[k] = (static_cast<Real>(DOFn) * temp) * (DT * 100.0) * (DT * 100.0); //Here \tau = 100 * dt 
    v[k] = 1.0;
    G[k] = 0.0;
    //cout<<"Q["<<k<<"] = "<<Q[k]<<endl;
    //cout<<"v["<<k<<"] = "<<v[k]<<endl;
    //cout<<"G["<<k<<"] = "<<G[k]<<endl;
    for (k = 1; k < NHL; k++){
        Q[k] = temp * (DT * 100.0) * (DT * 100.0) ;
        v[k] = 0.1;
        G[k] = 0.0;
        //cout<<"Q["<<k<<"] = "<<Q[k]<<endl;
        //cout<<"v["<<k<<"] = "<<v[k]<<endl;
        //cout<<"G["<<k<<"] = "<<G[k]<<endl;
    }
    //cout<<Q[2]<<endl;
    for (int k = 0; k < NHL; k++){
        ovQ[k] = 1.0 / Q[k];
        //cout<<"Q["<<k<<"] = "<<Q[k]<<" 1/Q = "<< ovQ[k]<<endl;
    }
    //cout<<"calculated HERE"<<endl;
    wsy[0] = 1.0 / (2.0 - pow (2.0, (1.0 / 3.0)));
    wsy[2] = 1.0 / (2.0 - pow (2.0, (1.0 / 3.0)));
    wsy[1] = 1.0 - wsy[0] - wsy[2];
    for (int isy = 0; isy < nsy; isy ++){
        tsy[isy] = DT2 * wsy[isy] / (static_cast<Real>(nres));
        //cout<<"wsy["<<isy<<"] = "<<wsy[isy]<<endl;
        //cout<<"tsy["<<isy<<"] = "<<tsy[isy]<<endl;
    }
    //cout<<"init done! HOUA"<< endl;

    return;
}

void INIT_Nose() { 
    int k;   // NHL counting parameters
    int j;   // DOFn counting parameters
    int indND;
    Q1.    resize(NHL * DOFn, 0);
    ovQ1.  resize(NHL * DOFn, 0);
    v1.    resize(NHL * DOFn, 0);
    Q2.    resize(NHL * DOFn, 0);
    ovQ2.  resize(NHL * DOFn, 0);
    v2.    resize(NHL * DOFn, 0);
    fric.  resize(NHL * DOFn, 0);
    wsy.  resize(NHL, 0);
    tsy.  resize(nsy, 0);
    
    //4th-order Suzuki-Yolanda Factorisation weight
>>>>>>> origin/main
    wsy[0] = 1.0 / (2.0 - pow (2.0, (1.0 / 3.0)));
    wsy[2] = 1.0 / (2.0 - pow (2.0, (1.0 / 3.0)));
    wsy[1] = 1.0 - wsy[0] - wsy[2];
    
<<<<<<< HEAD
=======
    for (int isy = 0; isy < nsy; isy++){
        tsy[isy] = DT2 * wsy[isy] / (static_cast<Real>(nres));
    }
    
>>>>>>> origin/main
    for (k = 0; k < NHL; k++){
        if (k == 0){
            for (j = 0; j < DOFn; j++){
                indND = k * DOFn + j;
<<<<<<< HEAD
                v1[indND] = 0.1;
                v2[indND] = 1;
                Q1[j] = temp * (DT * 100) *  (DT * 100);
                ovQ1[j] = 1.0 / Q1[j];
                Q2[j] = (DT * 100) * (DT * 100) * temp;
                ovQ2[j] = 1.0/ Q2[j];
                fric[j] = 0.1;
=======
                Q1[indND] = temp;
                ovQ1[indND] = 1.0 / Q1[indND];
                Q2[indND] = temp;
                ovQ2[indND] = 1.0/ Q2[indND];
                fric[indND] = 0.1;
>>>>>>> origin/main
            }
        }
        else {
            for (j = 0; j < DOFn; j++){
                indND = k * DOFn + j;
<<<<<<< HEAD
                v1[indND] = 0.1;
                v2[indND] = 1;
                Q1[j] =  DT * 100 * temp *  (DT * 100);
                ovQ1[j] = 1.0 / Q1[j];
                Q2[j] = temp *  (DT * 100) *  (DT * 100);
                ovQ2[j] = 1.0/ Q2[j];
=======
                Q1[indND] = (static_cast<Real>(DOFn)) * temp;
                ovQ1[indND] = 1.0 / Q1[indND];
                Q2[indND] = temp;
                ovQ2[indND] = 1.0/ Q2[indND];
>>>>>>> origin/main
                fric[indND] = 0.1;
            }
        }
    }
<<<<<<< HEAD
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
=======
    /*
    for (k = 0; k < NHL; k++) {
        cout<<k<<"  NHL"<<endl;
        for (j = 0; j < DOFn; j++){
            indND = k * DOFn + j;
            cout<<indND<<" = "<<k <<" * DOFn + "<< j <<endl;
            cout<<" Q1[indND] = "<<Q1[indND] ;
            cout<<" ovQ1[indND] = 1.0 / Q1[indND] = " <<ovQ1[indND];
            cout<<" Q2[indND] = temp  = "<<Q2[indND];
            cout<<" ovQ2[indND] = 1.0/ Q2[indND] = "<< ovQ2[indND] ;
            cout<<" fric[indND] = 0.1 or "<<fric[indND] <<endl;
        }
    }*/
    return;
}
void MOVE_Nose(Real dt) {
    Real g(0); 
    Real h(0); 
    Real h_sum(0); 
    Real dti(0); 
    Real ev2(0);
    int j, k, ires, indND, isy;
    for (j = 0; j < DOFn; j++){
        for (isy = 0; isy < nsy; isy++) {
            for (ires = 0; ires < nres; ires++){
                h_sum = 0;
                for (k = 0; k < NHL; k++){
                    indND = k * DOFn + j;
                    g = (Q1[indND] * v1[indND] * v1[indND] - temp ) * ovQ2[indND];
                    v2[indND] += tsy[isy] * g;
                    ev2 = exp(- 2.0 * tsy[isy] * v2[indND]);
                    h_sum += v1[indND] * v1[indND] * ev2 * ev2;
                }
                h_sum *=  (static_cast<Real>(DOFn)) /  (static_cast<Real>(DOFn + 1));
                h = sqrt(  (static_cast<Real>(DOFn)) * temp / (V[j] * V[j] * mass[j] + h_sum) );
                
                for (k = 0; k < NHL; k++){
                    indND = k * DOFn + j;
                    v1[indND] *= (h * exp( - 2.0 * tsy[isy] * v2[indND] ) );
                }
                for (k = 0; k < NHL; k++){
                    indND = k * DOFn + j;
                    g = (Q1[indND] * v1[indND] * v1[indND] - temp ) * ovQ2[indND];
                    v2[indND] += tsy[isy] * g;
>>>>>>> origin/main
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
<<<<<<< HEAD
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
=======
    Real g(0); 
    g = (static_cast<Real> (NHL)) * temp;
    int j, k, indND;
    for (j = 0; j < DOFn; j++) {
        a = F[j] * V[j] / g;
        b = F[j] * F[j] * ovmass[j] / g;
        rb = sqrt(b);
        arg = (dt * 2) * rb * 0.5;
        if (arg > 0.00001){
            for (k = 0; k < NHL; k++){
                indND =  k * DOFn + j;
                st = (1.0 / rb ) * sinh(arg) + (a / b )* (cosh(arg) - 1.0);
>>>>>>> origin/main
                sdot = cosh(arg) + (a / rb) * sinh(arg);
            }
        }
        else {
<<<<<<< HEAD
            for (int k = 0; k < NHL; k++){
=======
            for (k = 0; k < NHL; k++){
>>>>>>> origin/main
                indND =  k * DOFn + j;
                st = ((((b * a / 24.0) * DT2 + b / 6.0) * DT2 + 0.5 * a ) * DT2 + 1.0);
                sdot = (((b * a / 6.0) + b * 0.5) * DT2 + a) * DT2 + 1.0;
            }
        }
<<<<<<< HEAD
        V[j] = (V[j] + F[j] * ovmass[j] * st) / sdot;
        KEne += (0.5 * mass[j] * V[j] * V[j]);
        for (int k = 0; k < NHL; k++){
=======
        V[j] = (V[j] + F[j] * ovmass[j] * st) / sdot; 
        for (k = 0; k < NHL; k++){
>>>>>>> origin/main
            indND =  k * DOFn + j;
            v1[indND] /= sdot;
        }
    }
    return;
}
void MOVE_Nose_R(Real dt){
    PEne = 0;
    for (int j = 0; j < DOFn; j++){
<<<<<<< HEAD
        R[j] += V[j] * DT2;
=======
        R[j] += V[j] * dt;
>>>>>>> origin/main
        PEne += 0.5 * mass[j] * Omega[j] * Omega[j] * (R[j] - Req[j]) * (R[j] - Req[j]);
    }
    return;
}
void MOVE_NoseOU(Real dt){
    normal_distribution<Real> normal_dist(0.0,1.0);
    Real rand_gauss;
<<<<<<< HEAD

    for (int j = 0; j < DOFn; j++){
        for (int k = 0; k < NHL; k++){
            int indND = k * DOFn + j;
            rand_gauss = normal_dist(gen);
            v2[indND] = exp( - fric[indND] * DT) * v2[indND] + sqrt((1.0 - exp( - 2.0 * fric[indND] * DT)) * temp * 2.0 * fric[indND] / ( Q2[j] )) * rand_gauss;
        }
    }
    return;
}

=======
    for (int j = 0; j < DOFn; j++){
        
        for (int k = 0; k < NHL; k++){
            int indND = k * DOFn + j;
            rand_gauss = normal_dist(gen);
            v2[indND] = exp( - fric[indND] * dt ) * v2[indND] + sqrt((1.0 - exp( - 2.0 * fric[indND] * dt)) / ( 0.5 * fric[indND])) * rand_gauss; 
        }
    }
    return;
}
void MOVE_NHC(Real dt){
    int ires;
    int k;  //Counting index for the nhc (0 to NHL-1)
    int kk; //Counting index for the nhc (NHL - 1 to 0)
    int kind;//Counting index for the nhc-1 
    int isy;
    //cout<<"hEre Q"<<endl;
    for (ires = 0; ires < nres; ires++) {
        for (isy = 0; isy < nsy; isy++){
            //cout<<"isy["<<isy<<"] = "<<isy<<" ires = "<<ires<<endl;
            //[1]Update the thermostat force G_j = v_{j - 1} ** 2 / Q - d_j k_B T
            k = 0;
            //cout<<"P[0] = "<<P[0]<<endl;
            G[k] = 0.0;
            for (int j = 0; j < DOFn; j++){
                //cout<< P[j] * P[j] * ovmass[j] - temp<<endl;
                //cout<<"G["<<k<<"] = "<<G[k]<<endl;
                G[k] = G[k] + P[j] * P[j] * ovmass[j] - temp;
                //cout<<"G["<<k<<"] = "<<G[k]<<endl;
            }
            for (k = 0; k < (NHL - 1); k++){
                kk = k + 1;
                //cout<<"v["<<k<<"] = "<<v[k]<<endl; 
                G[kk] = v[k] * v[k] * ovQ[k] - temp;
                //cout<<"G["<<kk<<"] = "<<G[kk]<<endl;
                }
            k = (NHL - 1);
            //cout<<"v[3] = "<<v[k]<<endl; 

            //[2]Update imaginary momenta Last ones
            k = 0;
            kk = NHL - k - 1;
            v[kk] += 0.5 * G[kk] * tsy[isy];
            
            //[3]Update imaginary momenta First (NHL-1)
            for (k = 1; k < NHL; k++){
                kk = NHL - k - 1;
                kind = kk + 1;
                v[kk] *= exp(- 0.25 * tsy[isy] * v[kind] * ovQ[kind]);
                v[kk] += 0.5 * tsy[isy] * G[kk];
                v[kk] *= exp(- 0.25 * tsy[isy] * v[kind] * ovQ[kind]);
            }
            /*
            for (k = 0; k < NHL; k++){
                cout<<"v["<<k<<"] = "<<v[k]<<endl;
            } 
            */
            
            //[4] Coupling with phyical momenta
            for (int j = 0; j < DOFn; j++){
                //cout<<"P["<<j<<"]_b4 = "<<P[j]<<endl;
                P[j] *= exp(- tsy[isy] * v[0] * ovQ[0]);
                //cout<<"P["<<j<<"]  = "<<P[j]<<endl;
            }
            
            //[5]Update the thermostat force G_j = v_{j - 1} ** 2 / Q - d_j k_B T
            k = 0;
            G[k] = 0;
            for (int j = 0; j < DOFn; j++){
                G[k] +=  P[j] * P[j] * ovmass[j] -  temp;
            }
            //cout<<"G["<<k<<"] = "<<G[k]<<endl;
            for (k = 0; k < (NHL - 1); k++){
                kk = k + 1;
                G[kk] = v[k] * v[k] * ovQ[k] - 1.0 * temp;
                //cout<<"G["<<kk<<"] = "<<G[kk]<<endl;
            }
            /*
            for (k = 0; k < NHL; k++) { 
                cout<<"G["<<k<<"] = "<<G[k]<<endl;
            }*/
           
            //[6]Update the last imaginary momenta
            k = 0;
            kk = NHL - k - 1;
            v[kk] += 0.5 * G[kk] * tsy[isy];
            
            //[7]Update imaginary momenta First (NHL-1)
            for (k = 1; k < NHL; k++){
                kk = NHL - k - 1;
                kind = kk + 1;
                v[kk] *= exp(- 0.25 * tsy[isy] * v[kind] * ovQ[kind]);
                v[kk] += 0.5 * tsy[isy] * G[kk];
                v[kk] *= exp(- 0.25 * tsy[isy] * v[kind] * ovQ[kind]);
            }

            //cout<<v[kk]<<endl;
            /*
            for (k = 0; k < NHL; k++) {
                cout<<"v["<<k<<"] = "<<v[k]<<endl;
            }*/
           // cout<< "This is number "<< isy << " SY Factorisation Section"<<endl;
        }
        //cout<<"Das ist number "<<ires<<" RESPA Faktoris Sektion"<<endl;
    }
    //cout<<"One NVT"<<endl;
    return;
}














































>>>>>>> origin/main

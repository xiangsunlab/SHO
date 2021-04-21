#include "sinr.h"
int main(){
    std::cout;
    cout.precision(16);
    normal_distribution<Real> normal_dist(0.0,1.0);
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
    
    //cout<<"Propagated to "<<t<<" samping"<<endl;
    }
    for (int tt = 0; tt < 10000; tt++){
        cout<< position[tt] << "  "<< velocity[tt] << endl; 
    }
    return 0;
}
void INIT_HAMILTONIAN(){
    //normal_distribution<Real> normal_dist(0.0,1.0);

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
    position.resize(10000,0);  //0.01 as a grid and ranging from -50 to +50
    velocity.resize(10000,0);  //0.01 as a grid and ranging from -50 to +50
    
    DT  = 0.01;
    DT2 = DT * 0.5;
    beta = 1.0;
    temp = 1.0 / beta;

    //Initial mass for each particle
    for (int j = 0; j < DOFn; j++) {
        mass  [j] = 2.0;
        ovmass[j] = 1.0 / mass[j];
    }
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
        PEne = 0.5 * mass[j] * Omega[j] * Omega[j] * (R[j] - Req[j]) * (R[j] - Req[j]);
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
        V[j] = P[j] * ovmass[j];
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

    for (int step = 0; step < nsteps; step++){
        //int ind = step + nsteps * t;
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
        
        /*for (int j = 0; j < DOFn; j++) {
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
        //    cout<<R[j]<<"  "<<V[j]<<endl;}
	/*
        for (int j = 0; j < DOFn; j++) {
            cout<<"P["<<j<<"] = "<<P[j]<<endl;
        }
        cout<<"---------C'est la "<<i<<"th step--------------------"<<endl;*/
        MOVE_Nose(DT2);
        Energy = KEne + PEne;
	//cout<<ind<<endl;
        //E[ind] = Energy;
        //U[ind] = PEne;
        //KE[ind]= KEne;
  
        if (abs(R[0]) >= 50 || abs(V[0]) >= 50){
            
            //cout<<o<<endl;
        }
        else{
            int q = int((R[0] - Req[0] + 50 ) * 100 );
            int p = int((V[0] * mass[0] + 50) * 100);
            velocity[p] += 1;
            position[q] += 1;
        }
          
        //cout<<velocity[p]<<endl;
        //cout<<Energy<<"  "<<PEne<<"  "<<KEne<<endl;
        
    }
    return;
}
void PrepOneTrajVV(int t){
    INIT_Nose();
    INIT_SAMPLING_Cl();
    UPDATE_F();
    //cout<<"R[0] = "<<R[0]<<" V[0] = "<<V[0]<<endl;
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
    Q1.    resize(DOFn, 0);
    ovQ1.  resize(DOFn, 0);
    v1.    resize(NHL * DOFn, 0);
    Q2.    resize(DOFn, 0);
    ovQ2.  resize(DOFn, 0);
    v2.    resize(NHL * DOFn, 0);
    fric.  resize(NHL * DOFn, 0);
    wsy.  resize(NHL, 0);
    tsy.  resize(nsy, 0);
    
    //4th-order Suzuki-Yolanda Factorisation weight
    wsy[0] = 1.0 / (2.0 - pow (2.0, (1.0 / 3.0)));
    wsy[2] = 1.0 / (2.0 - pow (2.0, (1.0 / 3.0)));
    wsy[1] = 1.0 - wsy[0] - wsy[2];
    
    for (int isy = 0; isy < nsy; isy++){
        tsy[isy] = DT2 * wsy[isy] / (static_cast<Real>(nres));
    }
    
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
    /*cout<<"Inin Nose "<<endl;
    for (k = 0; k < NHL; k++) {
        cout<<k<<"  NHL"<<endl;
        for (j = 0; j < DOFn; j++){
            indND = k * DOFn + j;
            cout<<indND<<" = "<<k <<" * DOFn + "<< j <<endl;
            cout<<" Q1[j] = "<<Q1[j] ;
            cout<<" ovQ1[j] = 1.0 / Q1[j] = " <<ovQ1[j];
            cout<<" Q2[j] = temp  = "<<Q2[j];
            cout<<" ovQ2[j] = 1.0/ Q2[j] = "<< ovQ2[j] ;
            cout<<" fric["<<indND<<"] = 0.1 or "<<fric[indND] <<endl;
        }
    }
    cout<<"------------------INIT Nose ends here---------------------"<<endl;*/
    return;
}
void MOVE_Nose(Real dt){
    Real g(0); 
    Real h(0); 
    Real h_sum(0); 
    Real dti(0); 
    Real ev2(0);
    int j, k, ires, indND, isy;
    
    //cout<<"MOVE Nose begins here"<<endl;
    for (j = 0; j < DOFn; j++){
        for (isy = 0; isy < nsy; isy++) {
            //cout<<"This is No._"<<isy<<"Suzuki Yolanda Factorisation"<<endl;
            for (ires = 0; ires < nres; ires++){
                //cout<<"This is No._"<<ires<<"RESPA Factorisation"<<endl;
                h_sum = 0;
                for (k = 0; k < NHL; k++){
                    //cout<<"h_sum B4(indND = "<<indND<<") = "<<h_sum<<endl;
                    indND = k * DOFn + j;
                    //cout<<"v1["<<indND<<"] = "<< v1[indND]<<endl;
                    g = (Q1[j] * v1[indND] * v1[indND] - temp ) * ovQ2[j];
                    //cout<<" g = "<<g<<endl;
                    //cout<<"v2["<<indND<<"] = "<<v2[indND]<<endl;
                    //cout<<"tsy["<<isy<<"] = "<<tsy[isy]<<endl;
                    v2[indND] += (0.5 * tsy[isy] * g);
                    //cout<<"v2["<<indND<<"] = "<<v2[indND]<<endl;
                    
                    ev2 = exp(- 2.0 * tsy[isy] * v2[indND]);
                    //cout<<"ev2 = "<<ev2<<endl;
                    
                    h_sum = h_sum + (Q1[j] * v1[indND] * v1[indND] * ev2);
                    //cout<<"h_sum @ (indND = "<<indND<<") = "<<h_sum<<endl;
                    //cout<< (Q1[j] * v1[indND] * v1[indND] * ev2)<<endl; 
                }
                //h_sum *= ( (static_cast<Real>(NHL)) /  (static_cast<Real>(NHL + 1)));
                //cout<<"h_sum = "<<h_sum<<endl;
                h = sqrt(  (static_cast<Real>(NHL)) * temp / (V[j] * V[j] * mass[j] + h_sum) );
                //cout<<"h = "<<h<<endl;
                for (k = 0; k < NHL; k++){
                    indND = k * DOFn + j;
                    //cout<<"v1["<<indND<<"] = "<<v1[indND]<<endl;
                    //cout<<"v2["<<indND<<"] = "<<v2[indND]<<endl;
                    v1[indND] *= (h * exp( - tsy[isy] * v2[indND] ) );
                    //cout<<"v1["<<indND<<"] = "<<v1[indND]<<endl;
                }
                for (k = 0; k < NHL; k++){
                    indND = k * DOFn + j;
                    //cout<<"v1["<<indND<<"] = "<<v1[indND]<<endl;
                    g = (Q1[j] * v1[indND] * v1[indND] - temp ) * ovQ2[j];
                    //cout<<" g = "<<g<<endl;
                    v2[indND] += (0.5 * tsy[isy] * g);
                    //cout<<"v2["<<indND<<"] = "<<v2[indND]<<endl;
                }
            }
        }
    }
    KEne = 0;
    for (j = 0; j < DOFn; j++){
        KEne += 0.5 * mass[j] * V[j] * V[j];
    }
    //cout<<"ENDO OF MOVE NOSE"<<endl;
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
    g = (static_cast<Real> (NHL)) * temp;
    int j, k, indND;
    for (j = 0; j < DOFn; j++) {
        //cout<<"F["<<j<<"] = "<<F[j]<<endl;
        //cout<<"V["<<j<<"] = "<<V[j]<<endl;
        a = F[j] * V[j] / g;
        //cout<<"a = "<<a<<endl;
        b = F[j] * F[j] * ovmass[j] / g;
        //cout<<"b = "<<b<<endl;
        rb = sqrt(b);
        arg = (dt * 2.0) * rb * 0.5;
        //cout<<"arg = "<<endl;
        if (arg > 0.00001){
            for (k = 0; k < NHL; k++){
                indND =  k * DOFn + j;
                //cout<<"v1["<<indND<<"] = "<<v1[indND]<<endl;
                st = (1.0 / rb ) * sinh(arg) + ( a / b )* (cosh(arg) - 1.0);
                //cout<<" s = "<<st<<endl;
                sdot = cosh(arg) + (a / rb) * sinh(arg);
                //cout<<"sdot = "<<sdot<<endl;
            }
        }
        else {
            for (k = 0; k < NHL; k++){
                indND =  k * DOFn + j;
                st = ((((b * a / 24.0) * DT2 + b / 6.0) * DT2 + 0.5 * a ) * DT2 + 1.0);
                sdot = (((b * a / 6.0) + b * 0.5) * DT2 + a) * DT2 + 1.0;
            }
        }
        V[j] = (V[j] + F[j] * ovmass[j] * st) / sdot; 
        //cout<<"V["<<j<<"] = "<<V[j]<<endl;
        for (k = 0; k < NHL; k++){
            indND =  k * DOFn + j;
            //cout<<" k * DOFn + j = "<<indND<<endl;
            v1[indND] /= sdot;
            //cout<<"v1["<<indND<<"] = "<<v1[indND]<<endl;
        }
    }
    //cout<<"END OF MOVE Nose V"<<endl;
    return;
}
void MOVE_Nose_R(Real dt){
    PEne = 0;
    for (int j = 0; j < DOFn; j++){
        R[j] += V[j] * dt;
        //cout<<"R["<<j<<"] = "<<R[j]<<endl;
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
            //cout<<"random number = "<<rand_gauss<<endl;
            //cout<<"v2["<<indND<<"] = "<<v2[indND]<<endl;
            v2[indND] = exp( - fric[indND] * dt ) * v2[indND] + sqrt((1.0 - exp( - 2.0 * fric[indND] * dt)) * temp * 2.0 * fric[indND] / ( Q2[j] )) * rand_gauss;
            //cout<<"v2["<<indND<<"] = "<<v2[indND]<<endl; 
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















































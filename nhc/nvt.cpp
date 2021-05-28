#include "nvt.h"
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
    }
    for (int t = 0; t < 10000; t++){
        cout<< position[t] << "  "<< velocity[t] << endl; 
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
    INIT_SAMPLING_Cl();
    //INIT_Langevin();  //Langevin dynamics propagation
    INIT_NHC();
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
        MOVE_NHC(DT2);
        //cout<<"updated here1"<<endl;
        MOVE_VV_V(DT2);
        //cout<<"updated here2"<<endl;
        /*for (int j = 0; j < DOFn; j++) {
            cout<<"R["<<j<<"] = "<<R[j]<<endl;
        }
        for (int j = 0; j < DOFn; j++) {
            cout<<"P["<<j<<"] = "<<P[j]<<endl;
        }*/
        
        MOVE_VV_R(DT);
        //cout<<"updated here3"<<endl;
        /*
        for (int j = 0; j < DOFn; j++) {
            cout<<"R["<<j<<"] = "<<R[j]<<endl;
        }*/
        
      
        
        //MOVE_OU(DT);
        /*
        for (int j = 0; j < DOFn; j++) {
            cout<<"F["<<j<<"] = "<<F[j]<<endl;
        }*/
        
        UPDATE_F(); 
        MOVE_VV_V(DT2);
        //for (int j = 0; j < DOFn; j++) {
        //    cout<<R[j]<<"  "<<P[j]<<endl;}
	/*
        for (int j = 0; j < DOFn; j++) {
            cout<<"P["<<j<<"] = "<<P[j]<<endl;
        }
        cout<<"---------C'est la "<<i<<"th step--------------------"<<endl;*/
        MOVE_NHC(DT2);
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
    return;
}
void PrepOneTrajVV(int t){

    INIT_SAMPLING_Cl();
    INIT_NHC();
    UPDATE_F();

    for (int i = 0; i < nsteps; i++){
        int ind = i + nsteps * t;
        MOVE_NHC(DT2);
        MOVE_VV_V(DT2);
        MOVE_VV_R(DT);


        UPDATE_F();
        MOVE_VV_V(DT2);
        MOVE_NHC(DT2);

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




















































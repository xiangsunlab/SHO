#include "nve.h"
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
        //}
        PropOneTrajVV(t);
    }
    for (int t = 0; t < 10000; t++){
//        cout<< position[t] << "  "<< velocity[t] << endl; 
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

    UPDATE_F();
    //Check Code
    
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


    for (int i = 0; i < nsteps; i++){
        int ind = i + nsteps * t;
        MOVE_VV_V(DT2);
        /*
        for (int j = 0; j < DOFn; j++) {
            cout<<"R["<<j<<"] = "<<R[j]<<endl;
        }
        for (int j = 0; j < DOFn; j++) {
            cout<<"P["<<j<<"] = "<<P[j]<<endl;
        }*/

        MOVE_VV_R(DT);
        /*
        for (int j = 0; j < DOFn; j++) {
            cout<<"R["<<j<<"] = "<<R[j]<<endl;
        }*/
        
        UPDATE_F();
        /*
        for (int j = 0; j < DOFn; j++) {
            cout<<"F["<<j<<"] = "<<F[j]<<endl;
        }*/

        MOVE_VV_V(DT2); 
        for (int j = 0; j < DOFn; j++) {
            cout<<R[j]<<"  "<<P[j]<<endl;}
	/*
        for (int j = 0; j < DOFn; j++) {
            cout<<"P["<<j<<"] = "<<P[j]<<endl;
        }
        cout<<"---------C'est la "<<i<<"th step--------------------"<<endl;*/
            
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




























































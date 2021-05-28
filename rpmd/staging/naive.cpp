#include"naive.h"
int main(){
    std::cout;
    cout.precision(16);
    normal_distribution<Real> normal_dist(0.0,1.0);
    INIT_HAMILTONIAN();
    gen.seed(153);
    //cout<<"Hamiltonian done"<<endl;
    for (int t = 0; t < sampling; t++){
        INIT_SAMPLING_Wigner();
        
        PropOneTrajVV(t);
    }
    for (int t = 0; t < 10000; t++){
        cout<< position[t] << "  "<< velocity[t] << endl;
    }
    
    return 0;
}

void INIT_HAMILTONIAN(){
    R     .resize(DOFn * nbeads,0);    //Cartisian Co-ordinate
    P     .resize(DOFn * nbeads,0);    //Cartisian Momenta
    V     .resize(DOFn * nbeads,0);    //Cartisian Velocity
    F     .resize(DOFn * nbeads,0);    //Cartisian FF Force (FF represents Force Field)
    F1    .resize(DOFn * nbeads,0);    //Ring Polymer Force
    stgR  .resize(DOFn * nbeads,0);
    stgV  .resize(DOFn * nbeads,0);
    stgP  .resize(DOFn * nbeads,0);
    stgF  .resize(DOFn * nbeads,0);
     
    mass  .resize(DOFn,0);    //Real Mass
    vmass. resize(DOFn * nbeads,0);    //Staging Mass
    ovmass.resize(DOFn,0);    //one over Real Mass
    ovvmass.resize(DOFn * nbeads,0);   //one over Staging Mass
    
    Req   .resize(DOFn,0);    //Shifted SHO
    Omega .resize(DOFn,0);    //External FF frequencies
    E     .resize(sampling * nsteps,0);
    U     .resize(sampling * nsteps,0);
    KE    .resize(sampling * nsteps,0);
    w2    .resize(DOFn,0);


    position.resize(10000,0);  //0.01 as a grid and ranging from -50 to +50
    velocity.resize(10000,0);  //0.01 as a grid and ranging from -50 to +50

    for (int j = 0; j < DOFn; j++) {
        mass  [j] = 2.0;
        ovmass[j] = 1.0 / mass[j];
    }
    for (int j = 0; j < DOFn; j++) {
        Req  [j] = 2.33;
        Omega[j] = 1.0;
        w2[j]    = Omega[j] * Omega[j] / (static_cast<Real> (nbeads));
    }
    for (int j = 0; j < DOFn; j++){
        for (int i = 0; i < nbeads; i++){
             int ind = j * nbeads + i;
             if (i == 0){
                 vmass[ind] = 0.0;
                 ovvmass[ind] = mass[j];
             }
             else{
                 vmass[ind]= mass[j] * (static_cast<Real> (i + 1) / static_cast<Real> (i)); //* (static_cast<Real> (nbeads)) / (2.0 * pi * hbar * 2.0 * pi * hbar);
                 ovvmass[ind] = 1.0 / vmass[ind];
             }
             //cout<<vmass[ind]<<"  "<<ovvmass[ind]<<" "<<mass[j]<<endl;
        }
    }
    return;
}

void UPDATE_RP_F(){
    vector<Real> D;    //Displacement accumulation between particles
    D. resize(DOFn * nbeads, 0); 
    
    int ind(0);
    int vor_ind(0);
    int back_ind(0);
    int vor_i(0);
    int back_i(0);
    
    for (int j = 0; j < DOFn; j++) {
        for (int i = 0; i < nbeads; i++) {
            ind = j * nbeads + i;
            vor_i = (i + 1 + nbeads) % nbeads; 
            back_i = (i - 1 + nbeads) % nbeads;
            vor_ind = j * nbeads + vor_i;
            back_ind = j * nbeads + back_i;
            D[ind] = 2.0 * R[ind] - R[vor_i] - R[back_i];
            F1[ind] = -wp2 * D[ind] * mass[j];
            PEne += 0.25 * mass[j] * wp2 * D[ind] * D[ind];
            //cout<<F1[ind]<<endl;
        }
    }
    return;
}

void UPDATE_EX_F() {
    for (int j = 0; j < DOFn; j++) {
        for (int i = 0; i < nbeads; i++) {
            int ind =  j * nbeads + i;
            F[ind] = - mass[j] * w2[j] * (R[ind] - Req[j]);
            //cout<<F[ind]<<endl;
            PEne += 0.5 * mass[j] * w2[j] * (R[ind] - Req[j]) * (R[ind] - Req[j]);
        }
    }
    PEne /= sqrt(static_cast<Real> (nbeads));
    return;
}

void MOVE_VV_V(Real dt) {
    for (int j = 0; j < DOFn; j++) {
        for (int i = 0; i < nbeads; i++) {
            int ind =  j * nbeads + i;
            P[ind] += (F[ind] + F1[ind]) * dt;
            //cout<<P[ind]<<endl;
        }
    }
    return;
}

void MOVE_STGV(Real dt){
    KEne = 0.0;
    for (int j = 0; j < DOFn; j++) {
        for (int i = 0; i < nbeads; i++) {
            int ind =  j * nbeads + i;
            P[ind] += F[ind] * dt;
            KEne += P[ind] * P[ind] * ovmass[j] * 0.5;
        }
    }
    //KEne /= (static_cast<Real> (nbeads))
    return;
}

void MOVE_VV_R(Real dt) {
    for (int j = 0; j < DOFn; j++) {
        for (int i = 0; i < nbeads; i++) {
            int ind =  j * nbeads + i;
            KEne += 0.5 * P[ind] * P[ind] * ovmass[j];
            V[ind] = P[ind] * ovvmass[j];
            R[ind] += V[ind] * dt;
        }
    }
    return;
}

void PropOneTrajVV(int t) {
    for (int nstep = 0; nstep < nsteps; nstep++){
        //KEne = 0; 
        //MOVE_VV_V(DT2);
        //PEne = 0;;
        //MOVE_VV_R(DT);
        //UPDATE_RP_F();
        //UPDATE_EX_F();
        //MOVE_VV_V(DT2);
        MOVE_STGV(DT2);
        
        
        NAIVE2STG();
        STGRPF();
        STGRPP(DT2);
        STGRPR(DT);
        PEne = 0;
        STGRPF();
        STGRPP(DT2);
        STG2NAIVE();
        
        UPDATE_EX_F();
        MOVE_STGV(DT2);
        if (nstep % 1000 == 0){
            int q = int((R[0] - Req[0] + 50 ) * 100 );
            int p = int((P[0] + 50) * 100);
            //cout<<P[0]<<"  "<<R[0]<<endl;
            if (q > 10000 || p > 10000 ){}
            else {
                velocity[p] += 1;
                position[q] += 1;
            
            Energy = PEne + KEne;//static_cast<Real>(nbeads));
            //cout<<Energy<<"  "<<PEne<<"  "<<KEne<<endl;
            }
        }
        else{}
    }
    return;
}

void INIT_SAMPLING_Wigner(){
    normal_distribution<Real> normal_dist(0.0,1.0);
    vector<Real> sigma_R;    //Standard derivation of each R[j] (STD = STandard Derivation)
    vector<Real> sigma_P;    //STD of each P[j];

    sigma_R.resize(DOFn,0);
    sigma_P.resize(DOFn,0);

    for (int j = 0; j < DOFn; j++) {
        //for (int i = 0; i < nbeads; i++) {
            //int ind = j * nbeads + i;
            sigma_R[j] = sqrt(hbar/ (2.0 * Omega[j] * tanh(0.5 * beta * hbar * Omega[j])) );
            sigma_P[j] = sigma_R[j] * Omega[j];
            //cout<<sigma_R[j]<<"  "<<sigma_P[j]<<endl;
        //}
    }
    for (int j = 0; j < DOFn; j++) {
        for (int i = 0; i < 1; i++) {
            int ind = j * nbeads + i;
            R[ind] = normal_dist(gen) * sigma_R[j] + Req[j];
            P[ind] = normal_dist(gen) * sigma_P[j] / sqrt(static_cast<Real> (nbeads));
            //cout<<R[ind]<<"  "<<P[ind]<<endl;
        }
    }
    UPDATE_EX_F();
    for (int j = 0; j < DOFn; j++) {
        for (int i = 1; i < nbeads; i++) { 
            int ind = j * nbeads + i;
            R[ind] = R[(j*nbeads)] + DT * P[j*nbeads] * static_cast<Real> (i+1);// 0.1 * normal_dist(gen);
            P[ind] = P[j*nbeads] + DT * F[j*nbeads] * static_cast<Real> (i+1); //0.1 * normal_dist(gen);
            //cout<<"R["<<(j*nbeads)<<"] = "<<R[j*nbeads]<<endl;
            //cout<<"P["<<j*nbeads<<"] = "<<P[j*nbeads]<<endl;
            //cout<<"DT = "<<DT<<endl;
            //cout<<R[ind]<<"  "<<P[ind]<<endl;
        }
    }
    NAIVE2STG();
    STG2NAIVE();
    STGRPF(); 
    return;
}

void INIT_SAMPLING_Cl(){
    normal_distribution<Real> normal_dist(0.0,1.0);
    vector<Real> sigma_R;    //Standard derivation of each R[j] (STD = STandard Derivation)
    vector<Real> sigma_P;    //STD of each P[j];
    sigma_R.resize(DOFn,0);
    sigma_P.resize(DOFn,0);

    for (int j = 0; j < DOFn; j++) {
        for (int i = 0; i < nbeads; i++) {
            int ind = j * nbeads + i;
            sigma_R[j] = sqrt(1.0 / (beta * Omega[j] * Omega[j] * mass[j]) );
            sigma_P[j] = sqrt(temp * mass[j]);
        }
    }
    for (int j = 0; j < DOFn; j++) {
        for (int i = 0; i < nbeads; i++) {
            int ind = j * nbeads + i;
            R[j] = normal_dist(gen) * sigma_R[j] + Req[j];
            P[j] = normal_dist(gen) * sigma_P[j];
        }
    }
    return;
}

void NAIVE2STG(){
    for (int j = 0; j < DOFn; j++){
        for (int i = 0; i < nbeads; i++){
            int ind = j * nbeads + i;
            if (i == 0){
                stgR[ind] = R[ind];
                stgP[ind] = P[ind];
            }
            else {
                int ind0 = j * nbeads;
                int ind1 = j * nbeads + (i + 1) % nbeads;
                stgR[ind] = R[ind] - ((static_cast<Real>(i)) * R[ind1] + R[ind0]) / (static_cast<Real> (i+1));
                stgP[ind] = P[ind] - ((static_cast<Real>(i)) * P[ind1] + P[ind0]) / (static_cast<Real> (i+1));
                //cout<<ind<<"  "<<R[ind]<<"  "<<ind1<<"  "<<R[ind1]<<"  "<<(R[ind]-((static_cast<Real> (i))*R[ind1] + R[ind0])/(static_cast<Real> (i+1)))<<"  "<<R[ind0]<<endl;
            }
        }
    }
    /*
    for (int j = 0; j < DOFn; j++){
        for (int i = 0; i < nbeads; i++){
            int ind = j * nbeads + i;
            cout<<R[ind]<<"  "<<stgR[ind]<<" "<<P[ind]<<"  " <<stgP[ind]<<endl;
        }
    }
    */
    return;
}

void STG2NAIVE(){
    for (int j = 0; j < DOFn; j++){
        int ind0 = j * nbeads;
        for (int i = 0; i < nbeads; i++){
            int ind = j * nbeads + i;
            if (i == 0){
                R[ind] = stgR[ind];
                P[ind] = stgP[ind];
            }
            else {
                Real a(0);
                Real b(0);
                for (int i1 = i; i1 < nbeads; i1++){
                    a += stgR[i1] * (static_cast<Real> (i)) / (static_cast<Real>(i1));
                    b += stgP[i1] * (static_cast<Real> (i)) / (static_cast<Real>(i1));
                }
                R[ind] = stgR[ind0] + a;
                P[ind] = stgP[ind0] + b;
            }
        }
    }
    /*
    for (int j = 0; j < DOFn; j++){
        for (int i = 0; i < nbeads; i++){
            int ind = j * nbeads + i;
            cout<<R[ind]<<"  "<<stgR[ind]<<" "<<P[ind]<<"  " <<stgP[ind]<<endl;
        }
    }*/
    return;
}

void STGRPF(){
    for (int j = 0; j < DOFn; j++){
        for (int i = 0; i < nbeads; i++){
            int ind = j * nbeads + i;
            stgF[ind] = - vmass[ind] * wp2 * stgR[ind];
            PEne += 0.5 * vmass[ind] * wp2 * stgR[ind] * stgR[ind];
            //cout<<stgF[ind]<<"  "<<stgR[ind]<<endl;
        }
    //PEne += 0.5 * mass[j] * wp2 * stgR[j * nbeads] * stgR[j * nbeads];
    }
    //PEne /=  (static_cast<Real>(nbeads));
    return; 
}


void STGRPP(Real dt){
    for (int j = 0; j < DOFn; j++){
        for (int i = 0; i < nbeads; i++){
            int ind = j * nbeads + i;
            stgP[ind] += stgF[ind] * dt;
            //cout<<stgF[ind]<<"  "<<stgP[ind]<<endl;
        }
    }
    return;
}

void STGRPR(Real dt){
    for (int j = 0; j < DOFn; j++){
        for (int i = 0; i < nbeads; i++){
            int ind = j * nbeads + i;
            //cout<<stgR[ind]<<"  ";
            stgR[ind] += stgP[ind] * dt * ovvmass[ind];
            //cout<<stgP[ind]<<"  "<<stgR[ind]<<endl;
        }
    }
    return;
}

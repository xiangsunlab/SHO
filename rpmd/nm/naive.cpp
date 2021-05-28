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
        //cout<<"Initialise Done"<<endl;
        INIT_NM();
       
        NAIVE2NM();
        NM2NAIVE();
        NAIVE2NM();
        //for (int i = 0; i < nbeads; i++){
        //    cout<<nmR[i]<<"  "<<nmP[i]<<endl;
        //}
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

void PropOneTrajVV(int t) {
    for (int nstep = 0; nstep < nsteps; nstep++){
        //for (int i = 0; i < nbeads; i++){
            //cout<<P[i]<<"  "<<R[i]<<"  "<<F[i]<<endl;
            //ut<<P[i]<<" "<<nmP[i]<<"  "<<R[i]<<"  "<<nmR[i]<<endl;
        //}
        NM_Momenta(DT2);
        
        //for (int i = 0; i < nbeads; i++){
            //cout<<eva00[i]<<"  "<<eva01[i]<<"  "<<eva10[i]<<"  "<<eva11[i]<<"  "<<endl;
            //cout<<P[i]<<"  "<<R[i]<<"  "<<F[i]<<endl;
            //cout<<P[i]<<" "<<nmP[i]<<"  "<<R[i]<<"  "<<nmR[i]<<endl;
        //}
        NAIVE2NM();
        //for (int i = 0; i < nbeads; i++){
            //cout<<P[i]<<"  "<<R[i]<<"  "<<F[i]<<endl;
            //cout<<P[i]<<" "<<nmP[i]<<"  "<<R[i]<<"  "<<nmR[i]<<endl;
        //}
        NM_Prop(DT);
        //for (int i = 0; i < nbeads; i++){
            //ut<<P[i]<<"  "<<R[i]<<"  "<<F[i]<<endl;
            //cout<<P[i]<<" "<<nmP[i]<<"  "<<R[i]<<"  "<<nmR[i]<<endl;
        //}
        NM2NAIVE();
        //for (int i = 0; i < nbeads; i++){
            //cout<<P[i]<<"  "<<R[i]<<"  "<<F[i]<<endl;
        //}
        NM_Force();
        //for (int i = 0; i < nbeads; i++){
            //cout<<P[i]<<"  "<<R[i]<<"  "<<F[i]<<endl;
        //}
        NM_Momenta(DT2);
        //for (int i = 0; i < nbeads; i++){
            //cout<<P[i]<<"  "<<R[i]<<"  "<<F[i]<<endl;
        //}
        //cout<<R[0]<<"  "<<P[0]<<endl;
        if (nstep % 10 == 0){
            int q = int((R[0] - Req[0] + 50 ) * 100 );
            int p = int((P[0] + 50) * 100);
            //cout<<P[0]<<"  "<<R[0]<<endl;
            if (q > 10000 || p > 10000 ){}
            else {
                velocity[p] += 1;
                position[q] += 1;
            
            //Energy = PEne + KEne;//static_cast<Real>(nbeads));
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
            P[ind] = normal_dist(gen) * sigma_P[j]; // sqrt(static_cast<Real> (nbeads));
            //cout<<R[ind]<<"  "<<P[ind]<<endl;
        }
    }
    for (int j = 0; j < DOFn; j++) {
        for (int i = 1; i < nbeads; i++) { 
            int ind = j * nbeads + i;
            Real gauss_dist;
            gauss_dist = normal_dist(gen);
            R[ind] = R[(j*nbeads)] + 0.1 * gauss_dist;
            gauss_dist = normal_dist(gen);
            P[ind] = P[j*nbeads] +  0.1 * gauss_dist;//DT * F[j*nbeads] * static_cast<Real> (i+1); //0.1 * normal_dist(gen);
            //cout<<R[ind]<<"  "<<P[ind]<<endl;
        }
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


void INIT_NM(){
    nmR      .resize(DOFn * nbeads, 0);
    nmP      .resize(DOFn * nbeads, 0);
    nmV      .resize(DOFn * nbeads, 0);
    nmF      .resize(DOFn * nbeads, 0);

    omega    .resize(nbeads, 0);
    fric     .resize(nbeads, 0);
    c1       .resize(DOFn * nbeads, 0);
    c2       .resize(DOFn * nbeads, 0);
    eva00    .resize(DOFn * nbeads, 0);
    eva01    .resize(DOFn * nbeads, 0);
    eva10    .resize(DOFn * nbeads, 0);
    eva11    .resize(DOFn * nbeads, 0);
    Uni      .resize(nbeads * nbeads, 0);
    UT       .resize(nbeads * nbeads, 0);
    
    //Initialise normal mode frequency
    for (int i = 0; i < nbeads; i++){
        omega[i] =  2.0 * wp * sin((static_cast<Real>(i)) * pi / (static_cast<Real> (nbeads)));
        if (i == 0) {
            fric[i] = 1.0;
        }
        else {
            fric[i] = 2.0 * omega[i];
        }
    }
    for (int j = 0; j < DOFn; j++){
        for (int i = 0; i < nbeads; i++){
            int ind = j * nbeads + i;
            c1[ind] = exp(- DT2 * fric[i]);
            c2[ind] = sqrt(1.0 - c1[i] * c1[i]) * sqrt(mass[j] / beta);
            //cout<<fric[i]<<"  "<<c1[ind]<<"  "<<c2[ind]<<endl;
        }
    }
    //Evolution Matrices
    for (int j = 0; j < DOFn; j++){
        for (int i = 0; i < nbeads; i++){
            int ind =  (j * nbeads + i);
            if (i == 0){
                eva00[ind] = 1.0;
                eva01[ind] = 0.0;
                eva10[ind] = DT/mass[j];
                eva11[ind] = 1.0; 
            }
            else {
                eva00[ind] = cos(omega[i] * DT);
                eva01[ind] = - omega[i] * mass[j] * sin(omega[i] * DT);
                eva10[ind] = sin(omega[i] * DT) /(omega[i] * mass[j]);
                eva11[ind] = eva00[ind];
            }    
        }
    }
    for (int k = 0; k < nbeads; k++){
        for (int j = 0; j < nbeads; j++){
            int ind = k * nbeads + j;
            if (k < nbeads/2){
                if (k == 0){
                    Uni[ind] = sqrt(1.0 / (static_cast<Real> (nbeads))) ;
                }
                else {
                    Uni[ind] = sqrt(2.0 / (static_cast<Real> (nbeads))) * cos(2.0 * pi * (static_cast<Real>(k * (j + 1)))/ (static_cast<Real> (nbeads)));
                }
            }
            else {
                if (k == (nbeads / 2)){
                    Uni[ind] = sqrt(1.0 / (static_cast<Real> (nbeads))) * pow(-1, (static_cast<Real>(j + 1)));
                }
                else {
                    Uni[ind] = sqrt(2.0 / (static_cast<Real> (nbeads))) * sin(2.0 * pi * (static_cast<Real>(k * (j + 1)))/ (static_cast<Real> (nbeads)));
                }
            }
        }
    }
    for (int i = 0; i < nbeads; i++){
        for (int j = 0; j < nbeads; j++){
            int indt = (i * nbeads + j);
            int ind =  (j * nbeads + i);
            UT[indt] = Uni[ind];
            //cout<<Uni[indt]<<"  ";
        }
        //cout<<endl;
    }
    for (int i = 0; i < nbeads; i++){
        for (int j = 0; j < nbeads; j++){
            int ind =  (i * nbeads + j);
            //cout<<UT[ind]<<"  ";
        }
        //cout<<endl;
    }
    for (int j = 0; j < DOFn; j++){
        for (int i = 0; i < nbeads; i++){
            int ind =  (j * nbeads + i);
            //cout<<eva00[ind]<<"  "<<eva01[ind]<<"  "<<eva10[ind]<<"  "<<eva11[ind]<<endl;
        }
    }
    NM_Force();
    return;
}

void NAIVE2NM(){
    for (int j = 0; j < DOFn; j++){    //Degree of freedom terms 
        for (int k = 0; k < nbeads; k++){
            int dofn = j * nbeads + k;
            nmR[dofn] = 0;
            nmP[dofn] = 0;
            for (int i = 0; i < nbeads; i++){
                int ind = k * nbeads + i;
                int dof2 = j * nbeads + i;
                nmR[dofn] += Uni[ind] * R[dof2];
                nmP[dofn] += Uni[ind] * P[dof2];
            }
        }
    }
    return;
}

void NM2NAIVE(){
    for (int j = 0; j < DOFn; j++){
        for (int i = 0; i < nbeads; i++){
            int ind = j * nbeads + i;
            R[ind] = 0;
            P[ind] = 0;
            for (int k = 0; k < nbeads; k++){
                int indt = k * nbeads + i;
                int ind1 = j * nbeads + k;
                R[ind] += Uni[indt] * nmR[ind1];
                P[ind] += Uni[indt] * nmP[ind1];
            }
        }
    }
    return;
}

void NM_Prop(Real dt){
    Real R(0);
    Real P(0);
    for (int j = 0; j < DOFn; j++){
        for (int k = 0; k < nbeads; k++){
            int ind = j * nbeads + k;
            P = eva00[ind] * nmP[ind] + eva01[ind] * nmR[ind];
            //cout<<P<<"  ";
            R = eva10[ind] * nmP[ind] + eva11[ind] * nmR[ind];
            //cout<<R<<"  ";
            nmP[ind] = P;
            nmR[ind] = R;
            //cout<<nmP[ind]<<"  "<<nmR[ind]<<endl;
        }
    }
    return;
}
void NM_Force(){
    for (int j = 0; j < DOFn; j++){
        for (int i = 0; i < nbeads; i++){
            int ind =  (j * nbeads + i);
            F[ind] = - mass[j] * w2[j] * (R[ind] - Req[j]);
        }
    }
    return;
}
void NM_Momenta(Real dt){
    for (int j = 0; j < DOFn; j++){
        for (int i = 0; i < nbeads; i++){
            int ind =  (j * nbeads + i);
            P[ind] += F[ind] * dt;
        }
    }
    return;
} 

/*
void NAIVE2NM() {
    for 
    return;
}*/

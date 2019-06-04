#include "constants.h"
#include "functions.h"

//Array functions and things

PetscInt c_index(PetscInt x,PetscInt y,PetscInt z,PetscInt comp,PetscInt ion,PetscInt Nx,PetscInt Ny,PetscInt Nz)
{
    return Nc*Ni* (Nx *(Ny*z + y) + x) + comp*Ni+ion;
//    return Nc*Ni* (Ny *(Nz*x + z) + y) + comp*Ni+ion;
}
PetscInt phi_index(PetscInt x,PetscInt y,PetscInt z,PetscInt comp,PetscInt Nx,PetscInt Ny,PetscInt Nz)
{
    return Nc* (Nx *(Ny*z + y) + x) + comp;
//    return Nc* (Ny *(Nz*x + z) + y) + comp;
}
PetscInt al_index(PetscInt x,PetscInt y,PetscInt z,PetscInt comp,PetscInt Nx,PetscInt Ny,PetscInt Nz)
{
    return (Nc-1)* (Nx *(Ny*z + y) + x) + comp;
//    return (Nc-1)* (Ny *(Nz*x + z) + y) + comp;
}
PetscInt xy_index(PetscInt x,PetscInt y,PetscInt z,PetscInt Nx,PetscInt Ny,PetscInt Nz)
{
    return Nx *(Ny*z + y)+x;
//    return Ny *(Nz*x + z)+y;
}
//Index based on Nv, which can change to either include or exclude alpha
PetscInt Ind_1(PetscInt x,PetscInt y,PetscInt z,PetscInt ion,PetscInt comp,PetscInt Nx,PetscInt Ny,PetscInt Nz)
{
//    return Nv*(Nx *(Ny*z + y)+x)+ion*Nc+comp;
    return Nv*(Ny *(Nz*x + z)+y)+ion*Nc+comp;
}
// Index based on solving c,phi, and alpha.
PetscInt Ind_2(PetscInt x,PetscInt y,PetscInt z,PetscInt ion,PetscInt comp,PetscInt nx,PetscInt ny,PetscInt nz)
{
    return ((Ni+2)*Nc-1)*(nx *(ny*z + y)+x)+ion*Nc+comp;
//    return ((Ni+2)*Nc-1)*(ny *(nz*x + z)+y)+ion*Nc+comp;
}

PetscErrorCode init_simstate(Vec state,struct SimState *state_vars,struct AppCtx *user)
{
    PetscErrorCode ierr;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt Nz = user->Nz;
    //Setup indices
    int x,y,z,comp,ion;
    PetscInt *c_ind = (PetscInt *) malloc(sizeof(PetscInt)*Nx*Ny*Nz*Nc*Ni);
    PetscInt *phi_ind = (PetscInt *) malloc(sizeof(PetscInt)*Nx*Ny*Nz*Nc);
    for(z=0; z<Nz ; z++){
        for(y = 0; y < Ny; y++){
            for(x = 0; x < Nx; x++){
                for(comp = 0; comp < Nc; comp++){
                    for(ion = 0; ion < Ni; ion++){
                        c_ind[c_index(x,y,z,comp,ion,Nx,Ny,Nz)] = Ind_1(x,y,z,ion,comp,Nx,Ny,Nz);
                    }
                    phi_ind[phi_index(x,y,z,comp,Nx,Ny,Nz)] = Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz);
                }
            }
        }
    }
    ierr = ISCreateGeneral(PETSC_COMM_WORLD,Nx*Ny*Nz*Ni*Nc,c_ind,PETSC_COPY_VALUES,&state_vars->c_ind); CHKERRQ(ierr);
    ierr = ISCreateGeneral(PETSC_COMM_WORLD,Nx*Ny*Nz*Nc,phi_ind,PETSC_COPY_VALUES,&state_vars->phi_ind); CHKERRQ(ierr);

    free(phi_ind);free(c_ind);
    if(!separate_vol) {
        PetscInt *al_ind = (PetscInt *) malloc(sizeof(PetscInt)*Nx*Ny*Nz*(Nc-1));
        for (z = 0; z < Nz; z++) {
            for (y = 0; y < Ny; y++) {
                for (x = 0; x < Nx; x++) {
                    for (comp = 0; comp < Nc - 1; comp++){
                        al_ind[al_index(x,y,z,comp,Nx,Ny,Nz)] = Ind_1(x,y,z,Ni+1,comp,Nx,Ny,Nz);
                    }
                }
            }
        }
        ierr = ISCreateGeneral(PETSC_COMM_WORLD, Nz*Nx * Ny * (Nc - 1), al_ind, PETSC_COPY_VALUES, &state_vars->al_ind);
        CHKERRQ(ierr);
        free(al_ind);
    }
    else{
        state_vars->alpha = (PetscReal*)malloc(sizeof(PetscReal)*Nx*Ny*Nz*(Nc-1));
    }
    extract_subarray(state,state_vars);
    return ierr;
}

PetscErrorCode extract_subarray(Vec state,struct SimState *state_vars)
{
    if(Profiling_on) {
        PetscLogEventBegin(event[2], 0, 0, 0, 0);
    }
    PetscErrorCode ierr;
    ierr = VecGetSubVector(state,state_vars->c_ind,&state_vars->c_vec); CHKERRQ(ierr);
    ierr = VecGetArray(state_vars->c_vec,&state_vars->c); CHKERRQ(ierr);

    ierr = VecGetSubVector(state,state_vars->phi_ind,&state_vars->phi_vec); CHKERRQ(ierr);
    ierr = VecGetArray(state_vars->phi_vec,&state_vars->phi); CHKERRQ(ierr);
    if(!separate_vol) {
        ierr = VecGetSubVector(state, state_vars->al_ind, &state_vars->al_vec);
        CHKERRQ(ierr);
        ierr = VecGetArray(state_vars->al_vec, &state_vars->alpha);
        CHKERRQ(ierr);
    }
    if(Profiling_on) {
        PetscLogEventEnd(event[2], 0, 0, 0, 0);
    }

    return ierr;

}

PetscErrorCode restore_subarray(Vec state,struct SimState *state_vars)
{
    if(Profiling_on) {
        PetscLogEventBegin(event[3], 0, 0, 0, 0);
    }
    PetscErrorCode ierr;

    ierr = VecRestoreArray(state_vars->c_vec,&state_vars->c); CHKERRQ(ierr);
    ierr = VecRestoreSubVector(state,state_vars->c_ind,&state_vars->c_vec); CHKERRQ(ierr);

    ierr = VecRestoreArray(state_vars->phi_vec,&state_vars->phi); CHKERRQ(ierr);
    ierr = VecRestoreSubVector(state,state_vars->phi_ind,&state_vars->phi_vec); CHKERRQ(ierr);

    if(!separate_vol) {
        ierr = VecRestoreArray(state_vars->al_vec, &state_vars->alpha);
        CHKERRQ(ierr);
        ierr = VecRestoreSubVector(state, state_vars->al_ind, &state_vars->al_vec);
        CHKERRQ(ierr);
        state_vars->alpha = NULL;
    }

    state_vars->c = NULL;
    state_vars->phi = NULL;
    if(Profiling_on) {
        PetscLogEventEnd(event[3], 0, 0, 0, 0);
    }

    return ierr;

}
PetscErrorCode extract_subarray_Read(Vec state,struct SimState *state_vars)
{
    if(Profiling_on) {
        PetscLogEventBegin(event[2], 0, 0, 0, 0);
    }
    PetscErrorCode ierr;
    ierr = VecGetSubVector(state,state_vars->c_ind,&state_vars->c_vec); CHKERRQ(ierr);
    ierr = VecGetArrayRead(state_vars->c_vec,&state_vars->c); CHKERRQ(ierr);

    ierr = VecGetSubVector(state,state_vars->phi_ind,&state_vars->phi_vec); CHKERRQ(ierr);
    ierr = VecGetArrayRead(state_vars->phi_vec,&state_vars->phi); CHKERRQ(ierr);
    if(!separate_vol) {
        ierr = VecGetSubVector(state, state_vars->al_ind, &state_vars->al_vec);
        CHKERRQ(ierr);
        ierr = VecGetArrayRead(state_vars->al_vec, &state_vars->alpha);
        CHKERRQ(ierr);
    }
    if(Profiling_on) {
        PetscLogEventEnd(event[2], 0, 0, 0, 0);
    }

    return ierr;

}

PetscErrorCode restore_subarray_Read(Vec state,struct SimState *state_vars)
{
    if(Profiling_on) {
        PetscLogEventBegin(event[3], 0, 0, 0, 0);
    }
    PetscErrorCode ierr;

    ierr = VecRestoreArrayRead(state_vars->c_vec,&state_vars->c); CHKERRQ(ierr);
    ierr = VecRestoreSubVector(state,state_vars->c_ind,&state_vars->c_vec); CHKERRQ(ierr);

    ierr = VecRestoreArrayRead(state_vars->phi_vec,&state_vars->phi); CHKERRQ(ierr);
    ierr = VecRestoreSubVector(state,state_vars->phi_ind,&state_vars->phi_vec); CHKERRQ(ierr);

    if(!separate_vol) {
        ierr = VecRestoreArrayRead(state_vars->al_vec, &state_vars->alpha);
        CHKERRQ(ierr);
        ierr = VecRestoreSubVector(state, state_vars->al_ind, &state_vars->al_vec);
        CHKERRQ(ierr);
        state_vars->alpha = NULL;
    }

    state_vars->c = NULL;
    state_vars->phi = NULL;
    if(Profiling_on) {
        PetscLogEventEnd(event[3], 0, 0, 0, 0);
    }

    return ierr;

}
PetscErrorCode copy_simstate(Vec current_state,struct SimState *state_vars_past)
{
    PetscErrorCode ierr;
    ierr = VecCopy(current_state,state_vars_past->v); CHKERRQ(ierr);
    ierr = extract_subarray(state_vars_past->v,state_vars_past); CHKERRQ(ierr);
    return ierr;
}

void init_arrays(struct AppCtx*user)
{
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt Nz = user->Nz;
    PetscInt nx = 2*width_size+1;
    PetscInt ny = 2*width_size+1;
    PetscInt nz = user->Nz;;
    //Flux quantities
    user->flux->mflux = (PetscReal*) malloc(Nx*Ny*Nz*Ni*Nc*sizeof(PetscReal));
    user->flux->dfdci = (PetscReal*) malloc(Nx*Ny*Nz*Ni*Nc*sizeof(PetscReal));
    user->flux->dfdce = (PetscReal*) malloc(Nx*Ny*Nz*Ni*Nc*sizeof(PetscReal));
    user->flux->dfdphim = (PetscReal*) malloc(Nx*Ny*Nz*Ni*Nc*sizeof(PetscReal));
    user->flux->wflow = (PetscReal*) malloc(Nx*Ny*Nz*(Nc-1)*sizeof(PetscReal));
    user->flux->dwdpi = (PetscReal*) malloc(Nx*Ny*Nz*(Nc-1)*sizeof(PetscReal));
    user->flux->dwdal = (PetscReal*) malloc(Nx*Ny*Nz*(Nc-1)*sizeof(PetscReal));

    //Gating variables (present)
    user->gate_vars->mNaT = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars->hNaT = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars->gNaT = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars->mNaP = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars->hNaP = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars->gNaP = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars->mKDR = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars->gKDR = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars->mKA = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars->hKA = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars->gKA = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars->yNMDA = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars->zNMDA = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars->dNMDA = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars->gNMDA = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    //Gating variables (past)
    user->gate_vars_past->mNaT = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars_past->hNaT = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars_past->gNaT = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars_past->mNaP = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars_past->hNaP = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars_past->gNaP = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars_past->mKDR = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars_past->gKDR = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars_past->mKA = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars_past->hKA = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars_past->gKA = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars_past->yNMDA = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars_past->zNMDA = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars_past->dNMDA = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gate_vars_past->gNMDA = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));



    //Excitation
    user->gexct->pNa = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gexct->pK = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gexct->pCl = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    user->gexct->pGlu = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));

    //Constant params
    user->con_vars->ao = (PetscReal*) malloc(Nx*Ny*Nz*Nc*sizeof(PetscReal));
    user->con_vars->zo = (PetscReal*) malloc(Nx*Ny*Nz*Nc*sizeof(PetscReal));
    user->con_vars->zeta1 = (PetscReal*) malloc(Nx*Ny*Nz*(Nc-1)*sizeof(PetscReal));
    user->con_vars->zetaalpha = (PetscReal*) malloc((Nc-1)*sizeof(PetscReal));
    user->con_vars->cm = (PetscReal*) malloc((Nc-1)*Nx*Ny*Nz*sizeof(PetscReal));
    user->con_vars->ell = (PetscReal *) malloc(Nx*Ny*Nz*sizeof(PetscReal));

    //Diffusion in ctx
    user->Dcs = (PetscReal*) malloc(Nx*Ny*Nz*Ni*Nc*3*sizeof(PetscReal));
    user->Dcb = (PetscReal*) malloc(Nx*Ny*Nz*Ni*Nc*3*sizeof(PetscReal));

    //Small Grid variables

    // Past membrane voltage storage
    user->vm_past = (PetscReal*) malloc(Nx*Ny*Nz*sizeof(PetscReal));
    //Grid Gating variables
    user->grid_gate_vars->mNaT = (PetscReal*) malloc(nx*ny*nz*sizeof(PetscReal));
    user->grid_gate_vars->hNaT = (PetscReal*) malloc(nx*ny*nz*sizeof(PetscReal));
    user->grid_gate_vars->gNaT = (PetscReal*) malloc(nx*ny*nz*sizeof(PetscReal));
    user->grid_gate_vars->mNaP = (PetscReal*) malloc(nx*ny*nz*sizeof(PetscReal));
    user->grid_gate_vars->hNaP = (PetscReal*) malloc(nx*ny*nz*sizeof(PetscReal));
    user->grid_gate_vars->gNaP = (PetscReal*) malloc(nx*ny*nz*sizeof(PetscReal));
    user->grid_gate_vars->mKDR = (PetscReal*) malloc(nx*ny*nz*sizeof(PetscReal));
    user->grid_gate_vars->gKDR = (PetscReal*) malloc(nx*ny*nz*sizeof(PetscReal));
    user->grid_gate_vars->mKA = (PetscReal*) malloc(nx*ny*nz*sizeof(PetscReal));
    user->grid_gate_vars->hKA = (PetscReal*) malloc(nx*ny*nz*sizeof(PetscReal));
    user->grid_gate_vars->gKA = (PetscReal*) malloc(nx*ny*nz*sizeof(PetscReal));
    user->grid_gate_vars->yNMDA = (PetscReal*) malloc(nx*ny*nz*sizeof(PetscReal));
    user->grid_gate_vars->zNMDA = (PetscReal*) malloc(nx*ny*nz*sizeof(PetscReal));
    user->grid_gate_vars->dNMDA = (PetscReal*) malloc(nx*ny*nz*sizeof(PetscReal));
    user->grid_gate_vars->gNMDA = (PetscReal*) malloc(nx*ny*nz*sizeof(PetscReal));

    //Grid state_vars
    user->grid_vars->c = (PetscReal*) malloc(Nc*Ni*nx*ny*nz*sizeof(PetscReal));
    user->grid_vars->phi = (PetscReal*) malloc(Nc*nx*ny*nz*sizeof(PetscReal));
    user->grid_vars->alpha = (PetscReal*) malloc((Nc-1)*nx*ny*nz*sizeof(PetscReal));
    user->grid_vars->v = NULL;
    user->grid_vars->phi_ind = NULL;
    user->grid_vars->phi_vec = NULL;
    user->grid_vars->c_ind = NULL;
    user->grid_vars->c_vec = NULL;
    user->grid_vars->al_ind = NULL;
    user->grid_vars->al_vec = NULL;

    //Grid past state_Vars

    user->grid_vars_past->c = (PetscReal*) malloc(Nc*Ni*nx*ny*nz*sizeof(PetscReal));
    user->grid_vars_past->phi = (PetscReal*) malloc(Nc*nx*ny*nz*sizeof(PetscReal));
    user->grid_vars_past->alpha = (PetscReal*) malloc((Nc-1)*nx*ny*nz*sizeof(PetscReal));
    user->grid_vars_past->v = NULL;
    user->grid_vars_past->phi_ind = NULL;
    user->grid_vars_past->phi_vec = NULL;
    user->grid_vars_past->c_ind = NULL;
    user->grid_vars_past->c_vec = NULL;
    user->grid_vars_past->al_ind = NULL;
    user->grid_vars_past->al_vec = NULL;

    //dt saving
    user->dt_space = (PetscReal*) malloc(Nx*Ny*sizeof(PetscReal));
    for (int y = 0; y < Ny; y++) {
        for (int x = 0; x < Nx; x++) {
            user->dt_space[xy_index(x,y,0,Nx,0,0)] = user->dt;
        }
    }



}

void parameter_dependence(struct AppCtx *user)
{
    PetscReal soma,soma_width, dendrite,dendrite_width;

//    soma=3;soma_width=1;
//    dendrite=1; dendrite_width=2;
//    soma = 0.1499;soma_width = 0.0301;
    soma = 0.039; soma_width = 0.00301;
    PetscReal dz = user->dz;



    struct ConstVars *con_vars = user->con_vars;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt Nz = user->Nz;
    PetscInt x,y,z;

    PetscReal cmt,sa,voli,vole;

    //Gating variables
    con_vars->pNaT = (PetscReal*)malloc(sizeof(PetscReal)*Nx*Ny*Nz);
    con_vars->pNaP = (PetscReal*)malloc(sizeof(PetscReal)*Nx*Ny*Nz);
    con_vars->pKDR = (PetscReal*)malloc(sizeof(PetscReal)*Nx*Ny*Nz);
    con_vars->pKA = (PetscReal*)malloc(sizeof(PetscReal)*Nx*Ny*Nz);
    con_vars->pNMDA = (PetscReal*)malloc(sizeof(PetscReal)*Nx*Ny*Nz);

    con_vars->pKIR = (PetscReal*)malloc(sizeof(PetscReal)*Nx*Ny*Nz);
    //Glial diffusion scaling
    con_vars->DNeuronScale = (PetscReal*)malloc(sizeof(PetscReal)*3*Nx*Ny*Nz);
    con_vars->DGliaScale = (PetscReal*)malloc(sizeof(PetscReal)*3*Nx*Ny*Nz);
    con_vars->DExtracellScale = (PetscReal*)malloc(sizeof(PetscReal)*3*Nx*Ny*Nz);
    for(z=0;z<Nz;z++) {
        for (y = 0; y < Ny; y++) {
            for (x = 0; x < Nx; x++){

//                con_vars->pNaT[xy_index(x,y,z,Nx,Ny)] = basepNaT*(z==3);
//                con_vars->pKA[xy_index(x,y,z,Nx,Ny)] = basepKA*(z==3);
//                con_vars->pNMDA[xy_index(x,y,z,Nx,Ny)] = basepNMDA*(z<2);
//                con_vars->pNaP[xy_index(x,y,z,Nx,Ny)] = basepNaP*(z<4);
//                con_vars->pKDR[xy_index(x,y,z,Nx,Ny)] = basepKDR*(z<4);

//                con_vars->pNaT[xy_index(x,y,z,Nx,Ny)] = basepNaT;
//                con_vars->pKA[xy_index(x,y,z,Nx,Ny)] = basepKA;
//                con_vars->pNMDA[xy_index(x,y,z,Nx,Ny)] = basepNMDA*((double)x)/(Nx-1);
//                con_vars->pNaP[xy_index(x,y,z,Nx,Ny)] = basepNaP*((double)z)/(Nz-1);
//                con_vars->pKDR[xy_index(x,y,z,Nx,Ny)] = basepKDR;


                con_vars->pKIR[xy_index(x,y,z,Nx,Ny,Nz)] = basepKIR;

                con_vars->DNeuronScale[xy_index(x,y,z,Nx,Ny,Nz)*3] = DNeuronMult[0]; //x-direction Neurons
                con_vars->DNeuronScale[xy_index(x,y,z,Nx,Ny,Nz)*3+1] = DNeuronMult[1]; //y-direction Neurons
                con_vars->DNeuronScale[xy_index(x,y,z,Nx,Ny,Nz)*3+2] = DNeuronMult[2]; //z-direction Neurons
                con_vars->DGliaScale[xy_index(x,y,z,Nx,Ny,Nz)*3] = DGliaMult[0]; //x-direction scale Glia
                con_vars->DGliaScale[xy_index(x,y,z,Nx,Ny,Nz)*3+1] = DGliaMult[1]; // y-direction scale glia
                con_vars->DGliaScale[xy_index(x,y,z,Nx,Ny,Nz)*3+2] = DGliaMult[2]; // z-direction scale glia
                con_vars->DExtracellScale[xy_index(x,y,z,Nx,Ny,Nz)*3] = DExtraMult[0]; //x-direction scale extracell
                con_vars->DExtracellScale[xy_index(x,y,z,Nx,Ny,Nz)*3+1] = DExtraMult[1]; // y-direction scale Extracell
                con_vars->DExtracellScale[xy_index(x,y,z,Nx,Ny,Nz)*3+2] = DExtraMult[2]; // z-direction scale Extracell


                // Modification of neuron size
                     cmt = 0.75e-3;            //membrane capacitance in mF/cm^2

                if(z*dz>=soma && z*dz<soma+soma_width){
                    // In the soma
                    sa = 1.586e-5;          //membrane area in cm^2
                    voli = 2.16e-9;         //intracellular volume in cm^3
                    vole = (0.15*voli);
                    con_vars->pNaT[xy_index(x,y,z,Nx,Ny,Nz)] = basepNaT;
                    con_vars->pKA[xy_index(x,y,z,Nx,Ny,Nz)] = basepKA;
                    con_vars->pNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = basepNMDA*0;
                    con_vars->pNaP[xy_index(x,y,z,Nx,Ny,Nz)] = basepNaP;
                    con_vars->pKDR[xy_index(x,y,z,Nx,Ny,Nz)] = basepKDR;

                    con_vars->DNeuronScale[xy_index(x,y,z,Nx,Ny,Nz)*3] = DNeuronMult[0]*0; //x-direction Neurons
                    con_vars->DNeuronScale[xy_index(x,y,z,Nx,Ny,Nz)*3+1] = DNeuronMult[1]*0; //y-direction Neurons
                } else if(z*dz<soma){//else if(z>=dendrite && z<(dendrite+dendrite_width)){
                    // In the dendrite
                    sa = 16.408e-5;          //membrane area in cm^2
//                    sa = ((1.586e-5-16.408e-5)/soma)*z*dz+16.408e-5;          //membrane area in cm^2
                    voli = 3.852e-9;         //intracellular volume in cm^3
                    vole = (0.15*voli);
                    con_vars->pNaT[xy_index(x,y,z,Nx,Ny,Nz)] = basepNaT*0;//*pow((z*dz/Lz),3);
                    con_vars->pKA[xy_index(x,y,z,Nx,Ny,Nz)] = basepKA*0;//*(1+0.14*pow(1-(z*dz/Lz),1));
                    con_vars->pNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = basepNMDA;//*pow(1-(z*dz/Lz),1);
                    con_vars->pNaP[xy_index(x,y,z,Nx,Ny,Nz)] = basepNaP;
                    con_vars->pKDR[xy_index(x,y,z,Nx,Ny,Nz)] = basepKDR;

                    if(z*dz>=soma-0.0036){
                        con_vars->pNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = 0;
                    }

                    con_vars->DNeuronScale[xy_index(x,y,z,Nx,Ny,Nz)*3] = DNeuronMult[0]*0; //x-direction Neurons
                    con_vars->DNeuronScale[xy_index(x,y,z,Nx,Ny,Nz)*3+1] = DNeuronMult[1]*0; //y-direction Neurons

                    if(z==0){
                        //          con_vars->DNeuronScale[xy_index(x,y,0,Nx,Ny)*3] = DNeuronMult[0]*0.5; //x-direction Neurons
                        //         con_vars->DNeuronScale[xy_index(x,y,0,Nx,Ny)*3+1] = DNeuronMult[1]*0.5; //y-direction Neurons

                    }
                } else{
                    //Below the soma
                    sa = 10.324e-5;          //membrane area in cm^2
//                    sa = ((1.586e-5-10.324e-5)/(soma+soma_width-Lz))*(z*dz-Lz)+10.324e-5;          //membrane area in cm^2
                    voli = 1.762e-9;         //intracellular volume in cm^3
                    vole = (0.15*voli);
                    con_vars->pNaT[xy_index(x,y,z,Nx,Ny,Nz)] = basepNaT*0;
                    con_vars->pKA[xy_index(x,y,z,Nx,Ny,Nz)] = basepKA*0;
                    con_vars->pNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = basepNMDA*0;//*(z*dz-(soma+soma_width))/(Lz-(soma+soma_width));
                    con_vars->pNaP[xy_index(x,y,z,Nx,Ny,Nz)] = basepNaP*0;//*(1+.3*(z*dz-(soma+soma_width))/(Lz-(soma+soma_width)));
                    con_vars->pKDR[xy_index(x,y,z,Nx,Ny,Nz)] = basepKDR*0;//*(1+.5*(z*dz-(soma+soma_width))/(Lz-(soma+soma_width)));
                    //        con_vars->DNeuronScale[xy_index(x,y,z,Nx,Ny)*3] = DNeuronMult[0]*.1*(z*dz-(soma+soma_width))/(Lz-(soma+soma_width)); //x-direction Neurons
                    //        con_vars->DNeuronScale[xy_index(x,y,z,Nx,Ny)*3+1] = DNeuronMult[1]*.1*(z*dz-(soma+soma_width))/(Lz-(soma+soma_width)); //y-direction Neurons
                }

                if(x==0&&y==0){
                    printf("z:%d,sa:%.10e\n",z,sa);
                }


                    con_vars->ell[xy_index(x,y,z,Nx,Ny,Nz)] = ((voli+vole)/sa);    //average membrane separation in cm
                     con_vars->cm[al_index(x,y,z,0,Nx,Ny,Nz)] = cmt*RTFC/FC/con_vars->ell[xy_index(x,y,z,Nx,Ny,Nz)];
//                sa = 1.586e-5;          //membrane area in cm^2
//                voli = 2.16e-9;         //intracellular volume in cm^3
//                vole = (0.15*voli);
                     con_vars->cm[al_index(x,y,z,1,Nx,Ny,Nz)] = cmt*RTFC/FC/((voli+vole)/sa);//con_vars->ell[xy_index(x,y,z,Nx,Ny)];     //membrane capacitance in mF/cm^2 converted to mmol/cm^3
            }
        }
    }

    for(z=0;z<Nz;z++){
        printf("z:%d,NaT:%1.5e,NaP:%1.5e,",z,con_vars->pNaT[xy_index(0,0,z,Nx,Ny,Nz)],con_vars->pNaP[xy_index(0,0,z,Nx,Ny,Nz)]);
        printf("KA:%1.5e,KDR:%1.5e",con_vars->pKA[xy_index(0,0,z,Nx,Ny,Nz)],con_vars->pKDR[xy_index(0,0,z,Nx,Ny,Nz)]);
        printf("NMDA:%1.5e,D_N:%1.5e\n",con_vars->pNMDA[xy_index(0,0,z,Nx,Ny,Nz)],con_vars->DNeuronScale[3*
                                                                                                        xy_index(0,0,z,Nx,Ny,Nz)]);
    }



    //Variables that get set during set_params to steady state necessary values
    con_vars->pNaKCl = (PetscReal*)malloc(sizeof(PetscReal)*Nx*Ny*Nz);
    con_vars->Imax = (PetscReal*)malloc(sizeof(PetscReal)*Nx*Ny*Nz);
    con_vars->Imaxg = (PetscReal*)malloc(sizeof(PetscReal)*Nx*Ny*Nz);
    con_vars->pNaLeak = (PetscReal*)malloc(sizeof(PetscReal)*Nx*Ny*Nz);
    con_vars->pNaLeakg = (PetscReal*)malloc(sizeof(PetscReal)*Nx*Ny*Nz);

//    con_vars->zo = (PetscReal*)malloc(sizeof(PetscReal)*Nc*Nx*Ny*Nz);
//    con_vars->ao = (PetscReal*)malloc(sizeof(PetscReal)*Nc*Nx*Ny*Nz);
//    con_vars->zeta1 = (PetscReal*)malloc(sizeof(PetscReal)*(Nc-1)*Nx*Ny*Nz);


}


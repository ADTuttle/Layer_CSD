#include "constants.h"
#include "functions.h"

void Load_Grid(struct AppCtx *user,PetscInt xi,PetscInt yi){
    if(Profiling_on) {
        PetscLogEventBegin(event[12], 0, 0, 0, 0);
    }
    struct SimState *state_vars = user->state_vars_past;
    struct SimState * grid_vars = user->grid_vars;

    struct GateType *gate_vars = user->gate_vars_past;
    struct GateType *grid_gate = user->grid_gate_vars;

    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt Nz = user->Nz;
    PetscInt nx = 2*width_size+1;
    PetscInt ny = 2*width_size+1;
    PetscInt nz = Nz;
    PetscInt xind,yind,zind,ion,comp,x,y,z;

    for( z=0;z<nz;z++){
        for(y = 0; y < ny; y++){
            for(x = 0; x < nx; x++){
                xind = x-width_size+xi;
                yind = y-width_size+yi;
                zind = z;
                //If in interior just copy
                if(xind > -1 && yind > -1 && xind < Nx && yind < Ny){

                    grid_gate->mNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mNaT[xy_index(xind,yind,zind,Nx,Ny,Nz)];
                    grid_gate->hNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hNaT[xy_index(xind,yind,zind,Nx,Ny,Nz)];
                    grid_gate->gNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNaT[xy_index(xind,yind,zind,Nx,Ny,Nz)];
                    grid_gate->mNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mNaP[xy_index(xind,yind,zind,Nx,Ny,Nz)];
                    grid_gate->hNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hNaP[xy_index(xind,yind,zind,Nx,Ny,Nz)];
                    grid_gate->gNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNaP[xy_index(xind,yind,zind,Nx,Ny,Nz)];
                    grid_gate->mKDR[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mKDR[xy_index(xind,yind,zind,Nx,Ny,Nz)];
                    grid_gate->gKDR[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gKDR[xy_index(xind,yind,zind,Nx,Ny,Nz)];
                    grid_gate->mKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mKA[xy_index(xind,yind,zind,Nx,Ny,Nz)];
                    grid_gate->hKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hKA[xy_index(xind,yind,zind,Nx,Ny,Nz)];
                    grid_gate->gKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gKA[xy_index(xind,yind,zind,Nx,Ny,Nz)];
                    grid_gate->yNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->yNMDA[xy_index(xind,yind,zind,Nx,Ny,Nz)];
                    grid_gate->gNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNMDA[xy_index(xind,yind,zind,Nx,Ny,Nz)];

                    for(comp = 0; comp < Nc; comp++){
                        for(ion = 0; ion < Ni; ion++){
                            grid_vars->c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)] =
                                    state_vars->c[c_index(xind,yind,zind,comp,ion,Nx,Ny,Nz)];
                        }
                        grid_vars->phi[phi_index(x,y,z,comp,Nx,Ny,Nz)] =
                                state_vars->phi[phi_index(xind,yind,zind,comp,Nx,Ny,Nz)];
                    }
                    for(comp = 0; comp < Nc-1; comp++){
                        grid_vars->alpha[al_index(x,y,z,comp,Nx,Ny,Nz)] =
                                state_vars->alpha[al_index(xind,yind,zind,comp,Nx,Ny,Nz)];
                    }
                }else{
                    //If not the interior get closest point
                    //Right side
                    if(xind == Nx && yind > -1 && yind < Ny){
                        grid_gate->mNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mNaT[xy_index(xind-1,yind,zind,Nx,Ny,Nz)];
                        grid_gate->hNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hNaT[xy_index(xind-1,yind,zind,Nx,Ny,Nz)];
                        grid_gate->gNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNaT[xy_index(xind-1,yind,zind,Nx,Ny,Nz)];
                        grid_gate->mNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mNaP[xy_index(xind-1,yind,zind,Nx,Ny,Nz)];
                        grid_gate->hNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hNaP[xy_index(xind-1,yind,zind,Nx,Ny,Nz)];
                        grid_gate->gNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNaP[xy_index(xind-1,yind,zind,Nx,Ny,Nz)];
                        grid_gate->mKDR[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mKDR[xy_index(xind-1,yind,zind,Nx,Ny,Nz)];
                        grid_gate->gKDR[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gKDR[xy_index(xind-1,yind,zind,Nx,Ny,Nz)];
                        grid_gate->mKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mKA[xy_index(xind-1,yind,zind,Nx,Ny,Nz)];
                        grid_gate->hKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hKA[xy_index(xind-1,yind,zind,Nx,Ny,Nz)];
                        grid_gate->gKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gKA[xy_index(xind-1,yind,zind,Nx,Ny,Nz)];
                        grid_gate->yNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->yNMDA[xy_index(xind-1,yind,zind,Nx,Ny,Nz)];
                        grid_gate->gNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNMDA[xy_index(xind-1,yind,zind,Nx,Ny,Nz)];

                        for(comp = 0; comp < Nc; comp++){
                            for(ion = 0; ion < Ni; ion++){
                                grid_vars->c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)] =
                                        state_vars->c[c_index(xind-1,yind,zind,comp,ion,Nx,Ny,Nz)];
                            }
                            grid_vars->phi[phi_index(x,y,z,comp,Nx,Ny,Nz)] =
                                    state_vars->phi[phi_index(xind-1,yind,zind,comp,Nx,Ny,Nz)];
                        }
                        for(comp = 0; comp < Nc-1; comp++){
                            grid_vars->alpha[al_index(x,y,z,comp,Nx,Ny,Nz)] =
                                    state_vars->alpha[al_index(xind-1,yind,zind,comp,Nx,Ny,Nz)];
                        }
                    }
                    //Top side
                    if(yind == Ny && xind > -1 && xind < Nx){
                        grid_gate->mNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mNaT[xy_index(xind,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->hNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hNaT[xy_index(xind,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->gNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNaT[xy_index(xind,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->mNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mNaP[xy_index(xind,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->hNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hNaP[xy_index(xind,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->gNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNaP[xy_index(xind,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->mKDR[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mKDR[xy_index(xind,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->gKDR[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gKDR[xy_index(xind,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->mKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mKA[xy_index(xind,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->hKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hKA[xy_index(xind,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->gKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gKA[xy_index(xind,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->yNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->yNMDA[xy_index(xind,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->gNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNMDA[xy_index(xind,yind-1,z,Nx,Ny,Nz)];

                        for(comp = 0; comp < Nc; comp++){
                            for(ion = 0; ion < Ni; ion++){
                                grid_vars->c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)] =
                                        state_vars->c[c_index(xind,yind-1,z,comp,ion,Nx,Ny,Nz)];
                            }
                            grid_vars->phi[phi_index(x,y,z,comp,Nx,Ny,Nz)] =
                                    state_vars->phi[phi_index(xind,yind-1,z,comp,Nx,Ny,Nz)];
                        }
                        for(comp = 0; comp < Nc-1; comp++){
                            grid_vars->alpha[al_index(x,y,z,comp,Nx,Ny,Nz)] =
                                    state_vars->alpha[al_index(xind,yind-1,z,comp,Nx,Ny,Nz)];
                        }
                    }
                    //left side
                    if(xind == -1 && yind > -1 && yind < Ny){
                        grid_gate->mNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mNaT[xy_index(xind+1,yind,zind,Nx,Ny,Nz)];
                        grid_gate->hNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hNaT[xy_index(xind+1,yind,zind,Nx,Ny,Nz)];
                        grid_gate->gNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNaT[xy_index(xind+1,yind,zind,Nx,Ny,Nz)];
                        grid_gate->mNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mNaP[xy_index(xind+1,yind,zind,Nx,Ny,Nz)];
                        grid_gate->hNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hNaP[xy_index(xind+1,yind,zind,Nx,Ny,Nz)];
                        grid_gate->gNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNaP[xy_index(xind+1,yind,zind,Nx,Ny,Nz)];
                        grid_gate->mKDR[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mKDR[xy_index(xind+1,yind,zind,Nx,Ny,Nz)];
                        grid_gate->gKDR[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gKDR[xy_index(xind+1,yind,zind,Nx,Ny,Nz)];
                        grid_gate->mKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mKA[xy_index(xind+1,yind,zind,Nx,Ny,Nz)];
                        grid_gate->hKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hKA[xy_index(xind+1,yind,zind,Nx,Ny,Nz)];
                        grid_gate->gKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gKA[xy_index(xind+1,yind,zind,Nx,Ny,Nz)];
                        grid_gate->yNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->yNMDA[xy_index(xind+1,yind,zind,Nx,Ny,Nz)];
                        grid_gate->gNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNMDA[xy_index(xind+1,yind,zind,Nx,Ny,Nz)];

                        for(comp = 0; comp < Nc; comp++){
                            for(ion = 0; ion < Ni; ion++){
                                grid_vars->c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)] =
                                        state_vars->c[c_index(xind+1,yind,zind,comp,ion,Nx,Ny,Nz)];
                            }
                            grid_vars->phi[phi_index(x,y,z,comp,Nx,Ny,Nz)] =
                                    state_vars->phi[phi_index(xind+1,yind,zind,comp,Nx,Ny,Nz)];
                        }
                        for(comp = 0; comp < Nc-1; comp++){
                            grid_vars->alpha[al_index(x,y,z,comp,Nx,Ny,Nz)] =
                                    state_vars->alpha[al_index(xind+1,yind,zind,comp,Nx,Ny,Nz)];
                        }
                    }
                    //Bottom side
                    if(yind == -1 && xind > -1 && xind < Nx){
                        grid_gate->mNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mNaT[xy_index(xind,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->hNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hNaT[xy_index(xind,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->gNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNaT[xy_index(xind,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->mNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mNaP[xy_index(xind,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->hNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hNaP[xy_index(xind,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->gNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNaP[xy_index(xind,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->mKDR[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mKDR[xy_index(xind,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->gKDR[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gKDR[xy_index(xind,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->mKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mKA[xy_index(xind,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->hKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hKA[xy_index(xind,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->gKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gKA[xy_index(xind,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->yNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->yNMDA[xy_index(xind,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->gNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNMDA[xy_index(xind,yind+1,z,Nx,Ny,Nz)];

                        for(comp = 0; comp < Nc; comp++){
                            for(ion = 0; ion < Ni; ion++){
                                grid_vars->c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)] =
                                        state_vars->c[c_index(xind,yind+1,z,comp,ion,Nx,Ny,Nz)];
                            }
                            grid_vars->phi[phi_index(x,y,z,comp,Nx,Ny,Nz)] =
                                    state_vars->phi[phi_index(xind,yind+1,z,comp,Nx,Ny,Nz)];
                        }
                        for(comp = 0; comp < Nc-1; comp++){
                            grid_vars->alpha[al_index(x,y,z,comp,Nx,Ny,Nz)] =
                                    state_vars->alpha[al_index(xind,yind+1,z,comp,Nx,Ny,Nz)];
                        }
                    }
                    //Top Right corner
                    if(xind == Nx && yind == Ny){
                        grid_gate->mNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mNaT[xy_index(xind-1,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->hNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hNaT[xy_index(xind-1,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->gNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNaT[xy_index(xind-1,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->mNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mNaP[xy_index(xind-1,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->hNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hNaP[xy_index(xind-1,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->gNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNaP[xy_index(xind-1,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->mKDR[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mKDR[xy_index(xind-1,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->gKDR[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gKDR[xy_index(xind-1,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->mKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mKA[xy_index(xind-1,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->hKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hKA[xy_index(xind-1,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->gKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gKA[xy_index(xind-1,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->yNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->yNMDA[xy_index(xind-1,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->gNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNMDA[xy_index(xind-1,yind-1,z,Nx,Ny,Nz)];

                        for(comp = 0; comp < Nc; comp++){
                            for(ion = 0; ion < Ni; ion++){
                                grid_vars->c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)] =
                                        state_vars->c[c_index(xind-1,yind-1,z,comp,ion,Nx,Ny,Nz)];
                            }
                            grid_vars->phi[phi_index(x,y,z,comp,Nx,Ny,Nz)] =
                                    state_vars->phi[phi_index(xind-1,yind-1,z,comp,Nx,Ny,Nz)];
                        }
                        for(comp = 0; comp < Nc-1; comp++){
                            grid_vars->alpha[al_index(x,y,z,comp,Nx,Ny,Nz)] =
                                    state_vars->alpha[al_index(xind-1,yind-1,z,comp,Nx,Ny,Nz)];
                        }
                    }
                    // Top left corner
                    if(yind == Ny && xind == -1){
                        grid_gate->mNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mNaT[xy_index(xind+1,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->hNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hNaT[xy_index(xind+1,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->gNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNaT[xy_index(xind+1,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->mNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mNaP[xy_index(xind+1,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->hNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hNaP[xy_index(xind+1,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->gNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNaP[xy_index(xind+1,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->mKDR[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mKDR[xy_index(xind+1,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->gKDR[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gKDR[xy_index(xind+1,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->mKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mKA[xy_index(xind+1,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->hKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hKA[xy_index(xind+1,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->gKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gKA[xy_index(xind+1,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->yNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->yNMDA[xy_index(xind+1,yind-1,z,Nx,Ny,Nz)];
                        grid_gate->gNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNMDA[xy_index(xind+1,yind-1,z,Nx,Ny,Nz)];

                        for(comp = 0; comp < Nc; comp++){
                            for(ion = 0; ion < Ni; ion++){
                                grid_vars->c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)] =
                                        state_vars->c[c_index(xind+1,yind-1,z,comp,ion,Nx,Ny,Nz)];
                            }
                            grid_vars->phi[phi_index(x,y,z,comp,Nx,Ny,Nz)] =
                                    state_vars->phi[phi_index(xind+1,yind-1,z,comp,Nx,Ny,Nz)];
                        }
                        for(comp = 0; comp < Nc-1; comp++){
                            grid_vars->alpha[al_index(x,y,z,comp,Nx,Ny,Nz)] =
                                    state_vars->alpha[al_index(xind+1,yind-1,z,comp,Nx,Ny,Nz)];
                        }
                    }
                    //Bottom left
                    if(xind == -1 && yind == -1){
                        grid_gate->mNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mNaT[xy_index(xind+1,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->hNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hNaT[xy_index(xind+1,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->gNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNaT[xy_index(xind+1,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->mNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mNaP[xy_index(xind+1,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->hNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hNaP[xy_index(xind+1,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->gNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNaP[xy_index(xind+1,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->mKDR[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mKDR[xy_index(xind+1,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->gKDR[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gKDR[xy_index(xind+1,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->mKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mKA[xy_index(xind+1,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->hKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hKA[xy_index(xind+1,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->gKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gKA[xy_index(xind+1,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->yNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->yNMDA[xy_index(xind+1,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->gNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNMDA[xy_index(xind+1,yind+1,z,Nx,Ny,Nz)];

                        for(comp = 0; comp < Nc; comp++){
                            for(ion = 0; ion < Ni; ion++){
                                grid_vars->c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)] =
                                        state_vars->c[c_index(xind+1,yind+1,z,comp,ion,Nx,Ny,Nz)];
                            }
                            grid_vars->phi[phi_index(x,y,z,comp,Nx,Ny,Nz)] =
                                    state_vars->phi[phi_index(xind+1,yind+1,z,comp,Nx,Ny,Nz)];
                        }
                        for(comp = 0; comp < Nc-1; comp++){
                            grid_vars->alpha[al_index(x,y,z,comp,Nx,Ny,Nz)] =
                                    state_vars->alpha[al_index(xind+1,yind+1,z,comp,Nx,Ny,Nz)];
                        }
                    }
                    //Bottom right
                    if(xind == Nx && yind == -1){
                        grid_gate->mNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mNaT[xy_index(xind-1,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->hNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hNaT[xy_index(xind-1,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->gNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNaT[xy_index(xind-1,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->mNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mNaP[xy_index(xind-1,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->hNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hNaP[xy_index(xind-1,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->gNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNaP[xy_index(xind-1,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->mKDR[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mKDR[xy_index(xind-1,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->gKDR[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gKDR[xy_index(xind-1,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->mKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mKA[xy_index(xind-1,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->hKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hKA[xy_index(xind-1,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->gKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gKA[xy_index(xind-1,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->yNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->yNMDA[xy_index(xind-1,yind+1,z,Nx,Ny,Nz)];
                        grid_gate->gNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNMDA[xy_index(xind-1,yind+1,z,Nx,Ny,Nz)];

                        for(comp = 0; comp < Nc; comp++){
                            for(ion = 0; ion < Ni; ion++){
                                grid_vars->c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)] =
                                        state_vars->c[c_index(xind-1,yind+1,z,comp,ion,Nx,Ny,Nz)];
                            }
                            grid_vars->phi[phi_index(x,y,z,comp,Nx,Ny,Nz)] =
                                    state_vars->phi[phi_index(xind-1,yind+1,z,comp,Nx,Ny,Nz)];
                        }
                        for(comp = 0; comp < Nc-1; comp++){
                            grid_vars->alpha[al_index(x,y,z,comp,Nx,Ny,Nz)] =
                                    state_vars->alpha[al_index(xind-1,yind+1,z,comp,Nx,Ny,Nz)];
                        }
                    }


                }
            }
        }
    }
    if(Profiling_on) {
        PetscLogEventEnd(event[12], 0, 0, 0, 0);
    }
}

void Unload_Grid(struct AppCtx *user,PetscInt x,PetscInt y){
    if(Profiling_on) {
        PetscLogEventBegin(event[13], 0, 0, 0, 0);
    }
    PetscInt comp,ion;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt Nz = user->Nz;
    PetscInt nx = 2*width_size+1;
    PetscInt ny = 2*width_size+1;
    for(PetscInt z=0;z<Nz;z++){
        for(comp = 0; comp < Nc; comp++){
            for(ion = 0; ion < Ni; ion++){
                user->state_vars->c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)] =
                        user->grid_vars->c[c_index(width_size,width_size,z,comp,ion,Nx,Ny,Nz)];
            }
            user->state_vars->phi[phi_index(x,y,z,comp,Nx,Ny,Nz)] =
                    user->grid_vars->phi[phi_index(width_size,width_size,z,comp,Nx,Ny,Nz)];
        }
        //Save the held variable
        for(comp = 0; comp < Nc-1; comp++){
            user->state_vars->alpha[al_index(x,y,z,comp,Nx,Ny,Nz)] =
                    user->grid_vars->alpha[al_index(width_size,width_size,z,comp,Nx,Ny,Nz)];
        }
        //Save the gating variables

        user->gate_vars->mNaT[xy_index(x,y,z,Nx,Ny,Nz)] = user->grid_gate_vars->mNaT[xy_index(width_size,width_size,z,Nx,Ny,Nz)];
        user->gate_vars->hNaT[xy_index(x,y,z,Nx,Ny,Nz)] = user->grid_gate_vars->hNaT[xy_index(width_size,width_size,z,Nx,Ny,Nz)];
        user->gate_vars->gNaT[xy_index(x,y,z,Nx,Ny,Nz)] = user->grid_gate_vars->gNaT[xy_index(width_size,width_size,z,Nx,Ny,Nz)];
        user->gate_vars->mNaP[xy_index(x,y,z,Nx,Ny,Nz)] = user->grid_gate_vars->mNaP[xy_index(width_size,width_size,z,Nx,Ny,Nz)];
        user->gate_vars->hNaP[xy_index(x,y,z,Nx,Ny,Nz)] = user->grid_gate_vars->hNaP[xy_index(width_size,width_size,z,Nx,Ny,Nz)];
        user->gate_vars->gNaP[xy_index(x,y,z,Nx,Ny,Nz)] = user->grid_gate_vars->gNaP[xy_index(width_size,width_size,z,Nx,Ny,Nz)];
        user->gate_vars->mKDR[xy_index(x,y,z,Nx,Ny,Nz)] = user->grid_gate_vars->mKDR[xy_index(width_size,width_size,z,Nx,Ny,Nz)];
        user->gate_vars->gKDR[xy_index(x,y,z,Nx,Ny,Nz)] = user->grid_gate_vars->gKDR[xy_index(width_size,width_size,z,Nx,Ny,Nz)];
        user->gate_vars->mKA[xy_index(x,y,z,Nx,Ny,Nz)] = user->grid_gate_vars->mKA[xy_index(width_size,width_size,z,Nx,Ny,Nz)];
        user->gate_vars->hKA[xy_index(x,y,z,Nx,Ny,Nz)] = user->grid_gate_vars->hKA[xy_index(width_size,width_size,z,Nx,Ny,Nz)];
        user->gate_vars->gKA[xy_index(x,y,z,Nx,Ny,Nz)] = user->grid_gate_vars->gKA[xy_index(width_size,width_size,z,Nx,Ny,Nz)];
        user->gate_vars->yNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = user->grid_gate_vars->yNMDA[xy_index(width_size,width_size,z,Nx,Ny,Nz)];
        user->gate_vars->gNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = user->grid_gate_vars->gNMDA[xy_index(width_size,width_size,z,Nx,Ny,Nz)];

    }
    if(Profiling_on) {
        PetscLogEventEnd(event[13], 0, 0, 0, 0);
    }
}

PetscErrorCode Grid_Residual(Vec Res,PetscInt xi,PetscInt yi,void *ctx)
{
    //Residual equation using derivative of the charge-capacitance relation
    // Volume not solved for here
    struct AppCtx * user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if(Profiling_on) {
        PetscLogEventBegin(event[1], 0, 0, 0, 0);
    }
    //Compute membrane ionic flux relation quantitites
    grid_ionmflux(user,xi,yi);

    //Compute membrane water flow related quantities
    grid_wflowm(user,0,0);

    PetscReal *c = user->grid_vars->c;
    PetscReal *phi = user->grid_vars->phi;
    PetscReal *al = user->grid_vars->alpha;
    PetscReal *cp = user->grid_vars_past->c;
    PetscReal *alp = user->grid_vars_past->alpha;
    PetscReal *phip = user->grid_vars_past->phi;

    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct FluxData *flux = user->flux;
    PetscReal *cm = user->con_vars->cm;
    PetscReal dt = user->dt;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    PetscInt Nx = 2*width_size+1;
    PetscInt Ny = 2*width_size+1;

    //Residual for concentration equations
    PetscReal Rcvx,Rcvy,Resc;
    PetscReal RcvxRight,RcvyUp;

    //Residual for fluxes in voltage differential equations
    PetscReal Rphx[Nc], Rphy[Nc], RphxRight[Nc], RphyUp[Nc];
    PetscReal Resph,ResphN;

    PetscReal alNc,alpNc;
    PetscInt ion,comp,x,y,z;


    for(x=0;x<Nx;x++) {
        for(y=0;y<Ny;y++) {
            //Init voltage tracking to zero
            for(comp=0;comp<Nc;comp++) {
                Rphx[comp]=0;
                Rphy[comp]=0;
                RphxRight[comp]=0;
                RphyUp[comp]=0;
            }
            for(ion=0;ion<Ni;ion++) {
                for(comp=0;comp<Nc-1;comp++) {
                    Rcvx = 0;
                    RcvxRight = 0;
                    if(x>0) {
                        //First difference term
                        Rcvx = Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]*(cp[c_index(x-1,y,0,comp,ion,Nx,
                                                                                   0,0)]+cp[c_index(x,y,0,
                                                                                                    comp,ion,
                                                                                                    Nx,0,0)])/2;
                        Rcvx = Rcvx*(log(c[c_index(x,y,0,comp,ion,Nx,0,0)])-log(c[c_index(x-1,y,0,comp,
                                                                                            ion,Nx,0,0)])+z_charge[ion]*(phi[phi_index(
                                x,y,0,comp,Nx,0,0)]-phi[phi_index(
                                x-1,y,0,comp,Nx,0,0)]))/dx*dt/dx;
                    }
                    //Add Second right moving difference
                    if(x<Nx-1) {
                        RcvxRight = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]*(cp[c_index(x,y,0,comp,ion,Nx,
                                                                                      0,0)]+cp[c_index(x+1,y,
                                                                                                       0,comp,
                                                                                                       ion,Nx,
                                                                                                       0,0)])/2;
                        RcvxRight = RcvxRight*(log(c[c_index(x+1,y,0,comp,ion,Nx,0,0)])-log(c[c_index(x,y,0,
                                                                                                        comp,
                                                                                                        ion,Nx,
                                                                                                        0,0)])+z_charge[ion]*(phi[phi_index(
                                x+1,y,0,comp,Nx,0,0)]-phi[phi_index(
                                x,y,0,comp,Nx,0,0)]))/dx*dt/dx;
                    }
                    Rcvy = 0;
                    RcvyUp = 0;
                    //Up down difference
                    if(y>0) {
                        Rcvy = Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]*(cp[c_index(x,y-1,0,comp,ion,
                                                                                     Nx,0,0)]+cp[c_index(x,
                                                                                                         y,
                                                                                                         0,
                                                                                                         comp,
                                                                                                         ion,
                                                                                                         Nx,
                                                                                                         0,0)])/2;
                        Rcvy = Rcvy*(log(c[c_index(x,y,0,comp,ion,Nx,0,0)])-log(c[c_index(x,y-1,0,comp,
                                                                                            ion,Nx,0,0)])+z_charge[ion]*(phi[phi_index(
                                x,y,0,comp,Nx,0,0)]-phi[phi_index(
                                x,y-1,0,comp,Nx,0,0)]))/dy*dt/dy;
                    }
                    //Next upward difference
                    if(y<Ny-1) {
                        RcvyUp = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]*(cp[c_index(x,y,0,comp,ion,Nx,
                                                                                     0,0)]+cp[c_index(x,
                                                                                                      y+1,0,
                                                                                                      comp,
                                                                                                      ion,
                                                                                                      Nx,
                                                                                                      0,0)])/2;
                        RcvyUp = RcvyUp*(log(c[c_index(x,y+1,0,comp,ion,Nx,0,0)])-log(c[c_index(x,y,0,comp,
                                                                                                  ion,Nx,
                                                                                                  0,0)])+z_charge[ion]*(phi[phi_index(
                                x,y+1,0,comp,Nx,0,0)]-phi[phi_index(
                                x,y,0,comp,Nx,0,0)]))/dy*dt/dy;
                    }
                    Resc = al[al_index(x,y,0,comp,Nx,0,0)]*c[c_index(x,y,0,comp,ion,Nx,0,0)]-alp[al_index(x,
                                                                                                          y,
                                                                                                          0,
                                                                                                          comp,
                                                                                                          Nx,
                                                                                                          0,0)]*cp[c_index(
                            x,y,0,comp,ion,Nx,0,0)];
                    Resc += Rcvx - RcvxRight + Rcvy - RcvyUp +flux->mflux[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;

                    ierr = VecSetValue(Res,Ind_2(x,y,0,ion,comp,Nx,0,0),Resc,INSERT_VALUES);CHKERRQ(ierr);

                    //Save values for voltage
                    Rphx[comp]+=z_charge[ion]*Rcvx;
                    Rphy[comp]+=z_charge[ion]*Rcvy;
                    RphxRight[comp]+=z_charge[ion]*RcvxRight;
                    RphyUp[comp]+=z_charge[ion]*RcvyUp;

                }
                //Set Extracellular values
                alNc = 1-al[al_index(x,y,0,0,Nx,0,0)]-al[al_index(x,y,0,1,Nx,0,0)];
                alpNc = 1-alp[al_index(x,y,0,0,Nx,0,0)]-alp[al_index(x,y,0,1,Nx,0,0)];
                comp = Nc-1;
                Rcvx = 0;
                RcvxRight = 0;
                if(x>0) {
                    //First difference term
                    Rcvx = Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]*(cp[c_index(x-1,y,0,comp,ion,Nx,
                                                                               0,0)]+cp[c_index(x,y,0,comp,
                                                                                                ion,Nx,
                                                                                                0,0)])/2;
                    Rcvx = Rcvx*(log(c[c_index(x,y,0,comp,ion,Nx,0,0)])-log(c[c_index(x-1,y,0,comp,ion,
                                                                                        Nx,0,0)])+z_charge[ion]*(phi[phi_index(
                            x,y,0,comp,Nx,0,0)]-phi[phi_index(
                            x-1,y,0,comp,Nx,0,0)]))/dx*dt/dx;
                }
                //Add Second right moving difference
                if(x<Nx-1) {
                    RcvxRight = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]*(cp[c_index(x,y,0,comp,ion,Nx,0,0)]+cp[c_index(
                            x+1,y,0,
                            comp,ion,Nx,
                            0,0)])/2;
                    RcvxRight = RcvxRight*(log(c[c_index(x+1,y,0,comp,ion,Nx,0,0)])-log(c[c_index(x,y,0,
                                                                                                    comp,ion,
                                                                                                    Nx,0,0)])+z_charge[ion]*(phi[phi_index(
                            x+1,y,0,comp,Nx,0,0)]-phi[phi_index(
                            x,y,0,comp,Nx,0,0)]))/dx*dt/dx;
                }
                Rcvy = 0;
                RcvyUp = 0;
                //Up down difference
                if(y>0) {
                    Rcvy = Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]*(cp[c_index(x,y-1,0,comp,ion,Nx,
                                                                                 0,0)]+cp[c_index(x,y,0,
                                                                                                  comp,
                                                                                                  ion,Nx,
                                                                                                  0,0)])/2;
                    Rcvy = Rcvy*(log(c[c_index(x,y,0,comp,ion,Nx,0,0)])-log(c[c_index(x,y-1,0,comp,ion,
                                                                                        Nx,0,0)])+z_charge[ion]*(phi[phi_index(
                            x,y,0,comp,Nx,0,0)]-phi[phi_index(
                            x,y-1,0,comp,Nx,0,0)]))/dy*dt/dy;
                }
                //Next upward difference
                if(y<Ny-1) {
                    RcvyUp = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]*(cp[c_index(x,y,0,comp,ion,Nx,0,0)]+cp[c_index(
                            x,y+1,0,
                            comp,ion,Nx,
                            0,0)])/2;
                    RcvyUp = RcvyUp*(log(c[c_index(x,y+1,0,comp,ion,Nx,0,0)])-log(c[c_index(x,y,0,comp,
                                                                                              ion,Nx,0,0)])+z_charge[ion]*(phi[phi_index(
                            x,y+1,0,comp,Nx,0,0)]-phi[phi_index(
                            x,y,0,comp,Nx,0,0)]))/dy*dt/dy;
                }
                Resc = alNc*c[c_index(x,y,0,comp,ion,Nx,0,0)]-alpNc*cp[c_index(x,y,0,comp,ion,Nx,0,0)];
                Resc += Rcvx - RcvxRight + Rcvy - RcvyUp +flux->mflux[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                //Add bath variables

                Resc -= sqrt(pow(Dcb[c_index(x,y,0,comp,ion,Nx,0,0)*2],2)+pow(Dcb[c_index(x,y,0,comp,ion,
                                                                                          Nx,0,0)*2+1],2))*(cp[c_index(
                        x,y,0,comp,ion,Nx,0,0)]+cbath[ion])/2.0*(log(c[c_index(
                        x,y,0,comp,ion,Nx,0,0)])-log(cbath[ion])+z_charge[ion]*phi[phi_index(x,y,0,comp,Nx,0,0)]-z_charge[ion]*phibath)*dt;
                ierr = VecSetValue(Res,Ind_2(x,y,0,ion,comp,Nx,0,0),Resc,INSERT_VALUES);CHKERRQ(ierr);

                //Save values for voltage
                Rphx[comp]+=z_charge[ion]*Rcvx;
                Rphy[comp]+=z_charge[ion]*Rcvy;
                RphxRight[comp]+=z_charge[ion]*RcvxRight;
                RphyUp[comp]+=z_charge[ion]*RcvyUp;
            }

            //Voltage Equations
            ResphN = 0;
            for(comp=0;comp<Nc-1;comp++) {
                Resph = cm[al_index(x+xi,y+yi,0,comp,Nx,0,0)]*(phi[phi_index(x,y,0,comp,Nx,0,0)]-phi[phi_index(x,y,0,
                                                                                                             Nc-
                                                                                                             1,Nx,0,0)])-cm[al_index(
                        x+xi,y+yi,0,comp,Nx,0,0)]*(phip[phi_index(
                        x,
                        y,
                        0,
                        comp,
                        Nx,
                        0,0)]-phip[phi_index(
                        x,y,0,Nc-1,Nx,0,0)]);
                for(ion=0;ion<Ni;ion++){
                    //Ion channel
                    Resph += z_charge[ion]*flux->mflux[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                }
                //Add the terms shared with extracell
                ResphN -= Resph; // Subtract total capacitance, subtract total ion channel flux
                Resph += Rphx[comp] - RphxRight[comp] + Rphy[comp] - RphyUp[comp];
                ierr = VecSetValue(Res,Ind_2(x,y,0,Ni,comp,Nx,0,0),Resph,INSERT_VALUES); CHKERRQ(ierr);
            }

            //Finish adding extracell
            comp = Nc-1;
            //Add bath contribution
            for(ion=0;ion<Ni;ion++){

                ResphN -= z_charge[ion]*sqrt(pow(Dcb[c_index(x,y,0,comp,ion,Nx,0,0)*2],2)+pow(Dcb[c_index(x,y,0,
                                                                                                          comp,ion,
                                                                                                          Nx,0,0)*2+1],2))*(cp[c_index(
                        x,y,0,comp,ion,Nx,0,0)]+cbath[ion])/2.0*(log(c[c_index(
                        x,y,0,comp,ion,Nx,0,0)])-log(cbath[ion])+z_charge[ion]*phi[phi_index(x,y,0,comp,Nx,0,0)]-z_charge[ion]*phibath)*dt;
            }
            ResphN += Rphx[comp] - RphxRight[comp] + Rphy[comp] - RphyUp[comp];
            ierr = VecSetValue(Res,Ind_2(x,y,0,Ni,comp,Nx,0,0),ResphN,INSERT_VALUES); CHKERRQ(ierr);
        }
    }

    ierr = VecAssemblyBegin(Res);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Res);CHKERRQ(ierr);

    if(Profiling_on) {
        PetscLogEventEnd(event[1], 0, 0, 0, 0);
    }
    return ierr;
}

PetscErrorCode Grid_Jacobian(Mat Jac,PetscInt xi,PetscInt yi,void *ctx) {
    //Jacobian equation using derivative of the charge-capacitance relation
    // Alpha is not solved here

    struct AppCtx *user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if (Profiling_on) {
        PetscLogEventBegin(event[0], 0, 0, 0, 0);
    }
    PetscReal *c = user->grid_vars->c;
    PetscReal *al = user->grid_vars->alpha;
    PetscReal *cp = user->grid_vars_past->c;

    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct FluxData *flux = user->flux;
    PetscReal dt = user->dt;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    PetscInt Nx = 2*width_size+1;
    PetscInt Ny = 2*width_size+1;
    struct ConstVars *con_vars = user->con_vars;
    PetscReal * cm = user->con_vars->cm;

    PetscInt ind = 0;
    PetscInt x, y, ion, comp;

    PetscReal Ftmpx, Fc0x, Fc1x, Fph0x, Fph1x;
    PetscReal Ftmpy, Fc0y, Fc1y, Fph0y, Fph1y;
    PetscReal Ac, Aphi, Avolt, AvoltN;

    PetscReal Fphph0x[Nc], Fphph1x[Nc];
    PetscReal Fphph0y[Nc], Fphph1y[Nc];

    //Ionic concentration equations
    for (x = 0; x < Nx; x++) {
        for (y = 0; y < Ny; y++) {
            for (comp = 0; comp < Nc; comp++) {
                Fphph0x[comp] = 0;
                Fphph1x[comp] = 0;
                Fphph0y[comp] = 0;
                Fphph1y[comp] = 0;
            }
            for (ion = 0; ion < Ni; ion++) {
                for (comp = 0; comp < Nc - 1; comp++) {
                    //Electrodiffusion contributions
                    Ftmpx = 0;
                    Fc0x = 0;
                    Fc1x = 0;
                    Fph0x = 0;
                    Fph1x = 0;
                    Ftmpy = 0;
                    Fc0y = 0;
                    Fc1y = 0;
                    Fph0y = 0;
                    Fph1y = 0;
                    if (x < Nx - 1) {
                        Ftmpx = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]*
                                (cp[c_index(x,y,0,comp,ion,Nx,0,0)]+cp[c_index(x+1,y,0,comp,ion,Nx,0,0)])/2/dx*
                                dt/dx;
                        Fc0x = Ftmpx / c[c_index(x,y,0,comp,ion,Nx,0,0)];
                        Fph0x = z_charge[ion] * Ftmpx;
                        // Right c with left c (-Fc0x)

                        ierr = MatSetValue(Jac,Ind_2(x+1,y,0,ion,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),-Fc0x,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Right c with left phi (-Fph0x)
                        ierr = MatSetValue(Jac,Ind_2(x+1,y,0,ion,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),-Fph0x,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                        //Right phi with left c in voltage eqn
                        ierr = MatSetValue(Jac,Ind_2(x+1,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),
                                           -z_charge[ion] * Fc0x,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if (x > 0) {
                        Ftmpx = Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]*
                                (cp[c_index(x-1,y,0,comp,ion,Nx,0,0)]+cp[c_index(x,y,0,comp,ion,Nx,0,0)])/2/dx*
                                dt/dx;
                        Fc1x = Ftmpx / c[c_index(x,y,0,comp,ion,Nx,0,0)];
                        Fph1x = z_charge[ion] * Ftmpx;
                        //left c with right c (-Fc1x)
                        ierr = MatSetValue(Jac,Ind_2(x-1,y,0,ion,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),-Fc1x,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Left c with right phi (-Fph1x)
                        ierr = MatSetValue(Jac,Ind_2(x-1,y,0,ion,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),-Fph1x,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                        //Left phi with right c in voltage eqn
                        ierr = MatSetValue(Jac,Ind_2(x-1,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),
                                           -z_charge[ion] * Fc1x,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if (y < Ny - 1) {
                        Ftmpy = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]*
                                (cp[c_index(x,y,0,comp,ion,Nx,0,0)]+cp[c_index(x,y+1,0,comp,ion,Nx,0,0)])/2/dy*
                                dt/dy;
                        Fc0y = Ftmpy / c[c_index(x,y,0,comp,ion,Nx,0,0)];
                        Fph0y = z_charge[ion] * Ftmpy;
                        // Upper c with lower c (-Fc0y)
                        ierr = MatSetValue(Jac,Ind_2(x,y+1,0,ion,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),-Fc0y,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Upper c with lower phi (-Fph0y)
                        ierr = MatSetValue(Jac,Ind_2(x,y+1,0,ion,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),-Fph0y,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                        //Upper phi with lower c in voltage eqn
                        ierr = MatSetValue(Jac,Ind_2(x,y+1,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),
                                           -z_charge[ion] * Fc0y,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if (y > 0) {
                        Ftmpy = Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]*
                                (cp[c_index(x,y-1,0,comp,ion,Nx,0,0)]+cp[c_index(x,y,0,comp,ion,Nx,0,0)])/2/dy*
                                dt/dy;
                        Fc1y = Ftmpy / c[c_index(x,y,0,comp,ion,Nx,0,0)];
                        Fph1y = z_charge[ion] * Ftmpy;
                        //Lower c with Upper c (-Fc1y)
                        ierr = MatSetValue(Jac,Ind_2(x,y-1,0,ion,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),-Fc1y,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Lower c with Upper phi (-Fph1y)
                        ierr = MatSetValue(Jac,Ind_2(x,y-1,0,ion,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),-Fph1y,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                        //Lower phi with upper c in voltage eqn
                        ierr = MatSetValue(Jac,Ind_2(x,y-1,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),
                                           -z_charge[ion] * Fc1y,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    //Diagonal term contribution
                    Ac = al[al_index(x,y,0,comp,Nx,0,0)]+Fc0x+Fc1x+Fc0y+Fc1y;
                    Aphi = Fph0x + Fph1x + Fph0y + Fph1y;

                    //Add up terms for voltage eqns
                    Fphph0x[comp] += z_charge[ion] * Fph0x;
                    Fphph1x[comp] += z_charge[ion] * Fph1x;
                    Fphph0y[comp] += z_charge[ion] * Fph0y;
                    Fphph1y[comp] += z_charge[ion] * Fph1y;

                    //membrane current contributions
                    Ac += flux->dfdci[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                    Aphi += flux->dfdphim[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                    // Different Compartment Terms
                    // C Extracellular with C Inside
                    ierr = MatSetValue(Jac,Ind_2(x,y,0,ion,Nc-1,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),
                                       -flux->dfdci[c_index(x,y,0,comp,ion,Nx,0,0)]*dt,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    // C Intra with C Extra
                    ierr = MatSetValue(Jac,Ind_2(x,y,0,ion,comp,Nx,0,0),Ind_2(x,y,0,ion,Nc-1,Nx,0,0),
                                       flux->dfdce[c_index(x,y,0,comp,ion,Nx,0,0)]*dt,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    // C Extracellular with Phi Inside
                    ierr = MatSetValue(Jac,Ind_2(x,y,0,ion,Nc-1,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),
                                       -flux->dfdphim[c_index(x,y,0,comp,ion,Nx,0,0)]*dt,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    // C Intra with Phi Extra
                    ierr = MatSetValue(Jac,Ind_2(x,y,0,ion,comp,Nx,0,0),Ind_2(x,y,0,Ni,Nc-1,Nx,0,0),
                                       -flux->dfdphim[c_index(x,y,0,comp,ion,Nx,0,0)]*dt,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Same compartment terms
                    // c with c
                    ierr = MatSetValue(Jac,Ind_2(x,y,0,ion,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),Ac,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    // c with phi
                    ierr = MatSetValue(Jac,Ind_2(x,y,0,ion,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),Aphi,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;

                    //Intra-Phi with c (voltage eqn)
                    ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),z_charge[ion]*
                                                                                                    (Fc0x + Fc1x + Fc0y +
                                                                                                    Fc1y +
                                                                                                                flux->dfdci[c_index(
                                                                                                                        x,y,0,
                                                                                                                        comp,
                                                                                                                        ion,Nx,
                                                                                                                        0,0)]*
                                                                                                                dt),INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //IntraPhi with c extra(volt eqn)
                    ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,ion,Nc-1,Nx,0,0),
                                       z_charge[ion] * (flux->dfdce[c_index(x,y,0,comp,ion,Nx,0,0)]*dt),INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Extra-Phi with intra-c (voltage eqn)
                    ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,Nc-1,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),
                                       -z_charge[ion] * (flux->dfdci[c_index(x,y,0,comp,ion,Nx,0,0)]*dt),INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;

                }
                //Extracellular terms
                comp = Nc - 1;
                //Electrodiffusion contributions
                Ftmpx = 0;
                Fc0x = 0;
                Fc1x = 0;
                Fph0x = 0;
                Fph1x = 0;
                Ftmpy = 0;
                Fc0y = 0;
                Fc1y = 0;
                Fph0y = 0;
                Fph1y = 0;
                if (x < Nx - 1) {
                    Ftmpx = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]*
                            (cp[c_index(x,y,0,comp,ion,Nx,0,0)]+cp[c_index(x+1,y,0,comp,ion,Nx,0,0)])/2/dx*dt/
                            dx;
                    Fc0x = Ftmpx / c[c_index(x,y,0,comp,ion,Nx,0,0)];
                    Fph0x = z_charge[ion] * Ftmpx;
                    // Right c with left c (-Fc0x)
                    ierr = MatSetValue(Jac,Ind_2(x+1,y,0,ion,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),-Fc0x,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Right c with left phi (-Fph0x)
                    ierr = MatSetValue(Jac,Ind_2(x+1,y,0,ion,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),-Fph0x,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;

                    // Right Phi with left c (voltage eqn)
                    ierr = MatSetValue(Jac,Ind_2(x+1,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),-z_charge[ion]*Fc0x,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if (x > 0) {
                    Ftmpx = Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]*
                            (cp[c_index(x-1,y,0,comp,ion,Nx,0,0)]+cp[c_index(x,y,0,comp,ion,Nx,0,0)])/2/dx*dt/
                            dx;
                    Fc1x = Ftmpx / c[c_index(x,y,0,comp,ion,Nx,0,0)];
                    Fph1x = z_charge[ion] * Ftmpx;
                    //left c with right c (-Fc1x)
                    ierr = MatSetValue(Jac,Ind_2(x-1,y,0,ion,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),-Fc1x,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Left c with right phi (-Fph1x)
                    ierr = MatSetValue(Jac,Ind_2(x-1,y,0,ion,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),-Fph1x,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;

                    // left Phi with right c (voltage eqn)
                    ierr = MatSetValue(Jac,Ind_2(x-1,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),-z_charge[ion]*Fc1x,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if (y < Ny - 1) {
                    Ftmpy = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]*
                            (cp[c_index(x,y,0,comp,ion,Nx,0,0)]+cp[c_index(x,y+1,0,comp,ion,Nx,0,0)])/2/dy*dt/
                            dy;
                    Fc0y = Ftmpy / c[c_index(x,y,0,comp,ion,Nx,0,0)];
                    Fph0y = z_charge[ion] * Ftmpy;
                    // Upper c with lower c (-Fc0y)
                    ierr = MatSetValue(Jac,Ind_2(x,y+1,0,ion,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),-Fc0y,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Upper c with lower phi (-Fph0y)
                    ierr = MatSetValue(Jac,Ind_2(x,y+1,0,ion,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),-Fph0y,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;

                    // Upper Phi with lower c (voltage eqn)
                    ierr = MatSetValue(Jac,Ind_2(x,y+1,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),-z_charge[ion]*Fc0y,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if (y > 0) {
                    Ftmpy = Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]*
                            (cp[c_index(x,y-1,0,comp,ion,Nx,0,0)]+cp[c_index(x,y,0,comp,ion,Nx,0,0)])/2/dy*dt/
                            dy;
                    Fc1y = Ftmpy / c[c_index(x,y,0,comp,ion,Nx,0,0)];
                    Fph1y = z_charge[ion] * Ftmpy;
                    //Lower c with Upper c (-Fc1y)
                    ierr = MatSetValue(Jac,Ind_2(x,y-1,0,ion,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),-Fc1y,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Lower c with Upper phi (-Fph1y)
                    ierr = MatSetValue(Jac,Ind_2(x,y-1,0,ion,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),-Fph1y,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;

                    // Lower Phi with upper c (voltage eqn)
                    ierr = MatSetValue(Jac,Ind_2(x,y-1,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),-z_charge[ion]*Fc1y,
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }

                //Diagonal term contribution
                Ac = (1-al[al_index(x,y,0,0,Nx,0,0)]-al[al_index(x,y,0,1,Nx,0,0)])+Fc0x+Fc1x+Fc0y+Fc1y;
                Aphi = Fph0x + Fph1x + Fph0y + Fph1y;

                Avolt = z_charge[ion] * (Fc0x + Fc1x + Fc0y + Fc1y);

                //Add up terms for voltage eqns
                Fphph0x[comp] += z_charge[ion] * Fph0x;
                Fphph1x[comp] += z_charge[ion] * Fph1x;
                Fphph0y[comp] += z_charge[ion] * Fph0y;
                Fphph1y[comp] += z_charge[ion] * Fph1y;

                //Membrane current contribution
                for (comp = 0; comp < Nc - 1; comp++) {
                    Ac -= flux->dfdce[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                    Aphi += flux->dfdphim[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                    Avolt -= z_charge[ion]*flux->dfdce[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                }
                //Add bath contributions
                Ftmpx = sqrt(pow(Dcb[c_index(x,y,0,Nc-1,ion,Nx,0,0)*2],2)+
                             pow(Dcb[c_index(x,y,0,Nc-1,ion,Nx,0,0)*2+1],2));
                Ac -= Ftmpx*(cp[c_index(x,y,0,Nc-1,ion,Nx,0,0)]+cbath[ion])/
                      (2 * c[c_index(x,y,0,Nc-1,ion,Nx,0,0)])*dt;
                Aphi -= Ftmpx*(cp[c_index(x,y,0,Nc-1,ion,Nx,0,0)]+cbath[ion])*z_charge[ion]/2*dt;

                Avolt -= z_charge[ion]*Ftmpx*(cp[c_index(x,y,0,Nc-1,ion,Nx,0,0)]+cbath[ion])/
                         (2 * c[c_index(x,y,0,Nc-1,ion,Nx,0,0)])*dt;

                //Insert extracell to extracell parts
                // c with c
                ierr = MatSetValue(Jac,Ind_2(x,y,0,ion,Nc-1,Nx,0,0),Ind_2(x,y,0,ion,Nc-1,Nx,0,0),Ac,INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
                // c with phi
                ierr = MatSetValue(Jac,Ind_2(x,y,0,ion,Nc-1,Nx,0,0),Ind_2(x,y,0,Ni,Nc-1,Nx,0,0),Aphi,INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;

                //phi with c (voltage eqn)
                ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,Nc-1,Nx,0,0),Ind_2(x,y,0,ion,Nc-1,Nx,0,0),Avolt,
                                   INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
            }
            //Derivative of charge-capacitance
            for (comp = 0; comp < Nc - 1; comp++) {
                if (x < Nx - 1) {
                    //Right phi with left phi (-Fph0x)
                    ierr = MatSetValue(Jac,Ind_2(x+1,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),-Fphph0x[comp],
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if (x > 0) {
                    //Left phi with right phi (-Fph1x)
                    ierr = MatSetValue(Jac,Ind_2(x-1,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),-Fphph1x[comp],
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if (y < Ny - 1) {
                    //Upper phi with lower phi (-Fph0y)
                    ierr = MatSetValue(Jac,Ind_2(x,y+1,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),-Fphph0y[comp],
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if (y > 0) {
                    //Lower phi with upper phi (-Fph1y)
                    ierr = MatSetValue(Jac,Ind_2(x,y-1,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),-Fphph1y[comp],
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                Avolt = cm[al_index(x+xi,y+yi,0,comp,Nx,0,0)]+Fphph0x[comp]+Fphph1x[comp]+Fphph0y[comp]+Fphph1y[comp];
                AvoltN = -cm[al_index(x+xi,y+yi,0,comp,Nx,0,0)];
                for (ion = 0; ion < Ni; ion++) {
                    Avolt += z_charge[ion]*flux->dfdphim[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                    AvoltN -= z_charge[ion]*flux->dfdphim[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                }

                //Intra-phi with Intra-phi
                ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),Avolt,INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
                //Intra-phi with extra-phi
                ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,Nc-1,Nx,0,0),AvoltN,INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
            }
            //Extracellular terms
            comp = Nc - 1;
            if (x < Nx - 1) {
                //Right phi with left phi (-Fph0x)
                ierr = MatSetValue(Jac,Ind_2(x+1,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),-Fphph0x[comp],
                                   INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
            }
            if (x > 0) {
                //Left phi with right phi (-Fph1x)
                ierr = MatSetValue(Jac,Ind_2(x-1,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),-Fphph1x[comp],
                                   INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
            }
            if (y < Ny - 1) {
                //Upper phi with lower phi (-Fph0y)
                ierr = MatSetValue(Jac,Ind_2(x,y+1,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),-Fphph0y[comp],
                                   INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
            }
            if (y > 0) {
                //Lower phi with upper phi (-Fph1y)
                ierr = MatSetValue(Jac,Ind_2(x,y-1,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),-Fphph1y[comp],
                                   INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
            }
            AvoltN = 0;

            for (int k = 0; k < Nc - 1; k++) {
                AvoltN += cm[al_index(x+xi,y+yi,0,comp,Nx,0,0)];
                Avolt = -cm[al_index(x+xi,y+yi,0,comp,Nx,0,0)];
                for (ion = 0; ion < Ni; ion++) {
                    Avolt -= z_charge[ion]*flux->dfdphim[c_index(x,y,0,k,ion,Nx,0,0)]*dt;
                    AvoltN += z_charge[ion]*flux->dfdphim[c_index(x,y,0,k,ion,Nx,0,0)]*dt;
                }
                //Extra-phi with Intra-phi
                ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,k,Nx,0,0),Avolt,INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
            }

            AvoltN += Fphph0x[comp] + Fphph1x[comp] + Fphph0y[comp] + Fphph1y[comp];

            //Bath terms
            for (ion = 0; ion < Ni; ion++) {
                Ftmpx = sqrt(pow(Dcb[c_index(x,y,0,Nc-1,ion,Nx,0,0)*2],2)+
                             pow(Dcb[c_index(x,y,0,Nc-1,ion,Nx,0,0)*2+1],2));
                AvoltN -= z_charge[ion]*Ftmpx*(cp[c_index(x,y,0,Nc-1,ion,Nx,0,0)]+cbath[ion])*z_charge[ion]/2*dt;
            }
            //extra-phi with extra-phi
            ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),AvoltN,INSERT_VALUES);
            CHKERRQ(ierr);
            ind++;

        }
    }

    ierr = MatAssemblyBegin(Jac, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Jac, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    if (Profiling_on) {
        PetscLogEventEnd(event[0], 0, 0, 0, 0);
    }
    return ierr;

}

PetscErrorCode Grid_Residual_algebraic(Vec Res,PetscInt xi,PetscInt yi,void *ctx)
{
    //Residual equation using algebraic version of the charge-capacitance relation
    //Alpha is solved for here
    struct AppCtx * user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if(Profiling_on) {
        PetscLogEventBegin(event[10], 0, 0, 0, 0);
    }
    //Compute membrane ionic flux relation quantitites
    grid_ionmflux(user,xi,yi);

    //Compute membrane water flow related quantities
    grid_wflowm(user,0,0);

    PetscReal *c = user->grid_vars->c;
    PetscReal *phi = user->grid_vars->phi;
    PetscReal *al = user->grid_vars->alpha;
    PetscReal *cp = user->grid_vars_past->c;
    PetscReal *alp = user->grid_vars_past->alpha;

    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct FluxData *flux = user->flux;
    PetscReal *cm = user->con_vars->cm;
    PetscReal dt = user->dt;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    PetscReal dz = user->dz;
    PetscInt Nx = 2*width_size+1;
    PetscInt Ny = 2*width_size+1;
    PetscInt Nz = user->Nz;

    //Residual for concentration equations
    PetscReal Rcvx,Rcvy,Rcvz,Resc;
    PetscReal RcvxRight,RcvyUp,RcvzTop;

    PetscReal alNc,alpNc;
    PetscInt ion,comp,x,y,z;

    for(z=0;z<Nz;z++){
        for(y = 0; y < Ny; y++){
            for(x = 0; x < Nx; x++){
                for(ion = 0; ion < Ni; ion++){
                    for(comp = 0; comp < Nc-1; comp++){
                        Rcvx = 0;
                        RcvxRight = 0;
                        if(x > 0){
                            //First difference term
                            Rcvx = Dcs[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)*3]*(cp[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                             cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2;
                            Rcvx = Rcvx*(log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-
                                         log(c[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)])+z_charge[ion]*
                                                                                   (phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-phi[phi_index(
                                                                                           x-1,y,z,comp,Nx,Ny,Nz)]))/dx*dt/dx;
                        }
                        //Add Second right moving difference
                        if(x < Nx-1){
                            RcvxRight = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                                cp[c_index(x+1,y,z,comp,ion,Nx,Ny,Nz)])/2;
                            RcvxRight = RcvxRight*(log(c[c_index(x+1,y,z,comp,ion,Nx,Ny,Nz)])-
                                                   log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])+z_charge[ion]*
                                                                                           (phi[phi_index(
                                                                                                   x+1,y,z,comp,Nx,Ny,Nz)]-phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]))/dx*dt/dx;
                        }
                        Rcvy = 0;
                        RcvyUp = 0;
                        //Up down difference
                        if(y > 0){
                            Rcvy = Dcs[c_index(x,y-1,z,comp,ion,Nx,Ny,Nz)*3+1]*(cp[c_index(x,y-1,z,comp,ion,Nx,Ny,Nz)]+
                                                                               cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2;
                            Rcvy = Rcvy*(log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-log(c[c_index(x,y-1,z,comp,ion,Nx,Ny,Nz)])+
                                    z_charge[ion]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-
                                                   phi[phi_index(x,y-1,z,comp,Nx,Ny,Nz)]))/dy*dt/dy;
                        }
                        //Next upward difference
                        if(y < Ny-1){
                            RcvyUp = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+1]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                               cp[c_index(x,y+1,z,comp,ion,Nx,Ny,Nz)])/2;
                            RcvyUp = RcvyUp*(log(c[c_index(x,y+1,z,comp,ion,Nx,Ny,Nz)])-
                                             log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])+z_charge[ion]*
                                                                                     (phi[phi_index(x,
                                                                                                    y+1,z,comp,Nx,Ny,Nz)]-phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]))/dy*dt/dy;
                        }
                        Rcvz = 0;
                        RcvzTop = 0;
                        if(z > 0){
                            //First difference term
                            Rcvz = Dcs[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)]+
                                                                               cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2;
                            Rcvz = Rcvz*(log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-
                                         log(c[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)])+z_charge[ion]*
                                                                                   (phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-phi[phi_index(x,y,
                                                                                                                                     z-
                                                                                                                                     1,comp,Nx,Ny,Nz)]))/dz*dt/dz;
                        }
                        //Add Second right moving difference
                        if(z < Nx-1){
                            RcvzTop = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                                cp[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])/2;
                            RcvzTop = RcvzTop*(log(c[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])-
                                               log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])+z_charge[ion]*
                                                                                       (phi[phi_index(x,y,
                                                                                                      z+1,comp,Nx,Ny,Nz)]-phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]))/dz*dt/dz;
                        }
                        Resc = al[al_index(x,y,z,comp,Nx,Ny,Nz)]*c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]-
                               alp[al_index(x,y,z,comp,Nx,Ny,Nz)]*cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                        Resc += Rcvx-RcvxRight+Rcvy-RcvyUp+Rcvz-RcvzTop+flux->mflux[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;

                        ierr = VecSetValue(Res,Ind_2(x,y,z,ion,comp,Nx,Ny,Nz),Resc,INSERT_VALUES);
                        CHKERRQ(ierr);

                    }
                    //Set Extracellular values
                    alNc = 1-al[al_index(x,y,z,0,Nx,Ny,Nz)]-al[al_index(x,y,z,1,Nx,Ny,Nz)];
                    alpNc = 1-alp[al_index(x,y,z,0,Nx,Ny,Nz)]-alp[al_index(x,y,z,1,Nx,Ny,Nz)];
                    comp = Nc-1;
                    Rcvx = 0;
                    RcvxRight = 0;
                    if(x > 0){
                        //First difference term
                        Rcvx = Dcs[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)*3]*(cp[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                         cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2;
                        Rcvx = Rcvx*(log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-log(c[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)])+
                                z_charge[ion]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-phi[phi_index(
                                        x-1,y,z,comp,Nx,Ny,Nz)]))/dx*dt/dx;
                    }
                    //Add Second right moving difference
                    if(x < Nx-1){
                        RcvxRight = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                            cp[c_index(x+1,y,z,comp,ion,Nx,Ny,Nz)])/2;
                        RcvxRight = RcvxRight*(log(c[c_index(x+1,y,z,comp,ion,Nx,Ny,Nz)])-
                                               log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])+z_charge[ion]*
                                                                                       (phi[phi_index(
                                                                                               x+1,y,z,comp,Nx,Ny,Nz)]-phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]))/dx*dt/dx;
                    }
                    Rcvy = 0;
                    RcvyUp = 0;
                    //Up down difference
                    if(y > 0){
                        Rcvy = Dcs[c_index(x,y-1,z,comp,ion,Nx,Ny,Nz)*3+1]*(cp[c_index(x,y-1,z,comp,ion,Nx,Ny,Nz)]+
                                                                           cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2;
                        Rcvy = Rcvy*(log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-log(c[c_index(x,y-1,z,comp,ion,Nx,Ny,Nz)])+
                                z_charge[ion]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-phi[phi_index(x,
                                                                                                y-1,z,comp,Nx,Ny,Nz)]))/dy*dt/dy;
                    }
                    //Next upward difference
                    if(y < Ny-1){
                        RcvyUp = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+1]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                           cp[c_index(x,y+1,z,comp,ion,Nx,Ny,Nz)])/2;
                        RcvyUp = RcvyUp*(log(c[c_index(x,y+1,z,comp,ion,Nx,Ny,Nz)])-log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])+
                                z_charge[ion]*(phi[phi_index(x,y+1,z,comp,Nx,Ny,Nz)]-
                                               phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]))/dy*dt/dy;
                    }
                    Rcvz = 0;
                    RcvzTop = 0;
                    if(z > 0){
                        //First difference term
                        Rcvz = Dcs[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)]+
                                                                           cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2;
                        Rcvz = Rcvz*(log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-
                                     log(c[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)])+z_charge[ion]*
                                                                               (phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-phi[phi_index(x,y,
                                                                                                                                 z-
                                                                                                                                 1,comp,Nx,Ny,Nz)]))/dz*dt/dz;
                    }
                    //Add Second right moving difference
                    if(z < Nx-1){
                        RcvzTop = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                            cp[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])/2;
                        RcvzTop = RcvzTop*(log(c[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])-
                                           log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])+z_charge[ion]*
                                                                                   (phi[phi_index(x,y,
                                                                                                  z+1,comp,Nx,Ny,Nz)]-phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]))/dz*dt/dz;
                    }
                    Resc = alNc*c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]-alpNc*cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                    Resc += Rcvx-RcvxRight+Rcvy-RcvyUp+Rcvz-RcvzTop+flux->mflux[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;
                    //Add bath variables

                    Resc -= sqrt(pow(Dcb[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3],2)+
                                 pow(Dcb[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+1],2)+
                                 pow(Dcb[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2],2))*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                                   cbath[ion])/2.0*(log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-log(cbath[ion])+
                                            z_charge[ion]*phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-z_charge[ion]*phibath)*dt;
                    ierr = VecSetValue(Res,Ind_2(x,y,z,ion,comp,Nx,Ny,Nz),Resc,INSERT_VALUES);
                    CHKERRQ(ierr);

                }
            }
        }
    }

    for(z=0;z<Nz;z++){
        for(y = 0; y < Ny; y++){
            for(x = 0; x < Nx; x++){
                //Residual for electroneutrality condition
                for(comp = 0; comp < Nc-1; comp++){
                    Resc = al[al_index(x,y,z,comp,Nx,Ny,Nz)]*cz(c,z_charge,x,y,z,Nx,Ny,Nz,comp,user)+
                           user->con_vars->zo[phi_index(x,y,z,comp,Nx,Ny,Nz)]*
                           user->con_vars->ao[phi_index(x,y,z,comp,Nx,Ny,Nz)];
                    ierr = VecSetValue(Res,Ind_2(x,y,z,Ni,comp,Nx,Ny,Nz),Resc,INSERT_VALUES);
                    CHKERRQ(ierr);
                }
                //Extracellular term
                comp = Nc-1;
                Resc = (1-al[al_index(x,y,z,0,Nx,Ny,Nz)]-al[al_index(x,y,z,1,Nx,Ny,Nz)])*
                       cz(c,z_charge,x,y,z,Nx,Ny,Nz,comp,user)+
                       user->con_vars->zo[phi_index(x,y,z,comp,Nx,Ny,Nz)]*user->con_vars->ao[phi_index(x,y,z,comp,Nx,Ny,Nz)];
                ierr = VecSetValue(Res,Ind_2(x,y,z,Ni,comp,Nx,Ny,Nz),Resc,INSERT_VALUES);
                CHKERRQ(ierr);

                //Residual for water flow
                //Plus modification to electroneutrality for non-zero mem.compacitance
                for(comp = 0; comp < Nc-1; comp++){
                    //Water flow
                    ierr = VecSetValue(Res,Ind_2(x,y,z,Ni+1,comp,Nx,Ny,Nz),al[al_index(x,y,z,comp,Nx,Ny,Nz)]-
                                                                          alp[al_index(x,y,z,comp,Nx,Ny,Nz)]+flux->wflow[al_index(x,y,z,comp,Nx,Ny,Nz)]*dt,INSERT_VALUES);
                    CHKERRQ(ierr);

                }
            }
        }
    }
    //Assemble before we add values in on top to modify the electroneutral.
    ierr = VecAssemblyBegin(Res);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Res);CHKERRQ(ierr);
    for(z=0;z<Nz;z++){
        for(y = 0; y < Ny; y++){
            for(x = 0; x < Nx; x++){
                // Add Modification to electroneutrality for non-zero mem.compacitance
                for(comp = 0; comp < Nc-1; comp++){
                    //Extracell voltage
                    ierr = VecSetValue(Res,Ind_2(x,y,z,Ni,Nc-1,Nx,Ny,Nz),-cm[al_index(x+xi,y+yi,z,comp,Nx,Ny,Nz)]*(phi[phi_index(x,y,z,
                                                                                                                           Nc-
                                                                                                                           1,Nx,Ny,Nz)]-
                                                                                                                  phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]),ADD_VALUES);
                    CHKERRQ(ierr);
                    //Intracell voltage mod
                    ierr = VecSetValue(Res,Ind_2(x,y,z,Ni,comp,Nx,Ny,Nz),-cm[al_index(x+xi,y+yi,z,comp,Nx,Ny,Nz)]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-
                                                                                                                  phi[phi_index(x,y,z,
                                                                                                                           Nc-1,Nx,Ny,Nz)]),ADD_VALUES);
                    CHKERRQ(ierr);
                }
            }
        }
    }

    ierr = VecAssemblyBegin(Res);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Res);CHKERRQ(ierr);
    if(Profiling_on) {
        PetscLogEventEnd(event[10], 0, 0, 0, 0);
    }
    return ierr;
}

PetscErrorCode
Grid_Jacobian_algebraic(Mat Jac,PetscInt xi, PetscInt yi,void *ctx)
{
    //Jacobian equation using algebraic version of the charge-capacitance relation
    // Alpha is solved for here
    struct AppCtx * user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if(Profiling_on) {
        PetscLogEventBegin(event[9], 0, 0, 0, 0);
    }
    PetscReal *c = user->grid_vars->c;
    PetscReal *al = user->grid_vars->alpha;
    PetscReal *cp = user->grid_vars_past->c;

    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct FluxData *flux = user->flux;
    PetscReal dt = user->dt;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    PetscReal dz = user->dz;
    PetscInt Nx = 2*width_size+1;
    PetscInt Ny = 2*width_size+1;
    PetscInt Nz = user->Nz;
    struct ConstVars *con_vars = user->con_vars;
    PetscReal *cm = user->con_vars->cm;

    PetscInt ind = 0;
    PetscInt x,y,z,ion,comp;

    PetscReal Ftmpx,Fc0x,Fc1x,Fph0x,Fph1x;
    PetscReal Ftmpy,Fc0y,Fc1y,Fph0y,Fph1y;
    PetscReal Ftmpz,Fc0z,Fc1z,Fph0z,Fph1z;
    PetscReal Ac,Aphi;


    //Ionic concentration equations
    for(z=0;z<Nz;z++){
        for(y = 0; y < Ny; y++){
            for(x = 0; x < Nx; x++){
                for(ion = 0; ion < Ni; ion++){
                    for(comp = 0; comp < Nc-1; comp++){
                        //Electrodiffusion contributions
                        Ftmpx = 0;Fc0x = 0;Fc1x = 0;Fph0x = 0;Fph1x = 0;
                        Ftmpy = 0;Fc0y = 0;Fc1y = 0;Fph0y = 0;Fph1y = 0;
                        Ftmpz = 0;Fc0z = 0;Fc1z = 0;Fph0z = 0;Fph1z = 0;
                        if(x < Nx-1){
                            Ftmpx = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                            cp[c_index(x+1,y,z,comp,ion,Nx,Ny,Nz)])/2/dx*dt/dx;
                            Fc0x = Ftmpx/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                            Fph0x = z_charge[ion]*Ftmpx;
                            // Right c with left c (-Fc0x)

                            ierr = MatSetValue(Jac,Ind_2(x+1,y,z,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc0x,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Right c with left phi (-Fph0x)
                            ierr = MatSetValue(Jac,Ind_2(x+1,y,z,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,Ni,comp,Nx,Ny,Nz),
                                               -Fph0x,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;

                        }
                        if(x > 0){
                            Ftmpx = Dcs[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)*3]*(cp[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                              cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2/dx*dt/dx;
                            Fc1x = Ftmpx/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                            Fph1x = z_charge[ion]*Ftmpx;
                            //left c with right c (-Fc1x)
                            ierr = MatSetValue(Jac,Ind_2(x-1,y,z,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc1x,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Left c with right phi (-Fph1x)
                            ierr = MatSetValue(Jac,Ind_2(x-1,y,z,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,Ni,comp,Nx,Ny,Nz),
                                               -Fph1x,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                        if(y < Ny-1){
                            Ftmpy = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+1]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                              cp[c_index(x,y+1,z,comp,ion,Nx,Ny,Nz)])/2/dy*dt/dy;
                            Fc0y = Ftmpy/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                            Fph0y = z_charge[ion]*Ftmpy;
                            // Upper c with lower c (-Fc0y)
                            ierr = MatSetValue(Jac,Ind_2(x,y+1,z,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc0y,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Upper c with lower phi (-Fph0y)
                            ierr = MatSetValue(Jac,Ind_2(x,y+1,z,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,Ni,comp,Nx,Ny,Nz),
                                               -Fph0y,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                        if(y > 0){
                            Ftmpy = Dcs[c_index(x,y-1,z,comp,ion,Nx,Ny,Nz)*3+1]*(cp[c_index(x,y-1,z,comp,ion,Nx,Ny,Nz)]+
                                                                                cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2/dy*dt/dy;
                            Fc1y = Ftmpy/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                            Fph1y = z_charge[ion]*Ftmpy;
                            //Lower c with Upper c (-Fc1y)
                            ierr = MatSetValue(Jac,Ind_2(x,y-1,z,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc1y,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Lower c with Upper phi (-Fph1y)
                            ierr = MatSetValue(Jac,Ind_2(x,y-1,z,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,Ni,comp,Nx,Ny,Nz),
                                               -Fph1y,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                        if(z < Nz-1){
                            Ftmpz = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                              cp[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])/2/dz*dt/dz;
                            Fc0z = Ftmpz/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                            Fph0z = z_charge[ion]*Ftmpz;
                            // Upper c with lower c (-Fc0y)
                            ierr = MatSetValue(Jac,Ind_2(x,y,z+1,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc0z,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Upper c with lower phi (-Fph0y)
                            ierr = MatSetValue(Jac,Ind_2(x,y,z+1,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,Ni,comp,Nx,Ny,Nz),
                                               -Fph0z,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                        if(z > 0){
                            Ftmpz = Dcs[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)]+
                                                                                cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2/dz*dt/dz;
                            Fc1z = Ftmpz/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                            Fph1z = z_charge[ion]*Ftmpz;
                            //Lower c with Upper c (-Fc1y)
                            ierr = MatSetValue(Jac,Ind_2(x,y,z-1,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc1z,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Lower c with Upper phi (-Fph1y)
                            ierr = MatSetValue(Jac,Ind_2(x,y,z-1,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,Ni,comp,Nx,Ny,Nz),
                                               -Fph1z,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                        //Diagonal term contribution
                        Ac = al[al_index(x,y,z,comp,Nx,Ny,Nz)]+Fc0x+Fc1x+Fc0y+Fc1y+Fc0z+Fc1z;
                        Aphi = Fph0x+Fph1x+Fph0y+Fph1y+Fph0z+Fph1z;


                        //membrane current contributions
                        Ac += flux->dfdci[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;
                        Aphi += flux->dfdphim[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;
                        // Different Compartment Terms
                        // C Extracellular with C Inside
                        ierr = MatSetValue(Jac,Ind_2(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ind_2(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -flux->dfdci[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // C Intra with C Extra
                        ierr = MatSetValue(Jac,Ind_2(x,y,z,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,ion,Nc-1,Nx,Ny,Nz),
                                           flux->dfdce[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // C Extracellular with Phi Inside
                        ierr = MatSetValue(Jac,Ind_2(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ind_2(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -flux->dfdphim[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // C Intra with Phi Extra
                        ierr = MatSetValue(Jac,Ind_2(x,y,z,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,Ni,Nc-1,Nx,Ny,Nz),
                                           -flux->dfdphim[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Volume terms
                        //C extra with intra alpha
                        ierr = MatSetValue(Jac,Ind_2(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ind_2(x,y,z,Ni+1,comp,Nx,Ny,Nz),
                                           -c[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)],INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //C intra with intra alpha
                        ierr = MatSetValue(Jac,Ind_2(x,y,z,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,Ni+1,comp,Nx,Ny,Nz),
                                           c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)],INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Same compartment terms
                        // c with c
                        ierr = MatSetValue(Jac,Ind_2(x,y,z,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,ion,comp,Nx,Ny,Nz),Ac,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // c with phi
                        ierr = MatSetValue(Jac,Ind_2(x,y,z,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,Ni,comp,Nx,Ny,Nz),Aphi,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                    }
                    //Extracellular terms
                    comp = Nc-1;
                    //Electrodiffusion contributions
                    Ftmpx = 0;Fc0x = 0;Fc1x = 0;Fph0x = 0;Fph1x = 0;
                    Ftmpy = 0;Fc0y = 0;Fc1y = 0;Fph0y = 0;Fph1y = 0;
                    Ftmpz = 0;Fc0z = 0;Fc1z = 0;Fph0z = 0;Fph1z = 0;
                    if(x < Nx-1){
                        Ftmpx = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                        cp[c_index(x+1,y,z,comp,ion,Nx,Ny,Nz)])/2/dx*dt/dx;
                        Fc0x = Ftmpx/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                        Fph0x = z_charge[ion]*Ftmpx;
                        // Right c with left c (-Fc0x)
                        ierr = MatSetValue(Jac,Ind_2(x+1,y,z,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc0x,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Right c with left phi (-Fph0x)
                        ierr = MatSetValue(Jac,Ind_2(x+1,y,z,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fph0x,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(x > 0){
                        Ftmpx = Dcs[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)*3]*(cp[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                          cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2/dx*dt/dx;
                        Fc1x = Ftmpx/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                        Fph1x = z_charge[ion]*Ftmpx;
                        //left c with right c (-Fc1x)
                        ierr = MatSetValue(Jac,Ind_2(x-1,y,z,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc1x,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Left c with right phi (-Fph1x)
                        ierr = MatSetValue(Jac,Ind_2(x-1,y,z,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fph1x,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(y < Ny-1){
                        Ftmpy = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+1]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                          cp[c_index(x,y+1,z,comp,ion,Nx,Ny,Nz)])/2/dy*dt/dy;
                        Fc0y = Ftmpy/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                        Fph0y = z_charge[ion]*Ftmpy;
                        // Upper c with lower c (-Fc0y)
                        ierr = MatSetValue(Jac,Ind_2(x,y+1,z,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc0y,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Upper c with lower phi (-Fph0y)
                        ierr = MatSetValue(Jac,Ind_2(x,y+1,z,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fph0y,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(y > 0){
                        Ftmpy = Dcs[c_index(x,y-1,z,comp,ion,Nx,Ny,Nz)*3+1]*(cp[c_index(x,y-1,z,comp,ion,Nx,Ny,Nz)]+
                                                                            cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2/dy*dt/dy;
                        Fc1y = Ftmpy/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                        Fph1y = z_charge[ion]*Ftmpy;
                        //Lower c with Upper c (-Fc1y)
                        ierr = MatSetValue(Jac,Ind_2(x,y-1,z,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc1y,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Lower c with Upper phi (-Fph1y)
                        ierr = MatSetValue(Jac,Ind_2(x,y-1,z,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fph1y,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(z < Nz-1){
                        Ftmpz = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                          cp[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])/2/dz*dt/dz;
                        Fc0z = Ftmpz/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                        Fph0z = z_charge[ion]*Ftmpz;
                        // Upper c with lower c (-Fc0y)
                        ierr = MatSetValue(Jac,Ind_2(x,y,z+1,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc0z,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Upper c with lower phi (-Fph0y)
                        ierr = MatSetValue(Jac,Ind_2(x,y,z+1,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fph0z,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(z > 0){
                        Ftmpz = Dcs[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)]+
                                                                            cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2/dz*dt/dz;
                        Fc1z = Ftmpz/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                        Fph1z = z_charge[ion]*Ftmpz;
                        //Lower c with Upper c (-Fc1y)
                        ierr = MatSetValue(Jac,Ind_2(x,y,z-1,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc1z,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Lower c with Upper phi (-Fph1y)
                        ierr = MatSetValue(Jac,Ind_2(x,y,z-1,ion,comp,Nx,Ny,Nz),Ind_2(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fph1z,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }

                    //Diagonal term contribution
                    Ac = (1-al[al_index(x,y,z,0,Nx,Ny,Nz)]-al[al_index(x,y,z,1,Nx,Ny,Nz)])+Fc0x+Fc1x+Fc0y+Fc1y+Fc0z+Fc1z;
                    Aphi = Fph0x+Fph1x+Fph0y+Fph1y+Fph0z+Fph1z;

                    //Membrane current contribution
                    for(comp = 0; comp < Nc-1; comp++){
                        Ac -= flux->dfdce[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;
                        Aphi += flux->dfdphim[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;
                    }
                    //Add bath contributions
                    Ftmpx = sqrt(pow(Dcb[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)*3],2)+
                                 pow(Dcb[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)*3+1],2)+
                                 pow(Dcb[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)*3+2],2));
                    Ac -= Ftmpx*(cp[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)]+cbath[ion])/(2*c[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)])*dt;
                    Aphi -= Ftmpx*(cp[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)]+cbath[ion])*z_charge[ion]/2*dt;

                    //Insert extracell to extracell parts
                    // c with c
                    ierr = MatSetValue(Jac,Ind_2(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ind_2(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ac,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    // c with phi
                    ierr = MatSetValue(Jac,Ind_2(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ind_2(x,y,z,Ni,Nc-1,Nx,Ny,Nz),Aphi,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }

            }
        }
    }

        //Electroneutrality charge-capcitance condition
    for(z=0;z<Nz;z++){
        for(y = 0; y < Ny; y++){
            for(x = 0; x < Nx; x++){
                //electroneutral-concentration entries
                for(ion = 0; ion < Ni; ion++){
                    for(comp = 0; comp < Nc-1; comp++){
                        //Phi with C entries
                        ierr = MatSetValue(Jac,Ind_2(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_2(x,y,z,ion,comp,Nx,Ny,Nz),
                                           z_charge[ion]*al[al_index(x,y,z,comp,Nx,Ny,Nz)],INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    //Phi with C extracellular one
                    comp = Nc-1;
                    ierr = MatSetValue(Jac,Ind_2(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_2(x,y,z,ion,comp,Nx,Ny,Nz),
                                       z_charge[ion]*
                                       (1-al[al_index(x,y,z,0,Nx,Ny,Nz)]-al[al_index(x,y,z,1,Nx,Ny,Nz)]),INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;

                }
                //electroneutrality-voltage entries
                Aphi = 0;
                for(comp = 0; comp < Nc-1; comp++){
                    Aphi -= cm[comp];
                }
                //extraphi with extra phi
                ierr = MatSetValue(Jac,Ind_2(x,y,z,Ni,Nc-1,Nx,Ny,Nz),Ind_2(x,y,z,Ni,Nc-1,Nx,Ny,Nz),Aphi,INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
                for(comp = 0; comp < Nc-1; comp++){
                    //The next 3 are inserted in init jacobian for the grid (only if cm is constant)
                    //Extra phi with intra phi
                ierr = MatSetValue(Jac,Ind_2(x,y,z,Ni,Nc-1,Nx,Ny,Nz),Ind_2(x,y,z,Ni,comp,Nx,Ny,Nz),cm[al_index(
                        x+xi,y+yi,z,comp,Nx,Ny,Nz)],INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                    // Intra phi with Extraphi
                ierr = MatSetValue(Jac,Ind_2(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_2(x,y,z,Ni,Nc-1,Nx,Ny,Nz),cm[al_index(
                        x+xi,y+yi,z,comp,Nx,Ny,Nz)],INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                    //Intra phi with Intra phi
                ierr = MatSetValue(Jac,Ind_2(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_2(x,y,z,Ni,comp,Nx,Ny,Nz),-cm[al_index(
                        x+xi,y+yi,z,comp,Nx,Ny,Nz)],INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                    //Extra phi with intra-Volume
                    ierr = MatSetValue(Jac,Ind_2(x,y,z,Ni,Nc-1,Nx,Ny,Nz),Ind_2(x,y,z,Ni+1,comp,Nx,Ny,Nz),
                                       -cz(c,z_charge,x,y,z,Nx,Ny,Nz,Nc-1,user),INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Intra phi with Intra Vol
                    ierr = MatSetValue(Jac,Ind_2(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_2(x,y,z,Ni+1,comp,Nx,Ny,Nz),
                                       cz(c,z_charge,x,y,z,Nx,Ny,Nz,comp,user),INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
            }
        }
    }
    //water flow
    for(z=0;z<Nz;z++){
        for(y = 0; y < Ny; y++){
            for(x = 0; x < Nx; x++){
                for(comp = 0; comp < Nc-1; comp++){
                    //Water flow volume fraction entries
                    //Volume to Volume
                    Ac = 1+(flux->dwdpi[al_index(x,y,z,comp,Nx,Ny,Nz)]*(con_vars->ao[phi_index(x,y,z,Nc-1,Nx,Ny,Nz)]/
                                                                       (pow(1-al[al_index(x,y,z,0,Nx,Ny,Nz)]-al[al_index(x,y,z,1,Nx,Ny,Nz)],2))+
                                                                     con_vars->ao[phi_index(x,y,z,comp,Nx,Ny,Nz)]/pow(al[al_index(x,y,z,comp,Nx,Ny,Nz)],2))+
                            flux->dwdal[al_index(x,y,z,comp,Nx,Ny,Nz)])*dt;
                    ierr = MatSetValue(Jac,Ind_2(x,y,z,Ni+1,comp,Nx,Ny,Nz),Ind_2(x,y,z,
                                                                                Ni+1,comp,Nx,Ny,Nz),Ac,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Off diagonal (from aNc=1-sum(ak))
                    for(PetscInt l = 0; l < comp; l++){
                        ierr = MatSetValue(Jac,Ind_2(x,y,z,Ni+1,comp,Nx,Ny,Nz),Ind_2(x,y,z,Ni+1,l,Nx,Ny,Nz),
                                           flux->dwdpi[al_index(x,y,z,comp,Nx,Ny,Nz)]*con_vars->ao[phi_index(x,y,z,Nc-1,Nx,Ny,Nz)]/
                                           pow(1-al[al_index(x,y,z,0,Nx,Ny,Nz)]-al[al_index(x,y,z,1,Nx,Ny,Nz)],2)*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    for(PetscInt l = comp+1; l < Nc-1; l++){
                        ierr = MatSetValue(Jac,Ind_2(x,y,z,Ni+1,comp,Nx,Ny,Nz),Ind_2(x,y,z,Ni+1,l,Nx,Ny,Nz),
                                           flux->dwdpi[al_index(x,y,z,comp,Nx,Ny,Nz)]*
                                           con_vars->ao[phi_index(x,y,z,Nc-1,Nx,Ny,Nz)]/
                                           ((1-al[al_index(x,y,z,0,Nx,Ny,Nz)]-al[al_index(x,y,z,1,Nx,Ny,Nz)])*
                                            (1-al[al_index(x,y,z,0,Nx,Ny,Nz)]-al[al_index(x,y,z,1,Nx,Ny,Nz)]))*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    for(ion = 0; ion < Ni; ion++){
                        //Volume to extra c
                        ierr = MatSetValue(Jac,Ind_2(x,y,z,Ni+1,comp,Nx,Ny,Nz),Ind_2(x,y,z,ion,Nc-1,Nx,Ny,Nz),
                                           flux->dwdpi[al_index(x,y,z,comp,Nx,Ny,Nz)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Volume to intra c
                        ierr = MatSetValue(Jac,Ind_2(x,y,z,Ni+1,comp,Nx,Ny,Nz),Ind_2(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -flux->dwdpi[al_index(x,y,z,comp,Nx,Ny,Nz)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                }
            }
        }
    }

    ierr = MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    if(Profiling_on) {
        PetscLogEventEnd(event[9], 0, 0, 0, 0);
    }
    return ierr;
}

int Newton_Solve_Grid(PetscInt xi, PetscInt yi,struct AppCtx *user) {


    if(Profiling_on) {
        PetscLogEventBegin(event[11], 0, 0, 0, 0);
    }

    PetscReal rsd;
    PetscErrorCode ierr = 0;
    PetscReal const *temp;

    PetscInt x,y,z,comp,ion;
    PetscInt Nx = 2*width_size+1;
    PetscInt Ny = 2*width_size+1;
    PetscInt Nz = user->Nz;


    PetscReal tol = reltol*array_max(user->grid_vars_past->c,(size_t)Nz*Nx*Ny*Ni*Nc);
    rsd = tol + 1;

    for (PetscInt iter = 0; iter < itermax; iter++) {

        ierr = Grid_Residual_algebraic(user->grid_slvr->Res, xi, yi, user);CHKERRQ(ierr);

        ierr = VecNorm(user->grid_slvr->Res, NORM_MAX, &rsd);CHKERRQ(ierr);
//        printf("(%d,%d),%.10e\n",xi,yi,rsd);

        if (rsd < tol) {
            if(Profiling_on) {
                PetscLogEventEnd(event[11], 0, 0, 0, 0);
            }
            return iter;
        }
        ierr = Grid_Jacobian_algebraic(user->grid_slvr->A, xi, yi, user);CHKERRQ(ierr);

        //Set the new operator
        ierr = KSPSetOperators(user->grid_slvr->ksp, user->grid_slvr->A, user->grid_slvr->A);CHKERRQ(ierr);

        //Solve
        ierr = KSPSolve(user->grid_slvr->ksp, user->grid_slvr->Res, user->grid_slvr->Q);CHKERRQ(ierr);


        ierr = VecGetArrayRead(user->grid_slvr->Q, &temp);CHKERRQ(ierr);
        for (z = 0; z < Nz; z++){
            for(y = 0; y < Ny; y++){
                for(x = 0; x < Nx; x++){

                    for(comp = 0; comp < Nc; comp++){
                        for(ion = 0; ion < Ni; ion++){
                            user->grid_vars->c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)] -= temp[Ind_2(x,y,z,ion,comp,Nx,Ny,Nz)];
                        }
                        user->grid_vars->phi[phi_index(x,y,z,comp,Nx,Ny,Nz)] -= temp[Ind_2(x,y,z,Ni,comp,Nx,Ny,Nz)];
                    }
                    for(comp = 0; comp < Nc-1; comp++){
                        user->grid_vars->alpha[al_index(x,y,z,comp,Nx,Ny,Nz)] -= temp[Ind_2(x,y,z,Ni+1,comp,Nx,Ny,Nz)];
                    }
                }
            }
        }
        ierr = VecRestoreArrayRead(user->grid_slvr->Q, &temp);CHKERRQ(ierr);

    }

    if (rsd > tol) {
        fprintf(stderr, "Netwon Iteration did not converge! Stopping...\n");
        exit(EXIT_FAILURE); /* indicate failure.*/
    }
    if(Profiling_on) {
        PetscLogEventEnd(event[11], 0, 0, 0, 0);
    }
    return itermax;
}
PetscErrorCode Update_Grid(PetscInt xi, PetscInt yi,PetscReal t,struct AppCtx *user)
{
    PetscErrorCode ierr = 0;


    PetscReal dt = user->dt;

    user->dt = user->dt_space[xy_index(xi,yi,0,user->Nx,0,0)];


    PetscInt Nx = 2*width_size+1;
    PetscInt Ny = 2*width_size+1;
    PetscInt Nz = user->Nz;
    PetscInt ion,comp,x,y,iter;

    //Load current variable into past variable
    memcpy(user->grid_vars_past->c,user->grid_vars->c,sizeof(PetscReal)*Nx*Ny*Nz*Nc*Ni);
    memcpy(user->grid_vars_past->phi,user->grid_vars->phi,sizeof(PetscReal)*Nx*Ny*Nz*Nc);
    memcpy(user->grid_vars_past->alpha,user->grid_vars->alpha,sizeof(PetscReal)*Nx*Ny*Nz*(Nc-1));
    //Calculate diffusion
    //compute diffusion coefficients
    grid_diff_coef(user->Dcs, user->grid_vars_past->alpha, 1, user,xi,yi);
    //Bath diffusion
    grid_diff_coef(user->Dcb, user->grid_vars_past->alpha, Batheps, user,xi,yi);

    excitation_grid(user, t - dt, xi, yi);

    PetscInt steps = 0;
    PetscInt NSteps = (PetscInt)floor(dt/user->dt); //1;
    PetscInt accepted_step = 0;

    while(steps<NSteps) {

        //Perform Newton Solve
        iter = Newton_Solve_Grid(xi, yi, user);

        //Check if we accept the step
        if (iter < 3 || accepted_step || NSteps >= Max_Grid_Refine) {
            steps++;

            //Update Gating variable
            gatevars_update_grid(user->grid_gate_vars, user->grid_vars, user->dt * 1e3, user);

            //Update Excitation
            excitation_grid(user,t-dt+user->dt*steps,xi,yi);

            accepted_step = 1;

            //Load current variable into past variable
            memcpy(user->grid_vars_past->c,user->grid_vars->c,sizeof(PetscReal)*Nx*Ny*Nz*Nc*Ni);
            memcpy(user->grid_vars_past->phi,user->grid_vars->phi,sizeof(PetscReal)*Nx*Ny*Nz*Nc);
            memcpy(user->grid_vars_past->alpha,user->grid_vars->alpha,sizeof(PetscReal)*Nx*Ny*Nz*(Nc-1));
            //Calculate diffusion
            //compute diffusion coefficients
            grid_diff_coef(user->Dcs, user->grid_vars_past->alpha, 1, user,xi,yi);
            //Bath diffusion
            grid_diff_coef(user->Dcb, user->grid_vars_past->alpha, Batheps, user,xi,yi);

            if(xi==16&&yi==16) {
                write_point(user->fp,user,t-dt+user->dt*steps,16,16);
            }

        } else {
            //If we aren't below cutoff. Half the time step.
            user->dt = user->dt / 2;
            NSteps = 2 * NSteps;
//            printf("Reducing step at (%d,%d) to %f\n",xi,yi,user->dt);
            //Reset current vars
            //Load current variable into past variable
            memcpy(user->grid_vars->c,user->grid_vars_past->c,sizeof(PetscReal)*Nx*Ny*Nz*Nc*Ni);
            memcpy(user->grid_vars->phi,user->grid_vars_past->phi,sizeof(PetscReal)*Nx*Ny*Nz*Nc);
            memcpy(user->grid_vars->alpha,user->grid_vars_past->alpha,sizeof(PetscReal)*Nx*Ny*Nz*(Nc-1));
        }
    }

    if(user->dt<dt){
        user->dt_space[xy_index(xi,yi,0,user->Nx,0,0)] = 2*user->dt;
    } else{
        user->dt_space[xy_index(xi,yi,0,user->Nx,0,0)] = user->dt;
    }


    user->dt = dt;
    return ierr;

}

PetscErrorCode Update_Solution(Vec current_state,PetscReal t,struct AppCtx *user)
{

    PetscErrorCode ierr = 0;
    PetscInt x,y,z,ion,comp;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt Nz = user->Nz;
    PetscInt nx = 2*width_size+1;
    PetscInt ny = 2*width_size+1;
    PetscReal vm_new;
    PetscReal threshhold = 0.1;//0.1; //mV threshhold for update guess.


    for(y=0;y<Ny;y++){
        for(x = 0; x < Nx; x++){
            //Look over all z's.
            z=0;
            while(z<Nz){
                vm_new = (user->state_vars->phi[phi_index(x,y,z,0,Nx,Ny,Nz)]-
                          user->state_vars->phi[phi_index(x,y,z,Nc-1,Nx,Ny,Nz)])*RTFC;

                //If it's above the threshhold. Or it previously was adaptively refined.
                if(fabs(vm_new-user->vm_past[xy_index(x,y,z,Nx,Ny,Nz)]) > threshhold || user->dt_space[xy_index(x,y,0,Nx,0,0)] < user->dt){
//            printf("Updating: (%d,%d)\n",x,y);
                    // Load new gridpoint
                    Load_Grid(user,x,y);
                    //Update new grid
                    Update_Grid(x,y,t,user);
                    //Save the held variable
                    Unload_Grid(user,x,y);
                    z=Nz; // Set to break loop
                }
                z++;
            }
        }
    }


    return ierr;

}

PetscErrorCode initialize_grid_jacobian(Mat Jac,struct AppCtx *user,int grid) {
    printf("Initializing Jacobian Memory\n");
    PetscErrorCode ierr;
    PetscInt Nx;
    PetscInt Ny;
    PetscInt Nz;
    if(grid) {
        Nx = 2 * width_size + 1;
        Ny = 2 * width_size + 1;
    }else{
        Nx = user->Nx;
        Ny = user->Ny;
    }
    PetscInt ind = 0;
    PetscInt x,y,ion,comp;

    //Ionic concentration equations
    for(x=0;x<Nx;x++) {
        for(y=0;y<Ny;y++) {
            for(ion=0;ion<Ni;ion++) {
                for(comp=0;comp<Nc-1;comp++) {
                    //Electrodiffusion contributions

                    if(x<Nx-1) {
                        // Right c with left c (-Fc0x)
                        ierr = MatSetValue(Jac,Ind_2(
                                x+1,y,0,ion,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),0,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Right c with left phi (-Fph0x)
                        ierr = MatSetValue(Jac,Ind_2(
                                x+1,y,0,ion,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),0,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        if (use_en_deriv&&!grid) {
                            //Right phi with left c in voltage eqn
                            ierr = MatSetValue(Jac,Ind_2(
                                    x+1,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),0,INSERT_VALUES);CHKERRQ(ierr);
                            ind++;
                        }


                    }
                    if(x>0) {
                        //left c with right c (-Fc1x)
                        ierr = MatSetValue(Jac,Ind_2(
                                x-1,y,0,ion,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),0,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Left c with right phi (-Fph1x)
                        ierr = MatSetValue(Jac,Ind_2(
                                x-1,y,0,ion,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),0,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        if (use_en_deriv&&!grid) {
                            //Left phi with right c in voltage eqn
                            ierr = MatSetValue(Jac,Ind_2(x-1,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),0,
                                               INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                    }
                    if(y<Ny-1) {
                        // Upper c with lower c (-Fc0y)
                        ierr = MatSetValue(Jac,Ind_2(x,y+1,0,ion,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),0,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Upper c with lower phi (-Fph0y)
                        ierr = MatSetValue(Jac,Ind_2(x,y+1,0,ion,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),0,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        if (use_en_deriv&&!grid) {
                            //Upper phi with lower c in voltage eqn
                            ierr = MatSetValue(Jac,Ind_2(x,y+1,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),0,
                                               INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                    }
                    if(y>0) {
                        //Lower c with Upper c (-Fc1y)
                        ierr = MatSetValue(Jac,Ind_2(x,y-1,0,ion,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),0,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Lower c with Upper phi (-Fph1y)
                        ierr = MatSetValue(Jac,Ind_2(x,y-1,0,ion,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),0,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        if (use_en_deriv&&!grid) {
                            //Lower phi with upper c in voltage eqn
                            ierr = MatSetValue(Jac,Ind_2(x,y-1,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),0,
                                               INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                    }

                    // Different Compartment Terms
                    // C Extracellular with C Inside
                    ierr = MatSetValue(Jac,Ind_2(x,y,0,ion,Nc-1,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // C Intra with C Extra
                    ierr = MatSetValue(Jac,Ind_2(x,y,0,ion,comp,Nx,0,0),Ind_2(x,y,0,ion,Nc-1,Nx,0,0),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // C Extracellular with Phi Inside
                    ierr = MatSetValue(Jac,Ind_2(x,y,0,ion,Nc-1,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // C Intra with Phi Extra
                    ierr = MatSetValue(Jac,Ind_2(x,y,0,ion,comp,Nx,0,0),Ind_2(x,y,0,Ni,Nc-1,Nx,0,0),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    if(!separate_vol||grid) {
                        //Volume terms
                        //C extra with intra alpha
                        ierr = MatSetValue(Jac,Ind_2(x,y,0,ion,Nc-1,Nx,0,0),Ind_2(x,y,0,Ni+1,comp,Nx,0,0),0,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //C intra with intra alpha
                        ierr = MatSetValue(Jac,Ind_2(x,y,0,ion,comp,Nx,0,0),Ind_2(x,y,0,Ni+1,comp,Nx,0,0),0,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    //Same compartment terms
                    // c with c
                    ierr = MatSetValue(Jac,Ind_2(x,y,0,ion,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // c with phi
                    ierr = MatSetValue(Jac,Ind_2(x,y,0,ion,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    if (use_en_deriv&&!grid) {
                        //Intra-Phi with c (voltage eqn)
                        ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),0,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //IntraPhi with c extra(volt eqn)
                        ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,ion,Nc-1,Nx,0,0),0,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                        //Extra-Phi with intra-c (voltage eqn)
                        ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,Nc-1,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),0,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }

                }
                //Extracellular terms
                comp = Nc-1;
                //Electrodiffusion contributions
                if(x<Nx-1) {
                    // Right c with left c (-Fc0x)
                    ierr = MatSetValue(Jac,Ind_2(
                            x+1,y,0,ion,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Right c with left phi (-Fph0x)
                    ierr = MatSetValue(Jac,Ind_2(
                            x+1,y,0,ion,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    if (use_en_deriv&&!grid) {
                        // left Phi with right c (voltage eqn)
                        ierr = MatSetValue(Jac,Ind_2(
                                x+1,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),0,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                }
                if(x>0) {
                    //left c with right c (-Fc1x)
                    ierr = MatSetValue(Jac,Ind_2(
                            x-1,y,0,ion,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Left c with right phi (-Fph1x)
                    ierr = MatSetValue(Jac,Ind_2(
                            x-1,y,0,ion,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    if (use_en_deriv&&!grid) {
                        // left Phi with right c (voltage eqn)
                        ierr = MatSetValue(Jac,Ind_2(
                                x-1,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),0,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                }
                if(y<Ny-1) {
                    // Upper c with lower c (-Fc0y)
                    ierr = MatSetValue(Jac,Ind_2(x,y+1,0,ion,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Upper c with lower phi (-Fph0y)
                    ierr = MatSetValue(Jac,Ind_2(x,y+1,0,ion,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    if (use_en_deriv&&!grid) {
                        // Upper Phi with lower c (voltage eqn)
                        ierr = MatSetValue(Jac,Ind_2(x,y+1,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),0,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                }
                if(y>0) {
                    //Lower c with Upper c (-Fc1y)
                    ierr = MatSetValue(Jac,Ind_2(x,y-1,0,ion,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Lower c with Upper phi (-Fph1y)
                    ierr = MatSetValue(Jac,Ind_2(x,y-1,0,ion,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),0,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    if (use_en_deriv&&!grid) {
                        // Lower Phi with upper c (voltage eqn)
                        ierr = MatSetValue(Jac,Ind_2(x,y-1,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),0,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                }
                //Insert extracell to extracell parts
                // c with c
                ierr = MatSetValue(Jac,Ind_2(x,y,0,ion,Nc-1,Nx,0,0),Ind_2(x,y,0,ion,Nc-1,Nx,0,0),0,INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                // c with phi
                ierr = MatSetValue(Jac,Ind_2(x,y,0,ion,Nc-1,Nx,0,0),Ind_2(x,y,0,Ni,Nc-1,Nx,0,0),0,INSERT_VALUES);CHKERRQ(ierr);
                ind++;
                if (use_en_deriv&&!grid) {
                    //phi with c (voltage eqn)
                    ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,Nc-1,Nx,0,0),Ind_2(x,y,0,ion,Nc-1,Nx,0,0),0,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
            }
            if (use_en_deriv&&!grid) {
                //Derivative of charge-capacitance
                for (comp = 0; comp < Nc - 1; comp++) {
                    if (x < Nx - 1) {
                        //Right phi with left phi (-Fph0x)
                        ierr = MatSetValue(Jac,Ind_2(x+1,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),0,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if (x > 0) {
                        //Left phi with right phi (-Fph1x)
                        ierr = MatSetValue(Jac,Ind_2(x-1,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),0,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if (y < Ny - 1) {
                        //Upper phi with lower phi (-Fph0y)
                        ierr = MatSetValue(Jac,Ind_2(x,y+1,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),0,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if (y > 0) {
                        //Lower phi with upper phi (-Fph1y)
                        ierr = MatSetValue(Jac,Ind_2(x,y-1,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),0,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    //Intra-phi with Intra-phi
                    ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),0,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Intra-phi with extra-phi
                    ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,Nc-1,Nx,0,0),0,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                //Extracellular terms
                comp = Nc - 1;
                if (x < Nx - 1) {
                    //Right phi with left phi (-Fph0x)
                    ierr = MatSetValue(Jac,Ind_2(x+1,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),0,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if (x > 0) {
                    //Left phi with right phi (-Fph1x)
                    ierr = MatSetValue(Jac,Ind_2(x-1,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),0,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if (y < Ny - 1) {
                    //Upper phi with lower phi (-Fph0y)
                    ierr = MatSetValue(Jac,Ind_2(x,y+1,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),0,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if (y > 0) {
                    //Lower phi with upper phi (-Fph1y)
                    ierr = MatSetValue(Jac,Ind_2(x,y-1,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),0,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }

                for (int k = 0; k < Nc - 1; k++) {
                    //Extra-phi with Intra-phi
                    ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,k,Nx,0,0),0,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                //extra-phi with extra-phi
                ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),0,INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
            }

        }
    }
    if(!use_en_deriv||grid) {
        //Electroneutrality charge-capcitance condition
        for (x = 0; x < Nx; x++) {
            for (y = 0; y < Ny; y++) {
                //electroneutral-concentration entries
                for (ion = 0; ion < Ni; ion++) {
                    for (comp = 0; comp < Nc - 1; comp++) {
                        //Phi with C entries
                        ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),0,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    //Phi with C extracellular one
                    comp = Nc - 1;
                    ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),0,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;

                }
                //electroneutrality-voltage entries

                //extraphi with extra phi
                ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,Nc-1,Nx,0,0),Ind_2(x,y,0,Ni,Nc-1,Nx,0,0),0,INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
                for (comp = 0; comp < Nc - 1; comp++) {
                    //Extra phi with intra phi
                    ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,Nc-1,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),0,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    // Intra phi with Extraphi
                    ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,Nc-1,Nx,0,0),0,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Intra phi with Intra phi
                    ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),0,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    if(!separate_vol||grid) {
                        //Extra phi with intra-Volume
                        ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,Nc-1,Nx,0,0),Ind_2(x,y,0,Ni+1,comp,Nx,0,0),0,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Intra phi with Intra Vol
                        ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni+1,comp,Nx,0,0),0,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(grid) {
                        //Extra phi with intra phi
                        ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,Nc-1,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),0,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // Intra phi with Extraphi
                        ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,Nc-1,Nx,0,0),0,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Intra phi with Intra phi
                        ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni,comp,Nx,0,0),Ind_2(x,y,0,Ni,comp,Nx,0,0),0,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }

                }
            }
        }
    }
    if(!separate_vol||grid) {
        //water flow
        for (x = 0; x < Nx; x++) {
            for (y = 0; y < Ny; y++) {
                for (comp = 0; comp < Nc - 1; comp++) {
                    //Water flow volume fraction entries
                    //Volume to Volume
                    ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni+1,comp,Nx,0,0),Ind_2(x,y,0,Ni+1,comp,Nx,0,0),0,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Off diagonal (from aNc=1-sum(ak))
                    for (PetscInt l = 0; l < comp; l++) {
                        ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni+1,comp,Nx,0,0),Ind_2(x,y,0,Ni+1,l,Nx,0,0),0,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    for (PetscInt l = comp + 1; l < Nc - 1; l++) {
                        ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni+1,comp,Nx,0,0),Ind_2(x,y,0,Ni+1,l,Nx,0,0),0,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    for (ion = 0; ion < Ni; ion++) {
                        //Volume to extra c
                        ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni+1,comp,Nx,0,0),Ind_2(x,y,0,ion,Nc-1,Nx,0,0),0,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Volume to intra c
                        ierr = MatSetValue(Jac,Ind_2(x,y,0,Ni+1,comp,Nx,0,0),Ind_2(x,y,0,ion,comp,Nx,0,0),0,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                }
            }
        }
    }
    ierr = MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);


    return ierr;
}

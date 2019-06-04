#include "constants.h"
#include "functions.h"


PetscErrorCode newton_solve(Vec current_state,struct Solver *slvr,struct AppCtx *user)
{

    PetscReal rsd;
    PetscErrorCode ierr = 0;
    PetscReal const *temp;
    PetscInt num_iter,comp,ion,x,y,z;
    PetscReal rnorm;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt Nz = user->Nz;


    PetscLogDouble tic,toc;

    //Diffusion in each compartment
    //Has x and y components
    //x will be saved at first positions (0,3,6,...)
    //y at 2nd (1,4,7,...)
    //z at 3rd (2,5,8,...)
    //still use c_index(x,y,comp,ion,Nx), but with ind*3 or ind*3+1 or ind*3+2

    extract_subarray(current_state,user->state_vars);

    PetscReal tol = reltol*array_max(user->state_vars->c,(size_t)Nx*Ny*Nz*Ni*Nc);
    restore_subarray(current_state,user->state_vars);
    rsd = tol+1;

    for(PetscInt iter=0;iter<10;iter++)
    {
        if(separate_vol){
            if(use_en_deriv){
                ierr = calc_residual_no_vol(user->slvr->snes,current_state,slvr->Res,user);CHKERRQ(ierr);

            } else{
                ierr = calc_residual_algebraic_no_vol(user->slvr->snes,current_state,slvr->Res,user);CHKERRQ(ierr);
            }

        }else{
            if(use_en_deriv){
                ierr = calc_residual(user->slvr->snes,current_state,slvr->Res,user);CHKERRQ(ierr);

            } else{
                ierr = calc_residual_algebraic(user->slvr->snes,current_state,slvr->Res,user);CHKERRQ(ierr);
            }

        }

        ierr = VecNorm(slvr->Res,NORM_MAX,&rsd);CHKERRQ(ierr);
        if(rsd<tol)
        {
            if(details)
            {
                printf("Iteration: %d, Residual: %.10e\n",iter,rsd);
            }
            return ierr;
        }
        if(separate_vol){
            if(use_en_deriv){
                ierr = calc_jacobian_no_vol(user->slvr->snes,current_state,slvr->A,slvr->A, user);CHKERRQ(ierr);

            }else{
                ierr = calc_jacobian_algebraic_no_vol(user->slvr->snes,current_state,slvr->A,slvr->A, user);CHKERRQ(ierr);

            }

        }else{
            if(use_en_deriv){
                ierr = calc_jacobian(user->slvr->snes,current_state,slvr->A,slvr->A, user);CHKERRQ(ierr);

            }else{
                ierr = calc_jacobian_algebraic(user->slvr->snes,current_state,slvr->A,slvr->A, user);CHKERRQ(ierr);

            }
        }
        //Set the new operator
        ierr = KSPSetOperators(slvr->ksp,slvr->A,slvr->A);CHKERRQ(ierr);

        //Solve
        PetscTime(&tic);
        ierr = KSPSolve(slvr->ksp,slvr->Res,slvr->Q);CHKERRQ(ierr);
        PetscTime(&toc);

        ierr = KSPGetIterationNumber(user->slvr->ksp,&num_iter); CHKERRQ(ierr);
        ierr =  KSPGetResidualNorm(user->slvr->ksp,&rnorm); CHKERRQ(ierr);


        if(details) {
            printf("KSP Solve time: %f, iter num:%d, norm: %.10e\n",toc-tic,num_iter,rnorm);
        }

//        PetscTime(&tic);
        ierr = VecGetArrayRead(slvr->Q,&temp); CHKERRQ(ierr);
        extract_subarray(current_state,user->state_vars);
        for(z=0;z<Nz;z++){
            for(y = 0; y < Ny; y++){
                for(x = 0; x < Nx; x++){
                    for(comp = 0; comp < Nc; comp++){
                        for(ion = 0; ion < Ni; ion++){
                            user->state_vars->c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)] -= temp[Ind_1(x,y,z,ion,comp,Nx,Ny,Nz)];
                        }
                        user->state_vars->phi[phi_index(x,y,z,comp,Nx,Ny,Nz)] -= temp[Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz)];
                    }
                    if(!separate_vol){
                        for(comp = 0; comp < Nc-1; comp++){
                            user->state_vars->alpha[al_index(x,y,z,comp,Nx,Ny,Nz)] -= temp[Ind_1(x,y,z,Ni+1,comp,Nx,Ny,Nz)];
                        }
                    }
                }
            }
        }

        ierr = VecRestoreArrayRead(slvr->Q,&temp);
        restore_subarray(current_state,user->state_vars);

        if(details)
        {
            printf("Iteration: %d, Residual: %.10e\n",iter,rsd);
        }
    }

    if(rsd>tol)
    {
        fprintf(stderr, "Netwon Iteration did not converge! Findal residual: %.10e. Stopping...\n",rsd);
        exit(EXIT_FAILURE);
    }

    return ierr;
}


PetscErrorCode calc_residual(SNES snes,Vec current_state,Vec Res,void *ctx)
{
    //Residual equation using derivative of the charge-capacitance relation
    //Volume is solved for here
    struct AppCtx * user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if(Profiling_on) {
        PetscLogEventBegin(event[1], 0, 0, 0, 0);
    }
    ierr = extract_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    //Compute membrane ionic flux relation quantitites
    ionmflux(user);

    //Compute membrane water flow related quantities
    wflowm(user);


    PetscReal *c = user->state_vars->c;
    PetscReal *phi = user->state_vars->phi;
    PetscReal *al = user->state_vars->alpha;
    PetscReal *cp = user->state_vars_past->c;
    PetscReal *alp = user->state_vars_past->alpha;
    PetscReal *phip = user->state_vars_past->phi;
    PetscReal *cm = user->con_vars->cm;

    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct FluxData *flux = user->flux;
    PetscReal dt = user->dt;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    PetscReal dz = user->dz;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt Nz = user->Nz;

    //Residual for concentration equations
    PetscReal Rcvx,Rcvy,Rcvz,Resc;
    PetscReal RcvxRight,RcvyUp,RcvzTop;

    //Residual for fluxes in voltage differential equations
    PetscReal Rphx[Nc], Rphy[Nc],Rphz[Nc], RphxRight[Nc], RphyUp[Nc],RphzTop[Nc];
    PetscReal Resph,ResphN;

    PetscReal alNc,alpNc;
    PetscInt ion,comp,x,y,z;

    for(z=0;z<Nz;z++){
        for(y = 0; y < Ny; y++){
            for(x = 0; x < Nx; x++){

                //Init voltage tracking to zero
                for(comp = 0; comp < Nc; comp++){
                    Rphx[comp] = 0;
                    Rphy[comp] = 0;
                    RphxRight[comp] = 0;
                    RphyUp[comp] = 0;
                    Rphz[comp]=0;
                    RphzTop[comp]=0;
                }
                for(ion = 0; ion < Ni; ion++){
                    for(comp = 0; comp < Nc-1; comp++){
                        Rcvx = 0;
                        RcvxRight = 0;
                        if(x > 0){
                            //First difference term
                            Rcvx = Dcs[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)*3]*(cp[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)]+cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2;
                            Rcvx = Rcvx*(log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-log(c[c_index(
                                    x-1,y,z,comp,ion,Nx,Ny,Nz)])+
                                    z_charge[ion]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-
                                                   phi[phi_index(x-1,y,z,comp,Nx,Ny,Nz)]))/dx*dt/dx;
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
                                    z_charge[ion]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-phi[phi_index(x,y-1,z,comp,Nx,Ny,Nz)]))/dy*dt/dy;
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
                        //Top Bottom difference
                        if(z > 0){
                            Rcvz = Dcs[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)]+
                                                                               cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2;
                            Rcvz = Rcvz*(log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-log(c[c_index(x,y,
                                                                                               z-1,comp,ion,Nx,Ny,Nz)])+
                                         z_charge[ion]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-phi[phi_index(x,y,
                                                                                                         z-1,comp,Nx,Ny,Nz)]))/dz*dt/dz;
                        }
                        //Next upward difference
                        if(z < Nz-1){
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

                        ierr = VecSetValue(Res,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Resc,INSERT_VALUES);
                        CHKERRQ(ierr);

                        //Save values for voltage
                        Rphx[comp] += z_charge[ion]*Rcvx;
                        Rphy[comp] += z_charge[ion]*Rcvy;
                        Rphz[comp] += z_charge[ion]*Rcvz;
                        RphxRight[comp] += z_charge[ion]*RcvxRight;
                        RphyUp[comp] += z_charge[ion]*RcvyUp;
                        RphzTop[comp] += z_charge[ion]*RcvzTop;

                    }
                    //Set Extracellular values
                    alNc = 1-al[al_index(x,y,z,0,Nx,Ny,Nz)]-al[al_index(x,y,z,1,Nx,Ny,Nz)];
                    alpNc = 1-alp[al_index(x,y,z,0,Nx,Ny,Nz)]-alp[al_index(x,y,z,1,Nx,Ny,Nz)];
                    comp = Nc-1;
                    Rcvx = 0;
                    RcvxRight = 0;
                    if(x > 0){
                        //First difference term
                        Rcvx = Dcs[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)*3]*(cp[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)]+cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2;
                        Rcvx = Rcvx*(log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-log(c[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)])+
                                     z_charge[ion]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-
                                                    phi[phi_index(x-1,y,z,comp,Nx,Ny,Nz)]))/dx*dt/dx;
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
                                     z_charge[ion]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-phi[phi_index(x,y-1,z,comp,Nx,Ny,Nz)]))/dy*dt/dy;
                    }
                    //Next upward difference
                    if(y < Ny-1){
                        RcvyUp = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+1]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                           cp[c_index(x,y+1,z,comp,ion,Nx,Ny,Nz)])/2;
                        RcvyUp = RcvyUp*(log(c[c_index(x,y+1,z,comp,ion,Nx,Ny,Nz)])-
                                         log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])+z_charge[ion]*
                                                                                 (phi[phi_index(x,y+1,z,comp,Nx,Ny,Nz)]-phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]))/dy*dt/dy;
                    }
                    Rcvz = 0;
                    RcvzTop = 0;
                    //Top Bottom difference
                    if(z > 0){
                        Rcvz = Dcs[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)]+
                                                                           cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2;
                        Rcvz = Rcvz*(log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-log(c[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)])+
                                     z_charge[ion]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-phi[phi_index(x,y,z-1,comp,Nx,Ny,Nz)]))/dz*dt/dz;
                    }
                    //Next upward difference
                    if(z < Nz-1){
                        RcvzTop = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                            cp[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])/2;
                        RcvzTop = RcvzTop*(log(c[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])-
                                           log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])+z_charge[ion]*
                                                                                   (phi[phi_index(x,y,
                                                                                                  z+1,comp,Nx,Ny,Nz)]-phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]))/dz*dt/dz;
                    }

                    Resc = alNc*c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]-
                           alpNc*cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                    Resc += Rcvx-RcvxRight+Rcvy-RcvyUp+Rcvz-RcvzTop+flux->mflux[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;
                    //Add bath variables

                    Resc -= sqrt(pow(Dcb[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3],2)+
                                 pow(Dcb[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+1],2)+
                                 pow(Dcb[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2],2))*
                            (cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+cbath[ion])/2.0*
                            (log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-log(cbath[ion])+
                            z_charge[ion]*phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-z_charge[ion]*phibath)*dt;
                    ierr = VecSetValue(Res,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Resc,INSERT_VALUES);
                    CHKERRQ(ierr);

                    //Save values for voltage
                    Rphx[comp] += z_charge[ion]*Rcvx;
                    Rphy[comp] += z_charge[ion]*Rcvy;
                    Rphz[comp] += z_charge[ion]*Rcvz;
                    RphxRight[comp] += z_charge[ion]*RcvxRight;
                    RphyUp[comp] += z_charge[ion]*RcvyUp;
                    RphzTop[comp] += z_charge[ion]*RcvzTop;
                }

                //Voltage Equations
                ResphN = 0;
                for(comp = 0; comp < Nc-1; comp++){
                    Resph = cm[al_index(x,y,z,comp,Nx,Ny,Nz)]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-phi[phi_index(x,y,z,
                                                                                                             Nc-1,Nx,Ny,Nz)])-
                            cm[al_index(x,y,z,comp,Nx,Ny,Nz)]*(phip[phi_index(x,y,z,comp,Nx,Ny,Nz)]-phip[phi_index(x,y,z,Nc-
                                                                                                                       1,Nx,Ny,Nz)]);
                    for(ion = 0; ion < Ni; ion++){
                        //Ion channel
                        Resph += z_charge[ion]*flux->mflux[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;
                    }
                    //Add the terms shared with extracell
                    ResphN -= Resph; // Subtract total capacitance, subtract total ion channel flux
                    Resph += Rphx[comp]-RphxRight[comp]+Rphy[comp]-RphyUp[comp]+Rphz[comp]-RphzTop[comp];
                    ierr = VecSetValue(Res,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Resph,INSERT_VALUES);
                    CHKERRQ(ierr);
                }

                //Finish adding extracell
                comp = Nc-1;
                //Add bath contribution
                for(ion = 0; ion < Ni; ion++){

                    ResphN -= z_charge[ion]*sqrt(pow(Dcb[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3],2)+
                                                 pow(Dcb[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+1],2)+
                                                 pow(Dcb[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2],2))*
                              (cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+cbath[ion])/2.0*
                              (log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-log(cbath[ion])+
                                    z_charge[ion]*phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-z_charge[ion]*phibath)*dt;
                }
                ResphN += Rphx[comp]-RphxRight[comp]+Rphy[comp]-RphyUp[comp]+Rphz[comp]-RphzTop[comp];
                ierr = VecSetValue(Res,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),ResphN,INSERT_VALUES);
                CHKERRQ(ierr);
            }
        }
    }

    for(z=0;z<Nz;z++){
        for(y = 0; y < Ny; y++){
            for(x = 0; x < Nx; x++){

                //Residual for water flow
                //Plus modification to electroneutrality for non-zero mem.compacitance
                for(comp = 0; comp < Nc-1; comp++){
                    //Water flow
                    ierr = VecSetValue(Res,Ind_1(x,y,z,Ni+1,comp,Nx,Ny,Nz),al[al_index(x,y,z,comp,Nx,Ny,Nz)]
                                                                          -alp[al_index(x,y,z,comp,Nx,Ny,Nz)]+
                                                                        flux->wflow[al_index(x,y,z,comp,Nx,Ny,Nz)]*
                                                                        dt,INSERT_VALUES);
                    CHKERRQ(ierr);

                }
            }
        }
    }
    //Assemble before we add values in on top to modify the electroneutral.
    ierr = VecAssemblyBegin(Res);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Res);CHKERRQ(ierr);
    ierr = restore_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);

    if(Profiling_on) {
        PetscLogEventEnd(event[1], 0, 0, 0, 0);
    }
    return ierr;
}

PetscErrorCode
calc_jacobian(SNES snes,Vec current_state, Mat A, Mat Jac,void *ctx)
{
    //Jacobian equation using derivative of the charge-capacitance relation
    // Volume is solved for here
    struct AppCtx * user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if(Profiling_on) {
        PetscLogEventBegin(event[0], 0, 0, 0, 0);
    }
    ierr = extract_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    PetscReal *c = user->state_vars->c;
    PetscReal *al = user->state_vars->alpha;
    PetscReal *cp = user->state_vars_past->c;

    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct FluxData *flux = user->flux;
    PetscReal dt = user->dt;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    PetscReal dz = user->dz;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt Nz = user->Nz;
    struct ConstVars *con_vars = user->con_vars;
    PetscReal *cm = con_vars->cm;

    PetscInt ind = 0;
    PetscInt x,y,z,ion,comp;

    PetscReal Ftmpx,Fc0x,Fc1x,Fph0x,Fph1x;
    PetscReal Ftmpy,Fc0y,Fc1y,Fph0y,Fph1y;
    PetscReal Ftmpz,Fc0z,Fc1z,Fph0z,Fph1z;
    PetscReal Ac,Aphi,Avolt,AvoltN;

    PetscReal Fphph0x[Nc],Fphph1x[Nc];
    PetscReal Fphph0y[Nc],Fphph1y[Nc];
    PetscReal Fphph0z[Nc],Fphph1z[Nc];
    //Ionic concentration equations
    for(z=0;z<Nz;z++){
        for(y = 0; y < Ny; y++){
            for(x = 0; x < Nx; x++){

                for(comp = 0; comp < Nc; comp++){
                    Fphph0x[comp] = 0;
                    Fphph1x[comp] = 0;
                    Fphph0y[comp] = 0;
                    Fphph1y[comp] = 0;
                    Fphph0z[comp] = 0;
                    Fphph1z[comp] = 0;
                }
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

                            ierr = MatSetValue(Jac,Ind_1(x+1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc0x,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Right c with left phi (-Fph0x)
                            ierr = MatSetValue(Jac,Ind_1(x+1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                               -Fph0x,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;

                            //Right phi with left c in voltage eqn
                            ierr = MatSetValue(Jac,Ind_1(x+1,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -z_charge[ion]*Fc0x,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                        if(x > 0){
                            Ftmpx = Dcs[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)*3]*(cp[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                              cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2/dx*dt/dx;
                            Fc1x = Ftmpx/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                            Fph1x = z_charge[ion]*Ftmpx;
                            //left c with right c (-Fc1x)
                            ierr = MatSetValue(Jac,Ind_1(x-1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc1x,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Left c with right phi (-Fph1x)
                            ierr = MatSetValue(Jac,Ind_1(x-1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                               -Fph1x,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;

                            //Left phi with right c in voltage eqn
                            ierr = MatSetValue(Jac,Ind_1(x-1,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -z_charge[ion]*Fc1x,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                        if(y < Ny-1){
                            Ftmpy = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+1]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                              cp[c_index(x,y+1,z,comp,ion,Nx,Ny,Nz)])/2/dy*dt/dy;
                            Fc0y = Ftmpy/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                            Fph0y = z_charge[ion]*Ftmpy;
                            // Upper c with lower c (-Fc0y)
                            ierr = MatSetValue(Jac,Ind_1(x,y+1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc0y,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Upper c with lower phi (-Fph0y)
                            ierr = MatSetValue(Jac,Ind_1(x,y+1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                               -Fph0y,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;

                            //Upper phi with lower c in voltage eqn
                            ierr = MatSetValue(Jac,Ind_1(x,y+1,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -z_charge[ion]*Fc0y,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                        if(y > 0){
                            Ftmpy = Dcs[c_index(x,y-1,z,comp,ion,Nx,Ny,Nz)*3+1]*(cp[c_index(x,y-1,z,comp,ion,Nx,Ny,Nz)]+
                                                                                cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2/dy*dt/dy;
                            Fc1y = Ftmpy/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                            Fph1y = z_charge[ion]*Ftmpy;
                            //Lower c with Upper c (-Fc1y)
                            ierr = MatSetValue(Jac,Ind_1(x,y-1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc1y,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Lower c with Upper phi (-Fph1y)
                            ierr = MatSetValue(Jac,Ind_1(x,y-1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                               -Fph1y,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;

                            //Lower phi with upper c in voltage eqn
                            ierr = MatSetValue(Jac,Ind_1(x,y-1,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -z_charge[ion]*Fc1y,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                        if(z < Nz-1){
                            Ftmpz = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                              cp[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])/2/dz*dt/dz;
                            Fc0z = Ftmpz/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                            Fph0z = z_charge[ion]*Ftmpz;
                            // Bot c with Top c (-Fc0x)

                            ierr = MatSetValue(Jac,Ind_1(x,y,z+1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc0z,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Bot c with Top phi (-Fph0x)
                            ierr = MatSetValue(Jac,Ind_1(x,y,z+1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                               -Fph0z,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;

                            //Bot phi with Top c in voltage eqn
                            ierr = MatSetValue(Jac,Ind_1(x,y,z+1,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -z_charge[ion]*Fc0z,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                        if(z > 0){
                            Ftmpz = Dcs[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)]+
                                                                                cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2/dz*dt/dz;
                            Fc1z = Ftmpz/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                            Fph1z = z_charge[ion]*Ftmpz;
                            //Top c Bot right c (-Fc1x)
                            ierr = MatSetValue(Jac,Ind_1(x,y,z-1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc1z,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Top c Bot right phi (-Fph1x)
                            ierr = MatSetValue(Jac,Ind_1(x,y,z-1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                               -Fph1z,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;

                            //Left phi with right c in voltage eqn
                            ierr = MatSetValue(Jac,Ind_1(x,y,z-1,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -z_charge[ion]*Fc1z,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                        //Diagonal term contribution
                        Ac = al[al_index(x,y,z,comp,Nx,Ny,Nz)]+Fc0x+Fc1x+Fc0y+Fc1y+Fc0z+Fc1z;
                        Aphi = Fph0x+Fph1x+Fph0y+Fph1y+Fph0z+Fph1z;

                        //Add up terms for voltage eqns
                        Fphph0x[comp] += z_charge[ion]*Fph0x;
                        Fphph1x[comp] += z_charge[ion]*Fph1x;
                        Fphph0y[comp] += z_charge[ion]*Fph0y;
                        Fphph1y[comp] += z_charge[ion]*Fph1y;
                        Fphph0z[comp] += z_charge[ion]*Fph0z;
                        Fphph1z[comp] += z_charge[ion]*Fph1z;

                        //membrane current contributions
                        Ac += flux->dfdci[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;
                        Aphi += flux->dfdphim[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;
                        // Different Compartment Terms
                        // C Extracellular with C Inside
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -flux->dfdci[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // C Intra with C Extra
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),
                                           flux->dfdce[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // C Extracellular with Phi Inside
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -flux->dfdphim[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // C Intra with Phi Extra
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,Nc-1,Nx,Ny,Nz),
                                           -flux->dfdphim[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Volume terms
                        //C extra with intra alpha
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,Ni+1,comp,Nx,Ny,Nz),
                                           -c[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)],INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //C intra with intra alpha
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni+1,comp,Nx,Ny,Nz),
                                           c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)],INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Same compartment terms
                        // c with c
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Ac,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // c with phi
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Aphi,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                        //Intra-Phi with c (voltage eqn)
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           z_charge[ion]*(Fc0x+Fc1x+Fc0y+Fc1y+
                                                          flux->dfdci[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt),INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //IntraPhi with c extra(volt eqn)
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),
                                           z_charge[ion]*(flux->dfdce[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt),INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Extra-Phi with intra-c (voltage eqn)
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -z_charge[ion]*(flux->dfdci[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt),INSERT_VALUES);
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

                        ierr = MatSetValue(Jac,Ind_1(x+1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc0x,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Right c with left phi (-Fph0x)
                        ierr = MatSetValue(Jac,Ind_1(x+1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fph0x,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                        //Right phi with left c in voltage eqn
                        ierr = MatSetValue(Jac,Ind_1(x+1,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -z_charge[ion]*Fc0x,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(x > 0){
                        Ftmpx = Dcs[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)*3]*(cp[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                          cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2/dx*dt/dx;
                        Fc1x = Ftmpx/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                        Fph1x = z_charge[ion]*Ftmpx;
                        //left c with right c (-Fc1x)
                        ierr = MatSetValue(Jac,Ind_1(x-1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc1x,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Left c with right phi (-Fph1x)
                        ierr = MatSetValue(Jac,Ind_1(x-1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fph1x,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                        //Left phi with right c in voltage eqn
                        ierr = MatSetValue(Jac,Ind_1(x-1,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -z_charge[ion]*Fc1x,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(y < Ny-1){
                        Ftmpy = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+1]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                          cp[c_index(x,y+1,z,comp,ion,Nx,Ny,Nz)])/2/dy*dt/dy;
                        Fc0y = Ftmpy/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                        Fph0y = z_charge[ion]*Ftmpy;
                        // Upper c with lower c (-Fc0y)
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc0y,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Upper c with lower phi (-Fph0y)
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fph0y,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                        //Upper phi with lower c in voltage eqn
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -z_charge[ion]*Fc0y,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(y > 0){
                        Ftmpy = Dcs[c_index(x,y-1,z,comp,ion,Nx,Ny,Nz)*3+1]*(cp[c_index(x,y-1,z,comp,ion,Nx,Ny,Nz)]+
                                                                            cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2/dy*dt/dy;
                        Fc1y = Ftmpy/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                        Fph1y = z_charge[ion]*Ftmpy;
                        //Lower c with Upper c (-Fc1y)
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc1y,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Lower c with Upper phi (-Fph1y)
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fph1y,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                        //Lower phi with upper c in voltage eqn
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -z_charge[ion]*Fc1y,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(z < Nz-1){
                        Ftmpz = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                          cp[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])/2/dz*dt/dz;
                        Fc0z = Ftmpz/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                        Fph0z = z_charge[ion]*Ftmpz;
                        // Bot c with Top c (-Fc0x)

                        ierr = MatSetValue(Jac,Ind_1(x,y,z+1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc0z,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Bot c with Top phi (-Fph0x)
                        ierr = MatSetValue(Jac,Ind_1(x,y,z+1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fph0z,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                        //Bot phi with Top c in voltage eqn
                        ierr = MatSetValue(Jac,Ind_1(x,y,z+1,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -z_charge[ion]*Fc0z,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(z > 0){
                        Ftmpz = Dcs[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)]+
                                                                            cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2/dz*dt/dz;
                        Fc1z = Ftmpz/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                        Fph1z = z_charge[ion]*Ftmpz;
                        //Top c Bot right c (-Fc1x)
                        ierr = MatSetValue(Jac,Ind_1(x,y,z-1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc1z,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Top c Bot right phi (-Fph1x)
                        ierr = MatSetValue(Jac,Ind_1(x,y,z-1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fph1z,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                        //Left phi with right c in voltage eqn
                        ierr = MatSetValue(Jac,Ind_1(x,y,z-1,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -z_charge[ion]*Fc1z,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }

                    //Diagonal term contribution
                    Ac = (1-al[al_index(x,y,z,0,Nx,Ny,Nz)]-al[al_index(x,y,z,1,Nx,Ny,Nz)])+Fc0x+Fc1x+Fc0y+Fc1y+Fc0z+Fc1z;
                    Aphi = Fph0x+Fph1x+Fph0y+Fph1y+Fph0z+Fph1z;

                    Avolt = z_charge[ion]*(Fc0x+Fc1x+Fc0y+Fc1y+Fc0z+Fc1z);

                    //Add up terms for voltage eqns
                    Fphph0x[comp] += z_charge[ion]*Fph0x;
                    Fphph1x[comp] += z_charge[ion]*Fph1x;
                    Fphph0y[comp] += z_charge[ion]*Fph0y;
                    Fphph1y[comp] += z_charge[ion]*Fph1y;
                    Fphph0z[comp] += z_charge[ion]*Fph0z;
                    Fphph1z[comp] += z_charge[ion]*Fph1z;

                    //Membrane current contribution
                    for(comp = 0; comp < Nc-1; comp++){
                        Ac -= flux->dfdce[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;
                        Aphi += flux->dfdphim[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;
                        Avolt -= z_charge[ion]*flux->dfdce[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;
                    }
                    //Add bath contributions
                    Ftmpx = sqrt(pow(Dcb[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)*3],2)+
                                 pow(Dcb[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)*3+1],2)+
                                 pow(Dcb[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)*3+2],2));
                    Ac -= Ftmpx*(cp[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)]+cbath[ion])/(2*c[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)])*dt;
                    Aphi -= Ftmpx*(cp[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)]+cbath[ion])*z_charge[ion]/2*dt;

                    Avolt -= z_charge[ion]*Ftmpx*(cp[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)]+cbath[ion])/
                             (2*c[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)])*dt;

                    //Insert extracell to extracell parts
                    // c with c
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ac,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    // c with phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,Ni,Nc-1,Nx,Ny,Nz),Aphi,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;

                    //phi with c (voltage eqn)
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,ion,
                                                                              Nc-1,Nx,Ny,Nz),Avolt,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                //Derivative of charge-capacitance
                for(comp = 0; comp < Nc-1; comp++){
                    if(x < Nx-1){
                        //Right phi with left phi (-Fph0x)
                        ierr = MatSetValue(Jac,Ind_1(x+1,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fphph0x[comp],INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(x > 0){
                        //Left phi with right phi (-Fph1x)
                        ierr = MatSetValue(Jac,Ind_1(x-1,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fphph1x[comp],INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(y < Ny-1){
                        //Upper phi with lower phi (-Fph0y)
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fphph0y[comp],INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(y > 0){
                        //Lower phi with upper phi (-Fph1y)
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fphph1y[comp],INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(z < Nz-1){
                        //Bot phi with Top phi (-Fph0x)
                        ierr = MatSetValue(Jac,Ind_1(x,y,z+1,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fphph0z[comp],INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(z > 0){
                        //Top phi with Bot phi (-Fph1x)
                        ierr = MatSetValue(Jac,Ind_1(x,y,z-1,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fphph1z[comp],INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    Avolt = cm[al_index(x,y,z,comp,Nx,Ny,Nz)]+Fphph0x[comp]+Fphph1x[comp]+Fphph0y[comp]+Fphph1y[comp]+
                            Fphph0z[comp]+Fphph1z[comp];
                    AvoltN = -cm[al_index(x,y,z,comp,Nx,Ny,Nz)];
                    for(ion = 0; ion < Ni; ion++){
                        Avolt += z_charge[ion]*flux->dfdphim[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;
                        AvoltN -= z_charge[ion]*flux->dfdphim[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;
                    }

                    //Intra-phi with Intra-phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Avolt,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Intra-phi with extra-phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,
                                                                              Nc-1,Nx,Ny,Nz),AvoltN,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                //Extracellular terms
                comp = Nc-1;
                if(x < Nx-1){
                    //Right phi with left phi (-Fph0x)
                    ierr = MatSetValue(Jac,Ind_1(x+1,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                       -Fphph0x[comp],INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if(x > 0){
                    //Left phi with right phi (-Fph1x)
                    ierr = MatSetValue(Jac,Ind_1(x-1,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                       -Fphph1x[comp],INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if(y < Ny-1){
                    //Upper phi with lower phi (-Fph0y)
                    ierr = MatSetValue(Jac,Ind_1(x,y+1,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                       -Fphph0y[comp],INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if(y > 0){
                    //Lower phi with upper phi (-Fph1y)
                    ierr = MatSetValue(Jac,Ind_1(x,y-1,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                       -Fphph1y[comp],INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if(z < Nz-1){
                    //Bot phi with Top phi (-Fph0x)
                    ierr = MatSetValue(Jac,Ind_1(x,y,z+1,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                       -Fphph0z[comp],INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if(z > 0){
                    //Top phi with Bot phi (-Fph1x)
                    ierr = MatSetValue(Jac,Ind_1(x,y,z-1,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                       -Fphph1z[comp],INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                AvoltN = 0;

                for(int k = 0; k < Nc-1; k++){
                    AvoltN += cm[al_index(x,y,z,k,Nx,Ny,Nz)];
                    Avolt = -cm[al_index(x,y,z,k,Nx,Ny,Nz)];
                    for(ion = 0; ion < Ni; ion++){
                        Avolt -= z_charge[ion]*flux->dfdphim[c_index(x,y,z,k,ion,Nx,Ny,Nz)]*dt;
                        AvoltN += z_charge[ion]*flux->dfdphim[c_index(x,y,z,k,ion,Nx,Ny,Nz)]*dt;
                    }
                    //Extra-phi with Intra-phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,k,Nx,Ny,Nz),Avolt,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }

                AvoltN += Fphph0x[comp]+Fphph1x[comp]+Fphph0y[comp]+Fphph1y[comp]+Fphph0z[comp]+Fphph1z[comp];;

                //Bath terms
                for(ion = 0; ion < Ni; ion++){
                    Ftmpx = sqrt(pow(Dcb[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)*3],2)+
                                 pow(Dcb[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)*3+1],2)+
                                 pow(Dcb[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)*3+2],2));
                    AvoltN -= z_charge[ion]*Ftmpx*(cp[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)]+cbath[ion])*z_charge[ion]/2*dt;
                }
                //extra-phi with extra-phi
                ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),AvoltN,INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;

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
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni+1,comp,Nx,Ny,Nz),Ind_1(x,y,z,
                                                                                Ni+1,comp,Nx,Ny,Nz),Ac,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Off diagonal (from aNc=1-sum(ak))
                    for(PetscInt l = 0; l < comp; l++){
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni+1,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni+1,l,Nx,Ny,Nz),
                                           flux->dwdpi[al_index(x,y,z,comp,Nx,Ny,Nz)]*
                                           con_vars->ao[phi_index(x,y,z,Nc-1,Nx,Ny,Nz)]/
                                           pow(1-al[al_index(x,y,z,0,Nx,Ny,Nz)]-al[al_index(x,y,z,1,Nx,Ny,Nz)],2)*dt,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    for(PetscInt l = comp+1; l < Nc-1; l++){
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni+1,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni+1,l,Nx,Ny,Nz),
                                           flux->dwdpi[al_index(x,y,z,comp,Nx,Ny,Nz)]*
                                           con_vars->ao[phi_index(x,y,z,Nc-1,Nx,Ny,Nz)]/
                                           ((1-al[al_index(x,y,z,0,Nx,Ny,Nz)]-al[al_index(x,y,z,1,Nx,Ny,Nz)])*
                                            (1-al[al_index(x,y,z,0,Nx,Ny,Nz)]-al[al_index(x,y,z,1,Nx,Ny,Nz)]))*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    for(ion = 0; ion < Ni; ion++){
                        //Volume to extra c
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni+1,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),
                                           flux->dwdpi[al_index(x,y,z,comp,Nx,Ny,Nz)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Volume to intra c
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni+1,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
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

//    MatView(Jac,PETSC_VIEWER_STDOUT_SELF);


    if (A != Jac) {
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr); }

    ierr = restore_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    if(Profiling_on) {
        PetscLogEventEnd(event[0], 0, 0, 0, 0);
    }
    return ierr;
}

void volume_update(struct SimState *state_vars,struct SimState *state_vars_past, struct AppCtx *user)
{
    if(Profiling_on) {
        PetscLogEventBegin(event[7], 0, 0, 0, 0);
    }
    int x,y,z,comp;
    PetscReal dt = user->dt;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt Nz = user->Nz;
    for(int n=0;n<1;n++) {

//        memcpy(state_vars_past->alpha, state_vars->alpha, sizeof(PetscReal) * Nx * Ny * (Nc - 1));
        //Forward Euler update
/*
    wflowm(user->flux,user->state_vars_past,user->con_vars);
    for(x=0;x<Nx;x++){
        for(y=0;y<Ny;y++){
            for(comp=0;comp<Nc-1;comp++) {
                state_vars->alpha[al_index(x, y, comp,Nx)] = state_vars_past->alpha[al_index(x,y,comp,Nx)]+user->dt*user->flux->wflow[al_index(x,y,comp,Nx)];
            }
        }
    }
*/
//    /*
        //Backward Euler update
        PetscReal res, Func, Deriv, max_res;
        for (int iter = 0; iter < 10; iter++) {
            max_res = 0;
            wflowm(user);
            for (z = 0; z < Nz; z++){
                for(y = 0; y < Ny; y++){
                    for(x = 0; x < Nx; x++){
                        for(comp = 0; comp < Nc-1; comp++){

                            Func = state_vars->alpha[al_index(x,y,z,comp,Nx,Ny,Nz)]-
                                   state_vars_past->alpha[al_index(x,y,z,comp,Nx,Ny,Nz)]+
                                    dt*user->flux->wflow[al_index(x,y,z,comp,Nx,Ny,Nz)];

                            Deriv = 1+dt*user->flux->dwdal[al_index(x,y,z,comp,Nx,Ny,Nz)];

                            res = -Func/Deriv;

                            state_vars->alpha[al_index(x,y,z,comp,Nx,Ny,Nz)] += res;

                            if(fabs(res) > max_res){ max_res = fabs(res); }
                        }
                    }
                }
            }
            if (max_res < reltol) {
                if(Profiling_on) {
                    PetscLogEventEnd(event[7], 0, 0, 0, 0);
                }
                return; }
        }
    }
    fprintf(stderr,"Volume failed to update!\n");
    if(Profiling_on) {
        PetscLogEventEnd(event[7], 0, 0, 0, 0);
    }
    exit(EXIT_FAILURE); /* indicate failure.*/
}

PetscErrorCode calc_residual_no_vol(SNES snes,Vec current_state,Vec Res,void *ctx)
{
    //Residual equation using derivative of the charge-capacitance relation
    // Volume not solved for here
    struct AppCtx * user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if(Profiling_on) {
        PetscLogEventBegin(event[1], 0, 0, 0, 0);
    }
    ierr = extract_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    //Compute membrane ionic flux relation quantitites
    ionmflux(user);

    //Compute membrane water flow related quantities
    wflowm(user);

    PetscReal *c = user->state_vars->c;
    PetscReal *phi = user->state_vars->phi;
    PetscReal *al = user->state_vars->alpha;
    PetscReal *cp = user->state_vars_past->c;
    PetscReal *alp = user->state_vars_past->alpha;
    PetscReal *phip = user->state_vars_past->phi;

    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct FluxData *flux = user->flux;
    PetscReal *cm = user->con_vars->cm;
    PetscReal dt = user->dt;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    PetscReal dz = user->dz;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt Nz = user->Nz;

    //Residual for concentration equations
    PetscReal Rcvx,Rcvy,Rcvz,Resc;
    PetscReal RcvxRight,RcvyUp,RcvzTop;

    //Residual for fluxes in voltage differential equations
    PetscReal Rphx[Nc], Rphy[Nc], RphxRight[Nc], RphyUp[Nc],Rphz[Nc],RphzTop[Nc];
    PetscReal Resph,ResphN;

    PetscReal alNc,alpNc;
    PetscInt ion,comp,x,y,z;

    for(z=0;z<Nz;z++){
        for(y = 0; y < Ny; y++){
            for(x = 0; x < Nx; x++){
                //Init voltage tracking to zero
                for(comp = 0; comp < Nc; comp++){
                    Rphx[comp] = 0;
                    Rphy[comp] = 0;
                    RphxRight[comp] = 0;
                    RphyUp[comp] = 0;
                    Rphz[comp] = 0;
                    RphzTop[comp] = 0;
                }
                for(ion = 0; ion < Ni; ion++){
                    for(comp = 0; comp < Nc-1; comp++){
                        Rcvx = 0;
                        RcvxRight = 0;
                        if(x > 0){
                            //First difference term
                            Rcvx = Dcs[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)*3]*(cp[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                             cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2;
                            Rcvx = Rcvx*(log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-log(c[c_index(
                                    x-1,y,z,comp,ion,Nx,Ny,Nz)])+
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
                                    z_charge[ion]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-phi[phi_index(x,y-1,z,comp,Nx,Ny,Nz)]))/dy*dt/dy;
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
                        //Top Bottom difference
                        if(z > 0){
                            Rcvz = Dcs[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)]+
                                                                               cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2;
                            Rcvz = Rcvz*(log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-log(c[c_index(x,y,
                                                                                               z-1,comp,ion,Nx,Ny,Nz)])+
                                         z_charge[ion]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-phi[phi_index(x,y,
                                                                                                         z-1,comp,Nx,Ny,Nz)]))/dz*dt/dz;
                        }
                        //Next upward difference
                        if(z < Nz-1){
                            RcvzTop = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                                cp[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])/2;
                            RcvzTop = RcvzTop*(log(c[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])-
                                               log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])+
                                             z_charge[ion]*(phi[phi_index(x,y,z+1,comp,Nx,Ny,Nz)]-
                                                            phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]))/dz*dt/dz;
                        }
                        Resc = al[al_index(x,y,z,comp,Nx,Ny,Nz)]*c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]-
                               alp[al_index(x,y,z,comp,Nx,Ny,Nz)]*cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                        Resc += Rcvx-RcvxRight+Rcvy-RcvyUp+Rcvz-RcvzTop+flux->mflux[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;

                        ierr = VecSetValue(Res,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Resc,INSERT_VALUES);
                        CHKERRQ(ierr);

                        //Save values for voltage
                        Rphx[comp] += z_charge[ion]*Rcvx;
                        Rphy[comp] += z_charge[ion]*Rcvy;
                        Rphz[comp] += z_charge[ion]*Rcvz;
                        RphxRight[comp] += z_charge[ion]*RcvxRight;
                        RphyUp[comp] += z_charge[ion]*RcvyUp;
                        RphzTop[comp] += z_charge[ion]*RcvzTop;

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
                                               log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])+z_charge[ion]*(phi[phi_index(
                                x+1,y,z,comp,Nx,Ny,Nz)]-
                                                                                                      phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]))/dx*dt/dx;
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
                                z_charge[ion]*(phi[phi_index(x,y+1,z,comp,Nx,Ny,Nz)]-phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]))/dy*dt/dy;
                    }
                    Rcvz = 0;
                    RcvzTop = 0;
                    //Top bot difference
                    if(z > 0){
                        Rcvz = Dcs[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)]+
                                                                           cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2;
                        Rcvz = Rcvz*(log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-log(c[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)])+
                                     z_charge[ion]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-phi[phi_index(x,y,z-1,comp,Nx,Ny,Nz)]))/dz*dt/dz;
                    }
                    //Next upward difference
                    if(z < Nz-1){
                        RcvzTop = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                            cp[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])/2;
                        RcvzTop = RcvzTop*(log(c[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])-log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])+
                                         z_charge[ion]*(phi[phi_index(x,y,z+1,comp,Nx,Ny,Nz)]-
                                                        phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]))/dz*dt/dz;
                    }
                    Resc = alNc*c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]-alpNc*cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                    Resc += Rcvx-RcvxRight+Rcvy-RcvyUp+Rcvz-RcvzTop+flux->mflux[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;
                    //Add bath variables

                    Resc -= sqrt(pow(Dcb[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3],2)+
                                 pow(Dcb[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+1],2)+
                                 pow(Dcb[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2],2))*
                            (cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+cbath[ion])/2.0*
                            (log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-log(cbath[ion])+
                                    z_charge[ion]*phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-z_charge[ion]*phibath)*dt;
                    ierr = VecSetValue(Res,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Resc,INSERT_VALUES);
                    CHKERRQ(ierr);

                    //Save values for voltage
                    Rphx[comp] += z_charge[ion]*Rcvx;
                    Rphy[comp] += z_charge[ion]*Rcvy;
                    Rphz[comp] += z_charge[ion]*Rcvz;
                    RphxRight[comp] += z_charge[ion]*RcvxRight;
                    RphyUp[comp] += z_charge[ion]*RcvyUp;
                    RphzTop[comp] += z_charge[ion]*RcvzTop;
                }

                //Voltage Equations
                ResphN = 0;
                for(comp = 0; comp < Nc-1; comp++){
                    Resph = cm[al_index(x,y,z,comp,Nx,Ny,Nz)]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-phi[phi_index(x,y,z,
                                                                                                             Nc-1,Nx,Ny,Nz)])-
                            cm[al_index(x,y,z,comp,Nx,Ny,Nz)]*(phip[phi_index(x,y,z,comp,Nx,Ny,Nz)]-phip[phi_index(x,y,z,Nc-
                                                                                                                       1,Nx,Ny,Nz)]);
                    for(ion = 0; ion < Ni; ion++){
                        //Ion channel
                        Resph += z_charge[ion]*flux->mflux[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;
                    }
                    //Add the terms shared with extracell
                    ResphN -= Resph; // Subtract total capacitance, subtract total ion channel flux
                    Resph += Rphx[comp]-RphxRight[comp]+Rphy[comp]-RphyUp[comp]+Rphz[comp]-RphzTop[comp];
                    ierr = VecSetValue(Res,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Resph,INSERT_VALUES);
                    CHKERRQ(ierr);
                }

                //Finish adding extracell
                comp = Nc-1;
                //Add bath contribution
                for(ion = 0; ion < Ni; ion++){
                    ResphN -= z_charge[ion]*sqrt(pow(Dcb[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3],2)+
                                                 pow(Dcb[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+1],2)+
                                                 pow(Dcb[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2],2))*
                              (cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+cbath[ion])/2.0*
                              (log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-log(cbath[ion])+
                                            z_charge[ion]*phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-z_charge[ion]*phibath)*dt;
                }
                ResphN += Rphx[comp]-RphxRight[comp]+Rphy[comp]-RphyUp[comp]+Rphz[comp]-RphzTop[comp];
                ierr = VecSetValue(Res,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),ResphN,INSERT_VALUES);
                CHKERRQ(ierr);
            }
        }
    }

    ierr = VecAssemblyBegin(Res);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Res);CHKERRQ(ierr);
    ierr = restore_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    if(Profiling_on) {
        PetscLogEventEnd(event[1], 0, 0, 0, 0);
    }
    return ierr;
}

PetscErrorCode
calc_jacobian_no_vol(SNES snes,Vec current_state, Mat A, Mat Jac,void *ctx)
{
    //Jacobian equation using derivative of the charge-capacitance relation
    // Alpha is not solved here

    struct AppCtx * user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if(Profiling_on) {
        PetscLogEventBegin(event[0], 0, 0, 0, 0);
    }
    ierr = extract_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    PetscReal *c = user->state_vars->c;
    PetscReal *al = user->state_vars->alpha;
    PetscReal *cp = user->state_vars_past->c;

    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct FluxData *flux = user->flux;
    PetscReal dt = user->dt;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    PetscReal dz = user->dz;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt Nz = user->Nz;
    struct ConstVars *con_vars = user->con_vars;
    PetscReal *cm = con_vars->cm;

    PetscInt ind = 0;
    PetscInt x,y,z,ion,comp;

    PetscReal Ftmpx,Fc0x,Fc1x,Fph0x,Fph1x;
    PetscReal Ftmpy,Fc0y,Fc1y,Fph0y,Fph1y;
    PetscReal Ftmpz,Fc0z,Fc1z,Fph0z,Fph1z;
    PetscReal Ac,Aphi,Avolt,AvoltN;

    PetscReal Fphph0x[Nc],Fphph1x[Nc];
    PetscReal Fphph0y[Nc],Fphph1y[Nc];
    PetscReal Fphph0z[Nc],Fphph1z[Nc];

    //Ionic concentration equations
    for(z=0;z<Nz;z++){
        for(y = 0; y < Ny; y++){
            for(x = 0; x < Nx; x++){
                for(comp = 0; comp < Nc; comp++){
                    Fphph0x[comp] = 0;
                    Fphph1x[comp] = 0;
                    Fphph0y[comp] = 0;
                    Fphph1y[comp] = 0;
                    Fphph0z[comp] = 0;
                    Fphph1z[comp] = 0;
                }
                for(ion = 0; ion < Ni; ion++){
                    for(comp = 0; comp < Nc-1; comp++){
                        //Electrodiffusion contributions
                        Ftmpx = 0;Fc0x = 0;Fc1x = 0;Fph0x = 0;Fph1x = 0;
                        Ftmpy = 0;Fc0y = 0;Fc1y = 0;Fph0y = 0;Fph1y = 0;
                        Ftmpz = 0;Fc0z = 0;Fc1z = 0; Fph0z = 0;Fph1z = 0;
                        if(x < Nx-1){
                            Ftmpx = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                            cp[c_index(x+1,y,z,comp,ion,Nx,Ny,Nz)])/2/dx*dt/dx;
                            Fc0x = Ftmpx/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                            Fph0x = z_charge[ion]*Ftmpx;
                            // Right c with left c (-Fc0x)

                            ierr = MatSetValue(Jac,Ind_1(x+1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc0x,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Right c with left phi (-Fph0x)
                            ierr = MatSetValue(Jac,Ind_1(x+1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                               -Fph0x,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;

                            //Right phi with left c in voltage eqn
                            ierr = MatSetValue(Jac,Ind_1(x+1,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -z_charge[ion]*Fc0x,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                        if(x > 0){
                            Ftmpx = Dcs[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)*3]*(cp[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                              cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2/dx*dt/dx;
                            Fc1x = Ftmpx/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                            Fph1x = z_charge[ion]*Ftmpx;
                            //left c with right c (-Fc1x)
                            ierr = MatSetValue(Jac,Ind_1(x-1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc1x,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Left c with right phi (-Fph1x)
                            ierr = MatSetValue(Jac,Ind_1(x-1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                               -Fph1x,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;

                            //Left phi with right c in voltage eqn
                            ierr = MatSetValue(Jac,Ind_1(x-1,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -z_charge[ion]*Fc1x,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                        if(y < Ny-1){
                            Ftmpy = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+1]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                              cp[c_index(x,y+1,z,comp,ion,Nx,Ny,Nz)])/2/dy*dt/dy;
                            Fc0y = Ftmpy/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                            Fph0y = z_charge[ion]*Ftmpy;
                            // Upper c with lower c (-Fc0y)
                            ierr = MatSetValue(Jac,Ind_1(x,y+1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc0y,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Upper c with lower phi (-Fph0y)
                            ierr = MatSetValue(Jac,Ind_1(x,y+1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                               -Fph0y,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;

                            //Upper phi with lower c in voltage eqn
                            ierr = MatSetValue(Jac,Ind_1(x,y+1,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -z_charge[ion]*Fc0y,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                        if(y > 0){
                            Ftmpy = Dcs[c_index(x,y-1,z,comp,ion,Nx,Ny,Nz)*3+1]*(cp[c_index(x,y-1,z,comp,ion,Nx,Ny,Nz)]+
                                                                                cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2/dy*dt/dy;
                            Fc1y = Ftmpy/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                            Fph1y = z_charge[ion]*Ftmpy;
                            //Lower c with Upper c (-Fc1y)
                            ierr = MatSetValue(Jac,Ind_1(x,y-1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc1y,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Lower c with Upper phi (-Fph1y)
                            ierr = MatSetValue(Jac,Ind_1(x,y-1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                               -Fph1y,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;

                            //Lower phi with upper c in voltage eqn
                            ierr = MatSetValue(Jac,Ind_1(x,y-1,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -z_charge[ion]*Fc1y,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                        if(z < Nz-1){
                            Ftmpz = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                              cp[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])/2/dz*dt/dz;
                            Fc0z = Ftmpz/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                            Fph0z = z_charge[ion]*Ftmpz;
                            // Top c with Bot c (-Fc0x)

                            ierr = MatSetValue(Jac,Ind_1(x,y,z+1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc0z,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Top c with Bot phi (-Fph0x)
                            ierr = MatSetValue(Jac,Ind_1(x,y,z+1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                               -Fph0z,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;

                            //Bot phi with Top c in voltage eqn
                            ierr = MatSetValue(Jac,Ind_1(x,y,z+1,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -z_charge[ion]*Fc0z,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                        if(z > 0){
                            Ftmpz = Dcs[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)]+
                                                                                cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2/dz*dt/dz;
                            Fc1z = Ftmpz/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                            Fph1z = z_charge[ion]*Ftmpz;
                            //left c with right c (-Fc1x)
                            ierr = MatSetValue(Jac,Ind_1(x,y,z-1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc1z,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Left c with right phi (-Fph1x)
                            ierr = MatSetValue(Jac,Ind_1(x,y,z-1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                               -Fph1z,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;

                            //Left phi with right c in voltage eqn
                            ierr = MatSetValue(Jac,Ind_1(x,y,z-1,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -z_charge[ion]*Fc1z,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                        //Diagonal term contribution
                        Ac = al[al_index(x,y,z,comp,Nx,Ny,Nz)]+Fc0x+Fc1x+Fc0y+Fc1y+Fc0z+Fc1z;
                        Aphi = Fph0x+Fph1x+Fph0y+Fph1y+Fph0z+Fph1z;

                        //Add up terms for voltage eqns
                        Fphph0x[comp] += z_charge[ion]*Fph0x;
                        Fphph1x[comp] += z_charge[ion]*Fph1x;
                        Fphph0y[comp] += z_charge[ion]*Fph0y;
                        Fphph1y[comp] += z_charge[ion]*Fph1y;
                        Fphph0z[comp] += z_charge[ion]*Fph0z;
                        Fphph1z[comp] += z_charge[ion]*Fph1z;

                        //membrane current contributions
                        Ac += flux->dfdci[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;
                        Aphi += flux->dfdphim[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;
                        // Different Compartment Terms
                        // C Extracellular with C Inside
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -flux->dfdci[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // C Intra with C Extra
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),
                                           flux->dfdce[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // C Extracellular with Phi Inside
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -flux->dfdphim[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // C Intra with Phi Extra
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,Nc-1,Nx,Ny,Nz),
                                           -flux->dfdphim[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Same compartment terms
                        // c with c
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Ac,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // c with phi
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Aphi,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                        //Intra-Phi with c (voltage eqn)
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           z_charge[ion]*(Fc0x+Fc1x+Fc0y+Fc1y+Fc0z+Fc1z+
                                                          flux->dfdci[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt),INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //IntraPhi with c extra(volt eqn)
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),
                                           z_charge[ion]*(flux->dfdce[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt),INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Extra-Phi with intra-c (voltage eqn)
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -z_charge[ion]*(flux->dfdci[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt),INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                    }
                    //Extracellular terms
                    comp = Nc-1;
                    //Electrodiffusion contributions
                    Ftmpx = 0;Fc0x = 0;Fc1x = 0;Fph0x = 0;Fph1x = 0;
                    Ftmpy = 0;Fc0y = 0;Fc1y = 0;Fph0y = 0;Fph1y = 0;
                    Ftmpz = 0;Fc0z = 0;Fc1z = 0; Fph0z = 0;Fph1z = 0;
                    if(x < Nx-1){
                        Ftmpx = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                        cp[c_index(x+1,y,z,comp,ion,Nx,Ny,Nz)])/2/dx*dt/dx;
                        Fc0x = Ftmpx/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                        Fph0x = z_charge[ion]*Ftmpx;
                        // Right c with left c (-Fc0x)
                        ierr = MatSetValue(Jac,Ind_1(x+1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc0x,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Right c with left phi (-Fph0x)
                        ierr = MatSetValue(Jac,Ind_1(x+1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fph0x,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                        // Right Phi with left c (voltage eqn)
                        ierr = MatSetValue(Jac,Ind_1(x+1,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -z_charge[ion]*Fc0x,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(x > 0){
                        Ftmpx = Dcs[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)*3]*(cp[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                          cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2/dx*dt/dx;
                        Fc1x = Ftmpx/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                        Fph1x = z_charge[ion]*Ftmpx;
                        //left c with right c (-Fc1x)
                        ierr = MatSetValue(Jac,Ind_1(x-1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc1x,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Left c with right phi (-Fph1x)
                        ierr = MatSetValue(Jac,Ind_1(x-1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fph1x,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                        // left Phi with right c (voltage eqn)
                        ierr = MatSetValue(Jac,Ind_1(x-1,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -z_charge[ion]*Fc1x,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(y < Ny-1){
                        Ftmpy = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+1]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                          cp[c_index(x,y+1,z,comp,ion,Nx,Ny,Nz)])/2/dy*dt/dy;
                        Fc0y = Ftmpy/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                        Fph0y = z_charge[ion]*Ftmpy;
                        // Upper c with lower c (-Fc0y)
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc0y,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Upper c with lower phi (-Fph0y)
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fph0y,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                        // Upper Phi with lower c (voltage eqn)
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -z_charge[ion]*Fc0y,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(y > 0){
                        Ftmpy = Dcs[c_index(x,y-1,z,comp,ion,Nx,Ny,Nz)*3+1]*(cp[c_index(x,y-1,z,comp,ion,Nx,Ny,Nz)]+
                                                                            cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2/dy*dt/dy;
                        Fc1y = Ftmpy/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                        Fph1y = z_charge[ion]*Ftmpy;
                        //Lower c with Upper c (-Fc1y)
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc1y,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Lower c with Upper phi (-Fph1y)
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fph1y,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                        // Lower Phi with upper c (voltage eqn)
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -z_charge[ion]*Fc1y,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(z < Nz-1){
                        Ftmpz = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                          cp[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])/2/dz*dt/dz;
                        Fc0z = Ftmpz/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                        Fph0z = z_charge[ion]*Ftmpz;
                        // Top c with Bot c (-Fc0x)

                        ierr = MatSetValue(Jac,Ind_1(x,y,z+1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc0z,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Top c with Bot phi (-Fph0x)
                        ierr = MatSetValue(Jac,Ind_1(x,y,z+1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fph0z,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                        //Bot phi with Top c in voltage eqn
                        ierr = MatSetValue(Jac,Ind_1(x,y,z+1,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -z_charge[ion]*Fc0z,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(z > 0){
                        Ftmpz = Dcs[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)]+
                                                                            cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2/dz*dt/dz;
                        Fc1z = Ftmpz/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                        Fph1z = z_charge[ion]*Ftmpz;
                        //left c with right c (-Fc1x)
                        ierr = MatSetValue(Jac,Ind_1(x,y,z-1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc1z,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Left c with right phi (-Fph1x)
                        ierr = MatSetValue(Jac,Ind_1(x,y,z-1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fph1z,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                        //Left phi with right c in voltage eqn
                        ierr = MatSetValue(Jac,Ind_1(x,y,z-1,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -z_charge[ion]*Fc1z,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }

                    //Diagonal term contribution
                    Ac = (1-al[al_index(x,y,z,0,Nx,Ny,Nz)]-al[al_index(x,y,z,1,Nx,Ny,Nz)])+Fc0x+Fc1x+Fc0y+Fc1y+Fc0z+Fc1z;
                    Aphi = Fph0x+Fph1x+Fph0y+Fph1y+Fph0z+Fph1z;

                    Avolt = z_charge[ion]*(Fc0x+Fc1x+Fc0y+Fc1y+Fc0z+Fc1z);

                    //Add up terms for voltage eqns
                    Fphph0x[comp] += z_charge[ion]*Fph0x;
                    Fphph1x[comp] += z_charge[ion]*Fph1x;
                    Fphph0y[comp] += z_charge[ion]*Fph0y;
                    Fphph1y[comp] += z_charge[ion]*Fph1y;
                    Fphph0z[comp] += z_charge[ion]*Fph0z;
                    Fphph1z[comp] += z_charge[ion]*Fph1z;

                    //Membrane current contribution
                    for(comp = 0; comp < Nc-1; comp++){
                        Ac -= flux->dfdce[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;
                        Aphi += flux->dfdphim[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;
                        Avolt -= z_charge[ion]*flux->dfdce[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;
                    }
                    //Add bath contributions
                    Ftmpx = sqrt(pow(Dcb[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)*3],2)+
                                 pow(Dcb[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)*3+1],2)+
                                 pow(Dcb[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)*3+2],2));
                    Ac -= Ftmpx*(cp[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)]+cbath[ion])/(2*c[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)])*dt;
                    Aphi -= Ftmpx*(cp[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)]+cbath[ion])*z_charge[ion]/2*dt;

                    Avolt -= z_charge[ion]*Ftmpx*(cp[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)]+
                                                  cbath[ion])/(2*c[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)])*dt;

                    //Insert extracell to extracell parts
                    // c with c
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ac,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    // c with phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,Ni,Nc-1,Nx,Ny,Nz),Aphi,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;

                    //phi with c (voltage eqn)
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,ion,
                                                                              Nc-1,Nx,Ny,Nz),Avolt,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                //Derivative of charge-capacitance
                for(comp = 0; comp < Nc-1; comp++){
                    if(x < Nx-1){
                        //Right phi with left phi (-Fph0x)
                        ierr = MatSetValue(Jac,Ind_1(x+1,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fphph0x[comp],INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(x > 0){
                        //Left phi with right phi (-Fph1x)
                        ierr = MatSetValue(Jac,Ind_1(x-1,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fphph1x[comp],INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(y < Ny-1){
                        //Upper phi with lower phi (-Fph0y)
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fphph0y[comp],INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(y > 0){
                        //Lower phi with upper phi (-Fph1y)
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fphph1y[comp],INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(z < Nz-1){
                        //Top phi with Bot phi (-Fph0x)
                        ierr = MatSetValue(Jac,Ind_1(x,y,z+1,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fphph0z[comp],INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(z > 0){
                        //Bot phi with Top phi (-Fph1x)
                        ierr = MatSetValue(Jac,Ind_1(x,y,z-1,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fphph1z[comp],INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    Avolt = cm[al_index(x,y,z,comp,Nx,Ny,Nz)]+Fphph0x[comp]+Fphph1x[comp]+Fphph0y[comp]+Fphph1y[comp]+Fphph0z[comp]+Fphph1z[comp];
                    AvoltN = -cm[al_index(x,y,z,comp,Nx,Ny,Nz)];
                    for(ion = 0; ion < Ni; ion++){
                        Avolt += z_charge[ion]*flux->dfdphim[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;
                        AvoltN -= z_charge[ion]*flux->dfdphim[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;
                    }

                    //Intra-phi with Intra-phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Avolt,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Intra-phi with extra-phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,
                                                                              Nc-1,Nx,Ny,Nz),AvoltN,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                //Extracellular terms
                comp = Nc-1;
                if(x < Nx-1){
                    //Right phi with left phi (-Fph0x)
                    ierr = MatSetValue(Jac,Ind_1(x+1,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                       -Fphph0x[comp],INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if(x > 0){
                    //Left phi with right phi (-Fph1x)
                    ierr = MatSetValue(Jac,Ind_1(x-1,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                       -Fphph1x[comp],INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if(y < Ny-1){
                    //Upper phi with lower phi (-Fph0y)
                    ierr = MatSetValue(Jac,Ind_1(x,y+1,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                       -Fphph0y[comp],INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if(y > 0){
                    //Lower phi with upper phi (-Fph1y)
                    ierr = MatSetValue(Jac,Ind_1(x,y-1,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                       -Fphph1y[comp],INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if(z < Nz-1){
                    //Top phi with Bot phi (-Fph0x)
                    ierr = MatSetValue(Jac,Ind_1(x,y,z+1,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                       -Fphph0z[comp],INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                if(z > 0){
                    //Bot phi with Top phi (-Fph1x)
                    ierr = MatSetValue(Jac,Ind_1(x,y,z-1,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                       -Fphph1z[comp],INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                AvoltN = 0;

                for(int k = 0; k < Nc-1; k++){
                    AvoltN += cm[al_index(x,y,z,k,Nx,Ny,Nz)];
                    Avolt = -cm[al_index(x,y,z,k,Nx,Ny,Nz)];
                    for(ion = 0; ion < Ni; ion++){
                        Avolt -= z_charge[ion]*flux->dfdphim[c_index(x,y,z,k,ion,Nx,Ny,Nz)]*dt;
                        AvoltN += z_charge[ion]*flux->dfdphim[c_index(x,y,z,k,ion,Nx,Ny,Nz)]*dt;
                    }
                    //Extra-phi with Intra-phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,k,Nx,Ny,Nz),Avolt,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }

                AvoltN += Fphph0x[comp]+Fphph1x[comp]+Fphph0y[comp]+Fphph1y[comp]+Fphph0z[comp]+Fphph1z[comp];

                //Bath terms
                for(ion = 0; ion < Ni; ion++){
                    Ftmpx = sqrt(pow(Dcb[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)*3],2)+
                                 pow(Dcb[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)*3+1],2)+
                                 pow(Dcb[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)*3+2],2));
                    AvoltN -= z_charge[ion]*Ftmpx*(cp[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)]+cbath[ion])*z_charge[ion]/2*dt;
                }
                //extra-phi with extra-phi
                ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),AvoltN,INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;

            }
        }
    }

    ierr = MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    if (A != Jac) {
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr); }

    ierr = restore_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    if(Profiling_on) {
        PetscLogEventEnd(event[0], 0, 0, 0, 0);
    }
    return ierr;
}

PetscErrorCode calc_residual_algebraic(SNES snes,Vec current_state,Vec Res,void *ctx)
{
    //Residual equation using algebraic version of the charge-capacitance relation
    //Alpha is solved for here
    struct AppCtx * user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if(Profiling_on) {
        PetscLogEventBegin(event[1], 0, 0, 0, 0);
    }
    ierr = extract_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    //Compute membrane ionic flux relation quantitites
    ionmflux(user);


    //Compute membrane water flow related quantities
    wflowm(user);

    PetscReal *c = user->state_vars->c;
    PetscReal *phi = user->state_vars->phi;
    PetscReal *al = user->state_vars->alpha;
    PetscReal *cp = user->state_vars_past->c;
    PetscReal *alp = user->state_vars_past->alpha;
    PetscReal *phip = user->state_vars_past->phi;

    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct FluxData *flux = user->flux;
    PetscReal *cm = user->con_vars->cm;
    PetscReal dt = user->dt;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    PetscReal dz = user->dz;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
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
                            Rcvx = Rcvx*(log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-log(c[c_index(
                                    x-1,y,z,comp,ion,Nx,Ny,Nz)])+
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
                                    z_charge[ion]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-phi[phi_index(x,y-1,z,comp,Nx,Ny,Nz)]))/dy*dt/dy;
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
                            Rcvz = Dcs[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)]+
                                                                               cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2;
                            Rcvz = Rcvz*(log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-log(c[c_index(x,y,
                                                                                               z-1,comp,ion,Nx,Ny,Nz)])+
                                         z_charge[ion]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-
                                                        phi[phi_index(x,y,z-1,comp,Nx,Ny,Nz)]))/dz*dt/dz;
                        }
                        //Add Second right moving difference
                        if(z < Nz-1){
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

                        ierr = VecSetValue(Res,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Resc,INSERT_VALUES);
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
                                z_charge[ion]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-
                                               phi[phi_index(x-1,y,z,comp,Nx,Ny,Nz)]))/dx*dt/dx;
                    }
                    //Add Second right moving difference
                    if(x < Nx-1){
                        RcvxRight = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                            cp[c_index(x+1,y,z,comp,ion,Nx,Ny,Nz)])/2;
                        RcvxRight = RcvxRight*(log(c[c_index(x+1,y,z,comp,ion,Nx,Ny,Nz)])-
                                               log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])+z_charge[ion]*(phi[phi_index(
                                x+1,y,z,comp,Nx,Ny,Nz)]-
                                                                                                      phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]))/dx*dt/dx;
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
                        Rcvz = Dcs[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)]+
                                                                           cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2;
                        Rcvz = Rcvz*(log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-log(c[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)])+
                                     z_charge[ion]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-
                                                    phi[phi_index(x,y,z-1,comp,Nx,Ny,Nz)]))/dz*dt/dz;
                    }
                    //Add Second right moving difference
                    if(z < Nz-1){
                        RcvzTop = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                            cp[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])/2;
                        RcvzTop = RcvzTop*(log(c[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])-
                                           log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])+z_charge[ion]*
                                                                                   (phi[phi_index(x,y,
                                                                                                  z+1,comp,Nx,Ny,Nz)]-
                                                                                    phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]))/dz*dt/dz;
                    }
                    Resc = alNc*c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]-alpNc*cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                    Resc += Rcvx-RcvxRight+Rcvy-RcvyUp+Rcvz-RcvzTop+flux->mflux[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;
                    //Add bath variables

                    Resc -= sqrt(pow(Dcb[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3],2)+
                                 pow(Dcb[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+1],2)+
                                 pow(Dcb[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2],2))*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                                   cbath[ion])/2.0*(log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-log(cbath[ion])+
                                            z_charge[ion]*phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-z_charge[ion]*phibath)*dt;
                    ierr = VecSetValue(Res,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Resc,INSERT_VALUES);
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
                    ierr = VecSetValue(Res,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Resc,INSERT_VALUES);
                    CHKERRQ(ierr);
                }
                //Extracellular term
                comp = Nc-1;
                Resc = (1-al[al_index(x,y,z,0,Nx,Ny,Nz)]-al[al_index(x,y,z,1,Nx,Ny,Nz)])*
                       cz(c,z_charge,x,y,z,Nx,Ny,Nz,comp,user)+
                       user->con_vars->zo[phi_index(x,y,z,comp,Nx,Ny,Nz)]*user->con_vars->ao[phi_index(x,y,z,comp,Nx,Ny,Nz)];
                ierr = VecSetValue(Res,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Resc,INSERT_VALUES);
                CHKERRQ(ierr);

                //Residual for water flow
                //Plus modification to electroneutrality for non-zero mem.compacitance
                for(comp = 0; comp < Nc-1; comp++){
                    //Water flow
                    ierr = VecSetValue(Res,Ind_1(x,y,z,Ni+1,comp,Nx,Ny,Nz),al[al_index(x,y,z,comp,Nx,Ny,Nz)]-
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
                    ierr = VecSetValue(Res,Ind_1(x,y,z,Ni,Nc-1,Nx,Ny,Nz),-cm[comp]*(phi[phi_index(x,y,z,Nc-1,Nx,Ny,Nz)]-
                                                                                   phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]),ADD_VALUES);
                    CHKERRQ(ierr);
                    //Intracell voltage mod
                    ierr = VecSetValue(Res,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),-cm[comp]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-
                                                                                   phi[phi_index(x,y,z,Nc-1,Nx,Ny,Nz)]),ADD_VALUES);
                    CHKERRQ(ierr);
                }
            }
        }
    }

    ierr = VecAssemblyBegin(Res);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Res);CHKERRQ(ierr);


    ierr = restore_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    if(Profiling_on) {
        PetscLogEventEnd(event[1], 0, 0, 0, 0);
    }
    return ierr;
}

PetscErrorCode
calc_jacobian_algebraic(SNES snes,Vec current_state, Mat A, Mat Jac,void *ctx){
    //Jacobian equation using algebraic version of the charge-capacitance relation
    // Alpha is solved for here
    struct AppCtx *user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if(Profiling_on){
        PetscLogEventBegin(event[0],0,0,0,0);
    }
    ierr = extract_subarray_Read(current_state,user->state_vars);
    CHKERRQ(ierr);
    PetscReal *c = user->state_vars->c;
    PetscReal *al = user->state_vars->alpha;
    PetscReal *cp = user->state_vars_past->c;

    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct FluxData *flux = user->flux;
    PetscReal dt = user->dt;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    PetscReal dz = user->dz;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt Nz = user->Nz;
    struct ConstVars *con_vars = user->con_vars;
    PetscReal *cm = con_vars->cm;

    PetscInt ind = 0;
    PetscInt x,y,z,ion,comp;

    PetscReal Ftmpx,Fc0x,Fc1x,Fph0x,Fph1x;
    PetscReal Ftmpy,Fc0y,Fc1y,Fph0y,Fph1y;
    PetscReal Ftmpz,Fc0z,Fc1z,Fph0z,Fph1z;
    PetscReal Ac,Aphi;


    //Ionic concentration equations
    for(z = 0; z < Nz; z++){
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

                            ierr = MatSetValue(Jac,Ind_1(x+1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc0x,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Right c with left phi (-Fph0x)
                            ierr = MatSetValue(Jac,Ind_1(x+1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
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
                            ierr = MatSetValue(Jac,Ind_1(x-1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc1x,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Left c with right phi (-Fph1x)
                            ierr = MatSetValue(Jac,Ind_1(x-1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
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
                            ierr = MatSetValue(Jac,Ind_1(x,y+1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc0y,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Upper c with lower phi (-Fph0y)
                            ierr = MatSetValue(Jac,Ind_1(x,y+1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
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
                            ierr = MatSetValue(Jac,Ind_1(x,y-1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc1y,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Lower c with Upper phi (-Fph1y)
                            ierr = MatSetValue(Jac,Ind_1(x,y-1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                               -Fph1y,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                        if(z < Nz-1){
                            Ftmpz = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                              cp[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])/2/dz*dt/dz;
                            Fc0z = Ftmpz/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                            Fph0z = z_charge[ion]*Ftmpz;
                            // Top c with Bot c (-Fc0x)

                            ierr = MatSetValue(Jac,Ind_1(x,y,z+1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc0z,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Top c with Bot phi (-Fph0x)
                            ierr = MatSetValue(Jac,Ind_1(x,y,z+1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                               -Fph0z,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;

                        }
                        if(z > 0){
                            Ftmpz = Dcs[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)]+
                                                                                cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2/dz*dt/dz;
                            Fc1z = Ftmpz/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                            Fph1z = z_charge[ion]*Ftmpz;
                            //Top c with Bot c (-Fc1x)
                            ierr = MatSetValue(Jac,Ind_1(x,y,z-1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc1z,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Top c with Bot phi (-Fph1x)
                            ierr = MatSetValue(Jac,Ind_1(x,y,z-1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
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
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -flux->dfdci[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // C Intra with C Extra
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),
                                           flux->dfdce[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // C Extracellular with Phi Inside
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -flux->dfdphim[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // C Intra with Phi Extra
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,Nc-1,Nx,Ny,Nz),
                                           -flux->dfdphim[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Volume terms
                        //C extra with intra alpha
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,Ni+1,comp,Nx,Ny,Nz),
                                           -c[c_index(x,y,z,Nc-1,ion,Nx,Ny,Nz)],INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //C intra with intra alpha
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni+1,comp,Nx,Ny,Nz),
                                           c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)],INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Same compartment terms
                        // c with c
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Ac,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // c with phi
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Aphi,INSERT_VALUES);
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
                        ierr = MatSetValue(Jac,Ind_1(x+1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc0x,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Right c with left phi (-Fph0x)
                        ierr = MatSetValue(Jac,Ind_1(x+1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fph0x,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(x > 0){
                        Ftmpx = Dcs[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)*3]*(cp[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)]
                                                                          +cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2/dx*dt/dx;
                        Fc1x = Ftmpx/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                        Fph1x = z_charge[ion]*Ftmpx;
                        //left c with right c (-Fc1x)
                        ierr = MatSetValue(Jac,Ind_1(x-1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc1x,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Left c with right phi (-Fph1x)
                        ierr = MatSetValue(Jac,Ind_1(x-1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
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
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc0y,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Upper c with lower phi (-Fph0y)
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
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
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc1y,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Lower c with Upper phi (-Fph1y)
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fph1y,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(z < Nz-1){
                        Ftmpz = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                          cp[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])/2/dz*dt/dz;
                        Fc0z = Ftmpz/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                        Fph0z = z_charge[ion]*Ftmpz;
                        // Top c with Bot c (-Fc0x)

                        ierr = MatSetValue(Jac,Ind_1(x,y,z+1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc0z,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Top c with Bot phi (-Fph0x)
                        ierr = MatSetValue(Jac,Ind_1(x,y,z+1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fph0z,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                    }
                    if(z > 0){
                        Ftmpz = Dcs[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)]+
                                                                            cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2/dz*dt/dz;
                        Fc1z = Ftmpz/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                        Fph1z = z_charge[ion]*Ftmpz;
                        //Top c with Bot c (-Fc1x)
                        ierr = MatSetValue(Jac,Ind_1(x,y,z-1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc1z,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Top c with Bot phi (-Fph1x)
                        ierr = MatSetValue(Jac,Ind_1(x,y,z-1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
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
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ac,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    // c with phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,Ni,Nc-1,Nx,Ny,Nz),Aphi,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }

            }
        }
    }

    //Electroneutrality charge-capcitance condition
    for(z = 0; z < Nz; z++){
        for(y = 0; y < Ny; y++){
            for(x = 0; x < Nx; x++){
                //electroneutral-concentration entries
                for(ion = 0; ion < Ni; ion++){
                    for(comp = 0; comp < Nc-1; comp++){
                        //Phi with C entries
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           z_charge[ion]*al[al_index(x,y,z,comp,Nx,Ny,Nz)],INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    //Phi with C extracellular one
                    comp = Nc-1;
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                       z_charge[ion]*(1-al[al_index(x,y,z,0,Nx,Ny,Nz)]-al[al_index(x,y,z,1,Nx,Ny,Nz)]),
                                       INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;

                }
                //electroneutrality-voltage entries
                Aphi = 0;
                for(comp = 0; comp < Nc-1; comp++){
                    Aphi -= cm[comp];
                }
                //extraphi with extra phi
                ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,Ni,Nc-1,Nx,Ny,Nz),Aphi,INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
                for(comp = 0; comp < Nc-1; comp++){
                    //Extra phi with intra phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,
                                                 Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),cm[comp],INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    // Intra phi with Extraphi
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,
                                                                              Nc-1,Nx,Ny,Nz),cm[comp],INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Intra phi with Intra phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),-cm[comp],INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Extra phi with intra-Volume
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,Ni+1,comp,Nx,Ny,Nz),
                                       -cz(c,z_charge,x,y,z,Nx,Ny,Nz,Nc-1,user),INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Intra phi with Intra Vol
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni+1,comp,Nx,Ny,Nz),
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
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni+1,comp,Nx,Ny,Nz),Ind_1(x,y,z,
                                                                                Ni+1,comp,Nx,Ny,Nz),Ac,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Off diagonal (from aNc=1-sum(ak))
                    for(PetscInt l = 0; l < comp; l++){
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni+1,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni+1,l,Nx,Ny,Nz),
                                           flux->dwdpi[al_index(x,y,z,comp,Nx,Ny,Nz)]*
                                           con_vars->ao[phi_index(x,y,z,Nc-1,Nx,Ny,Nz)]/
                                           pow(1-al[al_index(x,y,z,0,Nx,Ny,Nz)]-al[al_index(x,y,z,1,Nx,Ny,Nz)],2)*dt,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    for(PetscInt l = comp+1; l < Nc-1; l++){
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni+1,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni+1,l,Nx,Ny,Nz),
                                           flux->dwdpi[al_index(x,y,z,comp,Nx,Ny,Nz)]*
                                           con_vars->ao[phi_index(x,y,z,Nc-1,Nx,Ny,Nz)]/
                                           ((1-al[al_index(x,y,z,0,Nx,Ny,Nz)]-al[al_index(x,y,z,1,Nx,Ny,Nz)])*
                                            (1-al[al_index(x,y,z,0,Nx,Ny,Nz)]-al[al_index(x,y,z,1,Nx,Ny,Nz)]))*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    for(ion = 0; ion < Ni; ion++){
                        //Volume to extra c
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni+1,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),
                                           flux->dwdpi[al_index(x,y,z,comp,Nx,Ny,Nz)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Volume to intra c
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni+1,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
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

    if (A != Jac) {
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr); }

    ierr = restore_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    if(Profiling_on) {
        PetscLogEventEnd(event[0], 0, 0, 0, 0);
    }
    return ierr;
}

PetscErrorCode calc_residual_algebraic_no_vol(SNES snes,Vec current_state,Vec Res,void *ctx)
{
    //Residual equation using algebraic version of the charge-capacitance relation
    //Alpha is not solved for here
    struct AppCtx * user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if(Profiling_on) {
        PetscLogEventBegin(event[1], 0, 0, 0, 0);
    }
    ierr = extract_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    //Compute membrane ionic flux relation quantitites
    ionmflux(user);


    //Compute membrane water flow related quantities
    wflowm(user);

    PetscReal *c = user->state_vars->c;
    PetscReal *phi = user->state_vars->phi;
    PetscReal *al = user->state_vars->alpha;
    PetscReal *cp = user->state_vars_past->c;
    PetscReal *alp = user->state_vars_past->alpha;
    PetscReal *phip = user->state_vars_past->phi;

    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct FluxData *flux = user->flux;
    PetscReal *cm = user->con_vars->cm;
    PetscReal dt = user->dt;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    PetscReal dz = user->dz;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
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
                            Rcvx = Rcvx*(log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-log(c[c_index(
                                    x-1,y,z,comp,ion,Nx,Ny,Nz)])+
                                    z_charge[ion]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-
                                                   phi[phi_index(x-1,y,z,comp,Nx,Ny,Nz)]))/dx*dt/dx;
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
                            RcvyUp = RcvyUp*(log(c[c_index(x,y+1,z,comp,ion,Nx,Ny,Nz)])-log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])+
                                    z_charge[ion]*(phi[phi_index(x,y+1,z,comp,Nx,Ny,Nz)]-
                                                   phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]))/dy*dt/dy;
                        }
                        Rcvz = 0;
                        RcvzTop = 0;
                        //Top Bot difference
                        if(z > 0){
                            Rcvz = Dcs[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)]+
                                                                               cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2;
                            Rcvz = Rcvz*(log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-log(c[c_index(x,y,
                                                                                               z-1,comp,ion,Nx,Ny,Nz)])+
                                         z_charge[ion]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-
                                                        phi[phi_index(x,y,z-1,comp,Nx,Ny,Nz)]))/dz*dt/dz;
                        }
                        //Next Topward difference
                        if(z < Nz-1){
                            RcvzTop = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                                cp[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])/2;
                            RcvzTop = RcvzTop*(log(c[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])-log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])+
                                             z_charge[ion]*(phi[phi_index(x,y,z+1,comp,Nx,Ny,Nz)]-
                                                            phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]))/dz*dt/dz;
                        }
                        Resc = al[al_index(x,y,z,comp,Nx,Ny,Nz)]*c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]-
                               alp[al_index(x,y,z,comp,Nx,Ny,Nz)]*cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                        Resc += Rcvx-RcvxRight+Rcvy-RcvyUp+Rcvz-RcvzTop+flux->mflux[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;

                        ierr = VecSetValue(Res,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Resc,INSERT_VALUES);
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
                                z_charge[ion]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-
                                               phi[phi_index(x-1,y,z,comp,Nx,Ny,Nz)]))/dx*dt/dx;
                    }
                    //Add Second right moving difference
                    if(x < Nx-1){
                        RcvxRight = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                            cp[c_index(x+1,y,z,comp,ion,Nx,Ny,Nz)])/2;
                        RcvxRight = RcvxRight*(log(c[c_index(x+1,y,z,comp,ion,Nx,Ny,Nz)])-
                                               log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])+z_charge[ion]*(phi[phi_index(
                                x+1,y,z,comp,Nx,Ny,Nz)]-
                                                                                                      phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]))/dx*dt/dx;
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
                        RcvyUp = RcvyUp*(log(c[c_index(x,y+1,z,comp,ion,Nx,Ny,Nz)])-log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])+
                                z_charge[ion]*(phi[phi_index(x,y+1,z,comp,Nx,Ny,Nz)]-
                                               phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]))/dy*dt/dy;
                    }
                    Rcvz = 0;
                    RcvzTop = 0;
                    //Top Bot difference
                    if(z > 0){
                        Rcvz = Dcs[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)]+
                                                                           cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2;
                        Rcvz = Rcvz*(log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-log(c[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)])+
                                     z_charge[ion]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-
                                                    phi[phi_index(x,y,z-1,comp,Nx,Ny,Nz)]))/dz*dt/dz;
                    }
                    //Next Topward difference
                    if(z < Nz-1){
                        RcvzTop = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                            cp[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])/2;
                        RcvzTop = RcvzTop*(log(c[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])-log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])+
                                           z_charge[ion]*(phi[phi_index(x,y,z+1,comp,Nx,Ny,Nz)]-
                                                          phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]))/dz*dt/dz;
                    }
                    Resc = alNc*c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]-alpNc*cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                    Resc += Rcvx-RcvxRight+Rcvy-RcvyUp+Rcvz-RcvzTop+flux->mflux[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt;
                    //Add bath variables

                    Resc -= sqrt(pow(Dcb[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3],2)+
                                 pow(Dcb[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+1],2)+
                                 pow(Dcb[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2],2))*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                                   cbath[ion])/2.0*(log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-log(cbath[ion])+
                                            z_charge[ion]*phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-z_charge[ion]*phibath)*dt;
                    ierr = VecSetValue(Res,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Resc,INSERT_VALUES);
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
                    ierr = VecSetValue(Res,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Resc,INSERT_VALUES);
                    CHKERRQ(ierr);
                }
                //Extracellular term
                comp = Nc-1;
                Resc = (1-al[al_index(x,y,z,0,Nx,Ny,Nz)]-al[al_index(x,y,z,1,Nx,Ny,Nz)])*
                       cz(c,z_charge,x,y,z,Nx,Ny,Nz,comp,user)+
                       user->con_vars->zo[phi_index(x,y,z,comp,Nx,Ny,Nz)]*user->con_vars->ao[phi_index(x,y,z,comp,Nx,Ny,Nz)];
                ierr = VecSetValue(Res,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Resc,INSERT_VALUES);
                CHKERRQ(ierr);
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
                    ierr = VecSetValue(Res,Ind_1(x,y,z,Ni,Nc-1,Nx,Ny,Nz),-cm[comp]*(phi[phi_index(x,y,z,Nc-1,Nx,Ny,Nz)]-
                                                                                   phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]),ADD_VALUES);
                    CHKERRQ(ierr);
                    //Intracell voltage mod
                    ierr = VecSetValue(Res,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),-cm[comp]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-
                                                                                   phi[phi_index(x,y,z,Nc-1,Nx,Ny,Nz)]),ADD_VALUES);
                    CHKERRQ(ierr);
                }
            }
        }
    }

    ierr = VecAssemblyBegin(Res);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Res);CHKERRQ(ierr);

    ierr = restore_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    if(Profiling_on) {
        PetscLogEventEnd(event[1], 0, 0, 0, 0);
    }
    return ierr;
}

PetscErrorCode
calc_jacobian_algebraic_no_vol(SNES snes,Vec current_state, Mat A, Mat Jac,void *ctx)
{
    //Jacobian equation using algebraic version of the charge-capacitance relation
    // Alpha is not solved for here
    struct AppCtx * user = (struct AppCtx *) ctx;
    PetscErrorCode ierr;
    if(Profiling_on) {
        PetscLogEventBegin(event[0], 0, 0, 0, 0);
    }
    ierr = extract_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    PetscReal *c = user->state_vars->c;
    PetscReal *al = user->state_vars->alpha;
    PetscReal *cp = user->state_vars_past->c;

    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct FluxData *flux = user->flux;
    PetscReal dt = user->dt;
    PetscReal dx = user->dx;
    PetscReal dy = user->dy;
    PetscReal dz = user->dz;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt Nz = user->Nz;
    struct ConstVars *con_vars = user->con_vars;
    PetscReal *cm = con_vars->cm;

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

                            ierr = MatSetValue(Jac,Ind_1(x+1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc0x,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Right c with left phi (-Fph0x)
                            ierr = MatSetValue(Jac,Ind_1(x+1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
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
                            ierr = MatSetValue(Jac,Ind_1(x-1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc1x,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Left c with right phi (-Fph1x)
                            ierr = MatSetValue(Jac,Ind_1(x-1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
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
                            ierr = MatSetValue(Jac,Ind_1(x,y+1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc0y,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Upper c with lower phi (-Fph0y)
                            ierr = MatSetValue(Jac,Ind_1(x,y+1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
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
                            ierr = MatSetValue(Jac,Ind_1(x,y-1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc1y,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Lower c with Upper phi (-Fph1y)
                            ierr = MatSetValue(Jac,Ind_1(x,y-1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                               -Fph1y,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                        if(z < Nz-1){
                            Ftmpz = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                              cp[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])/2/dz*dt/dz;
                            Fc0z = Ftmpz/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                            Fph0z = z_charge[ion]*Ftmpz;
                            // Right c with left c (-Fc0z)

                            ierr = MatSetValue(Jac,Ind_1(x,y,z+1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc0z,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Right c with left phi (-Fph0z)
                            ierr = MatSetValue(Jac,Ind_1(x,y,z+1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                               -Fph0z,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;

                        }
                        if(z > 0){
                            Ftmpz = Dcs[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)]+
                                                                                cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2/dz*dt/dz;
                            Fc1z = Ftmpz/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                            Fph1z = z_charge[ion]*Ftmpz;
                            //left c with right c (-Fc1x)
                            ierr = MatSetValue(Jac,Ind_1(x,y,z-1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                               -Fc1z,INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Left c with right phi (-Fph1x)
                            ierr = MatSetValue(Jac,Ind_1(x,y,z-1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
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
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -flux->dfdci[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // C Intra with C Extra
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),
                                           flux->dfdce[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // C Extracellular with Phi Inside
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -flux->dfdphim[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // C Intra with Phi Extra
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,Nc-1,Nx,Ny,Nz),
                                           -flux->dfdphim[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Same compartment terms
                        // c with c
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Ac,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // c with phi
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Aphi,INSERT_VALUES);
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
                        ierr = MatSetValue(Jac,Ind_1(x+1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc0x,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Right c with left phi (-Fph0x)
                        ierr = MatSetValue(Jac,Ind_1(x+1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
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
                        ierr = MatSetValue(Jac,Ind_1(x-1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc1x,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Left c with right phi (-Fph1x)
                        ierr = MatSetValue(Jac,Ind_1(x-1,y,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
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
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc0y,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Upper c with lower phi (-Fph0y)
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
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
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc1y,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Lower c with Upper phi (-Fph1y)
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,z,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fph1y,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if(z < Nz-1){
                        Ftmpz = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+
                                                                          cp[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])/2/dz*dt/dz;
                        Fc0z = Ftmpz/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                        Fph0z = z_charge[ion]*Ftmpz;
                        // Right c with left c (-Fc0z)

                        ierr = MatSetValue(Jac,Ind_1(x,y,z+1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc0z,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Right c with left phi (-Fph0z)
                        ierr = MatSetValue(Jac,Ind_1(x,y,z+1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
                                           -Fph0z,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                    }
                    if(z > 0){
                        Ftmpz = Dcs[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)*3+2]*(cp[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)]+
                                                                            cp[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2/dz*dt/dz;
                        Fc1z = Ftmpz/c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)];
                        Fph1z = z_charge[ion]*Ftmpz;
                        //left c with right c (-Fc1x)
                        ierr = MatSetValue(Jac,Ind_1(x,y,z-1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                           -Fc1z,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Left c with right phi (-Fph1x)
                        ierr = MatSetValue(Jac,Ind_1(x,y,z-1,ion,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),
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
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ac,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    // c with phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,ion,Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,Ni,Nc-1,Nx,Ny,Nz),Aphi,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }

            }
        }
    }

    //Electroneutrality charge-capacitence condition
    for(z=0;z<Nz;z++){
        for(y = 0; y < Ny; y++){
            for(x = 0; x < Nx; x++){
                //electroneutral-concentration entries
                for(ion = 0; ion < Ni; ion++){
                    for(comp = 0; comp < Nc-1; comp++){
                        //Phi with C entries
                        ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                                z_charge[ion]*al[al_index(x,y,z,comp,Nx,Ny,Nz)],INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    //Phi with C extracellular one
                    comp = Nc-1;
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,ion,comp,Nx,Ny,Nz),
                            z_charge[ion]*(1-al[al_index(x,y,z,0,Nx,Ny,Nz)]-al[al_index(x,y,z,1,Nx,Ny,Nz)]),INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                //electroneutrality-voltage entries
                Aphi = 0;
                for(comp = 0; comp < Nc-1; comp++){
                    Aphi -= cm[comp];
                }
                //extraphi with extra phi
                ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,Ni,Nc-1,Nx,Ny,Nz),Aphi,INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
                for(comp = 0; comp < Nc-1; comp++){
                    //Extra phi with intra phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,
                                                 Nc-1,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),cm[comp],INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    // Intra phi with Extraphi
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,
                                                                              Nc-1,Nx,Ny,Nz),cm[comp],INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Intra phi with Intra phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),Ind_1(x,y,z,Ni,comp,Nx,Ny,Nz),-cm[comp],INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
            }
            // Neuron-Glia glutamate exchange
//            ierr = MatSetValue(Jac,Ind_1(x,y,3,0,Nx),Ind_1(x,y,3,1,Nx),-glut_Bg*dt,INSERT_VALUES);CHKERRQ(ierr);
//            ierr = MatSetValue(Jac,Ind_1(x,y,3,1,Nx),Ind_1(x,y,3,0,Nx),((1-glut_gamma)*glut_Bn*glut_Re-glut_Bg*glut_Rg)*dt,INSERT_VALUES);CHKERRQ(ierr);
        }
    }

    ierr = MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    if (A != Jac) {
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr); }

    ierr = restore_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    if(Profiling_on) {
        PetscLogEventEnd(event[0], 0, 0, 0, 0);
    }
    return ierr;
}




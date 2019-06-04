#include "constants.h"
#include "functions.h"

PetscErrorCode calc_residual_linear_algebraic(SNES snes,Vec current_state,Vec Res,void *ctx)
{
    //Linear discretization
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
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt Nz = user->Nz;

    //Residual for concentration equations
    PetscReal Rcvx,Rcvy,Resc;
    PetscReal RcvxRight,RcvyUp;

    PetscReal alNc,alpNc;
    PetscInt ion,comp,x,y;


    for(x=0;x<Nx;x++) {
        for(y=0;y<Ny;y++) {
            for(ion=0;ion<Ni;ion++) {
                for(comp=0;comp<Nc-1;comp++) {
                    Rcvx = 0;
                    RcvxRight = 0;
                    if(x>0) {
                        //First difference term
                        // C term
                        Rcvx = Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]*(c[c_index(x,y,0,comp,ion,Nx,
                                                                                  0,0)]-c[c_index(x-1,y,0,
                                                                                                  comp,ion,
                                                                                                  Nx,
                                                                                                  0,0)]);
                        //Phi term
                        Rcvx += Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]*z_charge[ion]*(cp[c_index(x-1,y,0,
                                                                                                  comp,ion,Nx,
                                                                                                  0,0)]+cp[c_index(
                                x,
                                y,
                                0,
                                comp,
                                ion,
                                Nx,
                                0,0)])/2*(phi[phi_index(x,y,0,comp,Nx,0,0)]-phi[phi_index(x-1,y,0,comp,Nx,
                                                                                          0,0)]);
                        Rcvx = Rcvx/dx*dt/dx;
                    }
                    //Add Second right moving difference
                    if(x<Nx-1) {
                        // C term
                        RcvxRight = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]*(c[c_index(x+1,y,0,comp,ion,Nx,
                                                                                     0,0)]-c[c_index(x,y,0,
                                                                                                     comp,ion,
                                                                                                     Nx,0,0)]);
                        //Phi term
                        RcvxRight += Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]*z_charge[ion]*(cp[c_index(x+1,y,0,
                                                                                                     comp,ion,Nx,
                                                                                                     0,0)]+cp[c_index(
                                x,y,0,comp,ion,Nx,0,0)])/2*(phi[phi_index(x+1,y,0,comp,Nx,0,0)]-phi[phi_index(
                                x,y,0,comp,Nx,0,0)]);
                        RcvxRight = RcvxRight/dx*dt/dx;
                    }
                    Rcvy = 0;
                    RcvyUp = 0;
                    //Up down difference
                    if(y>0) {
                        // C term
                        Rcvy = Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]*(c[c_index(x,y,0,comp,ion,Nx,
                                                                                    0,0)]-c[c_index(x,
                                                                                                    y-1,0,
                                                                                                    comp,
                                                                                                    ion,Nx,
                                                                                                    0,0)]);
                        //Phi term
                        Rcvy += Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]*z_charge[ion]*(cp[c_index(x,
                                                                                                    y-1,0,
                                                                                                    comp,ion,
                                                                                                    Nx,0,0)]+cp[c_index(
                                x,
                                y,
                                0,
                                comp,
                                ion,
                                Nx,
                                0,0)])/2*(phi[phi_index(x,y,0,comp,Nx,0,0)]-phi[phi_index(x,y-1,0,comp,Nx,
                                                                                          0,0)]);
                        Rcvy = Rcvy/dy*dt/dy;
                    }
                    //Next upward difference
                    if(y<Ny-1) {
                        // C term
                        RcvyUp = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]*(c[c_index(x,y+1,0,comp,ion,
                                                                                    Nx,0,0)]-c[c_index(x,y,
                                                                                                       0,
                                                                                                       comp,
                                                                                                       ion,
                                                                                                       Nx,0,0)]);
                        //Phi term
                        RcvyUp += Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]*z_charge[ion]*(cp[c_index(x,y+1,0,
                                                                                                    comp,ion,Nx,
                                                                                                    0,0)]+cp[c_index(
                                x,
                                y,
                                0,
                                comp,
                                ion,
                                Nx,
                                0,0)])/2*(phi[phi_index(x,y+1,0,comp,Nx,0,0)]-phi[phi_index(x,y,0,comp,Nx,
                                                                                            0,0)]);
                        RcvyUp = RcvyUp/dy*dt/dy;
                    }
                    Resc = al[al_index(x,y,0,comp,Nx,0,0)]*c[c_index(x,y,0,comp,ion,Nx,0,0)]-alp[al_index(x,
                                                                                                          y,
                                                                                                          0,
                                                                                                          comp,
                                                                                                          Nx,
                                                                                                          0,0)]*cp[c_index(
                            x,y,0,comp,ion,Nx,0,0)];
                    Resc += Rcvx - RcvxRight + Rcvy - RcvyUp +flux->mflux[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;

                    ierr = VecSetValue(Res,Ind_1(x,y,0,ion,comp,Nx,0,0),Resc,INSERT_VALUES);CHKERRQ(ierr);

                }
                //Set Extracellular values
                alNc = 1-al[al_index(x,y,0,0,Nx,0,0)]-al[al_index(x,y,0,1,Nx,0,0)];
                alpNc = 1-alp[al_index(x,y,0,0,Nx,0,0)]-alp[al_index(x,y,0,1,Nx,0,0)];
                comp = Nc-1;
                Rcvx = 0;
                RcvxRight = 0;
                if(x>0) {
                    //First difference term
                    // C term
                    Rcvx = Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]*(c[c_index(x,y,0,comp,ion,Nx,0,0)]-c[c_index(
                            x-1,y,0,comp,
                            ion,Nx,0,0)]);
                    //Phi term
                    Rcvx += Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]*z_charge[ion]*(cp[c_index(x-1,y,0,comp,
                                                                                              ion,Nx,0,0)]+cp[c_index(
                            x,y,0,
                            comp,
                            ion,Nx,
                            0,0)])/2*(phi[phi_index(x,y,0,comp,Nx,0,0)]-phi[phi_index(x-1,y,0,comp,Nx,0,0)]);
                    Rcvx = Rcvx/dx*dt/dx;
                }
                //Add Second right moving difference
                if(x<Nx-1) {
                    // C term
                    RcvxRight = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]*(c[c_index(x+1,y,0,comp,ion,Nx,
                                                                                 0,0)]-c[c_index(x,y,0,comp,
                                                                                                 ion,Nx,0,0)]);
                    //Phi term
                    RcvxRight += Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]*z_charge[ion]*(cp[c_index(x+1,y,0,comp,
                                                                                                 ion,Nx,0,0)]+cp[c_index(
                            x,y,
                            0,
                            comp,
                            ion,
                            Nx,
                            0,0)])/2*(phi[phi_index(x+1,y,0,comp,Nx,0,0)]-phi[phi_index(x,y,0,comp,Nx,0,0)]);
                    RcvxRight = RcvxRight/dx*dt/dx;
                }
                Rcvy = 0;
                RcvyUp = 0;
                //Up down difference
                if(y>0) {
                    // C term
                    Rcvy = Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]*(c[c_index(x,y,0,comp,ion,Nx,
                                                                                0,0)]-c[c_index(x,y-1,0,
                                                                                                comp,
                                                                                                ion,Nx,0,0)]);
                    //Phi term
                    Rcvy += Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]*z_charge[ion]*(cp[c_index(x,y-1,0,
                                                                                                comp,ion,Nx,
                                                                                                0,0)]+cp[c_index(
                            x,y,
                            0,
                            comp,
                            ion,
                            Nx,0,0)])/2*(phi[phi_index(x,y,0,comp,Nx,0,0)]-phi[phi_index(x,y-1,0,comp,Nx,
                                                                                         0,0)]);
                    Rcvy = Rcvy/dy*dt/dy;
                }
                //Next upward difference
                if(y<Ny-1) {
                    // C term
                    RcvyUp = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]*(c[c_index(x,y+1,0,comp,ion,Nx,
                                                                                0,0)]-c[c_index(x,y,0,comp,
                                                                                                ion,Nx,0,0)]);
                    //Phi term
                    RcvyUp += Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]*z_charge[ion]*(cp[c_index(x,
                                                                                                y+1,0,comp,
                                                                                                ion,Nx,0,0)]+cp[c_index(
                            x,y,
                            0,
                            comp,
                            ion,
                            Nx,0,0)])/2*(phi[phi_index(x,y+1,0,comp,Nx,0,0)]-phi[phi_index(x,y,0,comp,Nx,
                                                                                           0,0)]);
                    RcvyUp = RcvyUp/dy*dt/dy;
                }
                Resc = alNc*c[c_index(x,y,0,comp,ion,Nx,0,0)]-alpNc*cp[c_index(x,y,0,comp,ion,Nx,0,0)];
                Resc += Rcvx - RcvxRight + Rcvy - RcvyUp +flux->mflux[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                //Add bath variables

                Resc -= sqrt(pow(Dcb[c_index(x,y,0,comp,ion,Nx,0,0)*2],2)+pow(Dcb[c_index(x,y,0,comp,ion,
                                                                                          Nx,0,0)*2+1],2))*(c[c_index(
                        x,y,0,comp,ion,Nx,0,0)]-cbath[ion]+z_charge[ion]*(cp[c_index(
                        x,y,0,comp,ion,Nx,0,0)]+cbath[ion])/2.0*(phi[phi_index(x,y,0,comp,Nx,0,0)]-phibath))*dt;
                ierr = VecSetValue(Res,Ind_1(x,y,0,ion,comp,Nx,0,0),Resc,INSERT_VALUES);CHKERRQ(ierr);

            }
        }
    }


    for(x=0;x<Nx;x++) {
        for(y=0;y<Ny;y++) {

            //Residual for electroneutrality condition
            for(comp=0;comp<Nc-1;comp++)
            {

                Resc = al[al_index(x,y,0,comp,Nx,0,0)]*cz(c,z_charge,x,y,0,Nx,0,0,comp,user)+user->con_vars->zo[phi_index(0,0,0,comp,
                                                                                                                          Nx,0,0)]*user->con_vars->ao[phi_index(
                        0,0,0,comp,Nx,0,0)];
                ierr = VecSetValue(Res,Ind_1(x,y,0,Ni,comp,Nx,0,0),Resc,INSERT_VALUES); CHKERRQ(ierr);
            }
            //Extracellular term
            comp=Nc-1;
            Resc = (1-al[al_index(x,y,0,0,Nx,0,0)]-al[al_index(x,y,0,1,Nx,0,0)])*
                           cz(c,z_charge,x,y,0,Nx,0,0,comp,user)+user->con_vars->zo[phi_index(
                    0,0,0,comp,Nx,0,0)]*user->con_vars->ao[phi_index(
                    0,0,0,comp,Nx,0,0)];
            ierr = VecSetValue(Res,Ind_1(x,y,0,Ni,comp,Nx,0,0),Resc,INSERT_VALUES); CHKERRQ(ierr);
        }
    }
    //Assemble before we add values in on top to modify the electroneutral.
    ierr = VecAssemblyBegin(Res);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Res);CHKERRQ(ierr);

    for(x=0;x<Nx;x++)
    {
        for(y=0;y<Ny;y++)
        {
            // Add Modification to electroneutrality for non-zero mem.compacitance
            for(comp=0;comp<Nc-1;comp++)
            {
                //Extracell voltage
                ierr = VecSetValue(Res,Ind_1(x,y,0,Ni,Nc-1,Nx,0,0),-cm[al_index(x,y,0,comp,Nx,Ny,Nz)]*(phi[phi_index(x,y,0,
                                                                                                                           Nc-
                                                                                                                           1,Nx,0,0)]-phi[phi_index(
                        x,y,
                        0,
                        comp,
                        Nx,0,0)]),ADD_VALUES);CHKERRQ(ierr);
                //Intracell voltage mod
                ierr = VecSetValue(Res,Ind_1(x,y,0,Ni,comp,Nx,0,0),-cm[al_index(x,y,0,comp,Nx,Ny,Nz)]*(phi[phi_index(x,y,0,comp,Nx,0,0)]-phi[phi_index(
                        x,y,
                        0,
                        Nc-
                        1,Nx,
                        0,0)]),ADD_VALUES);CHKERRQ(ierr);
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
calc_jacobian_linear_algebraic(SNES snes,Vec current_state, Mat A, Mat Jac,void *ctx)
{
    //Linear discretization
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
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt Nz = user->Nz;
    struct ConstVars *con_vars = user->con_vars;
    PetscReal *cm=con_vars->cm;

    PetscInt ind = 0;
    PetscInt x,y,ion,comp;

    PetscReal Fc0x,Fc1x,Fph0x,Fph1x;
    PetscReal Fc0y,Fc1y,Fph0y,Fph1y;
    PetscReal Ac,Aphi;

    PetscInt iter;
    ierr = SNESGetIterationNumber(snes,&iter); CHKERRQ(ierr);
    //Only calculate off diagonals on first iteration
    // For next iterations only diagonal changes
    if(iter==0) {
        //Ionic concentration equations
        for (x = 0; x < Nx; x++) {
            for (y = 0; y < Ny; y++) {
                for (ion = 0; ion < Ni; ion++) {
                    for (comp = 0; comp < Nc - 1; comp++) {
                        //Electrodiffusion contributions
                        Fc0x = 0;
                        Fc1x = 0;
                        Fph0x = 0;
                        Fph1x = 0;
                        Fc0y = 0;
                        Fc1y = 0;
                        Fph0y = 0;
                        Fph1y = 0;
                        if (x < Nx - 1) {
                            Fc0x = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]/dx*dt/dx;
                            Fph0x = z_charge[ion]*Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]*
                                    (cp[c_index(x,y,0,comp,ion,Nx,0,0)]+cp[c_index(x+1,y,0,comp,ion,Nx,
                                                                                   0,0)])/2/dx*dt/
                                    dx;
                            // Right c with left c (-Fc0x)
                            ierr = MatSetValue(Jac,Ind_1(x+1,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,
                                                                                        0,0),-Fc0x,
                                               INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Right c with left phi (-Fph0x)
                            ierr = MatSetValue(Jac,Ind_1(x+1,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-Fph0x,
                                               INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;

                        }
                        if (x > 0) {
                            Fc1x = Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]/dx*dt/dx;
                            Fph1x = z_charge[ion]*Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]*
                                    (cp[c_index(x-1,y,0,comp,ion,Nx,0,0)]+cp[c_index(x,y,0,comp,ion,Nx,
                                                                                     0,0)])/2/dx*dt/
                                    dx;
                            //left c with right c (-Fc1x)
                            ierr = MatSetValue(Jac,Ind_1(x-1,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,
                                                                                        0,0),-Fc1x,
                                               INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Left c with right phi (-Fph1x)
                            ierr = MatSetValue(Jac,Ind_1(x-1,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-Fph1x,
                                               INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                        if (y < Ny - 1) {
                            Fc0y = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]/dy*dt/dy;
                            Fph0y = z_charge[ion]*Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]*
                                    (cp[c_index(x,y,0,comp,ion,Nx,0,0)]+cp[c_index(x,y+1,0,comp,ion,Nx,
                                                                                   0,0)])/2/dy*dt/
                                    dy;
                            // Upper c with lower c (-Fc0y)
                            ierr = MatSetValue(Jac,Ind_1(x,y+1,0,ion,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,
                                                                                        0,0),-Fc0y,
                                               INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Upper c with lower phi (-Fph0y)
                            ierr = MatSetValue(Jac,Ind_1(x,y+1,0,ion,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-Fph0y,
                                               INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                        if (y > 0) {
                            Fc1y = Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]/dy*dt/dy;
                            Fph1y = z_charge[ion]*Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]*
                                    (cp[c_index(x,y-1,0,comp,ion,Nx,0,0)]+cp[c_index(x,y,0,comp,ion,Nx,
                                                                                     0,0)])/2/dy*dt/
                                    dy;
                            //Lower c with Upper c (-Fc1y)
                            ierr = MatSetValue(Jac,Ind_1(x,y-1,0,ion,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,
                                                                                        0,0),-Fc1y,
                                               INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                            //Lower c with Upper phi (-Fph1y)
                            ierr = MatSetValue(Jac,Ind_1(x,y-1,0,ion,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-Fph1y,
                                               INSERT_VALUES);
                            CHKERRQ(ierr);
                            ind++;
                        }
                        //Diagonal term contribution
                        Ac = al[al_index(x,y,0,comp,Nx,0,0)]+Fc0x+Fc1x+Fc0y+Fc1y;
                        Aphi = Fph0x + Fph1x + Fph0y + Fph1y;


                        //membrane current contributions
                        Ac += flux->dfdci[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                        Aphi += flux->dfdphim[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                        // Different Compartment Terms
                        // C Extracellular with C Inside
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,Nc-1,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),
                                           -flux->dfdci[c_index(x,y,0,comp,ion,Nx,0,0)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // C Intra with C Extra
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,ion,Nc-1,Nx,0,0),
                                           flux->dfdce[c_index(x,y,0,comp,ion,Nx,0,0)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // C Extracellular with Phi Inside
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,Nc-1,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),
                                           -flux->dfdphim[c_index(x,y,0,comp,ion,Nx,0,0)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // C Intra with Phi Extra
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,Ni,Nc-1,Nx,0,0),
                                           -flux->dfdphim[c_index(x,y,0,comp,ion,Nx,0,0)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Same compartment terms
                        // c with c
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),Ac,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // c with phi
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),Aphi,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                    }
                    //Extracellular terms
                    comp = Nc - 1;
                    //Electrodiffusion contributions
                    Fc0x = 0;
                    Fc1x = 0;
                    Fph0x = 0;
                    Fph1x = 0;
                    Fc0y = 0;
                    Fc1y = 0;
                    Fph0y = 0;
                    Fph1y = 0;
                    if (x < Nx - 1) {
                        Fc0x = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]/dx*dt/dx;
                        Fph0x = z_charge[ion]*Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]*
                                (cp[c_index(x,y,0,comp,ion,Nx,0,0)]+cp[c_index(x+1,y,0,comp,ion,Nx,0,0)])/2/dx*dt/dx;
                        // Right c with left c (-Fc0x)
                        ierr = MatSetValue(Jac,Ind_1(x+1,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),-Fc0x,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Right c with left phi (-Fph0x)
                        ierr = MatSetValue(Jac,Ind_1(x+1,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-Fph0x,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if (x > 0) {
                        Fc1x = Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]/dx*dt/dx;
                        Fph1x = z_charge[ion]*Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]*
                                (cp[c_index(x-1,y,0,comp,ion,Nx,0,0)]+cp[c_index(x,y,0,comp,ion,Nx,0,0)])/2/dx*dt/dx;
                        //left c with right c (-Fc1x)
                        ierr = MatSetValue(Jac,Ind_1(x-1,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),-Fc1x,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Left c with right phi (-Fph1x)
                        ierr = MatSetValue(Jac,Ind_1(x-1,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-Fph1x,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if (y < Ny - 1) {
                        Fc0y = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]/dy*dt/dy;
                        Fph0y = z_charge[ion]*Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]*
                                (cp[c_index(x,y,0,comp,ion,Nx,0,0)]+cp[c_index(x,y+1,0,comp,ion,Nx,0,0)])/2/dy*dt/dy;
                        // Upper c with lower c (-Fc0y)
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,0,ion,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),-Fc0y,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Upper c with lower phi (-Fph0y)
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,0,ion,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-Fph0y,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    if (y > 0) {
                        Fc1y = Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]/dy*dt/dy;
                        Fph1y = z_charge[ion]*Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]*
                                (cp[c_index(x,y-1,0,comp,ion,Nx,0,0)]+cp[c_index(x,y,0,comp,ion,Nx,0,0)])/2/dy*dt/dy;
                        //Lower c with Upper c (-Fc1y)
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,0,ion,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),-Fc1y,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Lower c with Upper phi (-Fph1y)
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,0,ion,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-Fph1y,
                                           INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }

                    //Diagonal term contribution
                    Ac = (1-al[al_index(x,y,0,0,Nx,0,0)]-al[al_index(x,y,0,1,Nx,0,0)])+Fc0x+Fc1x+Fc0y+Fc1y;
                    Aphi = Fph0x + Fph1x + Fph0y + Fph1y;

                    //Membrane current contribution
                    for (comp = 0; comp < Nc - 1; comp++) {
                        Ac -= flux->dfdce[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                        Aphi += flux->dfdphim[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                    }
                    //Add bath contributions
                    Ac -= sqrt(pow(Dcb[c_index(x,y,0,Nc-1,ion,Nx,0,0)*2],2)+
                               pow(Dcb[c_index(x,y,0,Nc-1,ion,Nx,0,0)*2+1],2))*dt;
                    Aphi -= sqrt(pow(Dcb[c_index(x,y,0,Nc-1,ion,Nx,0,0)*2],2)+
                                 pow(Dcb[c_index(x,y,0,Nc-1,ion,Nx,0,0)*2+1],2))*
                            (cp[c_index(x,y,0,Nc-1,ion,Nx,0,0)]+cbath[ion])*z_charge[ion]/2*dt;

                    //Insert extracell to extracell parts
                    // c with c
                    ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,Nc-1,Nx,0,0),Ind_1(x,y,0,ion,Nc-1,Nx,0,0),Ac,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    // c with phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,Nc-1,Nx,0,0),Ind_1(x,y,0,Ni,Nc-1,Nx,0,0),Aphi,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }

            }
        }

        //Electroneutrality charge-capacitence condition
        for (x = 0; x < Nx; x++) {
            for (y = 0; y < Ny; y++) {
                //electroneutral-concentration entries
                for (ion = 0; ion < Ni; ion++) {
                    for (comp = 0; comp < Nc - 1; comp++) {
                        //Phi with C entries
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),
                                           z_charge[ion] * al[al_index(x,y,0,comp,Nx,0,0)],INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                    }
                    //Phi with C extracellular one
                    comp = Nc - 1;
                    ierr = MatSetValue(Jac,Ind_1(x,y,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),
                                       z_charge[ion] * (1-al[al_index(x,y,0,0,Nx,0,0)]-al[al_index(x,y,0,1,Nx,0,0)]),INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;

                }
                //electroneutrality-voltage entries
                Aphi = 0;
                for (comp = 0; comp < Nc - 1; comp++) {
                    Aphi -= cm[al_index(x,y,0,comp,Nx,Ny,Nz)];
                }
                //extraphi with extra phi
                ierr = MatSetValue(Jac,Ind_1(x,y,0,Ni,Nc-1,Nx,0,0),Ind_1(x,y,0,Ni,Nc-1,Nx,0,0),Aphi,INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;
                for (comp = 0; comp < Nc - 1; comp++) {
                    //Extra phi with intra phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,0,Ni,Nc-1,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),cm[al_index(x,y,0,comp,Nx,Ny,Nz)],INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    // Intra phi with Extraphi
                    ierr = MatSetValue(Jac,Ind_1(x,y,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,Ni,Nc-1,Nx,0,0),cm[al_index(x,y,0,comp,Nx,Ny,Nz)],INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Intra phi with Intra phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-cm[al_index(x,y,0,comp,Nx,Ny,Nz)],INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
            }
        }
    } else {
        //Ionic concentration equations
        for (x = 0; x < Nx; x++) {
            for (y = 0; y < Ny; y++) {
                for (ion = 0; ion < Ni; ion++) {
                    for (comp = 0; comp < Nc - 1; comp++) {
                        //Electrodiffusion contributions
                        Fc0x = 0;
                        Fc1x = 0;
                        Fph0x = 0;
                        Fph1x = 0;
                        Fc0y = 0;
                        Fc1y = 0;
                        Fph0y = 0;
                        Fph1y = 0;
                        if (x < Nx - 1) {
                            Fc0x = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]/dx*dt/dx;
                            Fph0x = z_charge[ion]*Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]*(cp[c_index(x,y,0,comp,
                                                                                                    ion,Nx,0,0)]+cp[c_index(
                                    x+1,y,0,comp,ion,Nx,0,0)])/2/dx*dt/dx;
                        }
                        if (x > 0) {
                            Fc1x = Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]/dx*dt/dx;
                            Fph1x = z_charge[ion]*Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]*(cp[c_index(
                                    x-1,y,0,
                                    comp,ion,
                                    Nx,0,0)]+cp[c_index(
                                    x,y,0,comp,ion,Nx,0,0)])/2/dx*dt/dx;

                        }
                        if (y < Ny - 1) {
                            Fc0y = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]/dy*dt/dy;
                            Fph0y = z_charge[ion]*Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]*(cp[c_index(x,y,0,
                                                                                                      comp,ion,
                                                                                                      Nx,0,0)]+cp[c_index(
                                    x,y+1,0,comp,ion,Nx,0,0)])/2/dy*dt/dy;
                        }
                        if (y > 0) {
                            Fc1y = Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]/dy*dt/dy;
                            Fph1y = z_charge[ion]*Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]*(cp[c_index(x,
                                                                                                        y-
                                                                                                        1,
                                                                                                        0,
                                                                                                        comp,
                                                                                                        ion,Nx,
                                                                                                        0,0)]+cp[c_index(
                                    x,y,0,comp,ion,Nx,0,0)])/2/dy*dt/dy;
                        }
                        //Diagonal term contribution
                        Ac = al[al_index(x,y,0,comp,Nx,0,0)]+Fc0x+Fc1x+Fc0y+Fc1y;
                        Aphi = Fph0x + Fph1x + Fph0y + Fph1y;


                        //membrane current contributions
                        Ac += flux->dfdci[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                        Aphi += flux->dfdphim[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                        // Different Compartment Terms
                        // C Extracellular with C Inside
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,Nc-1,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),-flux->dfdci[c_index(
                                x,y,0,comp,ion,Nx,0,0)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // C Intra with C Extra
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,ion,Nc-1,Nx,0,0),
                                           flux->dfdce[c_index(x,y,0,comp,ion,Nx,0,0)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // C Extracellular with Phi Inside
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,Nc-1,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),
                                           -flux->dfdphim[c_index(x,y,0,comp,ion,Nx,0,0)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // C Intra with Phi Extra
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,Ni,Nc-1,Nx,0,0),
                                           -flux->dfdphim[c_index(x,y,0,comp,ion,Nx,0,0)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Same compartment terms
                        // c with c
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),Ac,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // c with phi
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),Aphi,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                    }
                    //Extracellular terms
                    comp = Nc - 1;
                    //Electrodiffusion contributions
                    Fc0x = 0;
                    Fc1x = 0;
                    Fph0x = 0;
                    Fph1x = 0;
                    Fc0y = 0;
                    Fc1y = 0;
                    Fph0y = 0;
                    Fph1y = 0;
                    if (x < Nx - 1) {
                        Fc0x = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]/dx*dt/dx;
                        Fph0x = z_charge[ion]*Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]*(cp[c_index(x,y,0,comp,ion,
                                                                                                Nx,0,0)]+cp[c_index(
                                x+1,y,0,comp,ion,Nx,0,0)])/2/dx*dt/dx;
                    }
                    if (x > 0) {
                        Fc1x = Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]/dx*dt/dx;
                        Fph1x = z_charge[ion]*Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]*(cp[c_index(x-1,y,0,
                                                                                                  comp,ion,
                                                                                                  Nx,0,0)]+cp[c_index(
                                x,y,0,comp,ion,Nx,0,0)])/2/dx*dt/dx;
                    }
                    if (y < Ny - 1) {
                        Fc0y = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]/dy*dt/dy;
                        Fph0y = z_charge[ion]*Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]*(cp[c_index(x,y,0,comp,
                                                                                                  ion,Nx,
                                                                                                  0,0)]+cp[c_index(
                                x,y+1,0,comp,ion,Nx,0,0)])/2/dy*dt/dy;
                    }
                    if (y > 0) {
                        Fc1y = Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]/dy*dt/dy;
                        Fph1y = z_charge[ion]*Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]*(cp[c_index(x,
                                                                                                    y-1,0,
                                                                                                    comp,
                                                                                                    ion,Nx,
                                                                                                    0,0)]+cp[c_index(
                                x,y,0,comp,ion,Nx,0,0)])/2/dy*dt/dy;
                    }

                    //Diagonal term contribution
                    Ac = (1-al[al_index(x,y,0,0,Nx,0,0)]-al[al_index(x,y,0,1,Nx,0,0)])+Fc0x+Fc1x+Fc0y+Fc1y;
                    Aphi = Fph0x + Fph1x + Fph0y + Fph1y;

                    //Membrane current contribution
                    for (comp = 0; comp < Nc - 1; comp++) {
                        Ac -= flux->dfdce[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                        Aphi += flux->dfdphim[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                    }
                    //Add bath contributions
                    Ac -= sqrt(pow(Dcb[c_index(x,y,0,Nc-1,ion,Nx,0,0)*2],2)+
                               pow(Dcb[c_index(x,y,0,Nc-1,ion,Nx,0,0)*2+1],2))*dt;
                    Aphi -= sqrt(pow(Dcb[c_index(x,y,0,Nc-1,ion,Nx,0,0)*2],2)+
                                 pow(Dcb[c_index(x,y,0,Nc-1,ion,Nx,0,0)*2+1],2))*
                            (cp[c_index(x,y,0,Nc-1,ion,Nx,0,0)]+cbath[ion])*z_charge[ion]/2*dt;

                    //Insert extracell to extracell parts
                    // c with c
                    ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,Nc-1,Nx,0,0),Ind_1(x,y,0,ion,Nc-1,Nx,0,0),Ac,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    // c with phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,Nc-1,Nx,0,0),Ind_1(x,y,0,Ni,Nc-1,Nx,0,0),Aphi,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }

            }
        }
    }
    ierr = MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

//    printf("Iter: %d, Number of inserts: %d\n",iter,ind);
    if (A != Jac) {
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr); }

    ierr = restore_subarray_Read(current_state,user->state_vars); CHKERRQ(ierr);
    if(Profiling_on) {
        PetscLogEventEnd(event[0], 0, 0, 0, 0);
    }
    return ierr;
}

PetscErrorCode calc_residual_linear_deriv(SNES snes,Vec current_state,Vec Res,void *ctx)
{
    //Linear discretization
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
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt Nz = user->Nz;

    //Residual for concentration equations
    PetscReal Rcvx,Rcvy,Resc;
    PetscReal RcvxRight,RcvyUp;

    //Residual for fluxes in voltage differential equations
    PetscReal Rphx[Nc], Rphy[Nc], RphxRight[Nc], RphyUp[Nc];
    PetscReal Resph,ResphN;

    PetscReal alNc,alpNc;
    PetscInt ion,comp,x,y;


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
                        // C term
                        Rcvx = Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]*(c[c_index(x,y,0,comp,ion,Nx,
                                                                                  0,0)]-c[c_index(x-1,y,0,
                                                                                                  comp,ion,
                                                                                                  Nx,
                                                                                                  0,0)]);
                        //Phi term
                        Rcvx += Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]*z_charge[ion]*(cp[c_index(x-1,y,0,
                                                                                                  comp,ion,Nx,
                                                                                                  0,0)]+cp[c_index(
                                x,
                                y,
                                0,
                                comp,
                                ion,
                                Nx,
                                0,0)])/2*(phi[phi_index(x,y,0,comp,Nx,0,0)]-phi[phi_index(x-1,y,0,comp,Nx,
                                                                                          0,0)]);
                        Rcvx = Rcvx/dx*dt/dx;
                    }
                    //Add Second right moving difference
                    if(x<Nx-1) {
                        // C term
                        RcvxRight = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]*(c[c_index(x+1,y,0,comp,ion,Nx,
                                                                                     0,0)]-c[c_index(x,y,0,
                                                                                                     comp,ion,
                                                                                                     Nx,0,0)]);
                        //Phi term
                        RcvxRight += Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]*z_charge[ion]*(cp[c_index(x+1,y,0,
                                                                                                     comp,ion,Nx,
                                                                                                     0,0)]+cp[c_index(
                                x,y,0,comp,ion,Nx,0,0)])/2*(phi[phi_index(x+1,y,0,comp,Nx,0,0)]-phi[phi_index(
                                x,y,0,comp,Nx,0,0)]);
                        RcvxRight = RcvxRight/dx*dt/dx;
                    }
                    Rcvy = 0;
                    RcvyUp = 0;
                    //Up down difference
                    if(y>0) {
                        // C term
                        Rcvy = Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]*(c[c_index(x,y,0,comp,ion,Nx,
                                                                                    0,0)]-c[c_index(x,
                                                                                                    y-1,0,
                                                                                                    comp,
                                                                                                    ion,Nx,
                                                                                                    0,0)]);
                        //Phi term
                        Rcvy += Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]*z_charge[ion]*(cp[c_index(x,
                                                                                                    y-1,0,
                                                                                                    comp,ion,
                                                                                                    Nx,0,0)]+cp[c_index(
                                x,
                                y,
                                0,
                                comp,
                                ion,
                                Nx,
                                0,0)])/2*(phi[phi_index(x,y,0,comp,Nx,0,0)]-phi[phi_index(x,y-1,0,comp,Nx,
                                                                                          0,0)]);
                        Rcvy = Rcvy/dy*dt/dy;
                    }
                    //Next upward difference
                    if(y<Ny-1) {
                        // C term
                        RcvyUp = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]*(c[c_index(x,y+1,0,comp,ion,
                                                                                    Nx,0,0)]-c[c_index(x,y,
                                                                                                       0,
                                                                                                       comp,
                                                                                                       ion,
                                                                                                       Nx,0,0)]);
                        //Phi term
                        RcvyUp += Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]*z_charge[ion]*(cp[c_index(x,y+1,0,
                                                                                                    comp,ion,Nx,
                                                                                                    0,0)]+cp[c_index(
                                x,
                                y,
                                0,
                                comp,
                                ion,
                                Nx,
                                0,0)])/2*(phi[phi_index(x,y+1,0,comp,Nx,0,0)]-phi[phi_index(x,y,0,comp,Nx,
                                                                                            0,0)]);
                        RcvyUp = RcvyUp/dy*dt/dy;
                    }
                    Resc = al[al_index(x,y,0,comp,Nx,0,0)]*c[c_index(x,y,0,comp,ion,Nx,0,0)]-alp[al_index(x,
                                                                                                          y,
                                                                                                          0,
                                                                                                          comp,
                                                                                                          Nx,
                                                                                                          0,0)]*cp[c_index(
                            x,y,0,comp,ion,Nx,0,0)];
                    Resc += Rcvx - RcvxRight + Rcvy - RcvyUp +flux->mflux[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;

                    ierr = VecSetValue(Res,Ind_1(x,y,0,ion,comp,Nx,0,0),Resc,INSERT_VALUES);CHKERRQ(ierr);

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
                    // C term
                    Rcvx = Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]*(c[c_index(x,y,0,comp,ion,Nx,0,0)]-c[c_index(
                            x-1,y,0,comp,
                            ion,Nx,0,0)]);
                    //Phi term
                    Rcvx += Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]*z_charge[ion]*(cp[c_index(x-1,y,0,comp,
                                                                                              ion,Nx,0,0)]+cp[c_index(
                            x,y,0,
                            comp,
                            ion,Nx,
                            0,0)])/2*(phi[phi_index(x,y,0,comp,Nx,0,0)]-phi[phi_index(x-1,y,0,comp,Nx,0,0)]);
                    Rcvx = Rcvx/dx*dt/dx;
                }
                //Add Second right moving difference
                if(x<Nx-1) {
                    // C term
                    RcvxRight = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]*(c[c_index(x+1,y,0,comp,ion,Nx,
                                                                                 0,0)]-c[c_index(x,y,0,comp,
                                                                                                 ion,Nx,0,0)]);
                    //Phi term
                    RcvxRight += Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]*z_charge[ion]*(cp[c_index(x+1,y,0,comp,
                                                                                                 ion,Nx,0,0)]+cp[c_index(
                            x,y,
                            0,
                            comp,
                            ion,
                            Nx,
                            0,0)])/2*(phi[phi_index(x+1,y,0,comp,Nx,0,0)]-phi[phi_index(x,y,0,comp,Nx,0,0)]);
                    RcvxRight = RcvxRight/dx*dt/dx;
                }
                Rcvy = 0;
                RcvyUp = 0;
                //Up down difference
                if(y>0) {
                    // C term
                    Rcvy = Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]*(c[c_index(x,y,0,comp,ion,Nx,
                                                                                0,0)]-c[c_index(x,y-1,0,
                                                                                                comp,
                                                                                                ion,Nx,0,0)]);
                    //Phi term
                    Rcvy += Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]*z_charge[ion]*(cp[c_index(x,y-1,0,
                                                                                                comp,ion,Nx,
                                                                                                0,0)]+cp[c_index(
                            x,y,
                            0,
                            comp,
                            ion,
                            Nx,0,0)])/2*(phi[phi_index(x,y,0,comp,Nx,0,0)]-phi[phi_index(x,y-1,0,comp,Nx,
                                                                                         0,0)]);
                    Rcvy = Rcvy/dy*dt/dy;
                }
                //Next upward difference
                if(y<Ny-1) {
                    // C term
                    RcvyUp = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]*(c[c_index(x,y+1,0,comp,ion,Nx,
                                                                                0,0)]-c[c_index(x,y,0,comp,
                                                                                                ion,Nx,0,0)]);
                    //Phi term
                    RcvyUp += Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]*z_charge[ion]*(cp[c_index(x,
                                                                                                y+1,0,comp,
                                                                                                ion,Nx,0,0)]+cp[c_index(
                            x,y,
                            0,
                            comp,
                            ion,
                            Nx,0,0)])/2*(phi[phi_index(x,y+1,0,comp,Nx,0,0)]-phi[phi_index(x,y,0,comp,Nx,
                                                                                           0,0)]);
                    RcvyUp = RcvyUp/dy*dt/dy;
                }
                Resc = alNc*c[c_index(x,y,0,comp,ion,Nx,0,0)]-alpNc*cp[c_index(x,y,0,comp,ion,Nx,0,0)];
                Resc += Rcvx - RcvxRight + Rcvy - RcvyUp +flux->mflux[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                //Add bath variables

                Resc -= sqrt(pow(Dcb[c_index(x,y,0,comp,ion,Nx,0,0)*2],2)+pow(Dcb[c_index(x,y,0,comp,ion,
                                                                                          Nx,0,0)*2+1],2))*(cp[c_index(
                        x,y,0,comp,ion,Nx,0,0)]+cbath[ion])/2.0*(log(c[c_index(
                        x,y,0,comp,ion,Nx,0,0)])-log(cbath[ion])+z_charge[ion]*phi[phi_index(x,y,0,comp,Nx,0,0)]-z_charge[ion]*phibath)*dt;
                ierr = VecSetValue(Res,Ind_1(x,y,0,ion,comp,Nx,0,0),Resc,INSERT_VALUES);CHKERRQ(ierr);

                //Save values for voltage
                Rphx[comp]+=z_charge[ion]*Rcvx;
                Rphy[comp]+=z_charge[ion]*Rcvy;
                RphxRight[comp]+=z_charge[ion]*RcvxRight;
                RphyUp[comp]+=z_charge[ion]*RcvyUp;
            }

            //Voltage Equations
            ResphN = 0;
            for(comp=0;comp<Nc-1;comp++) {
                Resph = cm[al_index(x,y,0,comp,Nx,Ny,Nz)]*(phi[phi_index(x,y,0,comp,Nx,0,0)]-phi[phi_index(x,y,0,Nc-
                                                                                                                1,Nx,0,0)])-cm[al_index(x,y,0,comp,Nx,Ny,Nz)]*(phip[phi_index(
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
                ierr = VecSetValue(Res,Ind_1(x,y,0,Ni,comp,Nx,0,0),Resph,INSERT_VALUES); CHKERRQ(ierr);
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
            ierr = VecSetValue(Res,Ind_1(x,y,0,Ni,comp,Nx,0,0),ResphN,INSERT_VALUES); CHKERRQ(ierr);
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
calc_jacobian_linear_deriv(SNES snes,Vec current_state, Mat A, Mat Jac,void *ctx)
{
    //Linear discretization
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
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt Nz = user->Nz;
    struct ConstVars *con_vars = user->con_vars;
    PetscReal *cm = con_vars->cm;

    PetscInt ind = 0;
    PetscInt x,y,ion,comp;

    PetscReal Ftmpx,Fc0x,Fc1x,Fph0x,Fph1x;
    PetscReal Fc0y,Fc1y,Fph0y,Fph1y;
    PetscReal Ac,Aphi,Avolt,AvoltN;

    PetscReal Fphph0x[Nc],Fphph1x[Nc];
    PetscReal Fphph0y[Nc],Fphph1y[Nc];

    PetscInt iter;
    ierr = SNESGetIterationNumber(snes,&iter); CHKERRQ(ierr);

    if(iter==0){
        for(x=0;x<Nx;x++) {
            for(y=0;y<Ny;y++) {
                for(comp=0;comp<Nc;comp++){
                    Fphph0x[comp]=0;
                    Fphph1x[comp]=0;
                    Fphph0y[comp]=0;
                    Fphph1y[comp]=0;
                }
                for(ion=0;ion<Ni;ion++) {
                    for(comp=0;comp<Nc-1;comp++) {
                        //Electrodiffusion contributions
                        Fc0x = 0;
                        Fc1x = 0;
                        Fph0x = 0;
                        Fph1x = 0;
                        Fc0y = 0;
                        Fc1y = 0;
                        Fph0y = 0;
                        Fph1y = 0;
                        if(x<Nx-1) {
                            Fc0x = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]/dx*dt/dx;
                            Fph0x = z_charge[ion]*Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]*
                                    (cp[c_index(x,y,0,comp,ion,Nx,0,0)]+cp[c_index(x+1,y,0,comp,ion,Nx,
                                                                                   0,0)])/2/dx*dt/dx;
                            // Right c with left c (-Fc0x)
                            ierr = MatSetValue(Jac,Ind_1(x+1,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,
                                                                                        0,0),-Fc0x,INSERT_VALUES);CHKERRQ(ierr);
                            ind++;
                            //Right c with left phi (-Fph0x)
                            ierr = MatSetValue(Jac,Ind_1(x+1,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-Fph0x,INSERT_VALUES);CHKERRQ(ierr);
                            ind++;

                            //Right phi with left c in voltage eqn
                            ierr = MatSetValue(Jac,Ind_1(
                                    x+1,y,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),-z_charge[ion]*Fc0x,INSERT_VALUES);CHKERRQ(ierr);
                            ind++;
                        }
                        if(x>0) {
                            Fc1x = Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]/dx*dt/dx;
                            Fph1x = z_charge[ion]*Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]*
                                    (cp[c_index(x-1,y,0,comp,ion,Nx,0,0)]+cp[c_index(x,y,0,comp,ion,Nx,
                                                                                     0,0)])/2/dx*dt/dx;
                            //left c with right c (-Fc1x)
                            ierr = MatSetValue(Jac,Ind_1(x-1,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,
                                                                                        0,0),-Fc1x,INSERT_VALUES);CHKERRQ(ierr);
                            ind++;
                            //Left c with right phi (-Fph1x)
                            ierr = MatSetValue(Jac,Ind_1(x-1,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-Fph1x,INSERT_VALUES);CHKERRQ(ierr);
                            ind++;

                            //Left phi with right c in voltage eqn
                            ierr = MatSetValue(Jac,Ind_1(
                                    x-1,y,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),-z_charge[ion]*Fc1x,INSERT_VALUES);CHKERRQ(ierr);
                            ind++;
                        }
                        if(y<Ny-1) {
                            Fc0y = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]/dy*dt/dy;
                            Fph0y = z_charge[ion]*Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]*
                                    (cp[c_index(x,y,0,comp,ion,Nx,0,0)]+cp[c_index(x,y+1,0,comp,ion,Nx,
                                                                                   0,0)])/2/dy*dt/dy;
                            // Upper c with lower c (-Fc0y)
                            ierr = MatSetValue(Jac,Ind_1(x,y+1,0,ion,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,
                                                                                        0,0),-Fc0y,INSERT_VALUES);CHKERRQ(ierr);
                            ind++;
                            //Upper c with lower phi (-Fph0y)
                            ierr = MatSetValue(Jac,Ind_1(x,y+1,0,ion,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-Fph0y,INSERT_VALUES);CHKERRQ(ierr);
                            ind++;

                            //Upper phi with lower c in voltage eqn
                            ierr = MatSetValue(Jac,Ind_1(x,y+1,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),-z_charge[ion]*Fc0y,INSERT_VALUES);CHKERRQ(ierr);
                            ind++;
                        }
                        if(y>0) {
                            Fc1y = Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]/dy*dt/dy;
                            Fph1y = z_charge[ion]*Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]*
                                    (cp[c_index(x,y-1,0,comp,ion,Nx,0,0)]+cp[c_index(x,y,0,comp,ion,Nx,
                                                                                     0,0)])/2/dy*dt/dy;
                            //Lower c with Upper c (-Fc1y)
                            ierr = MatSetValue(Jac,Ind_1(x,y-1,0,ion,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,
                                                                                        0,0),-Fc1y,INSERT_VALUES);CHKERRQ(ierr);
                            ind++;
                            //Lower c with Upper phi (-Fph1y)
                            ierr = MatSetValue(Jac,Ind_1(x,y-1,0,ion,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-Fph1y,INSERT_VALUES);CHKERRQ(ierr);
                            ind++;

                            //Lower phi with upper c in voltage eqn
                            ierr = MatSetValue(Jac,Ind_1(x,y-1,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),-z_charge[ion]*Fc1y,INSERT_VALUES);CHKERRQ(ierr);
                            ind++;
                        }
                        //Diagonal term contribution
                        Ac = al[al_index(x,y,0,comp,Nx,0,0)]+Fc0x+Fc1x+Fc0y+Fc1y;
                        Aphi = Fph0x + Fph1x + Fph0y + Fph1y;

                        //Add up terms for voltage eqns
                        Fphph0x[comp]+=z_charge[ion]*Fph0x;
                        Fphph1x[comp]+=z_charge[ion]*Fph1x;
                        Fphph0y[comp]+=z_charge[ion]*Fph0y;
                        Fphph1y[comp]+=z_charge[ion]*Fph1y;

                        //membrane current contributions
                        Ac+= flux->dfdci[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                        Aphi+= flux->dfdphim[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                        // Different Compartment Terms
                        // C Extracellular with C Inside
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,Nc-1,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),-flux->dfdci[c_index(x,
                                                                                                                              y,
                                                                                                                              0,
                                                                                                                              comp,
                                                                                                                              ion,
                                                                                                                              Nx,
                                                                                                                              0,0)]*dt,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        // C Intra with C Extra
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,ion,Nc-1,Nx,0,0),flux->dfdce[c_index(x,
                                                                                                                             y,
                                                                                                                             0,
                                                                                                                             comp,
                                                                                                                             ion,
                                                                                                                             Nx,
                                                                                                                             0,0)]*dt,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        // C Extracellular with Phi Inside
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,Nc-1,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-flux->dfdphim[c_index(x,
                                                                                                                               y,
                                                                                                                               0,
                                                                                                                               comp,
                                                                                                                               ion,
                                                                                                                               Nx,
                                                                                                                               0,0)]*dt,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        // C Intra with Phi Extra
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,Ni,Nc-1,Nx,0,0),-flux->dfdphim[c_index(x,
                                                                                                                               y,
                                                                                                                               0,
                                                                                                                               comp,
                                                                                                                               ion,
                                                                                                                               Nx,
                                                                                                                               0,0)]*dt,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Same compartment terms
                        // c with c
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),Ac,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        // c with phi
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),Aphi,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;

                        //Intra-Phi with c (voltage eqn)
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),z_charge[ion]*(Fc0x+Fc1x+Fc0y+Fc1y+flux->dfdci[c_index(
                                x,y,0,comp,ion,Nx,0,0)]*dt),INSERT_VALUES); CHKERRQ(ierr);
                        ind++;
                        //IntraPhi with c extra(volt eqn)
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,ion,Nc-1,Nx,0,0),z_charge[ion]*(flux->dfdce[c_index(
                                x,y,0,comp,ion,Nx,0,0)]*dt),INSERT_VALUES); CHKERRQ(ierr);
                        ind++;
                        //Extra-Phi with intra-c (voltage eqn)
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,Ni,Nc-1,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),-z_charge[ion]*(flux->dfdci[c_index(
                                x,y,0,comp,ion,Nx,0,0)]*dt),INSERT_VALUES); CHKERRQ(ierr);
                        ind++;

                    }
                    //Extracellular terms
                    comp = Nc-1;
                    //Electrodiffusion contributions
                    Fc0x = 0;
                    Fc1x = 0;
                    Fph0x = 0;
                    Fph1x = 0;
                    Fc0y = 0;
                    Fc1y = 0;
                    Fph0y = 0;
                    Fph1y = 0;
                    if(x<Nx-1) {
                        Fc0x = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]/dx*dt/dx;
                        Fph0x = z_charge[ion]*Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]*
                                (cp[c_index(x,y,0,comp,ion,Nx,0,0)]+cp[c_index(x+1,y,0,comp,ion,Nx,0,0)])/2/dx*dt/dx;
                        // Right c with left c (-Fc0x)
                        ierr = MatSetValue(Jac,Ind_1(
                                x+1,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),-Fc0x,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Right c with left phi (-Fph0x)
                        ierr = MatSetValue(Jac,Ind_1(
                                x+1,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-Fph0x,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;

                        //Right phi with left c in voltage eqn
                        ierr = MatSetValue(Jac,Ind_1(
                                x+1,y,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),-z_charge[ion]*Fc0x,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    if(x>0) {
                        Fc1x = Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]/dx*dt/dx;
                        Fph1x = z_charge[ion]*Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]*
                                (cp[c_index(x-1,y,0,comp,ion,Nx,0,0)]+cp[c_index(x,y,0,comp,ion,Nx,0,0)])/2/dx*dt/dx;
                        //left c with right c (-Fc1x)
                        ierr = MatSetValue(Jac,Ind_1(
                                x-1,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),-Fc1x,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Left c with right phi (-Fph1x)
                        ierr = MatSetValue(Jac,Ind_1(
                                x-1,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-Fph1x,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;

                        //Left phi with right c in voltage eqn
                        ierr = MatSetValue(Jac,Ind_1(
                                x-1,y,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),-z_charge[ion]*Fc1x,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    if(y<Ny-1) {
                        Fc0y = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]/dy*dt/dy;
                        Fph0y = z_charge[ion]*Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]*
                                (cp[c_index(x,y,0,comp,ion,Nx,0,0)]+cp[c_index(x,y+1,0,comp,ion,Nx,0,0)])/2/dy*dt/dy;
                        // Upper c with lower c (-Fc0y)
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,0,ion,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),-Fc0y,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Upper c with lower phi (-Fph0y)
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,0,ion,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-Fph0y,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;

                        //Upper phi with lower c in voltage eqn
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),-z_charge[ion]*Fc0y,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    if(y>0) {
                        Fc1y = Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]/dy*dt/dy;
                        Fph1y = z_charge[ion]*Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]*
                                (cp[c_index(x,y-1,0,comp,ion,Nx,0,0)]+cp[c_index(x,y,0,comp,ion,Nx,0,0)])/2/dy*dt/dy;
                        //Lower c with Upper c (-Fc1y)
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,0,ion,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),-Fc1y,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                        //Lower c with Upper phi (-Fph1y)
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,0,ion,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-Fph1y,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;

                        //Lower phi with upper c in voltage eqn
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),-z_charge[ion]*Fc1y,INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }

                    //Diagonal term contribution
                    Ac = (1-al[al_index(x,y,0,0,Nx,0,0)]-al[al_index(x,y,0,1,Nx,0,0)])+Fc0x+Fc1x+Fc0y+Fc1y;
                    Aphi = Fph0x + Fph1x + Fph0y + Fph1y;

                    Avolt = z_charge[ion]*(Fc0x+Fc1x+Fc0y+Fc1y);

                    //Add up terms for voltage eqns
                    Fphph0x[comp]+=z_charge[ion]*Fph0x;
                    Fphph1x[comp]+=z_charge[ion]*Fph1x;
                    Fphph0y[comp]+=z_charge[ion]*Fph0y;
                    Fphph1y[comp]+=z_charge[ion]*Fph1y;

                    //Membrane current contribution
                    for(comp=0;comp<Nc-1;comp++) {
                        Ac -= flux->dfdce[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                        Aphi += flux->dfdphim[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                        Avolt -= z_charge[ion]*flux->dfdce[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                    }
                    //Add bath contributions
                    Ftmpx=sqrt(pow(Dcb[c_index(x,y,0,Nc-1,ion,Nx,0,0)*2],2)+pow(Dcb[c_index(x,y,0,Nc-1,
                                                                                            ion,Nx,0,0)*2+1],2));
                    Ac -= Ftmpx*(cp[c_index(x,y,0,Nc-1,ion,Nx,0,0)]+cbath[ion])/(2*c[c_index(x,y,0,
                                                                                             Nc-1,ion,
                                                                                             Nx,0,0)])*dt;
                    Aphi -= Ftmpx*(cp[c_index(x,y,0,Nc-1,ion,Nx,0,0)]+cbath[ion])*z_charge[ion]/2*dt;

                    Avolt -= z_charge[ion]*Ftmpx*(cp[c_index(x,y,0,Nc-1,ion,Nx,0,0)]+cbath[ion])/(2*c[c_index(
                            x,y,0,Nc-1,ion,
                            Nx,0,0)])*dt;

                    //Insert extracell to extracell parts
                    // c with c
                    ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,Nc-1,Nx,0,0),Ind_1(x,y,0,ion,Nc-1,Nx,0,0),Ac,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    // c with phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,Nc-1,Nx,0,0),Ind_1(x,y,0,Ni,Nc-1,Nx,0,0),Aphi,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;

                    //phi with c (voltage eqn)
                    ierr = MatSetValue(Jac,Ind_1(x,y,0,Ni,Nc-1,Nx,0,0),Ind_1(x,y,0,ion,Nc-1,Nx,0,0),Avolt,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                //Derivative of charge-capacitance
                for(comp=0;comp<Nc-1;comp++) {
                    if(x<Nx-1) {
                        //Right phi with left phi (-Fph0x)
                        ierr = MatSetValue(Jac,Ind_1(x+1,y,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-Fphph0x[comp],INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    if(x>0) {
                        //Left phi with right phi (-Fph1x)
                        ierr = MatSetValue(Jac,Ind_1(x-1,y,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-Fphph1x[comp],INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    if(y<Ny-1) {
                        //Upper phi with lower phi (-Fph0y)
                        ierr = MatSetValue(Jac,Ind_1(x,y+1,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-Fphph0y[comp],INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    if(y>0) {
                        //Lower phi with upper phi (-Fph1y)
                        ierr = MatSetValue(Jac,Ind_1(x,y-1,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-Fphph1y[comp],INSERT_VALUES);CHKERRQ(ierr);
                        ind++;
                    }
                    Avolt = cm[al_index(x,y,0,comp,Nx,Ny,Nz)]+Fphph0x[comp]+Fphph1x[comp]+Fphph0y[comp]+Fphph1y[comp];
                    AvoltN = -cm[al_index(x,y,0,comp,Nx,Ny,Nz)];
                    for(ion=0;ion<Ni;ion++) {
                        Avolt+= z_charge[ion]*flux->dfdphim[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                        AvoltN-= z_charge[ion]*flux->dfdphim[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                    }

                    //Intra-phi with Intra-phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),Avolt,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                    //Intra-phi with extra-phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,Ni,Nc-1,Nx,0,0),AvoltN,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                //Extracellular terms
                comp = Nc-1;
                if(x<Nx-1) {
                    //Right phi with left phi (-Fph0x)
                    ierr = MatSetValue(Jac,Ind_1(x+1,y,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-Fphph0x[comp],INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(x>0) {
                    //Left phi with right phi (-Fph1x)
                    ierr = MatSetValue(Jac,Ind_1(x-1,y,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-Fphph1x[comp],INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(y<Ny-1) {
                    //Upper phi with lower phi (-Fph0y)
                    ierr = MatSetValue(Jac,Ind_1(x,y+1,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-Fphph0y[comp],INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                if(y>0) {
                    //Lower phi with upper phi (-Fph1y)
                    ierr = MatSetValue(Jac,Ind_1(x,y-1,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),-Fphph1y[comp],INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }
                AvoltN = 0;

                for(int k=0;k<Nc-1;k++) {
                    AvoltN += cm[al_index(x,y,0,comp,Nx,Ny,Nz)];
                    Avolt = -cm[al_index(x,y,0,comp,Nx,Ny,Nz)];
                    for(ion=0;ion<Ni;ion++) {
                        Avolt-= z_charge[ion]*flux->dfdphim[c_index(x,y,0,k,ion,Nx,0,0)]*dt;
                        AvoltN+= z_charge[ion]*flux->dfdphim[c_index(x,y,0,k,ion,Nx,0,0)]*dt;
                    }
                    //Extra-phi with Intra-phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,Ni,k,Nx,0,0),Avolt,INSERT_VALUES);CHKERRQ(ierr);
                    ind++;
                }

                AvoltN += Fphph0x[comp]+Fphph1x[comp]+Fphph0y[comp]+Fphph1y[comp];

                //Bath terms
                for(ion=0;ion<Ni;ion++) {
                    Ftmpx = sqrt(pow(Dcb[c_index(x,y,0,Nc-1,ion,Nx,0,0)*2],2)+pow(Dcb[c_index(x,y,0,
                                                                                              Nc-1,ion,
                                                                                              Nx,0,0)*2+1],2));
                    AvoltN -= z_charge[ion]*Ftmpx*(cp[c_index(x,y,0,Nc-1,ion,Nx,0,0)]+cbath[ion])*z_charge[ion]/2*dt;
                }
                //extra-phi with extra-phi
                ierr = MatSetValue(Jac,Ind_1(x,y,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),AvoltN,INSERT_VALUES);CHKERRQ(ierr);
                ind++;

            }
        }

    }
    else {
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
                        Fc0x = 0;
                        Fc1x = 0;
                        Fph0x = 0;
                        Fph1x = 0;
                        Fc0y = 0;
                        Fc1y = 0;
                        Fph0y = 0;
                        Fph1y = 0;
                        if (x < Nx - 1) {
                            Fc0x = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]/dx*dt/dx;
                            Fph0x = z_charge[ion]*Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]*
                                    (cp[c_index(x,y,0,comp,ion,Nx,0,0)]+cp[c_index(x+1,y,0,comp,ion,Nx,
                                                                                   0,0)])/2/dx*dt/dx;
                        }
                        if (x > 0) {
                            Fc1x = Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]/dx*dt/dx;
                            Fph1x = z_charge[ion]*Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]*
                                    (cp[c_index(x-1,y,0,comp,ion,Nx,0,0)]+cp[c_index(x,y,0,comp,ion,Nx,
                                                                                     0,0)])/2/dx*dt/dx;
                        }
                        if (y < Ny - 1) {
                            Fc0y = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]/dy*dt/dy;
                            Fph0y = z_charge[ion]*Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]*
                                    (cp[c_index(x,y,0,comp,ion,Nx,0,0)]+cp[c_index(x,y+1,0,comp,ion,Nx,
                                                                                   0,0)])/2/dy*dt/dy;
                        }
                        if (y > 0) {
                            Fc1y = Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]/dy*dt/dy;
                            Fph1y = z_charge[ion]*Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]*
                                    (cp[c_index(x,y-1,0,comp,ion,Nx,0,0)]+cp[c_index(x,y,0,comp,ion,Nx,
                                                                                     0,0)])/2/dy*dt/dy;
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
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,Nc-1,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),
                                           -flux->dfdci[c_index(x,y,0,comp,ion,Nx,0,0)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // C Intra with C Extra
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,ion,Nc-1,Nx,0,0),
                                           flux->dfdce[c_index(x,y,0,comp,ion,Nx,0,0)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // C Extracellular with Phi Inside
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,Nc-1,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),
                                           -flux->dfdphim[c_index(x,y,0,comp,ion,Nx,0,0)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // C Intra with Phi Extra
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,Ni,Nc-1,Nx,0,0),
                                           -flux->dfdphim[c_index(x,y,0,comp,ion,Nx,0,0)]*dt,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Same compartment terms
                        // c with c
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),Ac,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        // c with phi
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),Aphi,INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                        //Intra-Phi with c (voltage eqn)
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),z_charge[ion]*(Fc0x+Fc1x+Fc0y+Fc1y+flux->dfdci[c_index(
                                x,y,0,comp,ion,Nx,0,0)]*dt),INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //IntraPhi with c extra(volt eqn)
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,ion,Nc-1,Nx,0,0),z_charge[ion]*(flux->dfdce[c_index(
                                x,y,0,comp,ion,Nx,0,0)]*dt),INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;
                        //Extra-Phi with intra-c (voltage eqn)
                        ierr = MatSetValue(Jac,Ind_1(x,y,0,Ni,Nc-1,Nx,0,0),Ind_1(x,y,0,ion,comp,Nx,0,0),-z_charge[ion]*(flux->dfdci[c_index(
                                x,y,0,comp,ion,Nx,0,0)]*dt),INSERT_VALUES);
                        CHKERRQ(ierr);
                        ind++;

                    }
                    //Extracellular terms
                    comp = Nc - 1;
                    //Electrodiffusion contributions
                    Fc0x = 0;
                    Fc1x = 0;
                    Fph0x = 0;
                    Fph1x = 0;
                    Fc0y = 0;
                    Fc1y = 0;
                    Fph0y = 0;
                    Fph1y = 0;
                    if (x < Nx - 1) {
                        Fc0x = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]/dx*dt/dx;
                        Fph0x = z_charge[ion]*Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2]*
                                (cp[c_index(x,y,0,comp,ion,Nx,0,0)]+cp[c_index(x+1,y,0,comp,ion,Nx,0,0)])/2/dx*dt/dx;
                    }
                    if (x > 0) {
                        Fc1x = Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]/dx*dt/dx;
                        Fph1x = z_charge[ion]*Dcs[c_index(x-1,y,0,comp,ion,Nx,0,0)*2]*
                                (cp[c_index(x-1,y,0,comp,ion,Nx,0,0)]+cp[c_index(x,y,0,comp,ion,Nx,0,0)])/2/dx*dt/dx;
                    }
                    if (y < Ny - 1) {
                        Fc0y = Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]/dy*dt/dy;
                        Fph0y = z_charge[ion]*Dcs[c_index(x,y,0,comp,ion,Nx,0,0)*2+1]*
                                (cp[c_index(x,y,0,comp,ion,Nx,0,0)]+cp[c_index(x,y+1,0,comp,ion,Nx,0,0)])/2/dy*dt/dy;
                    }
                    if (y > 0) {
                        Fc1y = Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]/dy*dt/dy;
                        Fph1y = z_charge[ion]*Dcs[c_index(x,y-1,0,comp,ion,Nx,0,0)*2+1]*
                                (cp[c_index(x,y-1,0,comp,ion,Nx,0,0)]+cp[c_index(x,y,0,comp,ion,Nx,0,0)])/2/dy*dt/dy;
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
                    Ac -= Ftmpx*(cp[c_index(x,y,0,Nc-1,ion,Nx,0,0)]+cbath[ion])/(2*c[c_index(x,y,0,
                                                                                             Nc-1,ion,
                                                                                             Nx,0,0)])*dt;
                    Aphi -= Ftmpx*(cp[c_index(x,y,0,Nc-1,ion,Nx,0,0)]+cbath[ion])*z_charge[ion]/2*dt;

                    Avolt -= z_charge[ion]*Ftmpx*(cp[c_index(x,y,0,Nc-1,ion,Nx,0,0)]+cbath[ion])/(2*c[c_index(
                            x,y,0,Nc-1,
                            ion,Nx,0,0)])*dt;

                    //Insert extracell to extracell parts
                    // c with c
                    ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,Nc-1,Nx,0,0),Ind_1(x,y,0,ion,Nc-1,Nx,0,0),Ac,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    // c with phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,0,ion,Nc-1,Nx,0,0),Ind_1(x,y,0,Ni,Nc-1,Nx,0,0),Aphi,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;

                    //phi with c (voltage eqn)
                    ierr = MatSetValue(Jac,Ind_1(x,y,0,Ni,Nc-1,Nx,0,0),Ind_1(x,y,0,ion,Nc-1,Nx,0,0),Avolt,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                //Derivative of charge-capacitance
                for (comp = 0; comp < Nc - 1; comp++) {
                    Avolt = cm[comp] + Fphph0x[comp] + Fphph1x[comp] + Fphph0y[comp] + Fphph1y[comp];
                    AvoltN = -cm[comp];
                    for (ion = 0; ion < Ni; ion++) {
                        Avolt += z_charge[ion]*flux->dfdphim[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                        AvoltN -= z_charge[ion]*flux->dfdphim[c_index(x,y,0,comp,ion,Nx,0,0)]*dt;
                    }

                    //Intra-phi with Intra-phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),Avolt,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                    //Intra-phi with extra-phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,Ni,Nc-1,Nx,0,0),AvoltN,INSERT_VALUES);
                    CHKERRQ(ierr);
                    ind++;
                }
                //Extracellular terms
                comp = Nc - 1;
                AvoltN = 0;

                for (int k = 0; k < Nc - 1; k++) {
                    AvoltN += cm[k];
                    Avolt = -cm[k];
                    for (ion = 0; ion < Ni; ion++) {
                        Avolt -= z_charge[ion]*flux->dfdphim[c_index(x,y,0,k,ion,Nx,0,0)]*dt;
                        AvoltN += z_charge[ion]*flux->dfdphim[c_index(x,y,0,k,ion,Nx,0,0)]*dt;
                    }
                    //Extra-phi with Intra-phi
                    ierr = MatSetValue(Jac,Ind_1(x,y,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,Ni,k,Nx,0,0),Avolt,INSERT_VALUES);
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
                ierr = MatSetValue(Jac,Ind_1(x,y,0,Ni,comp,Nx,0,0),Ind_1(x,y,0,Ni,comp,Nx,0,0),AvoltN,INSERT_VALUES);
                CHKERRQ(ierr);
                ind++;

            }
        }
    }
//    printf("Iter:%d, inserts:%d\n",iter,ind);
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
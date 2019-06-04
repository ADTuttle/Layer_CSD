#include "functions.h"
#include "constants.h"
#include <string.h>
#include <math.h>


void print_all(struct AppCtx *user)
{
    PetscReal *Dcs = user->Dcs;
    PetscReal *Dcb = user->Dcb;
    struct ConstVars *con_vars = user->con_vars;
    struct FluxData *flux = user->flux;
    struct GateType *gvars = user->gate_vars;
    struct SimState *state_vars = user->state_vars;
    struct Solver *slvr = user->slvr;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt Nz = user->Nz;
    printf("ConstVars:\n");
    printf("%f,%f,%f\n",1e6*con_vars->pNaKCl[0],1e6*con_vars->Imax[0],1e6*con_vars->pNaLeak[0]);
    printf("%f,%f\n",1e6*con_vars->Imaxg[0],1e6*con_vars->pNaLeakg[0]);
    printf("%f,%f,%f\n",con_vars->ao[0],con_vars->ao[1],con_vars->ao[2]);
    printf("%f,%f,%f\n",con_vars->zo[0],con_vars->zo[1],con_vars->zo[2]);
    printf("%f,%f,%f\n",con_vars->kappa,con_vars->zeta1[0],con_vars->zeta1[1]);
    printf("%d,%f,%f\n",con_vars->S,1e6*con_vars->zetaalpha[0],1e6*con_vars->zetaalpha[1]);

    //Diffusion in each compartment
    //Has x and y components
    //x will be saved at even positions (0,2,4,...)
    //y at odd (1,3,5,...)
    //still use c_index(x,y,comp,ion,Nx), but with ind*2 or ind*2+1

    for(PetscInt ion=0;ion<Ni;ion++){
        for(PetscInt comp=0;comp<Nc;comp++){
            printf("Dcs: Ion %d, Comp %d ",ion,comp);
            printf("Dcs x: %.10e, Dcs y: %.10e,Dcs z: %.10e\n",Dcs[c_index(0,0,0,comp,ion,Nx,Ny,Nz)*3],
                   Dcs[c_index(0,0,0,comp,ion,Nx,Ny,Nz)*3+1],Dcs[c_index(0,0,0,comp,ion,Nx,Ny,Nz)*3+2]);
        }
    }
    printf("\n");

    //Bath diffusion

    for(PetscInt ion=0;ion<Ni;ion++){
        for(PetscInt comp=0;comp<Nc;comp++){
            printf("Dcb: Ion %d, Comp %d ",ion,comp);
            printf("Dcb x: %.10e, Dcb y: %.10e,Dcb z: %.10e\n",Dcb[c_index(0,0,0,comp,ion,Nx,Ny,Nz)*3],
                   Dcb[c_index(0,0,0,comp,ion,Nx,Ny,Nz)*3+1],Dcb[c_index(0,0,0,comp,ion,Nx,Ny,Nz)*3+2]);
        }
    }

    PetscInt x=0;PetscInt y=0;
    printf("\n");
    for(PetscInt ion=0;ion<Ni;ion++){
        for(PetscInt comp=0;comp<Nc;comp++){
            printf("Ion: %d, Comp %d, C: %f\n",ion,comp,state_vars->c[c_index(0,0,0,comp,ion,Nx,Ny,Nz)]);
        }
    }
    for(PetscInt comp=0;comp<Nc;comp++){
        printf("Comp %d, Phi: %f\n",comp,state_vars->phi[phi_index(0,0,0,comp,Nx,Ny,Nz)]*RTFC);
    }
    for(PetscInt comp=0;comp<Nc-1;comp++){
        printf("Comp %d, alpha: %.10e\n",comp,state_vars->alpha[al_index(0,0,0,comp,Nx,Ny,Nz)]);
    }
    printf("Gvars:\n");
    printf("NaT :%.10e,%.10e,%.10e\n",gvars->mNaT[0],gvars->hNaT[0],gvars->gNaT[0]);
    printf("NaP :%.10e,%.10e,%.10e\n",gvars->mNaP[0],gvars->hNaP[0],gvars->gNaP[0]);
    printf("KDR :%.10e,%.10e\n",gvars->mKDR[0],gvars->gKDR[0]);
    printf("KA :%.10e,%.10e,%.10e\n",gvars->mKA[0],gvars->hKA[0],gvars->gKA[0]);
    printf("NMDA :%.10e,%.10e\n",gvars->yNMDA[0],gvars->gNMDA[0]);
    printf("\n");
    //Compute membrane ionic flux relation quantitites

    for(PetscInt ion=0;ion<Ni;ion++){
        for(PetscInt comp=0;comp<Nc;comp++){
            printf("Ion: %d, Comp %d\n",ion,comp);
            printf("Flux: %.10e, dfdci: %.10e, dfdce: %.10e, dfdphim: %.10e\n",
                   flux->mflux[c_index(0,0,0,comp,ion,Nx,Ny,Nz)],
                   flux->dfdci[c_index(0,0,0,comp,ion,Nx,Ny,Nz)],
                   flux->dfdce[c_index(0,0,0,comp,ion,Nx,Ny,Nz)],flux->dfdphim[c_index(0,0,0,comp,ion,Nx,Ny,Nz)]);
        }
    }
    printf("\n");
    //Compute membrane water flow related quantities
    for(PetscInt comp=0;comp<Nc-1;comp++){
        printf("Comp: %d\n",comp);
        printf("wFlux: %.10e,%.10e,%.10e\n",flux->wflow[al_index(x,y,0,comp,Nx,Ny,Nz)],
               flux->dwdpi[al_index(x,y,0,comp,Nx,Ny,Nz)],flux->dwdal[al_index(x,y,0,comp,Nx,Ny,Nz)]);
    }
    printf("\n");
    // VecView(slvr->Res,PETSC_VIEWER_STDOUT_SELF);

    // VecView(slvr->Q,PETSC_VIEWER_STDOUT_SELF);

    return;

}

const char* getfield(char* line, int num)
{
    const char* tok;
    for (tok = strtok(line, ",");
         tok && *tok;
         tok = strtok(NULL, ",\n"))
    {
        if (!--num)
            return tok;
    }
    return NULL;
}

void find_print(int rowC, int colC, double valC, int iter)
{

    if(iter!=-1){
        return;
    }
//    printf("%d,%d,%.10e\n",rowC,colC,valC);
//    return;
    if(rowC>70){
        return;
    }

    int row,col;
    double val;

    //Open Julia Mat file.
    FILE *fp;
    fp = fopen("../../../Julia_work/2d_modules/matrix.txt","r");
    if(fp==NULL)
    {
        fprintf(stderr, "File not found....\n");
        exit(EXIT_FAILURE); /* indicate failure.*/
    }

    //Begin searching for the right row and col.
    char line[1024];
    while (fgets(line, 1024, fp))
    {
        char* tmp = strdup(line);
        row = -1; col=-1; val=0;

        row = atoi(getfield(tmp,1));


        tmp = strdup(line);
        col = atoi(getfield(tmp,2));

        if(row==rowC && col==colC)
        {
            tmp = strdup(line);
            val = atof(getfield(tmp,3));
            if(val!=0) {
                if (fabs(val - valC)/fabs(val) > 5e-16) {
                    printf("Row: %d, Col: %d, J: %f, C: %f, Diff: %fe-16\n", row, col, val, valC,
                           1e16 * (val - valC) / val);
                }
            }
            /*
            if(val==0){
                    printf("Row: %d, Col: %d, J: %f, C: %f, Diff: %fe-16\n",row,col,val,valC,1e16*(val-valC));

            }
            else{
                    printf("Row: %d, Col: %d, J: %f, C: %f, Diff: %fe-16,Abs: %f\n",row,col,val,valC,1e16*(val-valC)/val,1e16*(val-valC));
            }
             */
            return;
        }
        free(tmp);
    }


    fclose(fp);
    return;

}

void compare_res(double *Res, int iter)
{
    if(iter!=-1){
        return;
    }

    FILE *fp;
    fp = fopen("../../../Julia_work/2d_modules/Res.txt","r");
    if(fp==NULL)
    {
        fprintf(stderr, "File not found....\n");
        exit(EXIT_FAILURE); /* indicate failure.*/
    }

    //Begin searching for the right row and col.
    char line[1024];
    int row=0;
    double val;
    while (fgets(line, 1024, fp))
    {
        char* tmp = strdup(line);
        val=0;
        val = atof(getfield(tmp,1));

        if(val==0){
            printf("Row: %d, J: %.10e, C: %.10e, diff: %.10e\n",row,val,Res[row],val-Res[row]);
        }else{
            printf("Row: %d, J: %.10e, C: %.10e,abs: %.10e, diff: %.10e\n",row,val,Res[row],val-Res[row],(val-Res[row])/val);
        }
        row++;
        if(row>70){
            break;
        }
        free(tmp);
    }


    fclose(fp);
    return;
}

void write_data(FILE **fp,struct AppCtx*user,PetscInt numrecords,int start)
{
    if(Profiling_on) {
        PetscLogEventBegin(event[8], 0, 0, 0, 0);
    }
    struct SimState *state_vars = user->state_vars;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt Nz = user->Nz;
    if(!save_one_var) {
        if (start) {
            for(int z=0;z<Nz;z++){
                char name[16];
                if(z<10){
                    sprintf(name,"data_csd_0%d.txt",z);
                }else{
                    sprintf(name,"data_csd_%d.txt",z);
                }
                fp[z] = fopen(name,"w");
                fprintf(fp[z],"%d,%d,%d,%d,%d\n",Nx,Ny,numrecords,Nc,Ni);
            }
            write_data(fp, user, numrecords,0);
        } else{
            int ion,comp,x,y,z;
            for(z = 0; z < Nz; z++){
                for(ion = 0; ion < Ni; ion++){
                    for(comp = 0; comp < Nc; comp++){
                        for(y = 0; y < Ny; y++){
                            for(x = 0; x < Nx; x++){
                                if(x == Nx-1 & y == Ny-1){
                                    fprintf(fp[z],"%.10e\n",state_vars->c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]);
                                }else{
                                    fprintf(fp[z],"%.10e,",state_vars->c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]);
                                }
                            }
                        }
                    }
                }
                for(comp = 0; comp < Nc; comp++){
                    for(y = 0; y < Ny; y++){
                        for(x = 0; x < Nx; x++){
                            if(x == Nx-1 & y == Ny-1){
                                fprintf(fp[z],"%.10e\n",state_vars->phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]*RTFC);
                            }else{
                                fprintf(fp[z],"%.10e,",state_vars->phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]*RTFC);
                            }
                        }
                    }
                }
                for(comp = 0; comp < Nc-1; comp++){
                    for(y = 0; y < Ny; y++){
                        for(x = 0; x < Nx; x++){
                            if(x == Nx-1 & y == Ny-1){
                                fprintf(fp[z],"%.10e\n",state_vars->alpha[al_index(x,y,z,comp,Nx,Ny,Nz)]);
                            }else{
                                fprintf(fp[z],"%.10e,",state_vars->alpha[al_index(x,y,z,comp,Nx,Ny,Nz)]);
                            }
                        }
                    }
                }
            }
        }
    } else{
        if (start) {
            for(int z=0;z<Nz;z++){
                char name[16];
                if(z < 10){
                    sprintf(name,"data_csd_0%d.txt",z);
                }else{
                    sprintf(name,"data_csd_%d.txt",z);
                }
                fp[z] = fopen(name,"w");
                fprintf(fp[z],"%d,%d,%d,%d,%d\n",Nx,Ny,(int) floor(numrecords),0,0);
            }
            write_data(fp, user,numrecords, 0);
        } else {
            int ion, comp, x, y,z;
            comp = 0;
            for(z=0;z<Nz;z++){
                for(y = 0; y < Ny; y++){
                    for(x = 0; x < Nx; x++){
                        if(x == Nx-1 & y == Ny-1){
//                        fprintf(fp[z], "%.10e\n", state_vars->phi[phi_index(x, y, Nc-1,Nx)] * RTFC);
                        fprintf(fp[z],"%.10e\n",(state_vars->phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-
                                                 state_vars->phi[phi_index(x,y,z,Nc-1,Nx,Ny,Nz)])*RTFC);
                        }else{
//                        fprintf(fp[z], "%.10e,", state_vars->phi[phi_index(x, y,z, Nc-1,Nx,Ny)] * RTFC);
                            fprintf(fp[z],"%.10e,",(state_vars->phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-
                                                    state_vars->phi[phi_index(x,y,z,Nc-1,Nx,Ny,Nz)])*RTFC);
                        }
                    }
                }
            }
        }
    }
    if(Profiling_on) {
        PetscLogEventEnd(event[8], 0, 0, 0, 0);
    }
}
void write_point(FILE *fp,struct AppCtx* user,PetscReal t,PetscInt x,PetscInt y)
{
    if(Profiling_on) {
        PetscLogEventBegin(event[8], 0, 0, 0, 0);
    }
    struct SimState *state_vars = user->state_vars;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt Nz = user->Nz;
    int ion, comp;
    comp=0;
    fprintf(fp,"%f,",t);
    for(PetscInt z=0;z<Nz;z++){
        fprintf(fp,"%.10e,",(state_vars->phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-
                             state_vars->phi[phi_index(x,y,z,Nc-1,Nx,Ny,Nz)])*RTFC);
    }
    fprintf(fp,"\n");

}
void save_timestep(FILE *fp,struct AppCtx*user,PetscInt numrecords,int start)
{
    if(Profiling_on) {
        PetscLogEventBegin(event[8], 0, 0, 0, 0);
    }
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    if (start) {
        fprintf(fp, "%d,%d,%d,%d,%d\n", Nx, Ny, numrecords, 0, 0);
    } else {
        int x, y;

        for (y = 0; y < Ny; y++) {
            for (x = 0; x < Nx; x++) {
                if (x == Nx - 1 & y == Ny - 1) {
                    fprintf(fp, "%f\n", user->dt_space[xy_index(x,y,0,Nx,0,0)]);
                } else {
                    fprintf(fp, "%f,", user->dt_space[xy_index(x,y,0,Nx,0,0)]);
                }
            }
        }
    }



    if(Profiling_on) {
        PetscLogEventEnd(event[8], 0, 0, 0, 0);
    }
}

void record_measurements(FILE **fp_measures,struct AppCtx *user,PetscInt count,PetscInt numrecords,int start){
    if(start){
        if(start_at_steady) {
            fp_measures[0] = fopen("flux_csd.txt", "w");
            fp_measures[1] = fopen("grad_field.txt", "w");
            fp_measures[2] = fopen("measures.txt", "w");

            measure_flux(fp_measures[0],user,numrecords,start);
            velocity_field(fp_measures[1],user,numrecords,start);
//            calculate_measures(fp_measures[2],user,numrecords,start);
            calculate_energy(fp_measures[2],user,numrecords,start);
        }else{
            fp_measures[0] = fopen("flux_csd.txt", "a");
            fp_measures[1] = fopen("grad_field.txt", "a");
            fp_measures[2] = fopen("measures.txt", "a");
        }
    } else{
        measure_flux(fp_measures[0],user,numrecords,start);
        velocity_field(fp_measures[1],user,numrecords,start);
        //            calculate_measures(fp_measures[2],user,numrecords,start);
        calculate_energy(fp_measures[2],user,numrecords,start);

        if(count%100==0) {
            draw_csd(user);
        }
    }
}
void measure_flux(FILE *fp, struct AppCtx* user,PetscInt numrecords,int start)
{
    struct SimState *state_vars= user->state_vars;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt Nz = user->Nz;
    if (start) {
        fprintf(fp, "%d,%d,%d,%d\n", Nz, numrecords, Nc, Ni);
        measure_flux(fp, user, numrecords, 0);
    } else{
        //compute diffusion coefficients
        diff_coef(user->Dcs,state_vars->alpha,1,user);
        //Bath diffusion
        diff_coef(user->Dcb,state_vars->alpha,Batheps,user);

        PetscReal *c = state_vars->c;
        PetscReal *phi = state_vars->phi;
        PetscReal *al = state_vars->alpha;
        PetscReal *Dcs = user->Dcs;
        PetscReal *Dcb = user->Dcb;
        PetscReal dx = user->dx;
        PetscReal dy = user->dy;
        PetscReal dz = user->dz;

        PetscInt x,y,z,comp,ion;
        PetscReal Rcvz,RcvzTop;

        x = (PetscInt) floor(((double) Nx)/2);
        y = 3;//(PetscInt)floor(((double)Ny)/2);
        ionmflux(user);

        for(z = 0; z < Nz; z++){
            for(ion = 0; ion < Ni; ion++){
                for(comp = 0; comp < Nc; comp++){

                    Rcvz = 0;
                    RcvzTop = 0;
                    //Top Bottom difference
                    if(z > 0){
                        Rcvz = Dcs[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)*3+2]*
                                (c[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)]+c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])/2;
                        Rcvz = Rcvz*(log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])
                                -log(c[c_index(x,y,z-1,comp,ion,Nx,Ny,Nz)])+z_charge[ion]
                                *(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-phi[phi_index(x,y,z-1,comp,Nx,Ny,Nz)]))/dz/dz;
                    }
                    //Next upward difference
                    if(z < Nz-1){
                        RcvzTop = Dcs[c_index(x,y,z,comp,ion,Nx,Ny,Nz)*3+2]*
                                (c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]+c[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])/2;
                        RcvzTop = RcvzTop*(log(c[c_index(x,y,z+1,comp,ion,Nx,Ny,Nz)])
                                -log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])+z_charge[ion]
                                *(phi[phi_index(x,y,z+1,comp,Nx,Ny,Nz)]-phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]))/dz/dz;
                    }
                    fprintf(fp,"%.10e,%.10e,",Rcvz-RcvzTop,user->flux->mflux[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]);
                }
            }
            fprintf(fp,"\n");
        }
    }
}
void init_events(struct AppCtx *user)
{
    PetscLogDefaultBegin();


    PetscClassId id;
    PetscClassIdRegister("CSD",&id);

    PetscLogEventRegister("Jacobian",id,&event[0]);
    PetscLogEventRegister("Residual",id,&event[1]);
    PetscLogEventRegister("Extract Subarray",id,&event[2]);
    PetscLogEventRegister("Restore Subarray",id,&event[3]);
    PetscLogEventRegister("Gating Variables",id,&event[4]);
    PetscLogEventRegister("Ion Flux",id,&event[5]);
    PetscLogEventRegister("Water Flux",id,&event[6]);
    PetscLogEventRegister("Volume Update",id,&event[7]);
    PetscLogEventRegister("Write to File",id,&event[8]);
    PetscLogEventRegister("Predict Jacobian",id,&event[9]);
    PetscLogEventRegister("Predict Residual",id,&event[10]);
    PetscLogEventRegister("Predict Solve",id,&event[11]);
    PetscLogEventRegister("Load Grid",id,&event[12]);
    PetscLogEventRegister("Unload Grid",id,&event[13]);

    //Deactivate Petsc tracking
    PetscLogEvent deactivate;

    PetscLogEventDeactivateClass(MAT_CLASSID);
    PetscLogEventDeactivateClass(KSP_CLASSID); /* includes PC and KSP */
    PetscLogEventDeactivateClass(VEC_CLASSID);
//    PetscLogEventDeactivateClass(SNES_CLASSID);

    //Some events are leftover somehow

    PetscLogEventGetId("PCApply",&deactivate);
    PetscLogEventDeactivate(deactivate);
    PetscLogEventGetId("VecSet",&deactivate);
    PetscLogEventDeactivate(deactivate);
    PetscLogEventGetId("MatAssemblyBegin",&deactivate);
    PetscLogEventDeactivate(deactivate);
    PetscLogEventGetId("MatAssemblyEnd",&deactivate);
    PetscLogEventDeactivate(deactivate);
    PetscLogEventGetId("SNESLineSearch",&deactivate);
    PetscLogEventDeactivate(deactivate);
    PetscLogEventGetId("PCSetUp",&deactivate);
    PetscLogEventDeactivate(deactivate);

    PetscLogEventGetId("SNESFunctionEval",&deactivate);
    PetscLogEventDeactivate(deactivate);
    PetscLogEventGetId("SNESJacobianEval",&deactivate);
    PetscLogEventDeactivate(deactivate);

}

void save_file(struct AppCtx *user){

    FILE *fp;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt Nz = user->Nz;
    PetscInt x,y,ion,comp;

    struct SimState *state_vars = user->state_vars;
    struct GateType *gate_vars = user->gate_vars;

    fp = fopen("save_csd.txt","w");


    //Header
    fprintf(fp, "%d,%d,%d,%d,%d\n", Nx, Ny, 1, Nc, Ni);

    //Macro-variables
    for (ion = 0; ion < Ni; ion++) {
        for (comp = 0; comp < Nc; comp++) {
            for (y = 0; y < Ny; y++) {
                for (x = 0; x < Nx; x++) {
                    if (x == Nx - 1 & y == Ny - 1) {
                        fprintf(fp, "%.20e\n", state_vars->c[c_index(x,y,0,comp,ion,Nx,Ny,Nz)]);
                    } else {
                        fprintf(fp, "%.20e,", state_vars->c[c_index(x,y,0,comp,ion,Nx,Ny,Nz)]);
                    }
                }
            }
        }
    }
    for (comp = 0; comp < Nc; comp++) {
        for (y = 0; y < Ny; y++) {
            for (x = 0; x < Nx; x++) {
                if (x == Nx - 1 & y == Ny - 1) {
                    fprintf(fp, "%.20e\n", state_vars->phi[phi_index(x,y,0,comp,Nx,Ny,Nz)]);
                } else {
                    fprintf(fp, "%.20e,", state_vars->phi[phi_index(x,y,0,comp,Nx,Ny,Nz)]);
                }
            }
        }
    }
    for (comp = 0; comp < Nc - 1; comp++) {
        for (y = 0; y < Ny; y++) {
            for (x = 0; x < Nx; x++) {
                if (x == Nx - 1 & y == Ny - 1) {
                    fprintf(fp, "%.20e\n", state_vars->alpha[al_index(x,y,0,comp,Nx,Ny,Nz)]);
                } else {
                    fprintf(fp, "%.20e,", state_vars->alpha[al_index(x,y,0,comp,Nx,Ny,Nz)]);
                }
            }
        }
    }

    //Gating Variables
    for(y=0;y<Ny;y++){
        for(x=0;x<Nx;x++){
            if (x == Nx - 1 && y == Ny - 1) {
                fprintf(fp, "%.20e\n", gate_vars->mNaT[xy_index(x,y,0,Nx,Ny,Nz)]);
            } else {
                fprintf(fp, "%.20e,", gate_vars->mNaT[xy_index(x,y,0,Nx,Ny,Nz)]);
            }
        }
    }
    for(y=0;y<Ny;y++){
        for(x=0;x<Nx;x++){
            if (x == Nx - 1 && y == Ny - 1) {
                fprintf(fp, "%.20e\n", gate_vars->hNaT[xy_index(x,y,0,Nx,Ny,Nz)]);
            } else {
                fprintf(fp, "%.20e,", gate_vars->hNaT[xy_index(x,y,0,Nx,Ny,Nz)]);
            }
        }
    }
    for(y=0;y<Ny;y++){
        for(x=0;x<Nx;x++){
            if (x == Nx - 1 && y == Ny - 1) {
                fprintf(fp, "%.20e\n", gate_vars->mNaP[xy_index(x,y,0,Nx,Ny,Nz)]);
            } else {
                fprintf(fp, "%.20e,", gate_vars->mNaP[xy_index(x,y,0,Nx,Ny,Nz)]);
            }
        }
    }
    for(y=0;y<Ny;y++){
        for(x=0;x<Nx;x++){
            if (x == Nx - 1 && y == Ny - 1) {
                fprintf(fp, "%.20e\n", gate_vars->hNaP[xy_index(x,y,0,Nx,Ny,Nz)]);
            } else {
                fprintf(fp, "%.20e,", gate_vars->hNaP[xy_index(x,y,0,Nx,Ny,Nz)]);
            }
        }
    }
    for(y=0;y<Ny;y++){
        for(x=0;x<Nx;x++){
            if (x == Nx - 1 && y == Ny - 1) {
                fprintf(fp, "%.20e\n", gate_vars->mKDR[xy_index(x,y,0,Nx,Ny,Nz)]);
            } else {
                fprintf(fp, "%.20e,", gate_vars->mKDR[xy_index(x,y,0,Nx,Ny,Nz)]);
            }
        }
    }
    for(y=0;y<Ny;y++){
        for(x=0;x<Nx;x++){
            if (x == Nx - 1 && y == Ny - 1) {
                fprintf(fp, "%.20e\n", gate_vars->mKA[xy_index(x,y,0,Nx,Ny,Nz)]);
            } else {
                fprintf(fp, "%.20e,", gate_vars->mKA[xy_index(x,y,0,Nx,Ny,Nz)]);
            }
        }
    }
    for(y=0;y<Ny;y++){
        for(x=0;x<Nx;x++){
            if (x == Nx - 1 && y == Ny - 1) {
                fprintf(fp, "%.20e\n", gate_vars->hKA[xy_index(x,y,0,Nx,Ny,Nz)]);
            } else {
                fprintf(fp, "%.20e,", gate_vars->hKA[xy_index(x,y,0,Nx,Ny,Nz)]);
            }
        }
    }
    for(y=0;y<Ny;y++){
        for(x=0;x<Nx;x++){
            if (x == Nx - 1 && y == Ny - 1) {
                fprintf(fp, "%.20e\n", gate_vars->yNMDA[xy_index(x,y,0,Nx,Ny,Nz)]);
            } else {
                fprintf(fp, "%.20e,", gate_vars->yNMDA[xy_index(x,y,0,Nx,Ny,Nz)]);
            }
        }
    }
    for(y=0;y<Ny;y++){
        for(x=0;x<Nx;x++){
            if (x == Nx - 1 && y == Ny - 1) {
                fprintf(fp, "%.20e\n", gate_vars->zNMDA[xy_index(x,y,0,Nx,Ny,Nz)]);
            } else {
                fprintf(fp, "%.20e,", gate_vars->zNMDA[xy_index(x,y,0,Nx,Ny,Nz)]);
            }
        }
    }
    for(y=0;y<Ny;y++){
        for(x=0;x<Nx;x++){
            if (x == Nx - 1 && y == Ny - 1) {
                fprintf(fp, "%.20e\n", gate_vars->dNMDA[xy_index(x,y,0,Nx,Ny,Nz)]);
            } else {
                fprintf(fp, "%.20e,", gate_vars->dNMDA[xy_index(x,y,0,Nx,Ny,Nz)]);
            }
        }
    }


    fclose(fp);

}
void read_file(struct AppCtx *user)
{
    //Read data_csd.txt and copy it into the current state.
    int Nx,Ny,Nz,numrecords;
    struct SimState *state_vars = user->state_vars;
    struct GateType *gate_vars = user->gate_vars;
    PetscInt x,y,z,comp,ion;

    Nz=user->Nz;


    char *line = (char *) malloc(sizeof(char) * 1024 * 1024);//[1024];
    char *tmp;

    FILE *fp;
    fp = fopen("save_csd.txt","r");


    if(fp==NULL) {
        fp = fopen("data_csd.txt", "r");
        if (fp == NULL) {
            fprintf(stderr, "File not found....\n");
            exit(EXIT_FAILURE); /* indicate failure.*/
        }

        //Read top file details
        fgets(line, 1024 * 1024, fp);

        tmp = strdup(line);

        Nx = atoi(getfield(tmp, 1));

        tmp = strdup(line);
        Ny = atoi(getfield(tmp, 2));

        tmp = strdup(line);
        numrecords = atoi(getfield(tmp, 3));


        printf("%d,%d,%d\n", Nx, Ny, numrecords);

        //Zero out for safety
        for (x = 0; x < Nx; x++) {
            for (y = 0; y < Ny; y++) {
                for (comp = 0; comp < Nc; comp++) {
                    for (ion = 0; ion < Ni; ion++) {
                        state_vars->c[c_index(x,y,0,comp,ion,Nx,Ny,Nz)] = 0;
                    }
                    state_vars->phi[phi_index(x,y,0,comp,Nx,Ny,Nz)] = 0;
                }
                for (comp = 0; comp < Nc - 1; comp++) {
                    state_vars->alpha[al_index(x,y,0,comp,Nx,Ny,Nz)] = 0;
                }
            }
        }


        //Get to the last recorded value
        for (int count = 0; count < ((Ni + 2) * Nc - 1) * (numrecords - 1); count++) {
            fgets(line, 1024 * 1024, fp);
        }

        for (ion = 0; ion < Ni; ion++) {
            for (comp = 0; comp < Nc; comp++) {
                fgets(line, 1024 * 1024, fp);
                for (y = 0; y < Ny; y++) {
                    for (x = 0; x < Nx; x++) {
                        tmp = strdup(line);

                        state_vars->c[c_index(x,y,0,comp,ion,Nx,Ny,Nz)] = atof(getfield(tmp,
                                                                                       xy_index(x,y,0,Nx,Ny,Nz)+1));

                    }
                }
            }
        }
        for (comp = 0; comp < Nc; comp++) {
            fgets(line, 1024 * 1024, fp);
            for (y = 0; y < Ny; y++) {
                for (x = 0; x < Nx; x++) {
                    tmp = strdup(line);

                    state_vars->phi[phi_index(x,y,0,comp,Nx,Ny,Nz)] = atof(getfield(tmp,xy_index(x,y,0,Nx,Ny,Nz)+1));

                }
            }
        }
        for (comp = 0; comp < Nc - 1; comp++) {
            fgets(line, 1024 * 1024, fp);
            for (y = 0; y < Ny; y++) {
                for (x = 0; x < Nx; x++) {
                    tmp = strdup(line);

                    state_vars->alpha[al_index(x,y,0,comp,Nx,Ny,Nz)] = atof(getfield(tmp,xy_index(x,y,0,Nx,Ny,Nz)+1));

                }
            }
        }
        for(z=0;z<Nz;z++){
            for(y=0;y<Ny;y++){
                for(x=0;x<Nx;x++){
                    for(comp=0;comp<Nc;comp++){
                        for(ion=0;ion<Ni;ion++){
                            state_vars->c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)] = state_vars->c[c_index(x,y,0,comp,ion,Nx,Ny,Nz)];
                        }
                        state_vars->phi[phi_index(x,y,z,comp,Nx,Ny,Nz)] = state_vars->phi[phi_index(x,y,0,comp,Nx,Ny,Nz)];
                        if(comp<Nc-1){
                            state_vars->alpha[al_index(x,y,z,comp,Nx,Ny,Nz)] = state_vars->alpha[al_index(x,y,0,comp,Nx,Ny,Nz)];
                        }
                    }
                }
            }
        }
        gatevars_update(user->gate_vars,user->gate_vars,state_vars,0,user,1);
        gatevars_update(user->gate_vars_past,user->gate_vars_past,state_vars,0,user,1);
    } else{
        //Read top file details
        fgets(line, 1024 * 1024, fp);

        tmp = strdup(line);

        Nx = atoi(getfield(tmp, 1));

        tmp = strdup(line);
        Ny = atoi(getfield(tmp, 2));

        tmp = strdup(line);
        numrecords = atoi(getfield(tmp, 3));


        printf("%d,%d,%d\n", Nx, Ny, numrecords);

        //Zero out for safety
        for (x = 0; x < Nx; x++) {
            for (y = 0; y < Ny; y++) {
                for (comp = 0; comp < Nc; comp++) {
                    for (ion = 0; ion < Ni; ion++) {
                        state_vars->c[c_index(x,y,0,comp,ion,Nx,Ny,Nz)] = 0;
                    }
                    state_vars->phi[phi_index(x,y,0,comp,Nx,Ny,Nz)] = 0;
                }
                for (comp = 0; comp < Nc - 1; comp++) {
                    state_vars->alpha[al_index(x,y,0,comp,Nx,Ny,Nz)] = 0;
                }
            }
        }


        for (ion = 0; ion < Ni; ion++) {
            for (comp = 0; comp < Nc; comp++) {
                fgets(line, 1024 * 1024, fp);
                for (y = 0; y < Ny; y++) {
                    for (x = 0; x < Nx; x++) {
                        tmp = strdup(line);
                        state_vars->c[c_index(x,y,0,comp,ion,Nx,Ny,Nz)] = atof(getfield(tmp,
                                                                                       xy_index(x,y,0,Nx,Ny,Nz)+1));
                    }
                }
            }
        }
        for (comp = 0; comp < Nc; comp++) {
            fgets(line, 1024 * 1024, fp);
            for (y = 0; y < Ny; y++) {
                for (x = 0; x < Nx; x++) {
                    tmp = strdup(line);
                    state_vars->phi[phi_index(x,y,0,comp,Nx,Ny,Nz)] = atof(getfield(tmp,xy_index(x,y,0,Nx,Ny,Nz)+1));
                }
            }
        }
        for (comp = 0; comp < Nc - 1; comp++) {
            fgets(line, 1024 * 1024, fp);
            for (y = 0; y < Ny; y++) {
                for (x = 0; x < Nx; x++) {
                    tmp = strdup(line);
                    state_vars->alpha[al_index(x,y,0,comp,Nx,Ny,Nz)] = atof(getfield(tmp,xy_index(x,y,0,Nx,Ny,Nz)+1));
                }
            }
        }
        //Gating variables
        fgets(line, 1024 * 1024, fp);
        for (y=0;y<Ny;y++){
            for(x=0;x<Nx;x++){
                tmp = strdup(line);
                gate_vars->mNaT[xy_index(x,y,0,Nx,Ny,Nz)] = atof(getfield(tmp,xy_index(x,y,0,Nx,Ny,Nz)+1));
            }
        }
        fgets(line, 1024 * 1024, fp);
        for (y=0;y<Ny;y++){
            for(x=0;x<Nx;x++){
                tmp = strdup(line);
                gate_vars->hNaT[xy_index(x,y,0,Nx,Ny,Nz)] = atof(getfield(tmp,xy_index(x,y,0,Nx,Ny,Nz)+1));
            }
        }
        fgets(line, 1024 * 1024, fp);
        for (y=0;y<Ny;y++){
            for(x=0;x<Nx;x++){
                tmp = strdup(line);
                gate_vars->mNaP[xy_index(x,y,0,Nx,Ny,Nz)] = atof(getfield(tmp,xy_index(x,y,0,Nx,Ny,Nz)+1));
            }
        }
        fgets(line, 1024 * 1024, fp);
        for (y=0;y<Ny;y++){
            for(x=0;x<Nx;x++){
                tmp = strdup(line);
                gate_vars->hNaP[xy_index(x,y,0,Nx,Ny,Nz)] = atof(getfield(tmp,xy_index(x,y,0,Nx,Ny,Nz)+1));
            }
        }
        fgets(line, 1024 * 1024, fp);
        for (y=0;y<Ny;y++){
            for(x=0;x<Nx;x++){
                tmp = strdup(line);
                gate_vars->mKDR[xy_index(x,y,0,Nx,Ny,Nz)] = atof(getfield(tmp,xy_index(x,y,0,Nx,Ny,Nz)+1));
            }
        }
        fgets(line, 1024 * 1024, fp);
        for (y=0;y<Ny;y++){
            for(x=0;x<Nx;x++){
                tmp = strdup(line);
                gate_vars->mKA[xy_index(x,y,0,Nx,Ny,Nz)] = atof(getfield(tmp,xy_index(x,y,0,Nx,Ny,Nz)+1));
            }
        }
        fgets(line, 1024 * 1024, fp);
        for (y=0;y<Ny;y++){
            for(x=0;x<Nx;x++){
                tmp = strdup(line);
                gate_vars->hKA[xy_index(x,y,0,Nx,Ny,Nz)] = atof(getfield(tmp,xy_index(x,y,0,Nx,Ny,Nz)+1));
            }
        }
        fgets(line, 1024 * 1024, fp);
        for (y=0;y<Ny;y++){
            for(x=0;x<Nx;x++){
                tmp = strdup(line);
                gate_vars->yNMDA[xy_index(x,y,0,Nx,Ny,Nz)] = atof(getfield(tmp,xy_index(x,y,0,Nx,Ny,Nz)+1));
            }
        }
        fgets(line, 1024 * 1024, fp);
        for (y=0;y<Ny;y++){
            for(x=0;x<Nx;x++){
                tmp = strdup(line);
                gate_vars->zNMDA[xy_index(x,y,0,Nx,Ny,Nz)] = atof(getfield(tmp,xy_index(x,y,0,Nx,Ny,Nz)+1));
            }
        }
        fgets(line, 1024 * 1024, fp);
        for (y=0;y<Ny;y++){
            for(x=0;x<Nx;x++){
                tmp = strdup(line);
                gate_vars->dNMDA[xy_index(x,y,0,Nx,Ny,Nz)] = atof(getfield(tmp,xy_index(x,y,0,Nx,Ny,Nz)+1));
            }
        }

        //Copy over past vars and calculate g.
        PetscReal Gphi,v;
        PetscReal K_r = 2.3e-6;//34.9e-6;
        PetscReal npow = 1.5; //1.4;
        PetscReal Fglu;
        for(y=0;y<Ny;y++){
            for(x=0;x<Nx;x++){
                gate_vars->gNaT[xy_index(x,y,0,Nx,Ny,Nz)] =
                        pow(gate_vars->mNaT[xy_index(x,y,0,Nx,Ny,Nz)],3)*gate_vars->hNaT[xy_index(x,y,0,Nx,Ny,Nz)];
                gate_vars->gNaP[xy_index(x,y,0,Nx,Ny,Nz)] =
                        pow(gate_vars->mNaP[xy_index(x,y,0,Nx,Ny,Nz)],2)*gate_vars->hNaP[xy_index(x,y,0,Nx,Ny,Nz)];
                gate_vars->gKDR[xy_index(x,y,0,Nx,Ny,Nz)] = pow(gate_vars->mKDR[xy_index(x,y,0,Nx,Ny,Nz)],2);
                gate_vars->gKA[xy_index(x,y,0,Nx,Ny,Nz)] =
                        pow(gate_vars->mKA[xy_index(x,y,0,Nx,Ny,Nz)],2)*gate_vars->hKA[xy_index(x,y,0,Nx,Ny,Nz)];

                v = (state_vars->phi[phi_index(x,y,0,0,Nx,Ny,Nz)]-state_vars->phi[phi_index(x,y,0,Nc-1,Nx,Ny,Nz)])*RTFC;

                Fglu = pow(state_vars->c[c_index(x,y,0,Nc-1,3,Nx,Ny,Nz)],npow)
                       /(pow(state_vars->c[c_index(x,y,0,Nc-1,3,Nx,Ny,Nz)],npow)+pow(K_r,npow));
                Gphi = 1/(1+0.56*exp(-0.062*v));
                gate_vars->gNMDA[xy_index(x,y,0,Nx,Ny,Nz)] = (gate_vars->yNMDA[xy_index(x,y,0,Nx,Ny,Nz)]*Fglu
                                                             +Desensitize[0]*gate_vars->zNMDA[xy_index(x,y,0,Nx,Ny,Nz)]
                                                             +Desensitize[1]*gate_vars->dNMDA[xy_index(x,y,0,Nx,Ny,Nz)])*Gphi;

            }
        }
        for(z=0;z<Nz;z++){
            for(y=0;y<Ny;y++){
                for(x=0;x<Nx;x++){
                    for(comp=0;comp<Nc;comp++){
                        for(ion=0;ion<Ni;ion++){
                            state_vars->c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)] = state_vars->c[c_index(x,y,0,comp,ion,Nx,Ny,Nz)];
                        }
                        state_vars->phi[phi_index(x,y,z,comp,Nx,Ny,Nz)] = state_vars->phi[phi_index(x,y,0,comp,Nx,Ny,Nz)];
                        if(comp<Nc-1){
                            state_vars->alpha[al_index(x,y,z,comp,Nx,Ny,Nz)] = state_vars->alpha[al_index(x,y,0,comp,Nx,Ny,Nz)];
                        }
                    }
                    gate_vars->gNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNaT[xy_index(x,y,0,Nx,Ny,Nz)];
                    gate_vars->hNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hNaT[xy_index(x,y,0,Nx,Ny,Nz)];
                    gate_vars->mNaT[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mNaT[xy_index(x,y,0,Nx,Ny,Nz)];

                    gate_vars->gNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNaP[xy_index(x,y,0,Nx,Ny,Nz)];
                    gate_vars->hNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hNaP[xy_index(x,y,0,Nx,Ny,Nz)];
                    gate_vars->mNaP[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mNaP[xy_index(x,y,0,Nx,Ny,Nz)];

                    gate_vars->gKDR[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gKDR[xy_index(x,y,0,Nx,Ny,Nz)];
                    gate_vars->mKDR[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mKDR[xy_index(x,y,0,Nx,Ny,Nz)];

                    gate_vars->gKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gKA[xy_index(x,y,0,Nx,Ny,Nz)];
                    gate_vars->mKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->mKA[xy_index(x,y,0,Nx,Ny,Nz)];
                    gate_vars->hKA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->hKA[xy_index(x,y,0,Nx,Ny,Nz)];

                    gate_vars->gNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->gNMDA[xy_index(x,y,0,Nx,Ny,Nz)];
                    gate_vars->yNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->yNMDA[xy_index(x,y,0,Nx,Ny,Nz)];
                    gate_vars->zNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->zNMDA[xy_index(x,y,0,Nx,Ny,Nz)];
                    gate_vars->dNMDA[xy_index(x,y,z,Nx,Ny,Nz)] = gate_vars->dNMDA[xy_index(x,y,0,Nx,Ny,Nz)];

                }
            }
        }

        //Copy old gating variables
        //Save the gating variables
        memcpy(user->gate_vars_past->mNaT,user->gate_vars->mNaT,sizeof(PetscReal)*Nx*Ny*Nz);
        memcpy(user->gate_vars_past->hNaT,user->gate_vars->hNaT,sizeof(PetscReal)*Nx*Ny*Nz);
        memcpy(user->gate_vars_past->gNaT,user->gate_vars->gNaT,sizeof(PetscReal)*Nx*Ny*Nz);
        memcpy(user->gate_vars_past->mNaP,user->gate_vars->mNaP,sizeof(PetscReal)*Nx*Ny*Nz);
        memcpy(user->gate_vars_past->hNaP,user->gate_vars->hNaP,sizeof(PetscReal)*Nx*Ny*Nz);
        memcpy(user->gate_vars_past->gNaP,user->gate_vars->gNaP,sizeof(PetscReal)*Nx*Ny*Nz);
        memcpy(user->gate_vars_past->gKA,user->gate_vars->gKA,sizeof(PetscReal)*Nx*Ny*Nz);
        memcpy(user->gate_vars_past->hKA,user->gate_vars->hKA,sizeof(PetscReal)*Nx*Ny*Nz);
        memcpy(user->gate_vars_past->mKA,user->gate_vars->mKA,sizeof(PetscReal)*Nx*Ny*Nz);
        memcpy(user->gate_vars_past->mKDR,user->gate_vars->mKDR,sizeof(PetscReal)*Nx*Ny*Nz);
        memcpy(user->gate_vars_past->gKDR,user->gate_vars->gKDR,sizeof(PetscReal)*Nx*Ny*Nz);
        memcpy(user->gate_vars_past->yNMDA,user->gate_vars->yNMDA,sizeof(PetscReal)*Nx*Ny*Nz);
        memcpy(user->gate_vars_past->zNMDA,user->gate_vars->zNMDA,sizeof(PetscReal)*Nx*Ny*Nz);
        memcpy(user->gate_vars_past->dNMDA,user->gate_vars->dNMDA,sizeof(PetscReal)*Nx*Ny*Nz);
        memcpy(user->gate_vars_past->gNMDA,user->gate_vars->gNMDA,sizeof(PetscReal)*Nx*Ny*Nz);
    }



    free(tmp);
    free(line);
    fclose(fp);
    return;
    //Modify beginning of file
    fp = fopen("data_csd.txt", "r");
    if (fp == NULL) {
        fprintf(stderr, "File not found....\n");
        exit(EXIT_FAILURE); /* indicate failure.*/
    }
    //Read top file details
    fgets(line, 1024 * 1024, fp);

    tmp = strdup(line);

    Nx = atoi(getfield(tmp, 1));

    tmp = strdup(line);
    Ny = atoi(getfield(tmp, 2));

    tmp = strdup(line);
    numrecords = atoi(getfield(tmp, 3));
    fclose(fp);

    fp = fopen("data_csd.txt","r+");
    fseek( fp, 0, SEEK_SET );
    fprintf(fp, "%d,%d,%d,%d,%d\n", Nx, Ny, numrecords+(PetscInt)floor(Time/trecordstep)-1, Nc, Ni);

    fclose(fp);
}

void velocity_field(FILE *fp,struct AppCtx *user,PetscInt numrecords,int start) {


    if (start) {
            fprintf(fp, "%d,%d,%d,%d,%d,%d,%d,%d\n", user->Nx, user->Ny, numrecords, Nc, Ni, use_en_deriv, separate_vol,
                    Linear_Diffusion);

        velocity_field(fp, user, numrecords, 0);
    } else {
        PetscReal *c = user->state_vars->c;
        PetscReal *phi = user->state_vars->phi;
        PetscReal *al = user->state_vars->alpha;

        diff_coef(user->Dcs, al, 1, user);
        PetscReal *Dcs = user->Dcs;


        PetscReal dt = user->dt;
        PetscReal dx = user->dx;
        PetscReal dy = user->dy;
        PetscInt Nx = user->Nx;
        PetscInt Ny = user->Ny;
        PetscInt Nz = user->Nz;

        PetscInt x, y,z, ion, comp;

         PetscReal Gradleft, GradRight;
            PetscReal GradUp, GradDown;
            PetscInt count_x, count_y;
            PetscReal alNc, alNcRight, alNcUp;
            PetscReal vx, vy;

            PetscReal Gradx, Grady;

            z=0;
            for (ion = 0; ion < Ni; ion++) {
                for (comp = Nc - 1; comp < Nc; comp++) {
                    for (y = 0; y < Ny; y++) {
                        for (x = 0; x < Nx; x++){
                            count_x = 0;
                            count_y = 0;

                            Gradleft = 0;
                            GradRight = 0;
                            if(x > 0){
                                //First difference term
                                Gradleft = 1;//Dcs[c_index(x-1,y,comp,ion,Nx)*2]*(c[c_index(x-1,y,comp,ion,Nx)]+c[c_index(x,y,comp,ion,Nx)])/2;
                                Gradleft = Gradleft*(log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-
                                                     log(c[c_index(x-1,y,z,comp,ion,Nx,Ny,Nz)])+
                                                     z_charge[ion]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-
                                                                    phi[phi_index(x-1,y,z,comp,Nx,Ny,Nz)]))/dx;
                                count_x++;
                            }
                            //Add Second right moving difference
                            if(x < Nx-1){
                                GradRight = 1;//Dcs[c_index(x,y,comp,ion,Nx)*2]*(c[c_index(x,y,comp,ion,Nx)]+c[c_index(x+1,y,comp,ion,Nx)])/2;
                                GradRight = GradRight*(log(c[c_index(x+1,y,z,comp,ion,Nx,Ny,Nz)])-
                                                       log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])+
                                                       z_charge[ion]*(phi[phi_index(x+1,y,z,comp,Nx,Ny,Nz)]-
                                                                      phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]))/dx;
                                count_x++;
                            }
                            GradDown = 0;
                            GradUp = 0;
                            //Up down difference
                            if(y > 0){
                                GradDown = 1;//Dcs[c_index(x,y-1,comp,ion,Nx)*2+1]*(c[c_index(x,y-1,comp,ion,Nx)]+c[c_index(x,y,comp,ion,Nx)])/2;
                                GradDown = GradDown*(log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])-
                                                     log(c[c_index(x,y-1,z,comp,ion,Nx,Ny,Nz)])+
                                                     z_charge[ion]*(phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-
                                                                    phi[phi_index(x,y-1,z,comp,Nx,Ny,Nz)]))/dy;
                                count_y++;
                            }
                            //Next upward difference
                            if(y < Ny-1){
                                GradUp = 1;//Dcs[c_index(x,y,comp,ion,Nx)*2+1]*(c[c_index(x,y,comp,ion,Nx)]+c[c_index(x,y+1,comp,ion,Nx)])/2;
                                GradUp = GradUp*(log(c[c_index(x,y+1,z,comp,ion,Nx,Ny,Nz)])-
                                                 log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)])+
                                                 z_charge[ion]*(phi[phi_index(x,y+1,z,comp,Nx,Ny,Nz)]-
                                                                phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]))/dy;
                                count_y++;
                            }

                            Gradx = (Gradleft+GradRight)/count_x;
                            Grady = (GradUp+GradDown)/count_y;

                            if(x == Nx-1 & y == Ny-1){
                                fprintf(fp,"%f,%f\n",Gradx,Grady);
                            }else{
                                fprintf(fp,"%f,%f,",Gradx,Grady);
                            }

                        }
                    }
                }
            }
    }
}

void calculate_measures(FILE *fp, struct AppCtx *user,PetscInt numrecords,int start)
{
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt Nz = user->Nz;
    PetscInt x,y;

    PetscReal dx = Lx/Nx;
    PetscReal dy = Ly/Ny;

    PetscReal *c= user->state_vars->c;
    PetscReal *al = user->state_vars->alpha;

    PetscReal Neu_Gli_K_diff = 0;
    PetscReal Glia_K_per_amt = 0;
    PetscReal total_amt_K=0;
    PetscReal alN;
    for(x=0;x<Nx;x++){
        for(y=0;y<Ny;y++){
            alN = 1-al[al_index(x,y,0,0,Nx,Ny,Nz)]-al[al_index(x,y,0,1,Nx,Ny,Nz)];
            total_amt_K += al[al_index(x,y,0,0,Nx,Ny,Nz)]*c[c_index(x,y,0,0,1,Nx,0,0)]+al[al_index(x,y,0,1,
                                                                                                  Nx,0,0)]*c[c_index(
                    x,y,0,1,1,
                    Nx,0,0)]+alN*c[c_index(
                    x,y,0,Nc-1,1,Nx,0,0)];

            Neu_Gli_K_diff += al[al_index(x,y,0,0,Nx,0,0)]*c[c_index(x,y,0,0,1,Nx,0,0)]-al[al_index(x,y,0,
                                                                                                    1,Nx,
                                                                                                    0,0)]*c[c_index(
                    x,y,0,1,
                    1,Nx,0,0)];

            Glia_K_per_amt += al[al_index(x,y,0,1,Nx,0,0)]*c[c_index(x,y,0,1,1,Nx,0,0)];
        }
    }

    Neu_Gli_K_diff = Neu_Gli_K_diff*dx*dy;

    Glia_K_per_amt = (Glia_K_per_amt*dx*dy)/(total_amt_K*dx*dy);


    fprintf(fp,"%.10e,%.10e\n",Neu_Gli_K_diff,Glia_K_per_amt);
}
void calculate_energy(FILE *fp, struct AppCtx *user, PetscInt numrecords, int start){
    if (start) {
        fprintf(fp, "%d,%d,%d,%d,%d\n", user->Nx, user->Ny, numrecords, 0, 0);
        calculate_energy(fp, user, numrecords, 0);
    }else{
        PetscScalar *c = user->state_vars->c;
        PetscScalar *phi = user->state_vars->phi;
        PetscScalar  *al = user->state_vars->alpha;
        PetscInt Nx = user->Nx;
        PetscInt Ny = user->Ny;
        PetscInt Nz = user->Nz;
        PetscReal dz = user->dz;
        PetscReal *cm = user->con_vars->cm;

        PetscScalar Energy,alNc;
        PetscInt comp,ion,z;
        for(PetscInt y=0;y<Ny;y++){
            for(PetscInt x = 0; x < Nx; x++){
                Energy = 0;
                // Integrate over z
                for(z = 0; z < Nz; z++){
                    //Ionic contribution
                    for(comp = 0; comp < Nc-1; comp++){
                        //Immobile ion part
                        Energy += user->con_vars->ao[al_index(x,y,z,comp,Nx,Ny,Nz)]*log(user->con_vars->ao[al_index(x,y,z,comp,Nx,Ny,Nz)])/al[al_index(x,y,z,comp,Nx,Ny,Nz)];
                        //Mobile ions
                        for(ion = 0; ion < Ni; ion++){
                            Energy += al[al_index(x,y,z,comp,Nx,Ny,Nz)]*c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*
                                      log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]);
                        }
                        //ElectroPotential part
                        Energy += (cm[al_index(x,y,0,comp,Nx,Ny,Nz)]/2)*pow((phi[phi_index(x,y,z,comp,Nx,Ny,Nz)]-
                                                                            phi[phi_index(x,y,z,Nc-1,Nx,Ny,Nz)])*RTFC,2);
                    }
                    //Extracellular ion term
                    comp = Nc-1;
                    alNc = 1-al[al_index(x,y,z,0,Nx,Ny,Nz)]-al[al_index(x,y,z,1,Nx,Ny,Nz)];
                    //Immobile ion part
                    Energy += user->con_vars->ao[al_index(x,y,z,comp,Nx,Ny,Nz)]*log(user->con_vars->ao[al_index(x,y,z,comp,Nx,Ny,Nz)])/(alNc);
                    //Mobile ions
                    for(ion = 0; ion < Ni; ion++){
                        Energy += alNc*c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]*log(c[c_index(x,y,z,comp,ion,Nx,Ny,Nz)]);
                    }
                }
                Energy*=dz;
                //Write to file
                if(x == Nx-1 & y == Ny-1){
                    fprintf(fp,"%.10e\n",Energy);
                }else{
                    fprintf(fp,"%.10e,",Energy);
                }


            }
        }

    }
}

void draw_csd(struct AppCtx *user)
{

    PetscReal vm,threshhold;
    threshhold = -30;
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscInt Nz = user->Nz;
//    PetscInt z = 2;
    PetscInt x=0;
//    printf("Layer: %d\n",z);
    for(PetscInt y=0;y<Ny;y++){
        printf("_");
    }
    printf("\n");
//    for(PetscInt x=0;x<Nx;x++){
        for(PetscInt z=0;z<Nz;z++){
        printf("|");
        for(PetscInt y=0;y<Ny;y++){
            vm = user->state_vars->phi[phi_index(x,y,z,0,Nx,Ny,Nz)]-user->state_vars->phi[phi_index(x,y,z,Nc-1,Nx,Ny,Nz)];
            vm = vm*RTFC;
            if (vm > threshhold) {
                printf("x");
            } else {
                printf(" ");
            }
        }
        printf("|\n");
    }
    for(PetscInt y=0;y<Ny;y++){
        printf("_");
    }
    printf("\n");
}
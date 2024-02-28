//
// Created by zhangmh on 2023/4/30.
//
#include"symbol.h"

double** spatial_time_sequence( const grid_information * gi0){
    double dx=gi0->length/(gi0->xgridnumber);
    auto *pp=new double[(gi0->xgridnumber)+2];
    auto *tp=new double[(gi0->timegridnumber)+1];
    auto **ptr0=new double*[2];
    *ptr0=pp;*(ptr0+1)=tp;
    pp[0]=0;
    pp[(gi0->xgridnumber)+1]=gi0->length;
    for(int i=1;i<(gi0->xgridnumber)+1;i++){
        *(pp+i)=(i-0.5)*dx;
    }
    for(int j=0;j<(gi0->timegridnumber)+1;j++){
        *(tp+j)=j*(gi0->timestep);
    }
    return ptr0;
}

double **parameter_field(const grid_information *gi1, optical_thermal_parameter *ohp){
    double ** ppa= matrix0(5,gi1->xgridnumber+2);
    for(int i=0;i<(gi1->xgridnumber)+2;i++){
        *(*ppa+i)=ohp->absorptioncof;
        *(*(ppa+1)+i)=ohp->scattercof;
        *(*(ppa+2)+i)=ohp->absorptioncof+ohp->scattercof;  //extinction coefficient
        *(*(ppa+3)+i)=ohp->conductivity;
        *(*(ppa+4)+i)=ohp->heat_capacity;
    }
    return ppa;
}
double *temperature_initial_startmoment(const grid_information * gg,const optical_thermal_parameter *op){
    auto *tis0= new double [gg->xgridnumber+2]();
    for(int ii=0;ii<gg->xgridnumber+2;ii++){
        tis0[ii]=op->right_Tf;
    }
    return tis0;
}
double **temperature_initial(const optical_thermal_parameter *oh1,const grid_information *gg0){
    double **tis11= matrix0(2,gg0->xgridnumber+2);
    for(auto jj=0;jj<(gg0->xgridnumber)+2;jj++){
        *(*tis11+jj)=oh1->right_Tf;
        *(*(tis11+1)+jj)=oh1->right_Tf;
    }
    return tis11;
}
double ***radiative_intensity(const grid_information *gg2){
    double ** thh= matrix0(gg2->xgridnumber+2,gg2->anglegridnumber);
    double ** th11= matrix0(gg2->xgridnumber+2,gg2->anglegridnumber);
    auto ***t=new double **[2];
    *t=thh;*(t+1)=th11;
    return t;
}
double *** radiative_initial(const grid_information& g1,const optical_thermal_parameter& opp){
     double *** radiative_intensity0= radiative_intensity(&g1);
    for(auto kk=0;kk<2;kk++) {
        for (auto ii = 0; ii < g1.xgridnumber + 2; ii++) {
            for (auto jj = 0; jj < g1.anglegridnumber; jj++) {
                radiative_intensity0[kk][ii][jj] =opp.radiative_initial_field;
                }
            }
        }
    return radiative_intensity0;
}
//
// Created by zhangmh on 2023/4/30.
//
#include"symbol.h"
#include<iostream>
#include<cmath>
void radiative_boundary_condition(const grid_information& g0,const optical_thermal_parameter& otp,double *** radiative_intensity,const double * tem_old,const double *wv0){
    for(auto i=0;i<g0.anglegridnumber;i++){   //radiative boundary condition
         radiative_intensity[0][0][i]=0;
         radiative_intensity[0][g0.xgridnumber+1][i]=0;
         for(int j=0;j<g0.anglegridnumber;j++){
             if(wv0[j]<0){
                 radiative_intensity[0][0][i]+=radiative_intensity[1][0][j]*std::abs(wv0[j])*(1-otp.leftwallem)*2;
             }
             else if(wv0[j]>0) {
                 radiative_intensity[0][g0.xgridnumber+1][i]+=radiative_intensity[1][g0.xgridnumber+1][j]*std::abs(wv0[j])*(1-otp.rightwallem)*2;
             }
         }
         radiative_intensity[0][0][i]+=otp.leftwallem*sb*pow(REFRACTIVE_INDEX,2)*pow(tem_old[0],4)/PI;
         radiative_intensity[0][g0.xgridnumber+1][i]+=otp.rightwallem*sb*pow(REFRACTIVE_INDEX,2)*pow(tem_old[g0.xgridnumber+1],4)/PI;
    }
}
void radiative_intensity_solver(const grid_information& g0,const double *wv,double **ri,const double **tsi,const double * delta_om,const double* ex_co){
     double west_surface_intensity[g0.xgridnumber+1][g0.anglegridnumber];//temporary variables
     double east_surface_intensity[g0.xgridnumber+1][g0.anglegridnumber];
     double volume=g0.length/g0.xgridnumber;
     for(int i=0;i<g0.anglegridnumber;i++){
         if(wv[i]>0){
             west_surface_intensity[1][i]=ri[0][i];
             for(int j=1;j<g0.xgridnumber+1;j++){  //interior node   left side-> right side
                 ri[j][i]=(west_surface_intensity[j][i]*1*wv[i]/(g0.diff_factor)+volume*tsi[j][i]*delta_om[i])/(wv[i]*1/g0.diff_factor+ex_co[j]*volume*delta_om[i]);
                 east_surface_intensity[j][i]=(ri[j][i]+(g0.diff_factor-1)*west_surface_intensity[j][i])/g0.diff_factor;
                 if(j<g0.xgridnumber) {
                   west_surface_intensity[j + 1][i] = east_surface_intensity[j][i];
                 }
                 if(j==g0.xgridnumber){
                   ri[j+1][i]=east_surface_intensity[j][i];
                 }
            }
         }
         else if(wv[i]<0){
             east_surface_intensity[g0.xgridnumber][i]=ri[g0.xgridnumber+1][i];
             for(int jj=g0.xgridnumber;jj>0;jj--){//  right side->left side
                 ri[jj][i]=(east_surface_intensity[jj][i]*1*(-wv[i])/(g0.diff_factor)+volume*tsi[jj][i]*delta_om[i])/((-wv[i])*1/g0.diff_factor+ex_co[jj]*volume*delta_om[i]);
                 west_surface_intensity[jj][i]=(ri[jj][i]-g0.diff_factor*east_surface_intensity[jj][i])/(1-g0.diff_factor);
                 if(jj>1){
                     east_surface_intensity[jj-1][i]=west_surface_intensity[jj][i];
                 }
                 if(jj==1){
                     ri[jj-1][i]=west_surface_intensity[jj][i];
                 }
             }
         }
     }
}
void energy_boundary_condition1(const double &left_wall_tem,const double &right_wall_tem,double *tem,const grid_information& gi){
     tem[0]=left_wall_tem;
     tem[gi.xgridnumber+1]=right_wall_tem;
}

void energy_boundary_condition2(double **tem,const double *conductivity,const grid_information& gi,const optical_thermal_parameter& op){
    double delta_x=gi.length/gi.xgridnumber;
    tem[0][0]=(2*conductivity[0]*tem[1][1]+op.left_hf*op.left_Tf*delta_x)/(2*conductivity[0]+op.left_hf*delta_x);
    tem[0][gi.xgridnumber+1]=(2*conductivity[gi.xgridnumber+1]*tem[1][gi.xgridnumber]+op.right_hf*op.right_Tf*delta_x)/(2*conductivity[gi.xgridnumber+1]+op.right_hf*delta_x);
    std::cout<<tem[0][0]<<"  "<<tem[0][gi.xgridnumber+1]<<std::endl;
}
void energy_boundary_condition3(const grid_information& gi,const optical_thermal_parameter& op,double **tem,const double *heat_flux_b,const double *conductivity,const double &laser_power){
    double delta_x=gi.length/gi.xgridnumber;
    double a1=op.leftwallem*sb;
    double a2=op.rightwallem*sb;
    double b1=op.left_hf+2*conductivity[0]/delta_x;
    double b2=op.right_hf+2*conductivity[gi.xgridnumber+1]/delta_x;
    double c1=op.leftwallem*laser_power+op.left_hf*op.left_Tf+op.leftwallem*sb*pow(op.left_Tf,4)-heat_flux_b[0]+2*conductivity[0]/delta_x*tem[1][1];
    double c2=op.right_hf*op.right_Tf+op.rightwallem*sb*pow(op.right_Tf,4)-heat_flux_b[1]+2*conductivity[gi.xgridnumber+1]/delta_x*tem[1][gi.xgridnumber];
    tem[0][0]= equation_solver(300,a1,b1,c1);
    tem[0][gi.xgridnumber+1]= equation_solver(300,a2,b2,c2);
}
void energy_equation_solver(const double ** iteration_coe,double **tem,const grid_information& gi){
      for(int ii=1;ii<gi.xgridnumber;ii++){              //left_side->right_side
          tem[0][ii]=(iteration_coe[0][ii]*tem[1][ii+1]+iteration_coe[1][ii]*tem[0][ii-1]+iteration_coe[3][ii])/iteration_coe[4][ii];
      }
    tem[0][gi.xgridnumber]=(iteration_coe[0][gi.xgridnumber]*tem[0][gi.xgridnumber+1]+iteration_coe[1][gi.xgridnumber]*tem[0][gi.xgridnumber-1]+iteration_coe[3][gi.xgridnumber])/iteration_coe[4][gi.xgridnumber];
}


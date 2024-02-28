//
// Created by zhangmh on 2023/4/30.
//
#include"symbol.h"
#include<cmath>
#include<string>
#include<iomanip>
#include<fstream>
double ** matrix0(int rownumber,int colnumber){
    auto* ma=new double [rownumber*colnumber]();
    auto ** matrix=new double *[rownumber];
    for(int ii=0;ii<rownumber;ii++){
        *(matrix+ii)=ma+colnumber*ii;
    }
    return matrix;
}
double ** angle_division(const grid_information * gii){
          double dcta=PI/gii->anglegridnumber;
          double **app= matrix0(2,gii->anglegridnumber);
          for(auto k=0;k<gii->anglegridnumber;k++){
              *((*app)+k)=(k+1-0.5)*dcta;//discrete direction
              *(*(app+1)+k)=2*sin(*((*app)+k))*sin(dcta/2)*DFI;//size of solid angle
          }
          return app;
}
double* weight_value(const double *cta,const grid_information *gg){              //n is x positive axis
         auto* wv=new double[gg->anglegridnumber];
         for(int i=0;i<gg->anglegridnumber;i++){
             wv[i]=sin(*(cta+i))*cos(*(cta+i))*sin(PI/gg->anglegridnumber)*DFI;
         }
        return wv;
}
//linear anisotropic phase function(one dimension)
//double ** phase_function(const double *cta,const grid_information& gi,double factor){
//      double **pf0= matrix0(gi.anglegridnumber,gi.anglegridnumber);
//       for(int j=0;j<gi.anglegridnumber;j++){
//           for(int k=0;k<gi.anglegridnumber;k++){
//               pf0[j][k]=1+factor*cos(cta[j])*cos(cta[k]);
//           }
//       }
//     return pf0;
//}

void incident_radiation(double *ir, const double **ri,const double *delta_om,const grid_information& gii){
    for(int i=0;i<gii.xgridnumber+2;i++) {
         ir[i]=0;
        for (int j = 0; j < gii.anglegridnumber; j++) {
            ir[i]+=ri[i][j]*delta_om[j]*2*PI;
        }
    }
}
void heat_flux(double *hf,const double **ri,const double *wv,const grid_information& gi){
    for(int i=0;i<gi.xgridnumber+2;i++){
         hf[i]=0;
        for(int j=0;j<gi.anglegridnumber;j++){
            hf[i]+=ri[i][j]*wv[j]*2*PI;
        }
    }
}

void  heat_flux_boundary(double * heat_flux_boundary0,const double **ri,const double* tem_old,const grid_information& gi,const optical_thermal_parameter& op,double * wv){
      heat_flux_boundary0[0]=0;
      heat_flux_boundary0[1]=0;
      for(int i=0;i<gi.anglegridnumber;i++){
          if(wv[i]<0){
              heat_flux_boundary0[0]+=ri[0][i]*wv[i]*2*PI;
          }
          else if(wv[i]>0){
              heat_flux_boundary0[1]+=ri[gi.xgridnumber+1][i]*wv[i]*2*PI;
          }
      }
      heat_flux_boundary0[0]+=op.leftwallem*sb*pow(REFRACTIVE_INDEX,2)*pow(tem_old[0],4);
      heat_flux_boundary0[1]+=op.rightwallem*sb*pow(REFRACTIVE_INDEX,2)*pow(tem_old[gi.xgridnumber+1],4);
}
void transport_source_item(double **tsi0,const double* ir,const double* hf,const double* ad,const double* t,const double**parameter, const double &scattering_factor,const grid_information& gi ){
    for(int i=0;i<gi.xgridnumber+2;i++){
        for(int j=0;j<gi.anglegridnumber;j++){
             tsi0[i][j]=parameter[0][i]*pow(REFRACTIVE_INDEX,2)*sb*pow(t[i],4)/PI+parameter[1][i]/(4*PI)*(ir[i]+scattering_factor*cos(ad[j])*hf[i]);
         }
     }
}
void non_negative_produce(double ** ri,const grid_information& gi){
    for(int i=0;i<gi.xgridnumber+2;i++){
        for(int j=0;j<gi.anglegridnumber;j++){
            if(ri[i][j]<0){
                ri[i][j]=0;
            }
        }
    }
}
bool radiative_iteration_judge(const double **ri,const double **ri_old,const grid_information& gi){
      double ep_max=0;
      for(int ii=0;ii<gi.xgridnumber+2;ii++) {
          for (int jj = 0; jj < gi.anglegridnumber; jj++) {
              double ep1 = std::abs(ri[ii][jj] - ri_old[ii][jj]) / (ri[ii][jj] + 1e-40);
              if (ep1 >= ep_max) {
                  ep_max = ep1;
              }
          }
      }
      if(ep_max<=epsilon){
          return true;
      }
      else{
          return false;
      }
}
bool tem_iteration_judge(const double ** tem,const grid_information& gi){
    double ep_max0=0;
    for(int i=0;i<gi.xgridnumber+2;i++){
        double ep=std::abs(tem[0][i]-tem[1][i])/(tem[0][i]+1e-40);
        if(ep>ep_max0){
            ep_max0=ep;
        }
    }
    if(ep_max0<=epsilon){
        return true;
    }
    else{
        return false;
    }

}
void  renew_tem(double **tem,const grid_information& gi){
    for(int ii=0;ii<gi.xgridnumber+2;ii++){
        tem[1][ii]=tem[0][ii];
    }
}
void renew_tem_old_moment(const double * tem_new_moment, double * tem_old_moment,const grid_information& gi){
    for(int i=0;i<gi.xgridnumber+2;i++){
        tem_old_moment[i]=tem_new_moment[i];
    }

}

void renew_radiative_intensity(const double ** ri,double ** ri_old,const grid_information& gi){
    for(int i=0;i<gi.xgridnumber+2;i++){
        for(int j=0;j<gi.anglegridnumber;j++){
           ri_old[i][j]=ri[i][j];
        }
    }
}
void radiative_source_item(double ** rsi,const double *ir_new,const double *tem_old,const grid_information& gi,const double *ac) {    //no minus sign
    for(int i=0;i<gi.xgridnumber+2;i++){
        rsi[0][i]=-ac[i]*ir_new[i]; //Sc
        rsi[1][i]=4*ac[i]*pow(REFRACTIVE_INDEX,2)*sb*pow(tem_old[i],3);//Sp
    }
}
//iteration_coe matrix 5*(xgridnumber+1)
void iteration_co(double ** iteration_coe,const grid_information& gi,const optical_thermal_parameter& op,const double *old_moment_tem,const double ** rsi,const double *conductivity){
    iteration_coe[0][gi.xgridnumber]=2*conductivity[gi.xgridnumber+1]/(gi.length/gi.xgridnumber);
    iteration_coe[1][1]=2*conductivity[0]/(gi.length/gi.xgridnumber);
    for(int i=1;i<gi.xgridnumber+1;i++){
         if(i!=1&&i!=gi.xgridnumber) {
             iteration_coe[0][i] = 2 * conductivity[i + 1] * conductivity[i] / (gi.length / gi.xgridnumber) / (conductivity[i + 1] + conductivity[i]);//ae
             iteration_coe[1][i] = 2 * conductivity[i - 1] * conductivity[i] / (gi.length / gi.xgridnumber) / (conductivity[i - 1] +conductivity[i]);//aw
         }
         iteration_coe[2][i]=op.heat_capacity*(gi.length/gi.xgridnumber)/gi.timestep;//apo
         iteration_coe[3][i]=old_moment_tem[i]*iteration_coe[2][i]+(-rsi[0][i])*(gi.length/gi.xgridnumber);//b
         iteration_coe[4][i]=iteration_coe[0][i]+iteration_coe[1][i]+iteration_coe[2][i]+rsi[1][i]*(gi.length/gi.xgridnumber);//ap
    }
}
double equation_solver(double iv,const double& a,const double &b,const double &c ){
    double result=iv;
    do{
         iv=result;
        double f=a*pow(iv,4)+b*iv-c;
        double f_de=4*a*pow(iv,3)+b;
         result=iv-f/f_de;
    }while(std::abs(result-iv)>1e-5);
    return result;
}
void output(const std::string& str,const int &number,const double *spatial_time,const double *tem_info,const grid_information& gi,const optical_thermal_parameter& op){
    std::ofstream  file_out;
    std::string file_name=R"(C:\Users\zhangmh\Desktop\)"+str+std::to_string(number)+".txt";
    file_out.open(file_name);
    if(str=="transient"){
        for(int ii=0;ii<gi.timegridnumber+1;ii++){
            file_out<<spatial_time[ii]<<"   "<<tem_info[ii]-op.left_Tf<<std::fixed<<std::setprecision(4)<<std::endl;
        }
    }
    else if(str=="steady"){
        for(int jj=0;jj<gi.xgridnumber+2;jj++){
            file_out<<spatial_time[jj]/gi.length<<"   "<<tem_info[jj]/op.left_Tf<<std::fixed<<std::setprecision(4)<<std::endl;
        }
    }
}

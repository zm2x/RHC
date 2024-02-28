//
// Created by zhangmh on 2023/4/30.
//
//This program is used for computing radiation and conduction coupled problem in an emitting,absorbing,scattering participating media.
//one dimension condition(a very large slab)
#ifndef RAHC1_SYMBOL_H
#define RAHC1_SYMBOL_H
#include <string>
#define DFI   1.0
#define PI 3.14159
#define REFRACTIVE_INDEX 1.0
const double sb=5.67e-8;
const double epsilon=1e-6;
const double laser_time=1.0;
typedef struct {
    double length;//length of slab
    double timestep;
    double diff_factor; //linear difference factor
    int   xgridnumber;
    int  anglegridnumber;
    int  timegridnumber;
}grid_information;
typedef struct {
    double right_Tf;
    double left_Tf;
    double radiative_initial_field;
    double leftwallem;
    double rightwallem;
    double right_hf;
    double left_hf;
    double absorptioncof;
    double scattercof;
    double conductivity;
//    double laser_time;
//    double laser_power;
    double heat_capacity;
}optical_thermal_parameter;
double ** matrix0(int rownumber,int colnumber);
double ** spatial_time_sequence(const grid_information *gi);
double **parameter_field(const grid_information *gi, optical_thermal_parameter *ohp0);
double *temperature_initial_startmoment(const grid_information *g,const optical_thermal_parameter *op);
double **temperature_initial(const optical_thermal_parameter *ohp2,const grid_information *gin);
double ***radiative_intensity(const grid_information *gp);
double *** radiative_initial(const grid_information& g1,const optical_thermal_parameter& opp);
double **angle_division(const grid_information * gii);
double* weight_value(const double *cta,const grid_information *gg);
//double ** phase_function(const double *cta,const grid_information& gi,double factor);
void incident_radiation(double *ir, const double **ri,const double *delta_om,const grid_information& gii);
void heat_flux(double *hf,const double **ri,const double *wv,const grid_information& gi);
void  heat_flux_boundary(double * heat_flux_boundary0,const double **ri,const double* tem_old,const grid_information& gi,const optical_thermal_parameter& op,double * wv);
void transport_source_item(double **tsi0,const double* ir,const double* hf,const double* ad,const double* t,const double**parameter, const double& scattering_factor,const grid_information& gi );
void radiative_boundary_condition(const grid_information& g0,const optical_thermal_parameter& otp,double *** radiative_intensity,const double * tem_old,const double *wv0);
void radiative_intensity_solver(const grid_information& g0,const double *wv,double **ri,const double **tsi,const double * delta_om,const double* ex_co);
void energy_boundary_condition1(const double &left_wall_tem,const double &right_wall_tem,double *tem,const grid_information& gi);
void energy_boundary_condition2(double **tem,const double *conductivity,const grid_information& gi,const optical_thermal_parameter& op);
void energy_boundary_condition3(const grid_information& gi,const optical_thermal_parameter& op,double **tem,const double *heat_flux_b,const double *conductivity,const double &laser_power);
void non_negative_produce(double ** ri,const grid_information& gi);
bool radiative_iteration_judge(const double **ri,const double **ri_old,const grid_information& gi);
bool tem_iteration_judge(const double ** tem,const grid_information& gi);
void  renew_tem(double **tem,const grid_information& gi);
void renew_tem_old_moment(const double * tem_new_moment, double * tem_old_moment,const grid_information& gi);
void renew_radiative_intensity(const double ** ri,double ** ri_old,const grid_information& gi);
void radiative_source_item(double ** rsi,const double *ir_new,const double *tem_old,const grid_information& gi,const double *ac);
void iteration_co(double ** iteration_coe,const grid_information& gi,const optical_thermal_parameter& op,const double *ola_moment_tem,const double ** rsi,const double *ac);
void energy_equation_solver(const double ** iteration_coe,double **tem,const grid_information& gi);
void computation();
double equation_solver(double iv,const double& a,const double &b,const double &c );
void output(const std::string& str,const int &number,const double *spatial_time,const double *tem_info,const grid_information& gi,const optical_thermal_parameter& op);
#endif //RAHC1_SYMBOL_H

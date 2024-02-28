//
// Created by zhangmh on 2023/5/10.
//
#include"symbol.h"
#include<iostream>
void computation() {
    double scattering_factor = 0;
    double laser_power=50000;
    grid_information gi = {1.0, 0.001, 0.5, 200, 100,1800 };
    optical_thermal_parameter op = {500, 1000, 200, 1.0, 1.0, 2.0, 2.0, 1.0,0, 226.8,1.0e4};
    auto *tem_t_rightwall=new double [gi.timegridnumber+1]();
    double **spatial_time_sequence0 = spatial_time_sequence(&gi);
    double **parameter_fields = parameter_field(&gi, &op);
    double *tem_start = temperature_initial_startmoment( &gi,&op);
    double **tem = temperature_initial(&op, &gi);
    double ***radiative_intensity = radiative_initial(gi, op);
    double **angle_div = angle_division(&gi);
    double *wv = weight_value(angle_div[0], &gi);
    auto *ir0 = new double[gi.xgridnumber + 2]();
    auto *heat_flux0 = new double[gi.xgridnumber + 2]();
    double **tsi = matrix0(gi.xgridnumber + 2, gi.anglegridnumber);
    auto  * heat_flux_b=new double [2];
    double ** radiative_source_item0= matrix0(2,gi.xgridnumber+2);
    auto ** iteration_coe= matrix0(5,gi.xgridnumber+1);
    int radiative_iteration_number=0;
    int tem_iteration_number=0;
    for (int i = 0; i < gi.timegridnumber+1 ; i++) {
        tem_t_rightwall[i]=tem_start[gi.xgridnumber+1];
        do {
            tem_iteration_number++;
            renew_tem(tem,gi);
            do {
                radiative_iteration_number++;
                renew_radiative_intensity((const double **) radiative_intensity[0], radiative_intensity[1], gi);
                radiative_boundary_condition(gi, op, radiative_intensity, tem[1], wv);
                incident_radiation(ir0, (const double **) (radiative_intensity[1]), angle_div[1], gi);
                heat_flux(heat_flux0, (const double **) (radiative_intensity[1]), wv, gi);
                transport_source_item(tsi, ir0, heat_flux0, angle_div[0], tem[1], (const double **) parameter_fields,
                                      scattering_factor, gi);
                radiative_intensity_solver(gi, wv, radiative_intensity[0], (const double **) tsi, angle_div[1],
                                           parameter_fields[2]);
                non_negative_produce(radiative_intensity[0], gi);
            } while (!radiative_iteration_judge((const double **) radiative_intensity[0],
                                                (const double **) radiative_intensity[1], gi));
            heat_flux_boundary(heat_flux_b, (const double **) radiative_intensity[0], tem[1], gi, op, wv);
           //energy_boundary_condition2(tem, parameter_fields[3],  gi, op);
            energy_boundary_condition1(op.left_Tf,op.right_Tf,tem[0],gi);
//           if(spatial_time_sequence0[1][i]<laser_time){
//               energy_boundary_condition3(gi,op,tem,heat_flux_b,parameter_fields[3],laser_power);
//           }
//           else{
//               energy_boundary_condition3(gi,op,tem,heat_flux_b,parameter_fields[3],0);
//           }

            incident_radiation(ir0, (const double **) radiative_intensity[0], angle_div[1], gi);
            radiative_source_item(radiative_source_item0, ir0, tem[1], gi, parameter_fields[0]);
            iteration_co(iteration_coe, gi, op, tem_start, (const double **) radiative_source_item0,
                         parameter_fields[3]);
            energy_equation_solver((const double **) iteration_coe, tem, gi);
        }while(!tem_iteration_judge((const double **)tem,gi));
        if(i==43) {
            output("steady",0, spatial_time_sequence0[0], tem[0], gi, op);
        }
        else if(i==219){
                output("steady",1,spatial_time_sequence0[0],tem[0],gi,op);
        }
        else if(i==660){
            output("steady",2,spatial_time_sequence0[0],tem[0],gi,op);
        }
        else if(i==1763){
            output("steady",3,spatial_time_sequence0[0],tem[0],gi,op);
        }
       renew_tem_old_moment(tem[0],tem_start,gi);
    }
    //output("transient",spatial_time_sequence0[1],tem_t_rightwall,gi,op);
    std::cout<<"radiative_iteration_number"<<"   "<<radiative_iteration_number<<std::endl;
    std::cout<<"tem_iteration_number"<<"   "<<tem_iteration_number<<std::endl;
}
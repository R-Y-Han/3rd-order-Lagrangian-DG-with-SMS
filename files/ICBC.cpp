/**
 * @file ICBC.cpp
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief initial condition
 * @version 0.1
 * @date 2022-11-29
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#include <cmath>
#include "ICBC.h"

double ini_rho(double x, double y)
{
    double ans;
    switch(testcase)
    {
        case shockless_Noh:
            ans = 1;
            break;
        case Taylor_Green_vortex:
            ans = 1;
            break;
        
        default:
            ans = 0;
            break;
    }
    return ans;
}

double ini_rho_x(double x, double y)
{
    double ans;
    switch(testcase)
    {
        case shockless_Noh:
            ans = 0;
            break;
        case Taylor_Green_vortex:
            ans = 0;
            break;
        
        default:
            ans = 0;
            break;
    }
    return ans;
}

double ini_rho_y(double x, double y)
{
    double ans;
    switch(testcase)
    {
        case shockless_Noh:
            ans = 0;
            break;
        case Taylor_Green_vortex:
            ans = 0;
            break;
        
        default:
            ans = 0;
            break;
    }
    return ans;
}

double ini_rho_xx(double x, double y)
{
    double ans;
    switch(testcase)
    {
        case shockless_Noh:
            ans = 0;
            break;
        case Taylor_Green_vortex:
            ans = 0;
            break;
        
        default:
            ans = 0;
            break;
    }
    return ans;
}

double ini_rho_xy(double x, double y)
{
    double ans;
    switch(testcase)
    {
        case shockless_Noh:
            ans = 0;
            break;
        case Taylor_Green_vortex:
            ans = 0;
            break;
        
        default:
            ans = 0;
            break;
    }
    return ans;
}

double ini_rho_yy(double x, double y)
{
    double ans;
    switch(testcase)
    {
        case shockless_Noh:
            ans = 0;
            break;
        case Taylor_Green_vortex:
            ans = 0;
            break;
        
        default:
            ans = 0;
            break;
    }
    return ans;
}

double ini_ux(double x, double y)
{
    double ans;
    switch(testcase)
    {
        case shockless_Noh:
            ans = -x;
            break;
        case Taylor_Green_vortex:
            ans = sin(PI * x) * cos(PI * y);
            break;
        
        default:
        ans = 0;
        break;
    }
    return ans;
}

double ini_ux_x(double x, double y)
{
    double ans;
    switch(testcase)
    {
        case shockless_Noh:
            ans = -1;
            break;
        case Taylor_Green_vortex:
            ans = PI * cos(PI * x) * cos(PI * y);
            break;

        default:
            ans = 0;
            break;
    }
    return ans;
}

double ini_ux_y(double x, double y)
{
    double ans;
    switch(testcase)
    {
        case shockless_Noh:
            ans = 0;
            break;
        case Taylor_Green_vortex:
            ans = - PI * sin(PI * x) * sin(PI * y);
            break;
        
        default:
            ans = 0;
            break;
    }
    return ans;
}

double ini_ux_xx(double x, double y)
{
    double ans;
    switch(testcase)
    {
        case shockless_Noh:
            ans = 0;
            break;
        case Taylor_Green_vortex:
            ans = - PI * PI * sin(PI * x) * cos(PI * y);
            break;
        
        default:
            ans = -sin(x);
            break;
    }
    return ans;
}

double ini_ux_xy(double x, double y)
{
    double ans;
    switch(testcase)
    {
        case shockless_Noh:
            ans = 0;
            break;
        case Taylor_Green_vortex:
            ans = - PI * PI * cos(PI * x) * sin(PI * y);
            break;
        
        default:
            ans = 0;
            break;
    }
    return ans;
}

double ini_ux_yy(double x, double y)
{
    double ans;
    switch(testcase)
    {
        case shockless_Noh:
            ans = 0;
            break;
        case Taylor_Green_vortex:
            ans = - PI * PI * sin(PI * x) * cos(PI * y);
            break;
        
        default:
            ans = 0;
            break;
    }
    return ans;
}

double ini_uy(double x, double y)
{
    double ans;
    switch(testcase)
    {
        case shockless_Noh:
            ans = -y;
            break;
        case Taylor_Green_vortex:
            ans = - cos(PI * x) * sin(PI * y);
            break;
        
        default:
            ans = 0;
            break;
    }
    return ans;
}

double ini_uy_x(double x, double y)
{
    double ans;
    switch(testcase)
    {
        case shockless_Noh:
            ans = 0;
            break;
        case Taylor_Green_vortex:
            ans = PI * sin(PI * x) * sin(PI * y);
            break;
        
        default:
            ans = 0;
            break;
    }
    return ans;
}

double ini_uy_y(double x, double y)
{
    double ans;
    switch(testcase)
    {
        case shockless_Noh:
            ans = -1;
            break;
        case Taylor_Green_vortex:
            ans = - PI * cos(PI * x) * cos(PI * y);
            break;
        
        default:
            ans = 0;
            break;
    }
    return ans;
}

double ini_uy_xx(double x, double y)
{
    double ans;
    switch(testcase)
    {
        case shockless_Noh:
            ans = 0;
            break;
        case Taylor_Green_vortex:
            ans = PI * PI * cos(PI * x) * sin(PI * y);
            break;
        
        default:
            ans = 0;
            break;
    }
    return ans;
}

double ini_uy_xy(double x, double y)
{
    double ans;
    switch(testcase)
    {
        case shockless_Noh:
            ans = 0;
            break;
        case Taylor_Green_vortex:
            ans = PI * PI * sin(PI * x) * cos(PI * y);
            break;
        
        default:
            ans = 0;
            break;
    }
    return ans;
}

double ini_uy_yy(double x, double y)
{
    double ans;
    switch(testcase)
    {
        case shockless_Noh:
            ans = 0;
            break;
        case Taylor_Green_vortex:
            ans = PI * PI * cos(PI * x) * sin(PI * y);
        
        default:
            ans = 0;
            break;
    }
    return ans;
}

double ini_p(double x, double y)
{
    double ans;
    double gamma;
    switch(testcase)
    {
        case shockless_Noh:
            gamma = 5.0 / 3.0;
            ans = (gamma - 1)* 1 * 1;
            break;
        case Taylor_Green_vortex:
            ans = 0.25 * ( cos(2 * PI * x) + cos(2 * PI * y) ) + 1;
            break;
        
        default:
            ans = 0;
            break;
    }
    return ans;
}

double ini_p_x(double x, double y)
{
    double ans;
    switch(testcase)
    {
        case shockless_Noh:
            ans = 0;
            break;
        case Taylor_Green_vortex:
            ans = - PI * sin(2 * PI * x) / 2.0;
            break;
        
        default:
            ans = 0;
            break;
    }
    return ans;
}

double ini_p_y(double x, double y)
{
    double ans;
    switch(testcase)
    {
        case shockless_Noh:
            ans = 0;
            break;
        case Taylor_Green_vortex:
            ans = -PI * sin(2 * PI * y) / 2;
            break;
        
        default:
            ans = 0;
            break;
    }
    return ans;
}

double ini_p_xx(double x, double y)
{
    double ans;
    switch(testcase)
    {
        case shockless_Noh:
            ans = 0;
            break;
        case Taylor_Green_vortex:
            ans = - PI * PI * cos(2 * PI *x);
            break;
        
        default:
            ans = 0;
            break;
    }
    return ans;
}

double ini_p_xy(double x, double y)
{
    double ans;
    switch(testcase)
    {
        case shockless_Noh:
            ans = 0;
            break;
        case Taylor_Green_vortex:
            ans = 0;
            break;

        default:
            ans = 0;
            break;
    }
    return ans;
}

double ini_p_yy(double x, double y)
{
    double ans;
    switch(testcase)
    {
        case shockless_Noh:
            ans = 0;
            break;
        case Taylor_Green_vortex:
            ans = - PI * PI * cos(2 * PI * y);
            break;
        
        default:
            ans = 0;
            break;
    }
    return ans;
}

double ana_e(double x, double y, double t)
{
    double ans;
    switch(testcase)
    {
        case shockless_Noh:
            ans = 1.0 / pow((1-t),4.0 / 3.0);
            break;
        case Taylor_Green_vortex:
            ans = ini_p(x,y) / (7.0 / 5.0 -1);
        
        default:
            ans = 0;
            break;
    }
    return ans;
}
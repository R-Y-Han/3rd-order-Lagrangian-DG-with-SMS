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
        
        default:
        ans = sin(x);
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

        default:
            ans = cos(x);
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
        
        default:
            ans = 0;
            break;
    }
    return ans;
}
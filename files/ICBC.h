/**
 * @file ICBC.h
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief initial condition
 * @version 0.1
 * @date 2022-11-29
 * 
 * @copyright Copyright (c) 2022 R. Y. Han. All Rights Reserved.
 * 
 */

#ifndef ICBC_H
#define ICBC_H

#include "config.h"

/**
 * @name 初始条件
 * @brief 定义初值条件函数
 * @{
 */

/**
 * @brief 初始密度场
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_rho(double x, double y);

double ini_rho_x(double x, double y);
double ini_rho_y(double x, double y);
double ini_rho_xx(double x, double y);
double ini_rho_xy(double x, double y);
double ini_rho_yy(double x, double y);

/**
 * @brief 初始速度场的x方向分量
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_ux(double x, double y);

/**
 * @brief 初始速度场对x求导的x分量，后面几个都是求偏导
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_ux_x(double x, double y);

double ini_ux_y(double x, double y);

double ini_ux_xx(double x, double y);

double ini_ux_xy(double x, double y);

double ini_ux_yy(double x, double y);

/**
 * @brief 初始速度场y方向分量
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_uy(double x, double y);

/**
 * @brief 初始速度场对y求导的x分量
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_uy_x(double x, double y);

double ini_uy_y(double x, double y);

double ini_uy_xx(double x, double y);

double ini_uy_xy(double x, double y);

double ini_uy_yy(double x, double y);

/**
 * @brief 初始压力场
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_p(double x, double y);

double ini_p_x(double x, double y);

double ini_p_y(double x, double y);

double ini_p_xx(double x, double y);

double ini_p_xy(double x, double y);

double ini_p_yy(double x, double y);

/**
 * @brief 内能解析解
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @param t 时间
 * @return double 
 */
double ana_e(double x, double y, double t);
/** @} */


#endif
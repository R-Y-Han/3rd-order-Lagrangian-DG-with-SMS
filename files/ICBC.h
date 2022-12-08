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
 * @name ��ʼ����
 * @brief �����ֵ��������
 * @{
 */

/**
 * @brief ��ʼ�ܶȳ�
 * 
 * @param x ���������
 * @param y ����������
 * @return double 
 */
double ini_rho(double x, double y);

double ini_rho_x(double x, double y);
double ini_rho_y(double x, double y);
double ini_rho_xx(double x, double y);
double ini_rho_xy(double x, double y);
double ini_rho_yy(double x, double y);

/**
 * @brief ��ʼ�ٶȳ���x�������
 * 
 * @param x ���������
 * @param y ����������
 * @return double 
 */
double ini_ux(double x, double y);

/**
 * @brief ��ʼ�ٶȳ���x�󵼵�x���������漸��������ƫ��
 * 
 * @param x ���������
 * @param y ����������
 * @return double 
 */
double ini_ux_x(double x, double y);

double ini_ux_y(double x, double y);

double ini_ux_xx(double x, double y);

double ini_ux_xy(double x, double y);

double ini_ux_yy(double x, double y);

/**
 * @brief ��ʼ�ٶȳ�y�������
 * 
 * @param x ���������
 * @param y ����������
 * @return double 
 */
double ini_uy(double x, double y);

/**
 * @brief ��ʼ�ٶȳ���y�󵼵�x����
 * 
 * @param x ���������
 * @param y ����������
 * @return double 
 */
double ini_uy_x(double x, double y);

double ini_uy_y(double x, double y);

double ini_uy_xx(double x, double y);

double ini_uy_xy(double x, double y);

double ini_uy_yy(double x, double y);

/**
 * @brief ��ʼѹ����
 * 
 * @param x ���������
 * @param y ����������
 * @return double 
 */
double ini_p(double x, double y);

double ini_p_x(double x, double y);

double ini_p_y(double x, double y);

double ini_p_xx(double x, double y);

double ini_p_xy(double x, double y);

double ini_p_yy(double x, double y);

/**
 * @brief ���ܽ�����
 * 
 * @param x ���������
 * @param y ����������
 * @param t ʱ��
 * @return double 
 */
double ana_e(double x, double y, double t);
/** @} */


#endif
/**
 * @file nodal_solver.h
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief ����ڵ�ⷨ��
 * @version v1.01
 * @date 2022-12-01
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#ifndef NODAL_SOLVER_H
#define NODAL_SOLVER_H

const double vn = 0;    /**< �ⷨ���ٶȣ�������֪)*/
const double Pre = 0;   /**< �߽�ѹ��*/

/**
 * @brief �����(i,j)���ڵ㴦ÿ�����ڵ�Ԫ�ڸõ���ٶ��ع�ֵ��ƽ��
 * 
 * @param i �ڵ���б��
 * @param j �ڵ���б��
 * @return double* 
 */
double * u_average(int i, int j);

/**
 * @brief �����(i,j)����������SMS�ܶȵ�������rho_s - rho_nus
 * 
 * @param i ��������б��
 * @param j ��������б��
 * @return double 
 */
double density_SMS(int i, int j);

/**
 * @brief ������������ĵ���ٶȺ�corner force
 * 
 * @param i �ڵ��б��
 * @param j �ڵ��б��
 * @return double** 
 */
void ustar_f_center(int i, int j);

/**
 * @brief ����ڵ��ٶȺ�corner force
 * 
 * @param i �ڵ���б��
 * @param j �ڵ���б��
 * @return double **
 */
void nodal_velocity_corner_force(int i, int j);

/**
 * @brief ����߽��ٶ�
 * 
 * @param i �ڵ��б��
 * @param j �ڵ��б��
 * @param amu ���߽�������ѹ��ʱ���������ֵ
 */
void BCvelocity(int i, int j, double amu);

/**
 * @brief ��֪�ⷨ���ٶ�
 * 
 * @param i �ڵ��б��
 * @param j �ڵ��б��
 * @param vn �ⷨ���ٶȲ�Ϊ0ʱ��������
 */
void BCnormal_velocity(int i, int j, double vn);
/**
 * @brief ��֪�߽�ѹ��
 * 
 * @param i �ڵ��б��
 * @param j �ڵ��б��
 * @param pcorlast ��һ�ڵ��\sum(a mu)
 * @param pcornext ��һ�ڵ��\sum(a mu)
 */
void BCpressure(int i, int j, double pcorlast, double pcornext);

#endif
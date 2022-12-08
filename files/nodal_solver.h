/**
 * @file nodal_solver.h
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief 定义节点解法器
 * @version v1.01
 * @date 2022-12-01
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#ifndef NODAL_SOLVER_H
#define NODAL_SOLVER_H

const double vn = 0;    /**< 外法向速度（假设已知)*/
const double Pre = 0;   /**< 边界压力*/

/**
 * @brief 计算第(i,j)个节点处每个相邻单元在该点的速度重构值的平均
 * 
 * @param i 节点的行编号
 * @param j 节点的列编号
 * @return double* 
 */
double * u_average(int i, int j);

/**
 * @brief 计算第(i,j)个子网格中SMS密度的修正量rho_s - rho_nus
 * 
 * @param i 子网格的行编号
 * @param j 子网格的列编号
 * @return double 
 */
double density_SMS(int i, int j);

/**
 * @brief 计算大网格中心点的速度和corner force
 * 
 * @param i 节点行编号
 * @param j 节点列编号
 * @return double** 
 */
void ustar_f_center(int i, int j);

/**
 * @brief 计算节点速度和corner force
 * 
 * @param i 节点的行编号
 * @param j 节点的列编号
 * @return double **
 */
void nodal_velocity_corner_force(int i, int j);

/**
 * @brief 计算边界速度
 * 
 * @param i 节点行编号
 * @param j 节点列编号
 * @param amu 当边界条件是压力时输入的修正值
 */
void BCvelocity(int i, int j, double amu);

/**
 * @brief 已知外法向速度
 * 
 * @param i 节点行编号
 * @param j 节点列编号
 * @param vn 外法向速度不为0时的修正项
 */
void BCnormal_velocity(int i, int j, double vn);
/**
 * @brief 已知边界压力
 * 
 * @param i 节点行编号
 * @param j 节点列编号
 * @param pcorlast 上一节点的\sum(a mu)
 * @param pcornext 下一节点的\sum(a mu)
 */
void BCpressure(int i, int j, double pcorlast, double pcornext);

#endif
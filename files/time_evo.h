/**
 * @file time_evo.h
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief 时间更新相关函数
 * @version 0.1
 * @date 2022-12-03
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#ifndef TIME_EVO_H
#define TIME_EVO_H

/**
 * @brief RK第一步
 * 
 * @param vn n时刻的值
 * @param Rn n时刻右端项
 * @param Minv 系数矩阵的逆M·dv/dt = R
 * @param delta_t 时间步长
 * @param len 向量维数
 * @return double* 
 */
double * RK_step1(double * vn, double * Rn,
                  double ** Minv, double delta_t, int len);

/**
 * @brief RK第二步
 * 
 * @param vn n时刻值
 * @param vs1 s1步得出的值
 * @param Rs1 s1时刻右端项
 * @param Minv 系数矩阵的逆
 * @param delta_t 时间步长
 * @param len 向量维数
 * @return double* 
 */
double * RK_step2(double * vn, double * vs1, double * Rs1,
                  double ** Minv, double delta_t, int len);

/**
 * @brief RK第三步
 * 
 * @param vn n时刻值
 * @param vs2 s2时刻值
 * @param Rs2 s2时刻右端项
 * @param Minv 系数矩阵的逆 
 * @param delta_t 时间步长
 * @param len 向量维数
 * @return double* 
 */
double * RK_step3(double * vn, double * vs2, double * Rs2,
                  double ** Minv, double delta_t, int len);

/**
 * @brief 进行一次时间更新
 * 
 * @param dt 时间步长
 */
void one_time_step(double dt);

/**
 * @brief 更新几何量
 * 
 */
void update_geo();

/**
 * @brief 选择时间步长
 * 
 * @param dtn 原本的时间步长
 * 
 * @return double 
 */
double choosedt(double dtn);
#endif
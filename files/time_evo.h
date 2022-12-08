/**
 * @file time_evo.h
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief ʱ�������غ���
 * @version 0.1
 * @date 2022-12-03
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#ifndef TIME_EVO_H
#define TIME_EVO_H

/**
 * @brief RK��һ��
 * 
 * @param vn nʱ�̵�ֵ
 * @param Rn nʱ���Ҷ���
 * @param Minv ϵ���������M��dv/dt = R
 * @param delta_t ʱ�䲽��
 * @param len ����ά��
 * @return double* 
 */
double * RK_step1(double * vn, double * Rn,
                  double ** Minv, double delta_t, int len);

/**
 * @brief RK�ڶ���
 * 
 * @param vn nʱ��ֵ
 * @param vs1 s1���ó���ֵ
 * @param Rs1 s1ʱ���Ҷ���
 * @param Minv ϵ���������
 * @param delta_t ʱ�䲽��
 * @param len ����ά��
 * @return double* 
 */
double * RK_step2(double * vn, double * vs1, double * Rs1,
                  double ** Minv, double delta_t, int len);

/**
 * @brief RK������
 * 
 * @param vn nʱ��ֵ
 * @param vs2 s2ʱ��ֵ
 * @param Rs2 s2ʱ���Ҷ���
 * @param Minv ϵ��������� 
 * @param delta_t ʱ�䲽��
 * @param len ����ά��
 * @return double* 
 */
double * RK_step3(double * vn, double * vs2, double * Rs2,
                  double ** Minv, double delta_t, int len);

/**
 * @brief ����һ��ʱ�����
 * 
 * @param dt ʱ�䲽��
 */
void one_time_step(double dt);

/**
 * @brief ���¼�����
 * 
 */
void update_geo();

/**
 * @brief ѡ��ʱ�䲽��
 * 
 * @param dtn ԭ����ʱ�䲽��
 * 
 * @return double 
 */
double choosedt(double dtn);
#endif
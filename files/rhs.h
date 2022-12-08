/**
 * @file rhs.h
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief �����Ҷ˾���
 * @version 0.1
 * @date 2022-12-03
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#ifndef RHS_H
#define RHS_H

/**
 * @brief �õ����ڴ������еľֲ�����
 * 
 * @param pq ��ı��
 * @param oq ������ı��
 * @return double* 
 */
double * get_xieta(int pq, int oq);

/**
 * @brief Get the subcell in vertexneighbor
 * 
 * @param sq ��������
 * @param pq �ڵ���
 * @return int 
 */
int get_subcell_in_vertexneighbor(int sq, int pq);

/**
 * @brief �õ���Ԫ���ݵ��Ҷ�����
 * 
 * @param i ��������б��
 * @param j ��������б��
 * @return double* 
 */
double * specific_volume_matrix(int i, int j);

/**
 * @brief �õ���Ԫ�ٶȵ��Ҷ˾���
 * 
 * @param i �������б��
 * @param j �������б��
 * @return double** 
 */
double ** velocity_matrix(int i, int j);

/**
 * @brief �õ���Ԫ�������Ҷ�����
 * 
 * @param i ��������б��
 * @param j ��������б��
 * @return double* 
 */
double * energy_matrix(int i, int j);

#endif
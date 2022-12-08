/**
 * @file initial.h
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief �������񲢼�¼������Ϣ
 * @version 0.1
 * @date 2022-11-29
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#ifndef INITIAL_H
#define INITIAL_H

/**
 * @brief ���ɳ�ʼ���񻮷֣�
 * - �����ڴ�
 * - ��¼����ı��
 * - ��¼�ڵ�ı��
 * - ��¼�ڵ���������
 * - ��¼���񶥵��š�������������
 * - ��¼�������ţ����ڴ������ţ�������
 */
void generatemesh();

/**
 * @brief �Ե�(i,j)����ԪѰ�����ڵ�Ԫ��
 * ��(i,j-1)��ʼ��ʱ������һ�����ڵĵ�Ԫ��Ϊ��һ���
 * �洢һά��š�
 * 
 * @param i ����Ķ�ά�б��
 * @param j ����Ķ�ά�б��
 */
void element_findneighbor(int i, int j);

/**
 * @brief 
 * - �Ե�(i,j)���ڵ�Ѱ�����ڵ�����
 * ��(i,j)��������ʱ������һ�����ڵ�����ʼ��¼����ʱ��˳��
 * - �Ե�(i,j)���ڵ�Ѱ�����ڵĽڵ㣬
 * ��(i-1,j)��ʼ��ʱ�������һ�����ڵĽڵ㿪ʼ��¼����ʱ��˳��
 * - ������ߵģ����֣��߳����ⷨ����
 * 
 * @param i �ڵ���б��
 * @param j �ڵ���б��
 */
void node_findneighbor(int i, int j);

/**
 * @brief �Ե�(i,j)���ڵ��ʼ��ÿ����������ÿ�����ϵļ�Ȩ�ⷨ������corner force
 * 
 * @param i �ڵ��б��
 * @param j �ڵ��б��
 */
void node_initialsegment(int i, int j);

/**
 * @brief �����(i,j)����������ĺ�����
 * 
 * @param i �����б��
 * @param j �����б��
 * @return double* ���ĵĲο���Ԫ����
 * @note ���ص���ָ�룬�洢ʱ�Ǳ������ǵ�ɾ��ָ������ڴ�
 */
double * mass_center(int i, int j);

/**
 * @brief �����(i,j)���������������
 * 
 * @param i �����б��
 * @param j �����б��
 * @return double** ��������
 */
double ** mass_matrix(int i, int j);

/**
 * @brief �����(i,j)����������������������ͳ�ʼ���
 * 
 * @param i ��������б��
 * @param j ��������б��
 */
void subcell_mass_in_ij(int i, int j);

/**
 * @brief ��ʼ�������������������񲢼�¼��ʼ��Ϣ
 * 
 */
void initial();


#endif
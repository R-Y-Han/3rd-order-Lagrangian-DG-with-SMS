/**
 * @file rhs.h
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief 计算右端矩阵
 * @version 0.1
 * @date 2022-12-03
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#ifndef RHS_H
#define RHS_H

/**
 * @brief 得到点在大网格中的局部坐标
 * 
 * @param pq 点的编号
 * @param oq 大网格的编号
 * @return double* 
 */
double * get_xieta(int pq, int oq);

/**
 * @brief Get the subcell in vertexneighbor
 * 
 * @param sq 子网格编号
 * @param pq 节点编号
 * @return int 
 */
int get_subcell_in_vertexneighbor(int sq, int pq);

/**
 * @brief 得到单元比容的右端向量
 * 
 * @param i 大网格的行编号
 * @param j 大网格的列编号
 * @return double* 
 */
double * specific_volume_matrix(int i, int j);

/**
 * @brief 得到单元速度的右端矩阵
 * 
 * @param i 大网格行编号
 * @param j 大网格列编号
 * @return double** 
 */
double ** velocity_matrix(int i, int j);

/**
 * @brief 得到单元的总能右端向量
 * 
 * @param i 大网格的行编号
 * @param j 大网格的列编号
 * @return double* 
 */
double * energy_matrix(int i, int j);

#endif
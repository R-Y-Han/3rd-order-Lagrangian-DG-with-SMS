/**
 * @file config.h
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief Configuration for some definitions.
 * @version 1.01
 * @date 2022-11-28
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#ifndef CONFIG_H
#define CONFIG_H

#include <vector>
using namespace std;

/**
 * @name parameters for computation
 * @{
 */

const int pk = 6;   /**< 基函数个数*/

const int n_element = 5;   /**< 横向网格个数*/
const int m_element = 5;   /**< 纵向网格个数*/
const int n_point = 2 * n_element;   /**< 横向节点个数*/
const int m_point = 2 * m_element;   /**< 纵向节点个数*/

const double T = 0.2;   /**< 终止时间*/

enum test_cases{
    shockless_Noh = 1,
    Taylor_Green_vortex = 2
};  /**< 算例*/

const test_cases testcase = shockless_Noh;  /**< 算例*/

const double C_cfl = 0.1;   /**< CFL参数*/
const double C_m = 1.01;    /**< 时间步长相关参数*/

const double X = 1; /**< 计算区域横坐标大小*/
const double Y = 1; /**< 计算区域纵坐标大小*/

const double hx = X / (double ) n_point;  /**< 横向空间步长*/
const double hy = Y / (double ) m_point;  /**< 纵向空间步长*/



const double c_CFL = 0.1;   /**< CFL条件系数*/

const double nu_alpha = 0.6;  /**< 比容限制器参数*/
const double ux_alpha = 0.6;  /**< 横向速度限制器参数*/
const double uy_alpha = 0.6;  /**< 纵向速度限制器参数*/
const double tau_alpha =0.6; /**< 总能量限制参数*/

//先四个顶点，再四条边，再中点，从下往上，从左往右
const double ref_xi[9] = {-1, 1, 1, -1,
                          0,  1, 0, -1,
                          0};   /**< 参考空间九点横坐标，注意排序*/

const double ref_eta[9] = {-1, -1, 1, 1,
                           -1,  0, 1, 0,
                           0};  /**< 参考空间九点纵坐标，注意排序*/

/** @} 计算区域参数*/

/**
 * @name 流体参数和变量
 * @{
 */

/**
 * @brief EOS
 * 
 * @param rho 密度
 * @param e 内能
 * @return double 压力
 */
double EOS(double gamma, double rho, double e);   /**< EOS*/

/**
 * @struct Subcell
 * @brief 子网格定义
 * 
 */
struct SubElement{
    int q;  /**< 子网格编号，和左下角节点编号相同*/
    int father_element; /**< 所在大网格编号*/
    int vertex[4];  /**< 顶点编号*/
    double m_s; /**< 子网格质量*/
    double area_s;  /**< 子网格面积*/
    double avg_den; /**< 子网格平均密度*/
};

/**
 * @struct Omega
 * @brief 网格类定义
 * 
 */
struct Omega{
    int q;  /**< 网格编号，q = m*i + j*/
    int vertex[9];  /**< 顶点编号*/
    double vx[9];   /**< 顶点横坐标*/
    double vy[9];   /**< 顶点纵坐标*/
    double vx0[9];  /**< 顶点初始横坐标*/
    double vy0[9];  /**< 顶点初始纵坐标*/
    double xi_c;    /**< 质心参考横坐标*/
    double eta_c;   /**< 质心参考纵坐标*/
    double mass;    /**< 质量*/
    double gamma; /**< 定容比热*/

    double psi3avg, psi4avg, psi5avg;   /**< 基函数345的质量平均，减少运算量*/

    double * nu;    /**< 比容中心量*/
    double * ux;    /**< x速度分量中心量*/
    double * uy;    /**< y速度分量中心量*/
    double * tau;   /**< 总能中心量*/
    double c_center;    /**< 中心等熵声速*/

    double ** M;    /**< mass matrix*/
    double ** M_inv;    /**< inverse of mass matrix*/
    double rho_0;   /**< 初始中心密度*/

    vector<int> neighbor_element;   /**< 相邻单元编号*/

    int subcell[4]; /**< 子网格编号*/

    double phi_x(double xi, double eta);
    double phi_y(double xi, double eta);    /**< 参考空间到物理空间的映射*/
    double phi_x0(double xi, double eta);
    double phi_y0(double xi, double eta);   /**< 参考空间到初始物理空间的映射*/
    double Psi(int i, double xi, double eta);   /**< 基函数*/
    double Psi_xi(int i, double xi, double eta);    /**< 基函数对xi求导*/
    double Psi_eta(int i, double xi, double eta);    /**< 基函数对eta求导*/
    double ** getJacobiMatrix(double xi, double eta);   /**< Jacobi矩阵*/
    double ** getJacobiMatrix_0(double xi, double eta); /**< 初始Jacobi矩阵*/
    double detJacobi(double xi, double eta);    /**< Jacobi行列式*/
    double detJacobi_0(double xi, double eta);  /**< Jacobi初始行列式*/
    double edgedetJacobi(double xi, double eta, int pointa, int pointb);    /**< 边上三点插值的Jacobi行列式*/
};

/**
 * @struct node
 * @brief 网格节点定义
 * 
 *                      ·   last
 *                      |
 *           anlast  <==|   (subcell k)
 *                      |
 *          point (i,j) · —— —— —— · next
 *                            ||
 *                            \/
 *                          annext
 */
struct node{
    int q;  /**< 编号，q = (m_point + 1)*i + j*/
    int boundary;   /**< 是否为边界节点*/
    int ifcenter;   /**< 是否为大网格的中点*/
    int ifedge;     /**< 是否为边上的点*/
    double x, y;    /**< 节点物理坐标*/
    double x0, y0;  /**< 节点初始物理坐标*/
    double  upstarx, upstary; /**< 节点速度*/
    vector<int > neighbor_subcell;  /**< 相邻子网格编号*/
    vector<int > neighbor_node; /**< 相邻节点*/
    vector<double > weightnext;        /**< 第k个子网格内到下一点的积分权重*/
    vector<double > weightlast;        /**< 第k个子网格内到上一点的积分权重*/
    vector<int > segmentnext;   /**< 第k个相邻子网格内下一个顶点*/
    vector<int > segmentlast;   /**< 第k个子网格内上一个顶点*/
    vector<double > anext; /**< 第k个相邻子网格内到下一个顶点的加权外法向量*/
    vector<double > alast; /**< 第k个相邻子网格内到上一个顶点的加权外法向量*/
    vector<double* > Nnext;  /**< 参考单元上第k个相邻子网格内到下一个节点的单位外法向量*/
    vector<double* > Nlast;  /**< 参考单元上第k个相邻子网格内到上一个节点的单位外法向量*/
    vector<double* > nnext;  /**< 物理空间...的单位外法向量*/
    vector<double* > nlast;  /**< 物理空间...的单位外法向量*/
    vector<double* > flast; /**< 第k个相邻子网格内到下一个节点的corner force*/
    vector<double* > fnext; /**< 第k个相邻子网格内到上一个节点的corner force*/ 
};

extern SubElement subcell[n_point][m_point]; /**< 生成n_point*m_point个子网格*/
extern Omega o[n_element][m_element];   /**< 生成网格*/
extern node point[n_point+1][m_point+1];  /**< 生成节点*/

/** @}*/

/**
 * @name 数学求解工具
 * @{
 */

const double PI = 3.1415927;    /**< 圆周率*/

const int gpn = 16; /**< Gauss求积节点个数*/
extern double Gausspoint_xi[gpn];   /**< Gauss求积节点横坐标*/
extern double Gausspoint_eta[gpn];  /**< Gauss求积节点纵坐标*/
extern double Gaussweight[gpn]; /**< Gauss求积节点权重*/

/**
 * @brief 形函数
 * 
 * @param i i = 0,1,2,3,4,5,6,7,8,9
 * @param xi 参考空间横坐标
 * @param eta 参考空间纵坐标
 * @return 第i个形函数值
 * 
 */
double bp(int i, double xi, double eta);

/**
 * @brief 形函数对xi求导
 * 
 * @param i i = 0,1,2,3,4,5,6,7,8,9
 * @param xi 参考空间横坐标
 * @param eta 参考空间纵坐标
 * @return 第i个形函数导数值
 */
double bp_xi(int i, double xi, double eta);

/**
 * @brief 形函数对eta求导
 * 
 * @param i i = 0,1,2,3,4,5,6,7,8,9
 * @param xi 参考空间横坐标
 * @param eta 参考空间纵坐标
 * @return 第i个形函数倒数值
 */
double bp_eta(int i, double xi, double eta);

double bp_xixi(int i, double xi, double eta);

double bp_etaeta(int i, double xi, double eta);

double bp_xieta(int i, double xi, double eta);

/**
 * @brief 求解M的逆矩阵
 * 
 * @param M 需要求解的逆矩阵
 * @param len 矩阵大小
 * @return double** 
 */
double ** M_inverse(double ** M, int len);
/**@}*/


#endif
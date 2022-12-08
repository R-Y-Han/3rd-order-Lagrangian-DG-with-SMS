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

const int pk = 6;   /**< ����������*/

const int n_element = 5;   /**< �����������*/
const int m_element = 5;   /**< �����������*/
const int n_point = 2 * n_element;   /**< ����ڵ����*/
const int m_point = 2 * m_element;   /**< ����ڵ����*/

const double T = 0.2;   /**< ��ֹʱ��*/

enum test_cases{
    shockless_Noh = 1,
    Taylor_Green_vortex = 2
};  /**< ����*/

const test_cases testcase = shockless_Noh;  /**< ����*/

const double C_cfl = 0.1;   /**< CFL����*/
const double C_m = 1.01;    /**< ʱ�䲽����ز���*/

const double X = 1; /**< ��������������С*/
const double Y = 1; /**< ���������������С*/

const double hx = X / (double ) n_point;  /**< ����ռ䲽��*/
const double hy = Y / (double ) m_point;  /**< ����ռ䲽��*/



const double c_CFL = 0.1;   /**< CFL����ϵ��*/

const double nu_alpha = 0.6;  /**< ��������������*/
const double ux_alpha = 0.6;  /**< �����ٶ�����������*/
const double uy_alpha = 0.6;  /**< �����ٶ�����������*/
const double tau_alpha =0.6; /**< ���������Ʋ���*/

//���ĸ����㣬�������ߣ����е㣬�������ϣ���������
const double ref_xi[9] = {-1, 1, 1, -1,
                          0,  1, 0, -1,
                          0};   /**< �ο��ռ�ŵ�����꣬ע������*/

const double ref_eta[9] = {-1, -1, 1, 1,
                           -1,  0, 1, 0,
                           0};  /**< �ο��ռ�ŵ������꣬ע������*/

/** @} �����������*/

/**
 * @name ��������ͱ���
 * @{
 */

/**
 * @brief EOS
 * 
 * @param rho �ܶ�
 * @param e ����
 * @return double ѹ��
 */
double EOS(double gamma, double rho, double e);   /**< EOS*/

/**
 * @struct Subcell
 * @brief ��������
 * 
 */
struct SubElement{
    int q;  /**< �������ţ������½ǽڵ�����ͬ*/
    int father_element; /**< ���ڴ�������*/
    int vertex[4];  /**< ������*/
    double m_s; /**< ����������*/
    double area_s;  /**< ���������*/
    double avg_den; /**< ������ƽ���ܶ�*/
};

/**
 * @struct Omega
 * @brief �����ඨ��
 * 
 */
struct Omega{
    int q;  /**< �����ţ�q = m*i + j*/
    int vertex[9];  /**< ������*/
    double vx[9];   /**< ���������*/
    double vy[9];   /**< ����������*/
    double vx0[9];  /**< �����ʼ������*/
    double vy0[9];  /**< �����ʼ������*/
    double xi_c;    /**< ���Ĳο�������*/
    double eta_c;   /**< ���Ĳο�������*/
    double mass;    /**< ����*/
    double gamma; /**< ���ݱ���*/

    double psi3avg, psi4avg, psi5avg;   /**< ������345������ƽ��������������*/

    double * nu;    /**< ����������*/
    double * ux;    /**< x�ٶȷ���������*/
    double * uy;    /**< y�ٶȷ���������*/
    double * tau;   /**< ����������*/
    double c_center;    /**< ���ĵ�������*/

    double ** M;    /**< mass matrix*/
    double ** M_inv;    /**< inverse of mass matrix*/
    double rho_0;   /**< ��ʼ�����ܶ�*/

    vector<int> neighbor_element;   /**< ���ڵ�Ԫ���*/

    int subcell[4]; /**< ��������*/

    double phi_x(double xi, double eta);
    double phi_y(double xi, double eta);    /**< �ο��ռ䵽����ռ��ӳ��*/
    double phi_x0(double xi, double eta);
    double phi_y0(double xi, double eta);   /**< �ο��ռ䵽��ʼ����ռ��ӳ��*/
    double Psi(int i, double xi, double eta);   /**< ������*/
    double Psi_xi(int i, double xi, double eta);    /**< ��������xi��*/
    double Psi_eta(int i, double xi, double eta);    /**< ��������eta��*/
    double ** getJacobiMatrix(double xi, double eta);   /**< Jacobi����*/
    double ** getJacobiMatrix_0(double xi, double eta); /**< ��ʼJacobi����*/
    double detJacobi(double xi, double eta);    /**< Jacobi����ʽ*/
    double detJacobi_0(double xi, double eta);  /**< Jacobi��ʼ����ʽ*/
    double edgedetJacobi(double xi, double eta, int pointa, int pointb);    /**< ���������ֵ��Jacobi����ʽ*/
};

/**
 * @struct node
 * @brief ����ڵ㶨��
 * 
 *                      ��   last
 *                      |
 *           anlast  <==|   (subcell k)
 *                      |
 *          point (i,j) �� ���� ���� ���� �� next
 *                            ||
 *                            \/
 *                          annext
 */
struct node{
    int q;  /**< ��ţ�q = (m_point + 1)*i + j*/
    int boundary;   /**< �Ƿ�Ϊ�߽�ڵ�*/
    int ifcenter;   /**< �Ƿ�Ϊ��������е�*/
    int ifedge;     /**< �Ƿ�Ϊ���ϵĵ�*/
    double x, y;    /**< �ڵ���������*/
    double x0, y0;  /**< �ڵ��ʼ��������*/
    double  upstarx, upstary; /**< �ڵ��ٶ�*/
    vector<int > neighbor_subcell;  /**< ������������*/
    vector<int > neighbor_node; /**< ���ڽڵ�*/
    vector<double > weightnext;        /**< ��k���������ڵ���һ��Ļ���Ȩ��*/
    vector<double > weightlast;        /**< ��k���������ڵ���һ��Ļ���Ȩ��*/
    vector<int > segmentnext;   /**< ��k����������������һ������*/
    vector<int > segmentlast;   /**< ��k������������һ������*/
    vector<double > anext; /**< ��k�������������ڵ���һ������ļ�Ȩ�ⷨ����*/
    vector<double > alast; /**< ��k�������������ڵ���һ������ļ�Ȩ�ⷨ����*/
    vector<double* > Nnext;  /**< �ο���Ԫ�ϵ�k�������������ڵ���һ���ڵ�ĵ�λ�ⷨ����*/
    vector<double* > Nlast;  /**< �ο���Ԫ�ϵ�k�������������ڵ���һ���ڵ�ĵ�λ�ⷨ����*/
    vector<double* > nnext;  /**< ����ռ�...�ĵ�λ�ⷨ����*/
    vector<double* > nlast;  /**< ����ռ�...�ĵ�λ�ⷨ����*/
    vector<double* > flast; /**< ��k�������������ڵ���һ���ڵ��corner force*/
    vector<double* > fnext; /**< ��k�������������ڵ���һ���ڵ��corner force*/ 
};

extern SubElement subcell[n_point][m_point]; /**< ����n_point*m_point��������*/
extern Omega o[n_element][m_element];   /**< ��������*/
extern node point[n_point+1][m_point+1];  /**< ���ɽڵ�*/

/** @}*/

/**
 * @name ��ѧ��⹤��
 * @{
 */

const double PI = 3.1415927;    /**< Բ����*/

const int gpn = 16; /**< Gauss����ڵ����*/
extern double Gausspoint_xi[gpn];   /**< Gauss����ڵ������*/
extern double Gausspoint_eta[gpn];  /**< Gauss����ڵ�������*/
extern double Gaussweight[gpn]; /**< Gauss����ڵ�Ȩ��*/

/**
 * @brief �κ���
 * 
 * @param i i = 0,1,2,3,4,5,6,7,8,9
 * @param xi �ο��ռ������
 * @param eta �ο��ռ�������
 * @return ��i���κ���ֵ
 * 
 */
double bp(int i, double xi, double eta);

/**
 * @brief �κ�����xi��
 * 
 * @param i i = 0,1,2,3,4,5,6,7,8,9
 * @param xi �ο��ռ������
 * @param eta �ο��ռ�������
 * @return ��i���κ�������ֵ
 */
double bp_xi(int i, double xi, double eta);

/**
 * @brief �κ�����eta��
 * 
 * @param i i = 0,1,2,3,4,5,6,7,8,9
 * @param xi �ο��ռ������
 * @param eta �ο��ռ�������
 * @return ��i���κ�������ֵ
 */
double bp_eta(int i, double xi, double eta);

double bp_xixi(int i, double xi, double eta);

double bp_etaeta(int i, double xi, double eta);

double bp_xieta(int i, double xi, double eta);

/**
 * @brief ���M�������
 * 
 * @param M ��Ҫ���������
 * @param len �����С
 * @return double** 
 */
double ** M_inverse(double ** M, int len);
/**@}*/


#endif
/**
 * @file config.cpp
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief config.h中函数定义
 * @version 1.01
 * @date 2022-11-28
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#include <cmath>
#include "config.h"

/**
 * @brief EOS
 * 
 * @param rho 密度
 * @param e 内能
 * @return double 压力
 */
double EOS(double gamma, double rho, double e)
{
    double p;
    p = (gamma - 1) * rho * e;
    return p;
}

//************************//
//************************//
//************************//
//************************//
//************************//
//************************//

double Omega::phi_x(double xi, double eta)
{
    double xtemp = 0;
    for (int i=0; i<9; i++)
    {
        xtemp = xtemp + this->vx[i] * bp(i,xi,eta);
    }
    return xtemp;
}

double Omega::phi_y(double xi, double eta)
{
    double ytemp = 0;
    for (int i=0; i<9; i++)
    {
        ytemp = ytemp + this->vy[i] * bp(i,xi,eta);
    }
    return ytemp;
}

double Omega::phi_x0(double xi, double eta)
{
    double xtemp = 0;
    for (int i=0; i<9; i++)
    {
        xtemp = xtemp + this->vx0[i] * bp(i,xi,eta);
    }
    return xtemp;
}

double Omega::phi_y0(double xi, double eta)
{
    double ytemp = 0;
    for (int i=0; i<9; i++)
    {
        ytemp = ytemp + this->vy0[i] * bp(i,xi,eta);
    }
    return ytemp;
}

double Omega::Psi(int i, double xi, double eta)
{
    double ans;
    double xic, etac;
    xic = this->xi_c;
    etac = this->eta_c;
    switch(i)
    {
        case 0:
            ans = 1;
            break;
        case 1:
            ans = xi - xic;
            break;
        case 2:
            ans = eta - etac;
            break;
        case 3:
            ans = (xi - xic)*(xi - xic) / 2.0 - this->psi3avg;
            break;
        case 4:
            ans = (eta - etac)*(eta - etac) / 2.0 - this->psi4avg;
            break;
        case 5:
            ans = (xi - xic)*(eta - etac) - this->psi5avg;
            break;
        default:
            ans = 0;
            break;
    }
    return ans;
}

double Omega::Psi_xi(int i, double xi, double eta)
{
    double ans;
    double xic, etac;
    xic = this->xi_c;
    etac = this->eta_c;
    switch(i)
    {
        case 0:
            ans = 0;
            break;
        case 1:
            ans = 1;
            break;
        case 2:
            ans = 0;
            break;
        case 3:
            ans = xi - xic;
            break;
        case 4:
            ans = 0;
            break;
        case 5:
            ans = eta - etac;
            break;
        default:
            ans = 0;
            break;
    }
    return ans;
}

double Omega::Psi_eta(int i, double xi, double eta)
{
    double ans;
    double xic, etac;
    xic = this->xi_c;
    etac = this->eta_c;
    switch(i)
    {
        case 0:
            ans = 0;
            break;
        case 1:
            ans = 0;
            break;
        case 2:
            ans = 1;
            break;
        case 3:
            ans = 0;
            break;
        case 4:
            ans = eta - etac;
            break;
        case 5:
            ans = xi - xic;
            break;
        
        default:
            ans = 0;
            break;
    }
    return ans;
}

double ** Omega::getJacobiMatrix(double xi, double eta)
{
    double ** J = new double * [2];
    J[0] = new double [2];
    J[1] = new double [2];
    int i;
    J[0][0] = 0;
    J[0][1] = 0;
    J[1][0] = 0;
    J[1][1] = 0;
    for (i=0; i<9; i++)
    {
        J[0][0] = J[0][0] + this->vx[i] * bp_xi(i,xi,eta);
        J[0][1] = J[0][1] + this->vx[i] * bp_eta(i,xi,eta);
        J[1][0] = J[1][0] + this->vy[i] * bp_xi(i,xi,eta);
        J[1][1] = J[1][1] + this->vy[i] * bp_eta(i,xi,eta);
    }
    return J;
}

double ** Omega::getJacobiMatrix_0(double xi, double eta)
{
    double ** J_0 = new double * [2];
    J_0[0] = new double [2];
    J_0[1] = new double [2];
    int i;
    J_0[0][0] = 0;
    J_0[0][1] = 0;
    J_0[1][0] = 0;
    J_0[1][1] = 0;
    for (i=0; i<9; i++)
    {
        J_0[0][0] = J_0[0][0] + this->vx0[i] * bp_xi(i,xi,eta);
        J_0[0][1] = J_0[0][1] + this->vx0[i] * bp_eta(i,xi,eta);
        J_0[1][0] = J_0[1][0] + this->vy0[i] * bp_xi(i,xi,eta);
        J_0[1][1] = J_0[1][1] + this->vy0[i] * bp_eta(i,xi,eta);
    }
    return J_0;
}

double Omega::detJacobi(double xi, double eta)
{
    double ans;
    double a, b, c, d;
    
    a = 0;
    b = 0;
    c = 0;
    d = 0;
    for ( int i = 0; i < 9; i++)
    {
        a = a + this->vx[i] * bp_xi(i,xi,eta);
        b = b + this->vx[i] * bp_eta(i,xi,eta);
        c = c + this->vy[i] * bp_xi(i,xi,eta);
        d = d + this->vy[i] * bp_eta(i,xi,eta);
    }
    ans = a * d - b * c;
    return ans;
}

double Omega::detJacobi_0(double xi, double eta)
{
    double ans;
    double a, b, c, d;
    
    a = 0;
    b = 0;
    c = 0;
    d = 0;
    for ( int i = 0; i < 9; i++)
    {
        a = a + this->vx0[i] * bp_xi(i,xi,eta);
        b = b + this->vx0[i] * bp_eta(i,xi,eta);
        c = c + this->vy0[i] * bp_xi(i,xi,eta);
        d = d + this->vy0[i] * bp_eta(i,xi,eta);
    }
    ans = a * d - b * c;
    return ans;
}

/**
 * @brief 计算某条边上三点插值的Jacobi插值行列式
 * 
 * @param xi 所求点的拉氏横坐标
 * @param eta 所求点的拉氏纵坐标
 * @param pointa 边上某一点的*全局编号*
 * @param pointb 边上另一点的*全局编号*
 * @return double 
 */
double Omega::edgedetJacobi(double xi, double eta, int pointa, int pointb)
{
    double ans, common;
    int p1, p2, p3;
    int a1, a2, a3;
    for (p1=0; p1<9; p1++)
    {
        if (this->vertex[p1] == pointa)
        {
            break;
        }
    }
    for (p2=0; p2<9; p2++)
    {
        if (this->vertex[p2] == pointb)
        {
            break;
        }
    }

    double j0, j1;
    //找出第三个插值点
    ans = 0;
    if (ref_eta[p1] == ref_eta[p2])
    {
        common = ref_eta[p1];
        //在同一条横线上
        for (p3=0; p3<9; p3++)
        {
            if (p3 == p1 || p3 == p2)
            {
                continue;
            }
            if (ref_eta[p3] == common)
            {
                break;
            }
        }

        if (ref_xi[p1] == -1)
        {
            a1 = p1;
            if (ref_xi[p2] == 0)
            {
                a2 = p2;
                a3 = p3;
            }
            else{
                a2 = p3;
                a3 = p2;
            }
        }
        else if(ref_xi[p1] == 0)
        {
            a2 = p1;
            if (ref_xi[p2] == -1)
            {
                a1 = p2;
                a3 = p3;
            }
            else{
                a1 = p3;
                a3 = p2;
            }
        }
        else{
            a3 = p1;
            if (ref_xi[p2] == -1)
            {
                a1 = p2;
                a2 = p3;
            }
            else{
                a1 = p3;
                a2 = p2;
            }
        }

        j0 = this->vx[a1] * (xi - 0.5)
            + this->vx[a2] * (- 2 * xi)
            + this->vx[a3] * (xi + 0.5);
        j1 = this->vy[a1] * (xi - 0.5)
            + this->vy[a2] * (-2 * xi)
            + this->vy[a3] * (xi + 0.5);
        ans = j0 * j0 + j1 * j1;
        ans = sqrt(ans);
    }
    else if (ref_xi[p1] == ref_xi[p2])
    {
        common = ref_xi[p1];
        //在同一条竖线上
        for (p3=0; p3<9; p3++)
        {
            if (p3 == p1 || p3 == p2)
            {
                continue;
            }
            if (ref_xi[p3] == common)
            {
                break;
            }
        }

        if (ref_eta[p1] == -1)
        {
            a1 = p1;
            if (ref_eta[p2] == 0)
            {
                a2 = p2;
                a3 = p3;
            }
            else{
                a2 = p3;
                a3 = p2;
            }
        }
        else if(ref_eta[p1] == 0)
        {
            a2 = p1;
            if (ref_eta[p2] == -1)
            {
                a1 = p2;
                a3 = p3;
            }
            else{
                a1 = p3;
                a3 = p2;
            }
        }
        else{
            a3 = p1;
            if (ref_eta[p2] == -1)
            {
                a1 = p2;
                a2 = p3;
            }
            else{
                a1 = p3;
                a2 = p2;
            }
        }

        j0 = this->vx[a1] * (eta - 0.5)
           + this->vx[a2] * (- 2 * eta)
           + this->vx[a3] * (eta + 0.5);
        j1 = this->vy[a1] * (eta - 0.5)
           + this->vy[a2] * (- 2 * eta)
           + this->vy[a3] * (eta + 0.5);
        ans = j0 * j0 + j1 * j1;
        ans = sqrt(ans);
    }
    else{
        return 0;
    }

    return ans;
}


SubElement subcell[n_point][m_point]; /**< 生成n_point*m_point个子网格*/
Omega o[n_element][m_element];   /**< 生成网格*/
node point[n_point+1][m_point+1];  /**< 生成节点*/


//************************//
//************************//

double pos1 = -0.8611363115940526;
double pos2 = -0.3399810435848563;
double pos3 = 0.3399810435848563;
double pos4 = 0.8611363115940526;
double w1 = 0.3478548451374538;
double w2 = 0.6521451548625461;
double w3 = 0.6521451548625461;
double w4 = 0.3478548451374538;
double Gausspoint_xi[gpn] = {pos1, pos2, pos3, pos4,
                             pos1, pos2, pos3, pos4,
                             pos1, pos2, pos3, pos4,
                             pos1, pos2, pos3, pos4};
double Gausspoint_eta[gpn] = {pos1, pos1, pos1, pos1,
                              pos2, pos2, pos2, pos2,
                              pos3, pos3, pos3, pos3,
                              pos4, pos4, pos4, pos4};
double Gaussweight[gpn] = {w1*w1, w2*w1, w3*w1, w4*w1,
                           w1*w2, w2*w2, w3*w2, w4*w2,
                           w1*w3, w2*w3, w3*w3, w4*w3,
                           w1*w4, w2*w4, w3*w4, w4*w4};

//************************//
//************************//
//************************//

double bp(int i, double xi, double eta)
{
    double ans;
    double xil, xim, xir;
    double etal, etam, etar;
    xil = xi * (xi-1) / 2.0;
    xim = -(xi + 1) * (xi - 1);
    xir = xi * (xi + 1) / 2.0;
    etal = eta * (eta - 1) / 2.0;
    etam = - (eta + 1) * (eta - 1);
    etar = eta * (eta + 1) / 2.0;
    switch(i)
    {
        case 0:
            ans = xil * etal;
            break;
        case 1:
            ans = xir * etal;
            break;
        case 2:
            ans = xir * etar;
            break;
        case 3:
            ans = xil * etar;
            break;
        case 4:
            ans = xim * etal;
            break;
        case 5:
            ans = xir * etam;
            break;
        case 6:
            ans = xim * etar;
            break;
        case 7:
            ans = xil * etam;
            break;
        case 8:
            ans = xim * etam;
            break;
        default:
            ans = 0;
            break;
    }
    return ans;
}

double bp_xi(int i, double xi, double eta)
{
    double ans;
    double xil, xim, xir;
    double etal, etam, etar;
    xil = xi - 1.0 / 2.0;
    xim = - 2 * xi;
    xir = xi + 1.0 / 2.0;
    etal = eta * (eta - 1) / 2.0;
    etam = - (eta + 1) * (eta - 1);
    etar = eta * (eta + 1) / 2.0;
    switch(i)
    {
        case 0:
            ans = xil * etal;
            break;
        case 1:
            ans = xir * etal;
            break;
        case 2:
            ans = xir * etar;
            break;
        case 3:
            ans = xil * etar;
            break;
        case 4:
            ans = xim * etal;
            break;
        case 5:
            ans = xir * etam;
            break;
        case 6:
            ans = xim * etar;
            break;
        case 7:
            ans = xil * etam;
            break;
        case 8:
            ans = xim * etam;
            break;
        default:
            ans = 0;
            break;
    }
    return ans;
}

double bp_eta(int i, double xi, double eta)
{
    double ans;
    double xil, xim, xir;
    double etal, etam, etar;
    xil = xi * (xi-1) / 2.0;
    xim = -(xi + 1) * (xi - 1);
    xir = xi * (xi + 1) / 2.0;
    etal = eta - 1.0 / 2.0;
    etam = - 2 * eta;
    etar = eta + 1.0 / 2.0;
    switch(i)
    {
        case 0:
            ans = xil * etal;
            break;
        case 1:
            ans = xir * etal;
            break;
        case 2:
            ans = xir * etar;
            break;
        case 3:
            ans = xil * etar;
            break;
        case 4:
            ans = xim * etal;
            break;
        case 5:
            ans = xir * etam;
            break;
        case 6:
            ans = xim * etar;
            break;
        case 7:
            ans = xil * etam;
            break;
        case 8:
            ans = xim * etam;
            break;
        default:
            ans = 0;
            break;
    }
    return ans;
}

double bp_xixi(int i, double xi, double eta)
{
    double ans;
    double xil, xim, xir;
    double etal, etam, etar;
    xil = 1;
    xim = - 2;
    xir = 1;
    etal = eta * (eta - 1) / 2.0;
    etam = - (eta + 1) * (eta - 1);
    etar = eta * (eta + 1) / 2.0;
    switch(i)
    {
        case 0:
            ans = xil * etal;
            break;
        case 1:
            ans = xir * etal;
            break;
        case 2:
            ans = xir * etar;
            break;
        case 3:
            ans = xil * etar;
            break;
        case 4:
            ans = xim * etal;
            break;
        case 5:
            ans = xir * etam;
            break;
        case 6:
            ans = xim * etar;
            break;
        case 7:
            ans = xil * etam;
            break;
        case 8:
            ans = xim * etam;
            break;
        default:
            ans = 0;
            break;
    }
    return ans;
}

double bp_etaeta(int i, double xi, double eta)
{
    double ans;
    double xil, xim, xir;
    double etal, etam, etar;
    xil = xi * (xi-1) / 2.0;
    xim = -(xi + 1) * (xi - 1);
    xir = xi * (xi + 1) / 2.0;
    etal = 1;
    etam = - 2;
    etar = 1;
    switch(i)
    {
        case 0:
            ans = xil * etal;
            break;
        case 1:
            ans = xir * etal;
            break;
        case 2:
            ans = xir * etar;
            break;
        case 3:
            ans = xil * etar;
            break;
        case 4:
            ans = xim * etal;
            break;
        case 5:
            ans = xir * etam;
            break;
        case 6:
            ans = xim * etar;
            break;
        case 7:
            ans = xil * etam;
            break;
        case 8:
            ans = xim * etam;
            break;
        default:
            ans = 0;
            break;
    }
    return ans;
}

double bp_xieta(int i, double xi, double eta)
{
    double ans;
    double xil, xim, xir;
    double etal, etam, etar;
    xil = xi - 1.0 / 2.0;
    xim = - 2 * xi;
    xir = xi + 1.0 / 2.0;
    etal = eta - 1.0 / 2.0;
    etam = - 2 * eta;
    etar = eta + 1.0 / 2.0;
    switch(i)
    {
        case 0:
            ans = xil * etal;
            break;
        case 1:
            ans = xir * etal;
            break;
        case 2:
            ans = xir * etar;
            break;
        case 3:
            ans = xil * etar;
            break;
        case 4:
            ans = xim * etal;
            break;
        case 5:
            ans = xir * etam;
            break;
        case 6:
            ans = xim * etar;
            break;
        case 7:
            ans = xil * etam;
            break;
        case 8:
            ans = xim * etam;
            break;
        default:
            ans = 0;
            break;
    }
    return ans;
}

double ** M_inverse(double ** M, int len)
{
    int i, j, k;
    double temp;
    double ** M_inv = new double * [len];
    double ** Mtemp = new double * [len];
    for (i=0; i<len; i++)
    {
        Mtemp[i] = new double [len];
        M_inv[i] = new double [len];
    }

    for (i=0; i<len; i++)
    {
        for (j=0; j<len; j++)
        {
            Mtemp[i][j] = M[i][j];

            if (i == j)
            {
                M_inv[i][j] = 1;
            }
            else{
                M_inv[i][j] = 0;
            }
        }
    }

    for (i=0; i<len; i++)
    {
        temp =  Mtemp[i][i];
        for (j=0; j<len; j++)
        {
            Mtemp[i][j] = Mtemp[i][j] / temp;
            M_inv[i][j] = M_inv[i][j] / temp;
        }
        for (k=0; k<len; k++)
        {
            if ( k == i)
            {
                continue;
            }
            else{
                temp =  Mtemp[k][i];
                for (j=0; j<len; j++)
                {
                    Mtemp[k][j] = Mtemp[k][j] - temp * Mtemp[i][j];
                    M_inv[k][j] = M_inv[k][j] - temp * M_inv[i][j];
                }
            }
        }
    }

    for (i=0; i<len; i++)
    {
        delete[] Mtemp[i];
    }
    delete[] Mtemp;

    return M_inv;
}
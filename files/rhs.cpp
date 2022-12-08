/**
 * @file rhs.cpp
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief 
 * @version 0.1
 * @date 2022-12-03
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#include <iostream>
#include "rhs.h"
#include "config.h"
#include "ICBC.h"
#include <cmath>
using namespace std;

double * get_xieta(int pq, int oq)
{
    double * co = new double [2];

    int ok, ol, r;
    ok = oq / m_element;
    ol = oq % m_element;
    for (r=0; r<9; r++)
    {
        if (o[ok][ol].vertex[r] == pq)
        {
            break;
        }
    }
    if (o[ok][ol].vertex[r] != pq)
    {
        cout<<"point not in element"<<endl;
    }
    co[0] = ref_xi[r];
    co[1] = ref_eta[r];

    return co;
}

int get_subcell_in_vertexneighbor(int sq, int pq)
{
    int pi, pj, r, q;
    pi = pq / (m_point + 1);
    pj = pq % (m_point + 1);
    for (r=0; r<point[pi][pj].neighbor_subcell.size(); r++)
    {
        if (point[pi][pj].neighbor_subcell[r] == sq)
        {
            q = r;
            break;
        }
    }
    if (point[pi][pj].neighbor_subcell[r] != sq)
    {
        cout<<"subcell not arround point"<<endl;
    }
    return q;
}

double * specific_volume_matrix(int i, int j)
{
    double * R = new double [pk];
    int pi, pj, sk, sl, q, r, s, ba;
    double xit, etat;

    for (r=0; r<pk; r++)
    {
        R[r] = 0;
    }

    //下面计算边界通量
    for (r=0; r<4; r++)
    {
        //在每个子网格中计算
        sk = o[i][j].subcell[r] / m_point;
        sl = o[i][j].subcell[r] % m_point;
        
        for (s=0; s<4; s++)
        {
            //在每个子网格的顶点上计算
            pi = subcell[sk][sl].vertex[s] / (m_point + 1);
            pj = subcell[sk][sl].vertex[s] % (m_point + 1);
            //中心点不用计算
            if (point[pi][pj].ifcenter == 1)
            {
                continue;
            }

            //找出所求节点在大网格中对应的局部拉氏坐标
            double * co;
            co = get_xieta(point[pi][pj].q,o[i][j].q);
            xit = co[0];
            etat = co[1];
            delete[] co;
            
            //找出子网格在节点相邻网格的编号
            q = get_subcell_in_vertexneighbor(subcell[sk][sl].q,point[pi][pj].q);

            for (ba=0; ba<pk; ba++)
            {
                double tempx, tempy;
                tempx = point[pi][pj].alast[q] * point[pi][pj].nlast[q][0]
                      + point[pi][pj].anext[q] * point[pi][pj].nnext[q][0];
                
                tempy = point[pi][pj].alast[q] * point[pi][pj].nlast[q][1]
                      + point[pi][pj].anext[q] * point[pi][pj].nnext[q][1];
                
                R[ba] = R[ba] + o[i][j].Psi(ba,xit,etat)
                      * (tempx * point[pi][pj].upstarx
                       + tempy * point[pi][pj].upstary);
            }
        }
    }

    //下面计算体积分
    for (r=0; r<gpn; r++)
    {
        xit = Gausspoint_xi[r];
        etat = Gausspoint_eta[r];

        //下面计算在Gauss求积点的速度
        double ux, uy;
        ux = 0;
        uy = 0;
        for (s=0; s<pk; s++)
        {
            ux = ux + o[i][j].ux[s] * o[i][j].Psi(s,xit,etat);
            uy = uy + o[i][j].uy[s] * o[i][j].Psi(s,xit,etat);
        }

        //下面计算Gauss求积点处的Jacobi矩阵
        double ** J = o[i][j].getJacobiMatrix(xit,etat);

        //下面计算R的第ba个分量
        for (ba=0; ba<pk; ba++)
        {
            double psi_xi, psi_eta;
            psi_xi = o[i][j].Psi_xi(ba,xit,etat);
            psi_eta = o[i][j].Psi_eta(ba,xit,etat);

            R[ba] = R[ba] - Gaussweight[r] *
                    (   ux * (J[1][1] * psi_xi - J[1][0] * psi_eta)
                      + uy * (-J[0][1] * psi_xi + J[0][0] * psi_eta));
        }
//cout<<ux + o[i][j].phi_x0(xit,etat)<<endl;
        delete[] J[0];
        delete[] J[1];
        delete[] J;
    }

    return R;
}

double ** velocity_matrix(int i, int j)
{
    double ** R = new double * [pk];
    int pi, pj, sk, sl, q, r, s, ba;
    double xit, etat;
    for (r=0; r<pk; r++)
    {
        R[r] = new double [2];
    }
    for (r=0; r<pk; r++)
    {
        for (s=0; s<2; s++)
        {
            R[r][s] = 0;
        }
    }

    //下面计算边界通量
    for (r=0; r<4; r++)
    {
        //对每个子网格计算
        sk = o[i][j].subcell[r] / m_point;
        sl = o[i][j].subcell[r] % m_point;

        for (s=0; s<4; s++)
        {
            //在每个子网格的节点计算
            pi = subcell[sk][sl].vertex[s] / (m_point + 1);
            pj = subcell[sk][sl].vertex[s] % (m_point + 1);
            //中心点不用计算
            if (point[pi][pj].ifcenter == 1)
            {
                continue;
            }

            //找出所求节点在大网格中对应的局部坐标
            double * co;
            co = get_xieta(point[pi][pj].q,o[i][j].q);
            xit = co[0];
            etat = co[1];
            delete[] co;

            //找出子网格在节点相邻网格中的编号
            q = get_subcell_in_vertexneighbor(subcell[sk][sl].q,point[pi][pj].q);

            for (ba=0; ba<pk; ba++)
            {
                R[ba][0] = R[ba][0] + o[i][j].Psi(ba,xit,etat)
                                    * (point[pi][pj].flast[q][0]
                                     + point[pi][pj].fnext[q][0]);
                R[ba][1] = R[ba][1] + o[i][j].Psi(ba,xit,etat)
                                    * (point[pi][pj].flast[q][1]
                                     + point[pi][pj].fnext[q][1]);
            }
        }
    }

    //下面计算体积分
    for (r=0; r<gpn; r++)
    {
        xit = Gausspoint_xi[r];
        etat = Gausspoint_eta[r];

        //下面计算Gauss求积点的压力
        double rhot, ux, uy, e, sigma;
        rhot = 0;
        ux = 0;
        uy = 0;
        e = 0;
        for (s=0; s<pk; s++)
        {
            rhot = rhot + o[i][j].nu[s] * o[i][j].Psi(s,xit,etat);
            ux = ux + o[i][j].ux[s] * o[i][j].Psi(s,xit,etat);
            uy = uy + o[i][j].uy[s] * o[i][j].Psi(s,xit,etat);
            e = e + o[i][j].tau[s] * o[i][j].Psi(s,xit,etat);
        }
        rhot = 1.0 / rhot;
        e = e - 0.5 * (ux * ux + uy * uy);
        e = max(1e-9,e);
        sigma = - EOS(o[i][j].gamma,rhot,e);

        //下面计算Gauss求积点处Jacobi矩阵
        double ** J = o[i][j].getJacobiMatrix(xit,etat);

        //下面计算R的第ba个分量
        for (ba=0; ba<pk; ba++)
        {
            double psi_xi, psi_eta;
            psi_xi = o[i][j].Psi_xi(ba,xit,etat);
            psi_eta = o[i][j].Psi_eta(ba,xit,etat);

            R[ba][0] = R[ba][0] - Gaussweight[r] * sigma *
                       (J[1][1] * psi_xi - J[1][0] * psi_eta);
            R[ba][1] = R[ba][1] - Gaussweight[r] * sigma *
                       (-J[0][1] * psi_xi + J[0][0] * psi_eta);
        }

        delete[] J[0];
        delete[] J[1];
        delete[] J;
    }

    double ** Rt = new double * [2];
    Rt[0] = new double [pk];
    Rt[1] = new double [pk];
    for (r=0; r<pk; r++)
    {
        Rt[0][r] = R[r][0];
        Rt[1][r] = R[r][1];
    }
    
    for (r = 0; r<pk; r++)
    {
        delete[] R[r];
    }
    delete[] R;

    return Rt;
}

double * energy_matrix(int i, int j)
{
    double * R = new double [pk];
    int pi, pj, sk, sl, q, r, s, ba;
    double xit, etat;
    for (r=0; r<pk; r++)
    {
        R[r] = 0;
    }

    //下面计算边界通量
    for (r=0; r<4; r++)
    {
        //在每个子网格计算
        sk = o[i][j].subcell[r] / m_point;
        sl = o[i][j].subcell[r] % m_point;

        for (s=0; s<4; s++)
        {
            //在每个子网格的顶点计算
            pi = subcell[sk][sl].vertex[s] / (m_point + 1);
            pj = subcell[sk][sl].vertex[s] % (m_point + 1);
            //中心点不用计算
            if (point[pi][pj].ifcenter == 1)
            {
                continue;
            }

            //找出所求节点在大网格中对应的局部拉氏坐标
            double * co;
            co = get_xieta(point[pi][pj].q,o[i][j].q);
            xit = co[0];
            etat = co[1];
            delete[] co;
            
            //找出子网格在节点相邻网格的编号
            q = get_subcell_in_vertexneighbor(subcell[sk][sl].q,point[pi][pj].q);

            for (ba=0; ba<pk; ba++)
            {
                double tempx, tempy;
                tempx = point[pi][pj].flast[q][0] + point[pi][pj].fnext[q][0];
                tempy = point[pi][pj].flast[q][1] + point[pi][pj].fnext[q][1];
                
                R[ba] = R[ba] + o[i][j].Psi(ba,xit, etat)
                        * (tempx * point[pi][pj].upstarx
                           + tempy * point[pi][pj].upstary);
            }
        }
    }

    //下面计算体积分
    for (r=0; r<gpn; r++)
    {
        xit = Gausspoint_xi[r];
        etat = Gausspoint_eta[r];

        //下面计算Gauss求积点的速度和压力
        double rhot, ux, uy, e, sigma;
        rhot = 0;
        ux = 0;
        uy = 0;
        e = 0;
        for (s=0; s<pk; s++)
        {
            rhot = rhot + o[i][j].nu[s] * o[i][j].Psi(s,xit,etat);
            ux = ux + o[i][j].ux[s] * o[i][j].Psi(s,xit,etat);
            uy = uy + o[i][j].uy[s] * o[i][j].Psi(s,xit,etat);
            e = e + o[i][j].tau[s] * o[i][j].Psi(s,xit,etat);
        }
        rhot = 1.0 / rhot;
        e = e - 0.5 * (ux * ux + uy * uy);
        e = max(1e-9,e);
        sigma = - EOS(o[i][j].gamma,rhot,e);

        //下面计算Gauss求积点处Jacobi矩阵
        double ** J = o[i][j].getJacobiMatrix(xit,etat);

        //下面计算R的第ba个分量
        for (ba=0; ba<pk; ba++)
        {
            double psi_xi, psi_eta;
            psi_xi = o[i][j].Psi_xi(ba,xit,etat);
            psi_eta = o[i][j].Psi_eta(ba,xit,etat);

            R[ba] = R[ba] - Gaussweight[r] *  sigma
                     * (ux * (J[1][1] * psi_xi - J[0][1] * psi_eta)
                      + uy * (-J[0][1] * psi_xi + J[0][0] * psi_eta));
        }

        delete[] J[0];
        delete[] J[1];
        delete[] J;
    }

    return R;
}
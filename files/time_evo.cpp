/**
 * @file time_evo.cpp
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief 
 * @version 0.1
 * @date 2022-12-03
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved. 
 * 
 */

#include <iostream>
#include "time_evo.h"
#include "config.h"
#include "nodal_solver.h"
#include "rhs.h"
#include <cmath>
using namespace std;

double * RK_step1(double * vn, double * Rn,
                  double ** Minv, double delta_t, int len)
{
    double * vs1 = new double [len];
    double * temp = new double [len];
    int i, j;
    for (i=0; i<len; i++)
    {
        vs1[i] = 0;
        temp[i] = 0;
    }

    for (i=0; i<len; i++)
    {
        for (j=0; j<len; j++)
        {
            temp[i] = temp[i] + Minv[i][j] * Rn[j];
        }
    }

    for (i=0; i<len; i++)
    {
        vs1[i] = vn[i] + delta_t * temp[i];
    }

    delete[] temp;
    return vs1;
}

double * RK_step2(double * vn, double * vs1, double * Rs1,
                  double ** Minv, double delta_t, int len)
{
    double * vs2 = new double [len];
    double * temp = new double [len];
    int i, j;
    for (i=0; i<len; i++)
    {
        vs2[i] = 0;
        temp[i] = 0;
    }

    for (i=0; i<len; i++)
    {
        for (j=0; j<len; j++)
        {
            temp[i] = temp[i] + Minv[i][j] * Rs1[j];
        }
    }

    for (i=0; i<len; i++)
    {
        vs2[i] = 3.0 / 4.0 * vn[i] + 1.0 / 4.0 * vs1[i]
               + 1.0 / 4.0 * delta_t * temp[i];
    }
    
    delete[] temp;
    return vs2;
}

double * RK_step3(double * vn, double * vs2, double * Rs2,
                  double ** Minv, double delta_t, int len)
{
    double * vs3 = new double [len];
    double * temp = new double [len];
    int i, j;
    for (i=0; i<len; i++)
    {
        vs3[i] = 0;
        temp[i] = 0;
    }

    for (i=0; i<len; i++)
    {
        for (j=0; j<len; j++)
        {
            temp[i] = temp[i] + Minv[i][j] * Rs2[j];
        }
    }

    for (i=0; i<len; i++)
    {
        vs3[i] = 1.0 / 3.0 * vn[i] + 2.0 / 3.0 * vs2[i]
               + 2.0 / 3.0 * delta_t * temp[i];
    }
    
    delete[] temp;
    return vs3;
}

void one_time_step(double dt)
{
    int i,j, k, l, r, s;
    double ** xn = new double * [n_point+1];
    double ** yn = new double * [n_point+1];
    double ** xs1 = new double * [n_point+1];
    double ** ys1 = new double * [n_point+1];
    double ** xs2 = new double * [n_point+1];
    double ** ys2 = new double * [n_point+1];

    double *** nun = new double ** [n_element];
    double *** nus1 = new double ** [n_element];
    double *** nus2 = new double ** [n_element];
    double *** uxn = new double ** [n_element];
    double *** uxs1 = new double ** [n_element];
    double *** uxs2 = new double ** [n_element];
    double *** uyn = new double ** [n_element];
    double *** uys1 = new double ** [n_element];
    double *** uys2 = new double ** [n_element];
    double *** taun = new double ** [n_element];
    double *** taus1 = new double ** [n_element];
    double *** taus2 = new double ** [n_element];

    //初始化向量
    for (i=0; i<=n_point; i++)
    {
        xn[i] = new double [m_point+1];
        yn[i] = new double [m_point+1];
        xs1[i] = new double [m_point+1];
        ys1[i] = new double [m_point+1];
        xs2[i] = new double [m_point+1];
        ys2[i] = new double [m_point+1];
    }
    for (i=0; i<n_element; i++)
    {
        nun[i] = new double * [m_element];
        uxn[i] = new double * [m_element];
        uyn[i] = new double * [m_element];
        taun[i] = new double * [m_element];
        nus1[i] = new double * [m_element];
        uxs1[i] = new double * [m_element];
        uys1[i] = new double * [m_element];
        taus1[i] = new double * [m_element];
        nus2[i] = new double * [m_element];
        uxs2[i] = new double * [m_element];
        uys2[i] = new double * [m_element];
        taus2[i] = new double * [m_element];
        for (j=0; j<m_element; j++)
        {
            nun[i][j] = new double [pk];
            uxn[i][j] = new double [pk];
            uyn[i][j] = new double [pk];
            taun[i][j] = new double [pk];
            nus1[i][j] = new double [pk];
            uxs1[i][j] = new double [pk];
            uys1[i][j] = new double [pk];
            taus1[i][j] = new double [pk];
            nus2[i][j] = new double [pk];
            uxs2[i][j] = new double [pk];
            uys2[i][j] = new double [pk];
            taus2[i][j] = new double [pk];
        }
    }


    //*************第一次RK*************//
//cout<<"RK1"<<endl;
    //下面计算节点速度
    for (i=0; i<=n_point; i++)
    {
        for (j=0; j<=m_point; j++)
        {
            nodal_velocity_corner_force(i,j);
        }
    }

    //下面更新坐标
    for (i=0; i<=n_point; i++)
    {
        for (j=0; j<=m_point; j++)
        {
            //先记录n时刻x的值
            xn[i][j] = point[i][j].x;
            yn[i][j] = point[i][j].y;
            //对xy做第一步RK
            xs1[i][j] = xn[i][j] + dt * point[i][j].upstarx;
            ys1[i][j] = yn[i][j] + dt * point[i][j].upstary;
            //保存节点坐标
            point[i][j].x = xs1[i][j];
            point[i][j].y = ys1[i][j];
        }
    }

    //下面更新几何量
    update_geo();

    //下面更新物理量

    //下面更新大网格的物理量
    for (k=0; k<n_element; k++)
    {
        for (l=0; l<m_element; l++)
        {
            //先记录n时刻向量
            for (r=0; r<pk; r++)
            {
                nun[k][l][r] = o[k][l].nu[r];
                uxn[k][l][r] = o[k][l].ux[r];
                uyn[k][l][r] = o[k][l].uy[r];
                taun[k][l][r] = o[k][l].tau[r];
            }
            
            //用Minv带入计算防止保存的M_inv被篡改
            double ** Minv = new double * [pk];
            for (r=0; r<pk; r++)
            {
                Minv[r] = new double [pk];
            }

            for (r=0; r<pk; r++)
            {
                for (s=0; s<pk; s++)
                {
                    Minv[r][s] = o[k][l].M_inv[r][s];
                }
            }
            
            //更新几个中心量
            double * nuR, * nutemp, * tauR, * tautemp;
            double * uxtemp, * uytemp, ** uR;
            nuR = specific_volume_matrix(k,l);
            uR = velocity_matrix(k,l);
            tauR = energy_matrix(k,l);

            nutemp = RK_step1(nun[k][l],nuR,Minv,dt,pk);
            uxtemp = RK_step1(uxn[k][l],uR[0],Minv,dt,pk);
            uytemp = RK_step1(uyn[k][l],uR[1],Minv,dt,pk);
            tautemp = RK_step1(taun[k][l],tauR,Minv,dt,pk);

            for (r=0; r<pk; r++)
            {
                nus1[k][l][r] = nutemp[r];
                uxs1[k][l][r] = uxtemp[r];
                uys1[k][l][r] = uytemp[r];
                taus1[k][l][r] = tautemp[r];
                o[k][l].nu[r] = nutemp[r];
                o[k][l].ux[r] = uxtemp[r];
                o[k][l].uy[r] = uytemp[r];
                o[k][l].tau[r] = tautemp[r];
            }

            for (r=0; r<pk; r++)
            {
                delete[] Minv[r];
            }
            delete[] Minv;
            delete[] nuR;
            delete[] nutemp;
            delete[] tauR;
            delete[] tautemp;
            delete[] uxtemp;
            delete[] uytemp;
            delete[] uR[0];
            delete[] uR[1];
            delete[] uR;

            //更新声速
            double ac, pre_c, e, rho_c;
            e = o[k][l].tau[0] - 0.5 * (o[k][l].ux[0] * o[k][l].ux[0]
                                      + o[k][l].uy[0] * o[k][l].uy[0]);
            e = max(1e-9,e);
            rho_c = 1.0 / o[k][l].nu[0];
            pre_c = EOS(o[k][l].gamma,rho_c,e);
            ac = o[k][l].gamma * pre_c / rho_c;
            ac = sqrt(ac);
            o[k][l].c_center = ac;
        }
    }

    //下面更新每个子网格的平均密度
    for (i=0; i<n_point; i++)
    {
        for (j=0; j<m_point; j++)
        {
            subcell[i][j].avg_den = subcell[i][j].m_s / subcell[i][j].area_s;
        }
    }

//*/

    //第二次RK
//cout<<endl<<"RK2"<<endl;
    //下面计算节点速度
    for (i=0; i<=n_point; i++)
    {
        for (j=0; j<=m_point; j++)
        {
            nodal_velocity_corner_force(i,j);
        }
    }

    //下面更新坐标
    for (i=0; i<=n_point; i++)
    {
        for (j=0; j<=m_point; j++)
        {
            //对xy做第二步RK
            xs2[i][j] = 3.0 / 4.0 * xn[i][j] + 1.0 / 4.0 * xs1[i][j]
                      + 1.0 / 4.0 * dt * point[i][j].upstarx;
            ys2[i][j] = 3.0 / 4.0 * yn[i][j] + 1.0 / 4.0 * ys1[i][j]
                      + 1.0 / 4.0 * dt * point[i][j].upstary;
            //保存节点坐标
            point[i][j].x = xs2[i][j];
            point[i][j].y = ys2[i][j];
        }
    }

    //下面更新几何量
    update_geo();

    //下面更新物理量

    //下面更新每个大网格的物理量
    for (k=0; k<n_element; k++)
    {
        for (l=0; l<m_element; l++)
        {
            double ** Minv = new double * [pk];
            for (r=0; r<pk; r++)
            {
                Minv[r] = new double [pk];
            }

            for (r=0; r<pk; r++)
            {
                for (s=0; s<pk; s++)
                {
                    Minv[r][s] = o[k][l].M_inv[r][s];
                }
            }
            //更新几个中心量
            double * nuR, * nutemp, * tauR, * tautemp;
            double * uxtemp, * uytemp, ** uR;
            nuR = specific_volume_matrix(k,l);
            uR = velocity_matrix(k,l);
            tauR = energy_matrix(k,l);

            nutemp = RK_step2(nun[k][l],nus1[k][l],nuR,Minv,dt,pk);
            uxtemp = RK_step2(uxn[k][l],uxs1[k][l],uR[0],Minv,dt,pk);
            uytemp = RK_step2(uyn[k][l],uys1[k][l],uR[1],Minv,dt,pk);
            tautemp = RK_step2(taun[k][l],taus1[k][l],tauR,Minv,dt,pk);

            for (r=0; r<pk; r++)
            {
                nus2[k][l][r] = nutemp[r];
                uxs2[k][l][r] = uxtemp[r];
                uys2[k][l][r] = uytemp[r];
                taus2[k][l][r] = tautemp[r];
                o[k][l].nu[r] = nutemp[r];
                o[k][l].ux[r] = uxtemp[r];
                o[k][l].uy[r] = uytemp[r];
                o[k][l].tau[r] = tautemp[r];
            }

            for (r=0; r<pk; r++)
            {
                delete[] Minv[r];
            }
            delete[] Minv;
            delete[] nuR;
            delete[] nutemp;
            delete[] tauR;
            delete[] tautemp;
            delete[] uxtemp;
            delete[] uytemp;
            delete[] uR[0];
            delete[] uR[1];
            delete[] uR;

            //更新声速
            double ac, pre_c, e, rho_c;
            e = o[k][l].tau[0] - 0.5 * (o[k][l].ux[0] * o[k][l].ux[0]
                                      + o[k][l].uy[0] * o[k][l].uy[0]);
            e = max(1e-9,e);
            rho_c = 1.0 / o[k][l].nu[0];
            pre_c = EOS(o[k][l].gamma,rho_c,e);
            ac = o[k][l].gamma * pre_c / rho_c;
            ac = sqrt(ac);
            o[k][l].c_center = ac;
        }
    }

    //下面更新每个子网格的平均密度
    for (i=0; i<n_point; i++)
    {
        for (j=0; j<m_point; j++)
        {
            subcell[i][j].avg_den = subcell[i][j].m_s / subcell[i][j].area_s;
        }
    }

//*/
    //**************第三次RK********************//

    //下面计算节点速度
    for (i=0; i<=n_point; i++)
    {
        for (j=0; j<=m_point; j++)
        {
            nodal_velocity_corner_force(i,j);
        }
    }

    //下面更新坐标
    for (i=0; i<=n_point; i++)
    {
        for (j=0; j<=m_point; j++)
        {
            //对xy做第三步RK
            point[i][j].x = 1.0 / 3.0 * xn[i][j] + 2.0 / 3.0 * xs2[i][j]
                          + 2.0 / 3.0 * dt * point[i][j].upstarx;
            point[i][j].y = 1.0 / 3.0 * yn[i][j] + 2.0 / 3.0 * ys2[i][j]
                          + 2.0 / 3.0 * dt * point[i][j].upstary;
        }
    }

    //下面更新几何量
    update_geo();

    //下面更新物理量

    //下面更新每个大网格的物理量
    for (k=0; k<n_element; k++)
    {
        for (l=0; l<m_element; l++)
        {
            double ** Minv = new double * [pk];
            for (r=0; r<pk; r++)
            {
                Minv[r] = new double [pk];
            }

            for (r=0; r<pk; r++)
            {
                for (s=0; s<pk; s++)
                {
                    Minv[r][s] = o[k][l].M_inv[r][s];
                }
            }
            //更新几个中心量
            double * nuR, * nutemp, * tauR, * tautemp;
            double * uxtemp, * uytemp, ** uR;
            nuR = specific_volume_matrix(k,l);
            uR = velocity_matrix(k,l);
            tauR = energy_matrix(k,l);

            nutemp = RK_step3(nun[k][l],nus2[k][l],nuR,Minv,dt,pk);
            uxtemp = RK_step3(uxn[k][l],uxs2[k][l],uR[0],Minv,dt,pk);
            uytemp = RK_step3(uyn[k][l],uys2[k][l],uR[1],Minv,dt,pk);
            tautemp = RK_step3(taun[k][l],taus2[k][l],tauR,Minv,dt,pk);

            for (r=0; r<pk; r++)
            {
                o[k][l].nu[r] = nutemp[r];
                o[k][l].ux[r] = uxtemp[r];
                o[k][l].uy[r] = uytemp[r];
                o[k][l].tau[r] = tautemp[r];
            }

            for (r=0; r<pk; r++)
            {
                delete[] Minv[r];
            }
            delete[] Minv;
            delete[] nuR;
            delete[] nutemp;
            delete[] tauR;
            delete[] tautemp;
            delete[] uxtemp;
            delete[] uytemp;
            delete[] uR[0];
            delete[] uR[1];
            delete[] uR;

            //更新声速
            double ac, pre_c, e, rho_c;
            e = o[k][l].tau[0] - 0.5 * (o[k][l].ux[0] * o[k][l].ux[0]
                                      + o[k][l].uy[0] * o[k][l].uy[0]);
            e = max(1e-9,e);
            rho_c = 1.0 / o[k][l].nu[0];
            pre_c = EOS(o[k][l].gamma,rho_c,e);
            ac = o[k][l].gamma * pre_c / rho_c;
            ac = sqrt(ac);
            o[k][l].c_center = ac;
        }
    }

    //下面更新每个子网格的平均密度
    for (i=0; i<n_point; i++)
    {
        for (j=0; j<m_point; j++)
        {
            subcell[i][j].avg_den = subcell[i][j].m_s / subcell[i][j].area_s;
        }
    }

//*/

    for (i=0; i<=n_point; i++)
    {
        delete[] xn[i];
        delete[] yn[i];
        delete[] xs1[i];
        delete[] ys1[i];
        delete[] xs2[i];
        delete[] ys2[i];
    }
    for (i=0; i<n_element; i++)
    {
        for (j=0; j<m_element; j++)
        {
            delete[] nun[i][j];
            delete[] uxn[i][j];
            delete[] uyn[i][j];
            delete[] taun[i][j];
            delete[] nus1[i][j];
            delete[] uxs1[i][j];
            delete[] uys1[i][j];
            delete[] taus1[i][j];
            delete[] nus2[i][j];
            delete[] uxs2[i][j];
            delete[] uys2[i][j];
            delete[] taus2[i][j];
        }
        delete[] nun[i];
        delete[] uxn[i];
        delete[] uyn[i];
        delete[] taun[i];
        delete[] nus1[i];
        delete[] uxs1[i];
        delete[] uys1[i];
        delete[] taus1[i];
        delete[] nus2[i];
        delete[] uxs2[i];
        delete[] uys2[i];
        delete[] taus2[i];
    }
    delete[] xn;
    delete[] yn;
    delete[] xs1;
    delete[] ys1;
    delete[] xs2;
    delete[] ys2;

    return ;
}

void update_geo()
{
    int fk, fl, sk, sl, pi, pj ,r, s;
    
    //下面更新网格顶点位置
    for(fk=0; fk<n_element; fk++)
    {
        for (fl=0; fl<m_element; fl++)
        {
            for (r=0; r<9; r++)
            {
                pi = o[fk][fl].vertex[r] / (m_point + 1);
                pj = o[fk][fl].vertex[r] % (m_point + 1);
                o[fk][fl].vx[r] = point[pi][pj].x;
                o[fk][fl].vy[r] = point[pi][pj].y;
            }
        }
    }

    //下面更新子网格的面积
    for (sk=0; sk<n_point; sk++)
    {
        for (sl=0; sl<m_point; sl++)
        {
            double area, jt, xit, etat;
            double Al, Ar, Bb, Bt;
            //找到子网格在大网格中的位置
            fk = subcell[sk][sl].father_element / m_element;
            fl = subcell[sk][sl].father_element % m_element;
            for (r=0; r<4; r++)
            {
                if (o[fk][fl].subcell[r] == subcell[sk][sl].q)
                {
                    break;
                }
            }
            //找到子网格对应拉氏区域
            switch(r)
            {
                case 0:
                {
                    Al = -1;
                    Ar = 0;
                    Bb = -1;
                    Bt = 0;
                    break;
                }
                case 1:
                {
                    Al = 0;
                    Ar = 1;
                    Bb = -1;
                    Bt = 0;
                    break;
                }
                case 2:
                {
                    Al = 0;
                    Ar = 1;
                    Bb = 0;
                    Bt = 1;
                    break;
                }
                case 3:
                {
                    Al = -1;
                    Ar = 0;
                    Bb = 0;
                    Bt = 1;
                    break;
                }
                default:
                {
                    Al = 0;
                    Ar = 0;
                    Bb = 0;
                    Bt = 0;
                    break;
                }
            }
            //计算面积
            area = 0;
            for (s=0; s<gpn; s++)
            {
                xit = Gausspoint_xi[s];
                etat = Gausspoint_eta[s];
                if (xit >= Al && xit <=Ar)
                {
                    if (etat >= Bb && etat <= Bt)
                    {
                        //在子网格中
                        jt = o[fk][fl].detJacobi(xit,etat);
                    }
                }
                else{
                    jt = 0;
                }
                area = area + jt * Gaussweight[s];
            }
            subcell[sk][sl].area_s = area;
        }
    }

    //下面更新节点单位外法向量和加权边长
    for (pi=0; pi<=n_point; pi++)
    {
        for (pj=0; pj<=m_point; pj++)
        {
            for (r=0; r<point[pi][pj].neighbor_subcell.size(); r++)
            {
                sk = point[pi][pj].neighbor_subcell[r] / m_point;
                sl = point[pi][pj].neighbor_subcell[r] % m_point;

                fk = subcell[sk][sl].father_element / m_element;
                fl = subcell[sk][sl].father_element % m_element;

                //找到节点在大网格中的位置
                int loc;
                for (loc=0; loc<9; loc++)
                {
                    if (o[fk][fl].vertex[loc] == point[pi][pj].q)
                    {
                        break;
                    }
                }
                double xit, etat;
                xit = ref_xi[loc];
                etat = ref_eta[loc];

                //下面计算物理空间单位外法向量
                double * nl = new double [2];
                double * nn = new double [2];
                double ** J;
                double norm;

                J = o[fk][fl].getJacobiMatrix(xit,etat);
                //last
                nl[0] = J[1][1] * point[pi][pj].Nlast[r][0]
                      - J[1][0] * point[pi][pj].Nlast[r][1];
                nl[1] = - J[0][1] * point[pi][pj].Nlast[r][0]
                        + J[0][0] * point[pi][pj].Nlast[r][1];
                norm = nl[0] * nl[0] + nl[1] * nl[1];
                norm = sqrt(norm);
                nl[0] = nl[0] / norm;
                nl[1] = nl[1] / norm;
                point[pi][pj].nlast[r][0] = nl[0];
                point[pi][pj].nlast[r][1] = nl[1];

                //next
                nn[0] = J[1][1] * point[pi][pj].Nnext[r][0]
                      - J[1][0] * point[pi][pj].Nnext[r][1];
                nn[1] = - J[0][1] * point[pi][pj].Nnext[r][0]
                        + J[0][0] * point[pi][pj].Nnext[r][1];
                norm = nn[0] * nn[0] + nn[1] * nn[1];
                norm = sqrt(norm);
                nn[0] = nn[0] / norm;
                nn[1] = nn[1] / norm;
                point[pi][pj].nnext[r][0] = nn[0];
                point[pi][pj].nnext[r][1] = nn[1];

                //下面计算加权边长
                double al, an, detJ;
                //last
                detJ = o[fk][fl].edgedetJacobi(xit,etat,
                                    point[pi][pj].q,point[pi][pj].segmentlast[r]);
                al = point[pi][pj].weightlast[r] * detJ * 2;
                point[pi][pj].alast[r] = al;

                //next
                detJ = o[fk][fl].edgedetJacobi(xit,etat,
                                    point[pi][pj].q,point[pi][pj].segmentnext[r]);
                an = point[pi][pj].weightnext[r] * detJ * 2;
                point[pi][pj].anext[r] = an;

                delete[] J[0];
                delete[] J[1];
                delete[] J;
                delete[] nl;
                delete[] nn;
            }
        }
    }

    return ;
}

double choosedt(double dtn)
{
    double dt, ac, temp;
    int sk, sl, fk, fl, pi, pj, pk, pl, r, s;
    double ax, ay, bx, by, len;

    //下面计算每个子网格中最短边长
    //反正都是近似取值，两点间距就直接用直线计算
    dt = 1000;
    for (sk=0; sk<n_point; sk++)
    {
        for (sl=0; sl<m_point; sl++)
        {
            //先记录声速
            fk = subcell[sk][sl].father_element / m_element;
            fl = subcell[sk][sl].father_element % m_element;
            ac = o[fk][fl].c_center;

            //对每个子网格的顶点循环
            for (r=0; r<4; r++)
            {
                pi = subcell[sk][sl].vertex[r] / (m_point + 1);
                pj = subcell[sk][sl].vertex[r] % (m_point + 1);

                ax = point[pi][pj].x;
                ay = point[pi][pj].y;

                //下面计算到其他节点的时间
                for (s=r+1; s<4; s++)
                {
                    pk = subcell[sk][sl].vertex[s] / (m_point+1);
                    pl = subcell[sk][sl].vertex[s] % (m_point+1);

                    bx = point[pk][pl].x;
                    by = point[pk][pl].y;

                    //下面计算到相邻点的直线距离
                    len = (bx - ax) * (bx - ax) + (by - ay) * (by - ay);
                    len = sqrt(len);

                    //下面判断是否为最短时间
                    temp = len / ac;
                    if (temp < dt)
                    {
                        dt = temp;
                    }
                }
            }
        }
    }
    dt = C_cfl * dt;
    dt = min(dt,C_m * dtn);

    return dt;
}
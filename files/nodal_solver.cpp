/**
 * @file nodal_solver.cpp
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief �ڵ�ⷨ��
 * @version 0.1
 * @date 2022-12-01
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#include <iostream>
#include "nodal_solver.h"
#include "config.h"
#include "ICBC.h"
#include <cmath>
using namespace std;

double * u_average(int i, int j)
{
    double * uavg = new double [2];
    int sk, sl, fk, fl, r, q, loc, temp;
    double xit, etat;

    uavg[0] = 0;
    uavg[1] = 0;
    for (r=0; r<point[i][j].neighbor_subcell.size(); r++)
    {
        //(sk,sl)Ϊ����������ı��
        sk = point[i][j].neighbor_subcell[r] / m_point;
        sl = point[i][j].neighbor_subcell[r] % m_point;

        //(fk,fl)��ʱ�Ǵ�����ı��
        fk = subcell[sk][sl].father_element / m_element;
        fl = subcell[sk][sl].father_element % m_element;
        

        //�����жϴ˽ڵ��ڴ������еľֲ�����
        for (loc=0; loc<9; loc++)
        {
            if (o[fk][fl].vertex[loc] == point[i][j].q)
            {
                break;
            }
        }

        xit = ref_xi[loc];
        etat = ref_eta[loc];

        //�����ʱ�ٶ��ع�ֵ
        double uxt, uyt;
        uxt = 0;
        uyt = 0;
        for (temp=0; temp<pk; temp++)
        {
            uxt = uxt + o[fk][fl].ux[temp] * o[fk][fl].Psi(temp,xit,etat);
            uyt = uyt + o[fk][fl].uy[temp] * o[fk][fl].Psi(temp,xit,etat);
        }
        uavg[0] = uavg[0] + uxt;
        uavg[1] = uavg[1] + uyt;
    }
    uavg[0] = uavg[0] / (double ) point[i][j].neighbor_subcell.size();
    uavg[1] = uavg[1] / (double ) point[i][j].neighbor_subcell.size();
    return uavg;
}

double density_SMS(int i, int j)
{
    double rhoc;
    
    //ע��i,j��������ı��
    
    double rho_s, rho_nus;
    rho_s = subcell[i][j].avg_den;
    
    //k,l�Ǵ�����ı��
    int k, l, r, s;
    k = subcell[i][j].father_element / m_element;
    l = subcell[i][j].father_element % m_element;
    
    //�����ҵ������������
    double Al, Ar, Bb, Bt;
    for (r=0; r<4; r++)
    {
        if (o[k][l].subcell[r] == subcell[i][j].q)
        {
            break;
        }
    }
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

    //�����������������ƽ��
    double mt, xit, etat, rhot, jt;
    mt=0;
    for (s=0; s<gpn; s++)
    {
        xit = Gausspoint_xi[s];
        etat = Gausspoint_eta[s];

        rhot = 0;
        jt =0;
        if (xit >= Al && xit <= Ar)
        {
            if (etat >= Bb && etat <= Bt)
            {
                //����������
                rhot = 0;
                for (int tt=0; tt<pk; tt++)
                {
                    rhot = rhot + o[k][l].nu[tt] * o[k][l].Psi(tt,xit,etat);
                }
                rhot = 1.0 / rhot;
                jt = o[k][l].detJacobi(xit,etat);
            }
        }
        mt = mt + rhot * jt * Gaussweight[s];
    }

    rho_nus = mt / subcell[i][j].area_s;

    rhoc = (rho_s - rho_nus);
    
    return rhoc;
}

void ustar_f_center(int i, int j)
{
    double * ustar = new double [2];

    //���ĵ��ٶȾ��������ڸõ���ٶ�
    //��Ȼ�����ڸõ���������ƽ���ٶ�
    ustar = u_average(i,j);
    point[i][j].upstarx = ustar[0];
    point[i][j].upstary = ustar[1];

    //���ĵ㲻��Ҫ����f
    for (int k=0; k<point[i][j].flast.size(); k++)
    {
        point[i][j].flast[k][0] = 0;
        point[i][j].flast[k][1] = 0;
    }
    for (int k=0; k<point[i][j].fnext.size(); k++)
    {
        point[i][j].fnext[k][0] = 0;
        point[i][j].fnext[k][1] = 0;
    }
    delete[] ustar;
    
    return ;
}


void nodal_velocity_corner_force(int i, int j)
{
    //���ĵ�
    if (point[i][j].ifcenter == 1)
    {
        ustar_f_center(i,j);
        return ;
    }

    //���Ǳ߽��
    int ncell;
    ncell = point[i][j].neighbor_subcell.size();

    int sk, sl, fatherk, fatherl, q, r, s;
    int loc, nodel, noden;
    double xit, etat, mu, rho_sms, p_sms, e, sigma, temp;

    double * amulast = new double [ncell];  //��¼ÿ������������һ�ڵ��a_i mu_i
    double * amunext = new double [ncell];  //��¼ÿ������������һ�ڵ��a_i mu_i
    double * asigmanxlast = new double [ncell]; //��¼ÿ����������һ�ڵ��a_i sigma_c n_i x����
    double * asigmanylast = new double [ncell]; //��¼ÿ����������һ�ڵ��a_i sigma_c n_i y����
    double * asigmanxnext = new double [ncell]; //��¼ÿ����������һ�ڵ��a_i sigma_c n_i x����
    double * asigmanynext = new double [ncell]; //��¼ÿ����������һ�ڵ��a_i sigma_c n_i y����
    
    double * uavg = u_average(i,j); //(i,j)��ƽ���ٶȣ���Ϊ�ڵ��ٶȵĽ���
    double ** subuc = new double * [ncell]; //��¼ÿ���������е��ع��ٶ�uc
    for (r=0; r<ncell; r++)
    {
        subuc[r] = new double [2];
    }
    
    //�Խڵ��������������ѭ������¼ÿ�����������������Ϣ
    for (r=0; r<ncell; r++)
    {
        //�ȼ�¼��������
        sk = point[i][j].neighbor_subcell[r] / m_point;
        sl = point[i][j].neighbor_subcell[r] % m_point;

        //�ټ�¼���ڴ�����ı��
        fatherk = subcell[sk][sl].father_element / m_element;
        fatherl = subcell[sk][sl].father_element % m_element;

        //loc��¼�ڵ��ڴ������еı��
        for(loc=0; loc<9; loc++)
        {
            if (o[fatherk][fatherl].vertex[loc] == point[i][j].q )
            {
                break;
            }
        }
        //xit,etatΪ�ο��ռ�����
        xit = ref_xi[loc];
        etat = ref_eta[loc];

        //�������ڵ㴦�ڴ������е��ٶ��ع�
        double * uc = new double [2];
        uc[0] = 0;
        uc[1] = 0;
        for (s=0; s<pk; s++)
        {
            uc[0] = uc[0] + o[fatherk][fatherl].ux[s]
                          * o[fatherk][fatherl].Psi(s,xit,etat);
            uc[1] = uc[1] + o[fatherk][fatherl].uy[s]
                          * o[fatherk][fatherl].Psi(s,xit,etat);
        }
        subuc[r][0] = uc[0];
        subuc[r][1] = uc[1];
/*if (abs(uc[0])>1)
{
    cout<<"point "<<point[i][j].q<<"\t"<<xit<<"\t"<<etat<<endl;
    for (s=0; s<pk; s++)
    {
        cout<<o[fatherk][fatherl].ux[s]
            <<"\t"<<o[fatherk][fatherl].Psi(s,xit,etat)<<endl;
    }
    cout<<endl<<endl;
}*/
        //�жϼ��������������������������û��
        double ex, ey;
        ex = uavg[0] - uc[0];
        ey = uavg[1] - uc[1];
        double norm;
        norm = ex * ex + ey * ey;
        norm = sqrt(norm);
        if (norm < 1e-6)
        {
            ex = 1e-12;
            ey = 1e-12;
        }
        else{
            ex = ex / norm;
            ey = ey / norm;
        }

        //�������amu�����нڵ��ٶ����ٶ�ƽ������
        //��Ԫ�����ܶ�
        double rho_z;
        rho_z = 1.0 / o[fatherk][fatherl].nu[0];

        //��һ�ڵ�
        temp = (uavg[0] - uc[0]) * point[i][j].nlast[r][0]
             + (uavg[1] - uc[1]) * point[i][j].nlast[r][1];
        
        mu = rho_z * (o[fatherk][fatherl].c_center
                      + (o[fatherk][fatherl].gamma + 1) * temp / 2.0);
        
        //mu = mu * abs(ex * point[i][j].nlast[r][0]
        //              + ey * point[i][j].nlast[r][1]);

        amulast[r] = point[i][j].alast[r] * mu;

        //��һ�ڵ�
        temp = (uavg[0] - uc[0]) * point[i][j].nnext[r][0]
             + (uavg[1] - uc[1]) * point[i][j].nnext[r][1];
        
        mu = rho_z * (o[fatherk][fatherl].c_center
                      + (o[fatherk][fatherl].gamma + 1) * temp / 2.0);

        //mu = mu * abs(ex * point[i][j].nnext[r][0]
        //           + ey * point[i][j].nnext[r][1]);

        amunext[r] = point[i][j].anext[r] * mu;
        //�������asigma
        //�������ÿ���������ڽڵ��ع������ܺ��ܶ�
        e = 0;
        rho_sms = 0;
        for (s=0; s<pk; s++)
        {
            rho_sms = rho_sms + o[fatherk][fatherl].nu[s]
                              * o[fatherk][fatherl].Psi(s,xit,etat);
            e = e + o[fatherk][fatherl].tau[s]
                  * o[fatherk][fatherl].Psi(s,xit,etat);
        }
        //��ʱe������
        e = e - 0.5 * (uc[0] * uc[0] + uc[1] * uc[1]);
        e = max(1e-12,e);
        rho_sms = 1.0 / rho_sms;
        rho_sms = rho_sms + density_SMS(sk,sl);

        p_sms = EOS(o[fatherk][fatherl].gamma,rho_sms,e);
        sigma = - p_sms;
if (abs(e - 0.998333)<1e-5)
{
    cout<<i<<"\t"<<j<<"\t"<<point[i][j].ifedge<<"\t"<<point[i][j].boundary<<endl;
    cout<<e<<"\t"<<p_sms<<endl;
}

        asigmanxlast[r] = point[i][j].alast[r] * sigma * point[i][j].nlast[r][0];
        asigmanylast[r] = point[i][j].alast[r] * sigma * point[i][j].nlast[r][1];

        asigmanxnext[r] = point[i][j].anext[r] * sigma * point[i][j].nnext[r][0];
        asigmanynext[r] = point[i][j].anext[r] * sigma * point[i][j].nnext[r][1];

        delete[] uc;
    }

    //������½ڵ��ٶ�
    double topx, topy, bottom;
    topx = 0;
    topy = 0;
    bottom = 0;
    for (r=0; r<ncell; r++)
    {
        topx = topx + amulast[r] * subuc[r][0]
                    + amunext[r] * subuc[r][0]
                    - asigmanxlast[r]
                    - asigmanxnext[r];
        topy = topy + amulast[r] * subuc[r][1]
                    + amunext[r] * subuc[r][1]
                    - asigmanylast[r]
                    - asigmanynext[r];
        bottom = bottom + amulast[r] + amunext[r];
    }

    if (abs(bottom) < 1e-9)
    {
        point[i][j].upstarx = uavg[0];
        point[i][j].upstary = uavg[1];
    }
    else{
        point[i][j].upstarx = topx / bottom;
        point[i][j].upstary = topy / bottom;
    }

    if(point[i][j].boundary == 1)
    {
        //�߽��
        BCvelocity(i,j,bottom);
    }

//point[i][j].upstarx = ini_ux(point[i][j].x0,point[i][j].y0);
//point[i][j].upstary = ini_uy(point[i][j].x0,point[i][j].y0);

    //�������ÿ�������������ϵ�corner force
    double flx, fly, fnx, fny;
    double test1, test2;
    test1 = 0;
    test2 = 0;
    for (r=0; r<ncell; r++)
    {
        temp = point[i][j].upstarx - subuc[r][0];

        flx = asigmanxlast[r] + amulast[r] * temp;
        fnx = asigmanxnext[r] + amunext[r] * temp;
        
        temp = point[i][j].upstary - subuc[r][1];

        fly = asigmanylast[r] + amulast[r] * temp;
        fny = asigmanynext[r] + amunext[r] * temp;

        point[i][j].flast[r][0] = flx;
        point[i][j].flast[r][1] = fly;

        point[i][j].fnext[r][0] = fnx;
        point[i][j].fnext[r][1] = fny;

        test1 = test1 + flx + fnx;
        test2 = test2 + fly + fny;
    }
    if (point[i][j].boundary != 1)
    {
        if (abs(test1)>1e-5 || abs(test2)>1e-5)
        {
            cout<<"point "<<point[i][j].q<<"\t"<<i<<"\t"<<j<<endl;
            cout<<"if edge "<<point[i][j].ifedge<<endl;
            cout<<"sum f "<<test1<<"\t"<<test2<<endl;
            for (r=0; r<ncell; r++)
            {
                cout<<"amu r "<<amulast[r]<<"\t"<<amunext[r]<<endl;
            }
            cout<<endl;
        }
        
    }

    delete[] amulast;
    delete[] amunext;
    delete[] asigmanxlast;
    delete[] asigmanylast;
    delete[] asigmanxnext;
    delete[] asigmanynext;
    delete[] uavg;
    for (r=0; r<ncell; r++)
    {
        delete[] subuc[r];
    }
    delete[] subuc;

    return ;
}

void BCvelocity(int i, int j, double amu)
{
    double xt, yt, uxt, uyt;
    double Pl, Pn;
    double nstarx, nstary;
    switch(testcase)
    {
        case shockless_Noh:
            xt = point[i][j].x0;
            yt = point[i][j].y0;
            point[i][j].upstarx = ini_ux(xt,yt);
            point[i][j].upstary = ini_uy(xt,yt);
            break;
        case Taylor_Green_vortex:
            xt = point[i][j].x0;
            yt = point[i][j].y0;
            point[i][j].upstarx = ini_ux(xt,yt);
            point[i][j].upstary = ini_uy(xt,yt);
            break;
        
        default:
            BCnormal_velocity(i,j,vn);
            
            Pl = Pre / amu;
            Pn = Pre / amu;
            BCpressure(i,j,Pl,Pn);
            break;
    }
    return ;
}

void BCnormal_velocity(int i, int j, double vnt)
{
    double ux, uy, Pre;
    double nstarx1, nstary1, nstarx2, nstary2;
    int r, node1, node2, k1, l1, k2, l2;
    
    //���ҵ���һ�����ڱ߽����ⷨ��
    for (r=0; r<point[i][j].neighbor_subcell.size(); r++)
    {
        node1 = point[i][j].segmentlast[r];
        k1 = node1 / (m_point + 1);
        l1 = node1 % (m_point + 1);
        if (point[k1][l1].boundary == 1)
        {
            nstarx1 = point[i][j].nlast[r][0];
            nstary1 = point[i][j].nlast[r][1];
            break;
        }
        node1 = point[i][j].segmentnext[r];
        k1 = node1 / (m_point + 1);
        l1 = node1 % (m_point + 1);
        if (point[k1][l1].boundary == 1)
        {
            nstarx1 = point[i][j].nnext[r][0];
            nstary1 = point[i][j].nnext[r][1];
            break;
        }
    }
    //��ʱcode1�ǵ�һ�����ڱ߽�㣬nstarx1, nstary1�ǵ�һ���߽��ⷨ����

    //���ҵڶ����߽��
    for (r=0; r<point[i][j].neighbor_subcell.size(); r++)
    {
        node2 = point[i][j].segmentlast[r];
        k2 = node2 / (m_point + 1);
        l2 = node2 % (m_point + 1);
        if (point[k2][l2].boundary == 1)
        {
            if (node2 != node1)
            {
                nstarx2 = point[i][j].nlast[r][0];
                nstary2 = point[i][j].nlast[r][1];
                break;
            }
        }
        node2 = point[i][j].segmentnext[r];
        k2 = node2 / (m_point + 1);
        l2 = node2 % (m_point + 1);
        if (point[k2][l2].boundary == 1)
        {
            if (node2 != node1)
            {
                nstarx2 = point[i][j].nnext[r][0];
                nstary2 = point[i][j].nnext[r][1];
                break;
            }
        }
    }
    
    //�ж������ⷨ���Ƿ�ƽ��
    double para;
    para = nstarx1 * nstary2 - nstary1 * nstarx2;
    para = abs(para);

    if (para >= 1e-6)
    {
        //��ƽ��
        point[i][j].upstarx = vnt * (nstary2 - nstary1) / para;
        point[i][j].upstary = vnt * (- nstarx2 + nstarx1) / para;
        return ;
    }
    
    //ƽ��
    double nx, ny;
    nx = 0.5 * (nstarx1 + nstarx2);
    ny = 0.5 * (nstary1 + nstary2);

    Pre = ux * nx + uy * ny;

    point[i][j].upstarx = ux + Pre * nx + vnt;
    point[i][j].upstary = uy + Pre * ny + vnt;

    return ;
}

void BCpressure(int i, int j, double pcorlast, double pcornext)
{
    double ux, uy;
    double nstarx1, nstary1, nstarx2, nstary2;
    int r, node1, node2, k1, l1, k2, l2;
    
    //���ҵ���һ�����ڱ߽����ⷨ��
    for (r=0; r<point[i][j].neighbor_subcell.size(); r++)
    {
        node1 = point[i][j].segmentlast[r];
        k1 = node1 / (m_point + 1);
        l1 = node1 % (m_point + 1);
        if (point[k1][l1].boundary == 1)
        {
            nstarx1 = point[i][j].nlast[r][0];
            nstary1 = point[i][j].nlast[r][1];
            break;
        }
        node1 = point[i][j].segmentnext[r];
        k1 = node1 / (m_point + 1);
        l1 = node1 % (m_point + 1);
        if (point[k1][l1].boundary == 1)
        {
            nstarx1 = point[i][j].nnext[r][0];
            nstary1 = point[i][j].nnext[r][1];
            break;
        }
    }
    //��ʱcode1�ǵ�һ�����ڱ߽�㣬nstarx1, nstary1�ǵ�һ���߽��ⷨ����

    //���ҵڶ����߽��
    for (r=0; r<point[i][j].neighbor_subcell.size(); r++)
    {
        node2 = point[i][j].segmentlast[r];
        k2 = node2 / (m_point + 1);
        l2 = node2 % (m_point + 1);
        if (point[k2][l2].boundary == 1)
        {
            if (node2 != node1)
            {
                nstarx2 = point[i][j].nlast[r][0];
                nstary2 = point[i][j].nlast[r][1];
                break;
            }
        }
        node2 = point[i][j].segmentnext[r];
        k2 = node2 / (m_point + 1);
        l2 = node2 % (m_point + 1);
        if (point[k2][l2].boundary == 1)
        {
            if (node2 != node1)
            {
                nstarx2 = point[i][j].nnext[r][0];
                nstary2 = point[i][j].nnext[r][1];
                break;
            }
        }
    }
    //��node1, node2��˳ʱ������
    int temp, plast, pnext;
    double nstarxl, nstarxn, nstaryl, nstaryn;
    for (r=0; r<point[i][j].neighbor_node.size(); r++)
    {
        if (point[i][j].neighbor_node[r] == node1)
        {
            //node1�ȳ���
            plast = node1;
            pnext = node2;
            nstarxl = nstarx1;
            nstaryl = nstary1;
            nstarxn = nstarx2;
            nstaryn = nstary2;
            break;
        }
        if (point[i][j].neighbor_node[r] == node2)
        {
            //node2�ȳ���
            plast = node2;
            pnext = node1;
            nstarxl = nstarx2;
            nstaryl = nstary2;
            nstarxn = nstarx1;
            nstaryn = nstary1;
            break;
        }
    }

    point[i][j].upstarx = ux + pcorlast * nstarxl
                             + pcornext * nstarxn;
    point[i][j].upstary = uy + pcorlast * nstaryl
                             + pcornext * nstaryn;

    return ;
}
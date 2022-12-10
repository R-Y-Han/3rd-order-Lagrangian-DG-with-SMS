/**
 * @file initial.cpp
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief ��������ͼ�¼��Ӧ��Ϣ
 * @version 0.1
 * @date 2022-11-29
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#include <iostream>
#include "initial.h"
#include "config.h"
#include "ICBC.h"
#include <cmath>
using namespace std;

void generatemesh()
{
    int i, j;

    //�����¼�ڵ��ź����꣬�ж��Ƿ�Ϊ�߽�ڵ�
    for (i=0; i<=n_point; i++)
    {
        for (j=0; j<=m_point; j++)
        {
            //��¼���
            point[i][j].q = (m_point + 1) * i + j;

            //�ж��Ƿ�Ϊ�߽�ڵ�
            point[i][j].boundary = 0;
            if (i==0 || i == n_point || j == 0 || j == m_point)
            {
                point[i][j].boundary = 1;
            }
            //�ж��Ƿ�Ϊ�����������
            point[i][j].ifcenter = 0;
            if (i%2 == 1 && j%2 == 1)
            {
                point[i][j].ifcenter = 1;
            }
            //�ж��Ƿ�Ϊ���ϵĵ�
            point[i][j].ifedge = 0;
            if (i%2 == 0 && j%2 == 1)
            {
                point[i][j].ifedge = 1;
            }
            if (i%2 == 1 && j%2 == 0)
            {
                point[i][j].ifedge = 1;
            }

            //��¼�ڵ���������
            point[i][j].x = hx * i;
            point[i][j].y = hy * j;
            point[i][j].x0 = hx * i;
            point[i][j].y0 = hy * j;
            point[i][j].upstarx = 0;
            point[i][j].upstary = 0;
        }
    }

    //�����¼������
    for (i=0; i<n_element; i++)
    {
        for (j=0; j<m_element; j++)
        {
            //��¼������
            o[i][j].q = m_element * i + j;

            //��¼���񶥵�
            int tempi[9] = {2*i, 2*i+2, 2*i+2, 2*i,
                            2*i+1, 2*i+2, 2*i+1, 2*i,
                            2*i+1};
            int tempj[9] = {2*j, 2*j, 2*j+2, 2*j+2,
                            2*j, 2*j+1, 2*j+2, 2*j+1,
                            2*j+1};
            int it, jt;
            for (int k=0; k<9; k++)
            {
                it = tempi[k];
                jt = tempj[k];
                o[i][j].vertex[k] = point[it][jt].q;
                o[i][j].vx[k] = point[it][jt].x;
                o[i][j].vy[k] = point[it][jt].y;
                o[i][j].vx0[k] = point[it][jt].x0;
                o[i][j].vy0[k] = point[it][jt].y0;
            }
        }
    }

    //�����¼������ı�źͶ��㣬���ڴ�������
    for (i=0; i<n_point; i++)
    {
        for (j=0; j<m_point; j++)
        {
            subcell[i][j].q = m_point * i + j;
            subcell[i][j].father_element = m_element * (i/2) + j/2;
            int tempi[4] = {i, i+1, i+1, i};
            int tempj[4] = {j, j, j+1, j+1};
            int it, jt;
            for (int k=0; k<4; k++)
            {
                it = tempi[k];
                jt = tempj[k];
                subcell[i][j].vertex[k] = point[it][jt].q;
            }
        }
    }

    //�����¼����������������
    for (i=0; i<n_element; i++)
    {
        for (j=0; j<m_element; j++)
        {
            int tempi[4] = {2*i, 2*i+1, 2*i+1, 2*i};
            int tempj[4] = {2*j, 2*j, 2*j+1, 2*j+1};
            int it, jt;
            for (int k=0; k<4; k++)
            {
                it = tempi[k];
                jt = tempj[k];
                o[i][j].subcell[k] = subcell[it][jt].q;
            }
        }
    }
}

void element_findneighbor(int i, int j)
{
    o[i][j].neighbor_element.clear();

    int it, jt;
    int num;

    //�ж�(i,j-1)��Ԫ�Ƿ����
    if ( j == 0) 
    {
        //������
        ;
    }
    else{
        //����
        it = i;
        jt = j-1;
        num = o[it][jt].q;
        o[i][j].neighbor_element.push_back(num);
    }

    //�ж�(i+1,j)��Ԫ�Ƿ����
    if ( i == n_element - 1 ) 
    {
        //������
        ;
    }
    else{
        //����
        it = i+1;
        jt = j;
        num = o[it][jt].q;
        o[i][j].neighbor_element.push_back(num);
    }

    //�ж�(i,j+1)��Ԫ�Ƿ����
    if ( j == m_element-1) 
    {
        //������
        ;
    }
    else{
        //����
        it = i;
        jt = j+1;
        num = o[it][jt].q;
        o[i][j].neighbor_element.push_back(num);
    }

    //�ж�(i-1,j)��Ԫ�Ƿ����
    if ( i == 0) 
    {
        //������
        ;
    }
    else{
        //����
        it = i-1;
        jt = j;
        num = o[it][jt].q;
        o[i][j].neighbor_element.push_back(num);
    }

    return ;
}

void node_findneighbor(int i, int j)
{
    int ti, tj; //�������������ڵ���б�ź��б��
    int num;    //�������������ڵ��һά���
    
    //**************�ȼ�¼��������****************//
    point[i][j].neighbor_subcell.clear();
    
    //�ж�����(i,j)�Ƿ����
    if ( i == n_point)
    {
        //������
        ;
    }
    else{
        ti = i;
        if (j == m_point)
        {
            //�ڵ�Ϊ���Ͻ�
            ;
        }
        else{
            tj = j;
            num = subcell[ti][tj].q;
            point[i][j].neighbor_subcell.push_back(num);
        }
    }

    //�ж�����(i-1,j)�Ƿ����
    if ( i == 0)
    {
        //������
        ;
    }
    else{
        ti = i-1;
        if (j == m_point)
        {
            //������
            ;
        }
        else{
            tj = j;
            num = subcell[ti][tj].q;
            point[i][j].neighbor_subcell.push_back(num);
        }
    }

    //�ж�����(i-1,j-1)�Ƿ����
    if ( i == 0)
    {
        //������
        ;
    }
    else{
        ti = i-1;
        if (j == 0)
        {
            //������
            ;
        }
        else{
            tj = j-1;
            num = subcell[ti][tj].q;
            point[i][j].neighbor_subcell.push_back(num);
        }
    }

    //�ж�����(i,j-1)�Ƿ����
    if ( i == n_point)
    {
        //������
        ;
    }
    else{
        ti = i;
        if (j == 0)
        {
            //������
            ;
        }
        else{
            tj = j-1;
            num = subcell[ti][tj].q;
            point[i][j].neighbor_subcell.push_back(num);
        }
    }
    //*****************���������¼���******************//

    //*****************�����¼���ڽڵ�******************//
    point[i][j].neighbor_node.clear();

    //�ж�(i,j-1)�ڵ��Ƿ����
    if ( j == 0) 
    {
        //������
        ;
    }
    else{
        //����
        ti = i;
        tj = j-1;
        num = point[ti][tj].q;
        point[i][j].neighbor_node.push_back(num);
    }

    //�ж�(i+1,j)�ڵ��Ƿ����
    if ( i == n_point) 
    {
        //������
        ;
    }
    else{
        //����
        ti = i+1;
        tj = j;
        num = point[ti][tj].q;
        point[i][j].neighbor_node.push_back(num);
    }

    //�ж�(i,j+1)�ڵ��Ƿ����
    if ( j == m_point) 
    {
        //������
        ;
    }
    else{
        //����
        ti = i;
        tj = j+1;
        num = point[ti][tj].q;
        point[i][j].neighbor_node.push_back(num);
    }

    //�ж�(i-1,j)�ڵ��Ƿ����
    if ( i == 0) 
    {
        //������
        ;
    }
    else{
        //����
        ti = i-1;
        tj = j;
        num = point[ti][tj].q;
        point[i][j].neighbor_node.push_back(num);
    }
    //**********���ڽڵ��¼���******************//

    return ;
}

void node_initialsegment(int i, int j)
{
    //**********������㵽���ڽڵ�ĳ��Ⱥ��ⷨ��**********//
    point[i][j].weightnext.clear();
    point[i][j].weightlast.clear();
    point[i][j].segmentnext.clear();
    point[i][j].segmentlast.clear();
    point[i][j].anext.clear();
    point[i][j].alast.clear();
    point[i][j].Nlast.clear();
    point[i][j].Nnext.clear();
    point[i][j].nnext.clear();
    point[i][j].nlast.clear();
    point[i][j].flast.clear();
    point[i][j].fnext.clear();
    for (int r=0; r < point[i][j].neighbor_subcell.size(); r++)
    {
        int k, l;
        k = point[i][j].neighbor_subcell[r] / m_point;
        l = point[i][j].neighbor_subcell[r] % m_point;

        int loc;
        for (loc=0; loc<4; loc++)
        {
            if (subcell[k][l].vertex[loc] == point[i][j].q)
            {
                break;
            }
        }
        //��ʱloc��ʾpoint(i,j)��*������*(k,l)�еľֲ����

        int noden, nodel;
        if (loc == 3)
        {
            noden = 0;
        }
        else{
            noden = loc + 1;
        }

        if (loc == 0)
        {
            nodel = 3;
        }
        else{
            nodel = loc - 1;
        }
        //noder, nodelΪ(i,j)��һ������һ���ֲ�����
        point[i][j].segmentnext.push_back(subcell[k][l].vertex[noden]);
        point[i][j].segmentlast.push_back(subcell[k][l].vertex[nodel]);

        double * ftempnext = new double [2];
        double * ftemplast = new double [2];
        ftempnext[0] = 0;
        ftempnext[1] = 0;
        ftemplast[0] = r;
        ftemplast[1] = 0;
        point[i][j].fnext.push_back(ftempnext);
        point[i][j].flast.push_back(ftemplast);

        //��������Ȩ�ⷨ����
        int inext, ilast, jnext, jlast;
        inext = subcell[k][l].vertex[noden] / (m_point + 1);
        jnext = subcell[k][l].vertex[noden] % (m_point + 1);
        ilast = subcell[k][l].vertex[nodel] / (m_point + 1);
        jlast = subcell[k][l].vertex[nodel] % (m_point + 1);
        
        //����洢ÿ���ߵ�Ȩ��
        if (point[i][j].ifcenter == 1)
        {
            //�Ǵ������е�
            point[i][j].weightnext.push_back(1);
            point[i][j].weightlast.push_back(1);
        }
        else if (point[i][j].ifedge == 0)
        {
            //�Ƕ���
            point[i][j].weightnext.push_back(1.0/6.0);
            point[i][j].weightlast.push_back(1.0/6.0);
        }
        else{
            //�Ǳ��ϵĵ�
            if (point[inext][jnext].ifcenter == 1)
            {
                point[i][j].weightnext.push_back(4.0/6.0);
            }
            else{
                point[i][j].weightnext.push_back(2.0/6.0);
            }

            if (point[ilast][jlast].ifcenter == 1)
            {
                point[i][j].weightlast.push_back(4.0/6.0);
            }
            else{
                point[i][j].weightlast.push_back(2.0/6.0);
            }
        }

        //�������ÿ���ߵ�λ�ⷨ����
        //�ȼ���ο���Ԫ��λ������
        double * ntempnext = new double [2];
        double * ntemplast = new double [2];
        
        if (loc == 0)
        {
            //next
            ntempnext[0] = 0;
            ntempnext[1] = -1;
            point[i][j].Nnext.push_back(ntempnext);
            //last
            ntemplast[0] = -1;
            ntemplast[1] = 0;
            point[i][j].Nlast.push_back(ntemplast);
        }
        else if(loc == 1)
        {
            //next
            ntempnext[0] = 1;
            ntempnext[1] = 0;
            point[i][j].Nnext.push_back(ntempnext);
            //last
            ntemplast[0] = 0;
            ntemplast[1] = -1;
            point[i][j].Nlast.push_back(ntemplast);
        }
        else if(loc == 2)
        {
            //next
            ntempnext[0] = 0;
            ntempnext[1] = 1;
            point[i][j].Nnext.push_back(ntempnext);
            //last
            ntemplast[0] = 1;
            ntemplast[1] = 0;
            point[i][j].Nlast.push_back(ntemplast);
        }
        else{
            //next
            ntempnext[0] = -1;
            ntempnext[1] = 0;
            point[i][j].Nnext.push_back(ntempnext);
            //last
            ntemplast[0] = 0;
            ntemplast[1] = 1;
            point[i][j].Nlast.push_back(ntemplast);
        }

        //�����������ռ䵥λ�ⷨ����(n = |J|J^{-T} N)
        //���Ҽ����Ȩ�ⷨ����(weight * s * n)
        double detJ, xit, etat;
        double ** J;
        int Fcell, ifather, jfather, fatherloc;
        Fcell = subcell[k][l].father_element;
        ifather = Fcell / m_element;
        jfather = Fcell % m_element;

        //fatherloc��¼�ڴ������еľֲ����
        for (fatherloc=0; fatherloc<9; fatherloc++)
        {
            if(o[ifather][jfather].vertex[fatherloc] == point[i][j].q)
            {
                break;
            }
        }
        xit = ref_xi[fatherloc];
        etat = ref_eta[fatherloc];

        J = o[ifather][jfather].getJacobiMatrix(xit,etat);

        double norm;

        double * antempnext = new double [2];
        //next
        antempnext[0] = J[1][1] * ntempnext[0] - J[1][0] * ntempnext[1];
        antempnext[1] = - J[0][1] * ntempnext[0] + J[0][0] * ntempnext[1];
        norm = antempnext[0] * antempnext[0] + antempnext[1] * antempnext[1];
        norm = sqrt(norm);
        antempnext[0] = antempnext[0] / norm;
        antempnext[1] = antempnext[1] / norm;
        //��ʱantempnextΪ����ռ��ⷨ����
        point[i][j].nnext.push_back(antempnext);
        //�������ڱ߲�ֵ��|J|
        double atnext;
        detJ = o[ifather][jfather].edgedetJacobi(xit,etat,point[i][j].q,point[i][j].segmentnext[r]);
        atnext = point[i][j].weightnext[r] * detJ * 2;
        //��ʱatnextΪ����ռ��ϵļ�Ȩ�߳������г�2����Ϊ�ο��ռ�߳�Ϊ2
        point[i][j].anext.push_back(atnext);
        
        double * antemplast = new double [2];
        //last
        antemplast[0] = J[1][1] * ntemplast[0] - J[1][0] * ntemplast[1];
        antemplast[1] = - J[0][1] * ntemplast[0] + J[0][0] * ntemplast[1];
        norm = antemplast[0] * antemplast[0] + antemplast[1] * antemplast[1];
        norm = sqrt(norm);
        antemplast[0] = antemplast[0] / norm;
        antemplast[1] = antemplast[1] / norm;
        //��ʱantemplastΪ����ռ��ϵ�λ�ⷨ����
        point[i][j].nlast.push_back(antemplast);
        //�������ڱ߲�ֵ��|J|
        double atlast;
        detJ = o[ifather][jfather].edgedetJacobi(xit,etat,point[i][j].q,point[i][j].segmentlast[r]);
        atlast = point[i][j].weightlast[r] * detJ * 2;
        //��ʱatlastΪ����ռ��ϼ�Ȩ�߳�
        point[i][j].alast.push_back(atlast);

        delete[] J[0];
        delete[] J[1];
        delete[] J;
    }
    
    return ;
}

double * mass_center(int i, int j)
{
    double * ce = new double [2];

    int k;

    double masst;    //��Ԫ����
    double mf_xi, mf_eta;   //������

    //********��������Gauss���ּ���������������************//
    double rhot, jt;    //��¼�ܶȺ�Jacobi����ʽ�ĵ�ֵ
    double xt, yt;  //Gauss�����������
    masst = 0;
    mf_xi = 0;
    mf_eta = 0;
    for (k = 0; k < gpn; k++)
    {
        xt = o[i][j].phi_x0(Gausspoint_xi[k],Gausspoint_eta[k]);
        yt = o[i][j].phi_y0(Gausspoint_xi[k],Gausspoint_eta[k]);

        rhot = ini_rho(xt,yt);
        jt = o[i][j].detJacobi_0(Gausspoint_xi[k],Gausspoint_eta[k]);

        masst = masst + rhot * jt * Gaussweight[k];
        mf_xi = mf_xi + rhot * Gausspoint_xi[k] * jt * Gaussweight[k];
        mf_eta = mf_eta + rhot * Gausspoint_eta[k] * jt * Gaussweight[k];
    }
    o[i][j].mass = masst;

    //**********���������������***********//

    ce[0] = mf_xi / masst;  //���Ĳο��ռ������
    ce[1] = mf_eta / masst; //���Ĳο��ռ�������

    //������������345������ƽ��
    double a1, a2, a3;
    a1 = 0;
    a2 = 0;
    a3 = 0;
    for (k=0; k<gpn; k++)
    {
        xt = o[i][j].phi_x0(Gausspoint_xi[k],Gausspoint_eta[k]);
        yt = o[i][j].phi_y0(Gausspoint_xi[k],Gausspoint_eta[k]);

        rhot = ini_rho(xt,yt);
        jt = o[i][j].detJacobi_0(Gausspoint_xi[k],Gausspoint_eta[k]);

        a1 = a1 + rhot * (Gausspoint_xi[k] - ce[0])
                       * (Gausspoint_xi[k] - ce[0]) / 2.0 * jt * Gaussweight[k];
        a2 = a2 + rhot * (Gausspoint_eta[k] - ce[0])
                       * (Gausspoint_eta[k] - ce[0]) / 2.0 * jt * Gaussweight[k];
        a3 = a3 + rhot * (Gausspoint_xi[k] - ce[0])
                       * (Gausspoint_eta[k] - ce[0]) * jt * Gaussweight[k];
    }
    a1 = a1 / masst;
    a2 = a2 / masst;
    a3 = a3 / masst;

    o[i][j].psi3avg = a1;
    o[i][j].psi4avg = a2;
    o[i][j].psi5avg = a3;

    return ce;
}

double ** mass_matrix(int i, int j)
{
    int k, l;
    double rhot, jt, xt, yt, temp;
    double ** mma = new double * [pk];
    for ( k = 0; k < pk; k++)
    {
        mma[k] = new double [pk];
    }

    for (k = 0; k<pk; k++)
    {
        for (l=0; l<pk; l++)
        {
            mma[k][l] = 0;
            for (int r = 0; r <gpn; r++)
            {
                xt = o[i][j].phi_x0(Gausspoint_xi[r], Gausspoint_eta[r]);
                yt = o[i][j].phi_y0(Gausspoint_xi[r], Gausspoint_eta[r]);

                rhot = ini_rho(xt,yt);

                jt = o[i][j].detJacobi_0(Gausspoint_xi[r], Gausspoint_eta[r]);
                
                temp = o[i][j].Psi(k,Gausspoint_xi[r],Gausspoint_eta[r]);
                temp = temp * o[i][j].Psi(l,Gausspoint_xi[r],Gausspoint_eta[r]);
                
                mma[k][l] = mma[k][l] + rhot * temp * jt * Gaussweight[r];
            }
        }
    }

    return mma;
}

void subcell_mass_in_ij(int i, int j)
{
    double mass, area;
    int k, l, r, s;
    double xit, etat, xt, yt, rhot, jt, Al, Ar, Bb, Bt;
    for (r=0; r<4; r++)
    {
        k = o[i][j].subcell[r] / m_point;
        l = o[i][j].subcell[r] % m_point;
        //��������������
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

        mass = 0;
        area = 0;
        for (s=0; s<gpn; s++)
        {
            xit = Gausspoint_xi[s];
            etat = Gausspoint_eta[s];
            xt = o[i][j].phi_x0(xit,etat);
            yt = o[i][j].phi_y0(xit,etat);

            rhot = 0;
            jt =0;
            if (xit >= Al && xit <= Ar)
            {
                if (etat >= Bb && etat <= Bt)
                {
                    //����������
                    rhot = ini_rho(xt,yt);
                    jt = o[i][j].detJacobi_0(xit,etat);

                }
            }
            else{
                rhot = 0;
                jt = 0;
            }

            mass = mass + rhot * jt * Gaussweight[s];
            area = area + jt * Gaussweight[s];
        }

        subcell[k][l].m_s = mass;
        subcell[k][l].area_s = area;
        subcell[k][l].avg_den = mass / area;
    }
    return;
}

void initial()
{
    int i, j, k;
    //***********�����������񻮷�**********//
    generatemesh();
    
        //***********�����ÿ������Ѱ����������**********//
    for (i=0; i<n_element; i++)
    {
        for (j=0; j<m_element; j++)
        {
            element_findneighbor(i,j);
        }
    }

        //***********�����ÿ���ڵ�Ѱ�����ڽڵ�***********//
    for (i=0; i<=n_point; i++)
    {
        for (j=0; j<=m_point; j++)
        {
            node_findneighbor(i,j);
            node_initialsegment(i,j);
        }
    }


    //************�������ÿ����������������ģ���������������*********//
    for (i=0; i<n_element; i++)
    {
        for (j=0; j<m_element; j++)
        {
            //�����������
            double * ce;
            ce = mass_center(i,j);
            o[i][j].xi_c = ce[0];
            o[i][j].eta_c = ce[1];

            delete[] ce;

            //�������������������
            double ** Mat;
            Mat = mass_matrix(i,j);
            o[i][j].M = Mat;
            o[i][j].M_inv = M_inverse(Mat,pk);
        }
    }

    //�ڵ�(i,j)���������м���ÿ�����������������ʼ�����ƽ���ܶ�
    for (i=0; i<n_element; i++)
    {
        for (j=0; j<m_element; j++)
        {
            subcell_mass_in_ij(i,j);
        }
    }

//***********�����ÿ����Ԫ����ֵ***************//
    double xt, yt; //���ĵ���������
    double x_xi, x_eta, y_xi, y_eta, x_xixi, x_etaeta, x_xieta, y_xixi, y_etaeta, y_xieta;    //Jacobi������ĸ�����
    double taux, tauy, tauxx, tauyy, tauxy;
    for (i=0; i<n_element; i++)
    {
        for (j=0; j<m_element; j++)
        {
            //�ȸ����������ö��ݱ���
            switch(testcase)
            {
                case shockless_Noh:
                    o[i][j].gamma = 5.0 / 3.0;
                    break;
                case Taylor_Green_vortex:
                    o[i][j].gamma = 7.0 / 5.0;
                    break;
                default:
                    o[i][j].gamma = 5.0 / 3.0;
                    break;
            }

            o[i][j].nu = new double [pk];
            o[i][j].ux = new double [pk];
            o[i][j].uy = new double [pk];
            o[i][j].tau = new double [pk];

            xt = o[i][j].phi_x(o[i][j].xi_c,o[i][j].eta_c);
            yt = o[i][j].phi_y(o[i][j].xi_c,o[i][j].eta_c);
           
            x_xi = 0;
            x_eta = 0;
            y_xi = 0;
            y_eta = 0;
            x_xixi = 0;
            x_etaeta = 0;
            x_xieta = 0;
            y_xixi = 0;
            y_etaeta = 0;
            y_xieta = 0;
            for (k = 0; k< 9; k++)
            {
                x_xi = x_xi + o[i][j].vx[k] * bp_xi(k,o[i][j].xi_c,o[i][j].eta_c);
                x_eta = x_eta + o[i][j].vx[k] * bp_eta(k,o[i][j].xi_c,o[i][j].eta_c);
                y_xi = y_xi + o[i][j].vy[k] * bp_xi(k,o[i][j].xi_c,o[i][j].eta_c);
                y_eta = y_eta + o[i][j].vy[k] * bp_eta(k,o[i][j].xi_c,o[i][j].eta_c);
                x_xixi = x_xixi + o[i][j].vx[k] * bp_xixi(k,o[i][j].xi_c,o[i][j].eta_c);
                x_etaeta = x_etaeta + o[i][j].vx[k] * bp_etaeta(k,o[i][j].xi_c,o[i][j].eta_c);
                x_xieta = x_xieta + o[i][j].vx[k] * bp_xieta(k,o[i][j].xi_c,o[i][j].eta_c);
                y_xixi = y_xixi + o[i][j].vy[k] * bp_xixi(k,o[i][j].xi_c,o[i][j].eta_c);
                y_etaeta = y_etaeta + o[i][j].vy[k] * bp_etaeta(k,o[i][j].xi_c,o[i][j].eta_c);
                y_xieta = y_xieta + o[i][j].vy[k] * bp_xieta(k,o[i][j].xi_c,o[i][j].eta_c);
            }

            //********�����ʼ���ٶȳ�************//
            double rhot = ini_rho(xt,yt);
            o[i][j].nu[0] = 1.0 / rhot;
            o[i][j].nu[1] = ini_rho_x(xt,yt) * x_xi / (- rhot * rhot)
                          + ini_rho_y(xt,yt) * y_xi / (- rhot * rhot);
            o[i][j].nu[2] = ini_rho_x(xt,yt) * x_eta / (-rhot * rhot)
                          + ini_rho_y(xt,yt) * y_eta / (-rhot * rhot);
            o[i][j].nu[3] = (2 * ini_rho_x(xt,yt) * ini_rho_x(xt,yt) / (rhot*rhot*rhot) - ini_rho_xx(xt,yt) / (rhot*rhot)) * x_xi * x_xi
                          + (2 * ini_rho_x(xt,yt) * ini_rho_y(xt,yt) / (rhot*rhot*rhot) - ini_rho_xy(xt,yt) / (rhot*rhot)) * x_xi * y_xi
                          + ini_rho_x(xt,yt) * x_xixi / (-rhot * rhot)
                          + (2 * ini_rho_x(xt,yt) * ini_rho_y(xt,yt) / (rhot*rhot*rhot) - ini_rho_xy(xt,yt) / (rhot*rhot)) * x_xi * y_xi
                          + (2 * ini_rho_y(xt,yt) * ini_rho_y(xt,yt) / (rhot*rhot*rhot) - ini_rho_yy(xt,yt) / (rhot*rhot)) * y_xi * y_xi
                          + ini_rho_y(xt,yt) * y_xixi / (-rhot * rhot);
            o[i][j].nu[4] = (2 * ini_rho_x(xt,yt) * ini_rho_x(xt,yt) / (rhot*rhot*rhot) - ini_rho_xx(xt,yt) / (rhot*rhot)) * x_eta * x_eta
                          + (2 * ini_rho_x(xt,yt) * ini_rho_y(xt,yt) / (rhot*rhot*rhot) - ini_rho_xy(xt,yt) / (rhot*rhot)) * x_eta * y_eta
                          + ini_rho_x(xt,yt) * x_etaeta / (-rhot * rhot)
                          + (2 * ini_rho_x(xt,yt) * ini_rho_y(xt,yt) / (rhot*rhot*rhot) - ini_rho_xy(xt,yt) / (rhot*rhot)) * x_eta * y_eta
                          + (2 * ini_rho_y(xt,yt) * ini_rho_y(xt,yt) / (rhot*rhot*rhot) - ini_rho_yy(xt,yt) / (rhot*rhot)) * y_eta * y_eta
                          + ini_rho_y(xt,yt) * y_etaeta / (-rhot * rhot);
            o[i][j].nu[5] = (2 * ini_rho_x(xt,yt) * ini_rho_x(xt,yt) / (rhot*rhot*rhot) - ini_rho_xx(xt,yt) / (rhot*rhot)) * x_xi * x_eta
                          + (2 * ini_rho_x(xt,yt) * ini_rho_y(xt,yt) / (rhot*rhot*rhot) - ini_rho_xy(xt,yt) / (rhot*rhot)) * x_eta * y_xi
                          + ini_rho_x(xt,yt) * x_xieta / (-rhot * rhot)
                          + (2 * ini_rho_x(xt,yt) * ini_rho_y(xt,yt) / (rhot*rhot*rhot) - ini_rho_xy(xt,yt) / (rhot*rhot)) * x_xi * y_eta
                          + (2 * ini_rho_y(xt,yt) * ini_rho_y(xt,yt) / (rhot*rhot*rhot) - ini_rho_yy(xt,yt) / (rhot*rhot)) * y_xi * y_eta
                          + ini_rho_y(xt,yt) * y_xieta / (-rhot * rhot);

            o[i][j].ux[0] = ini_ux(xt,yt);
            o[i][j].ux[1] = ini_ux_x(xt,yt) * x_xi + ini_ux_y(xt,yt) * y_xi;
            o[i][j].ux[2] = ini_ux_x(xt,yt) * x_eta + ini_ux_y(xt,yt) * y_eta;
            o[i][j].ux[3] = ini_ux_xx(xt,yt) * x_xi * x_xi + ini_ux_xy(xt,yt) * x_xi * y_xi
                          + ini_ux_x(xt,yt) * x_xixi + ini_ux_xy(xt,yt) * x_xi * y_xi
                          + ini_ux_yy(xt,yt) * y_xi * y_xi + ini_ux_y(xt,yt) * y_xixi;
            o[i][j].ux[4] = ini_ux_xx(xt,yt) * x_eta * x_eta + ini_ux_xy(xt,yt) * x_eta * y_eta
                          + ini_ux_x(xt,yt) * x_etaeta + ini_ux_xy(xt,yt) * x_eta * y_eta
                          + ini_ux_yy(xt,yt) * y_eta * y_eta + ini_ux_y(xt,yt) * y_etaeta;
            o[i][j].ux[5] = ini_ux_xx(xt,yt) * x_xi * x_eta + ini_ux_xy(xt,yt) * x_eta * y_xi
                          + ini_ux_x(xt,yt) * x_xieta + ini_ux_xy(xt,yt) * x_xi * y_eta
                          + ini_ux_yy(xt,yt) * y_xi * y_eta + ini_ux_y(xt,yt) * y_xieta;
            
            o[i][j].uy[0] = ini_uy(xt,yt);
            o[i][j].uy[1] = ini_uy_x(xt,yt) * x_xi + ini_uy_y(xt,yt) * y_xi;
            o[i][j].uy[2] = ini_uy_x(xt,yt) * x_eta + ini_uy_y(xt,yt) * y_eta;
            o[i][j].uy[3] = ini_uy_xx(xt,yt) * x_xi * x_xi + ini_uy_xy(xt,yt) * x_xi * y_xi
                          + ini_uy_x(xt,yt) * x_xixi + ini_uy_xy(xt,yt) * x_xi * y_xi
                          + ini_uy_yy(xt,yt) * y_xi * y_xi + ini_uy_y(xt,yt) * y_xixi;
            o[i][j].uy[4] = ini_uy_xx(xt,yt) * x_eta * x_eta + ini_uy_xy(xt,yt) * x_eta * y_eta
                          + ini_uy_x(xt,yt) * x_etaeta + ini_uy_xy(xt,yt) * x_eta * y_eta
                          + ini_uy_yy(xt,yt) * y_eta * y_eta + ini_uy_y(xt,yt) * y_etaeta;
            o[i][j].uy[5] = ini_uy_xx(xt,yt) * x_xi * x_eta + ini_uy_xy(xt,yt) * x_eta * y_xi
                          + ini_uy_x(xt,yt) * x_xieta + ini_uy_xy(xt,yt) * x_xi * y_eta
                          + ini_uy_yy(xt,yt) * y_xi * y_eta + ini_uy_y(xt,yt) * y_xieta;
            
            double uxt, uyt, gammat;
            uxt = ini_ux(xt,yt);
            uyt = ini_uy(xt,yt);
            gammat = o[i][j].gamma;
            taux = ini_p_x(xt,yt) / ((gammat - 1) * rhot)
                 - ini_p(xt,yt) * ini_rho_x(xt,yt) / ((gammat - 1) * rhot * rhot)
                 + uxt * ini_ux_x(xt,yt) + uyt * ini_uy_x(xt,yt);
            tauy = ini_p_y(xt,yt) / ((gammat - 1) * rhot)
                 - ini_p(xt,yt) * ini_rho_y(xt,yt) / ((gammat - 1) * rhot * rhot)
                 + uxt * ini_ux_y(xt,yt) + uyt * ini_uy_y(xt,yt);
            
            tauxx = ini_ux_x(xt,yt) * ini_ux_x(xt,yt) + uxt * ini_ux_xx(xt,yt)
                  + ini_uy_x(xt,yt) * ini_uy_x(xt,yt) + uyt * ini_uy_xx(xt,yt)
                  + ini_p_xx(xt,yt) / ((gammat - 1) * rhot) - ini_p_x(xt,yt) * ini_rho_x(xt,yt) / ((gammat - 1) * rhot * rhot)
                  - (ini_p_x(xt,yt) * ini_rho_x(xt,yt) + ini_p(xt,yt) * ini_rho_xx(xt,yt)) / ((gammat - 1) * rhot * rhot)
                  + ini_p(xt,yt) * ini_rho_x(xt,yt) * 2 * rhot * ini_rho_x(xt,yt) / ((gammat - 1) * pow(rhot,4));
            
            tauyy = ini_ux_y(xt,yt) * ini_ux_y(xt,yt) + uxt * ini_ux_yy(xt,yt)
                  + ini_uy_y(xt,yt) * ini_uy_y(xt,yt) + uyt * ini_uy_yy(xt,yt)
                  + ini_p_yy(xt,yt) / ((gammat - 1) * rhot) - ini_p_y(xt,yt) * ini_rho_y(xt,yt) / ((gammat - 1) * rhot * rhot)
                  - (ini_p_y(xt,yt) * ini_rho_y(xt,yt) + ini_p(xt,yt) * ini_rho_yy(xt,yt)) / ((gammat - 1) * rhot * rhot)
                  + ini_p(xt,yt) * ini_rho_y(xt,yt) * 2 * rhot * ini_rho_y(xt,yt) / ((gammat - 1) * pow(rhot,4));

            tauxy = ini_ux_x(xt,yt) * ini_ux_y(xt,yt) + uxt * ini_ux_xy(xt,yt)
                  + ini_uy_x(xt,yt) * ini_uy_y(xt,yt) + uyt * ini_uy_xy(xt,yt)
                  + ini_p_xy(xt,yt) / ((gammat - 1) * rhot) - ini_p_x(xt,yt) * ini_rho_y(xt,yt) / ((gammat - 1) * rhot * rhot)
                  - (ini_p_y(xt,yt) * ini_rho_x(xt,yt) + ini_p(xt,yt) * ini_rho_xy(xt,yt)) / ((gammat - 1) * rhot * rhot)
                  + ini_p(xt,yt) * ini_rho_x(xt,yt) * 2 * rhot * ini_rho_y(xt,yt) / ((gammat - 1) * pow(rhot,4));

            o[i][j].tau[0] = ini_p(xt,yt) / ((gammat - 1) * rhot)
                           + 0.5 * (uxt * uxt + uyt * uyt);
            o[i][j].tau[1] = taux * x_xi + tauy * y_xi;
            o[i][j].tau[2] = taux * x_eta + tauy * y_eta;
            o[i][j].tau[3] = tauxx * x_xi * x_xi + tauxy * x_xi * y_xi
                           + taux * x_xixi + tauxy * x_xi * y_xi
                           + tauyy * y_xi * y_xi + tauy * y_xixi;
            
            o[i][j].tau[4] = tauxx * x_eta * x_eta + tauxy * x_eta * y_eta
                           + taux * x_etaeta + tauxy * x_eta * y_eta
                           + tauyy * y_eta * y_eta + tauy * y_etaeta;
            
            o[i][j].tau[5] = tauxx * x_xi * x_eta + tauxy * x_eta * y_xi
                           + taux * x_xieta + tauxy * x_xi * y_eta
                           + tauyy * y_xi * y_eta + tauy * y_xieta;

            o[i][j].rho_0 = rhot;

            o[i][j].c_center = gammat * ini_p(xt,yt) / rhot;
            o[i][j].c_center = sqrt(o[i][j].c_center);
        }
    }

    return ;
}


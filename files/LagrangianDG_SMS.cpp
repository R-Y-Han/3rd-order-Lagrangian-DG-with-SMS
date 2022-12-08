/**
* @mainpage  A high-order Lagrangian discontinuous Galerkin hydrodynamic method for quadratic cells using a subcell mesh stabilization
* <table>
* <tr><th>Project  <td> Code implementation for A high-order Lagrangian discontinuous Galerkin hydrodynamic method for quadratic cells using a subcell mesh stabilization
* <tr><th>Author   <td> R.Y. Han
* <tr><th>Institution   <td> IAPCM
* </table>
* @section   Project details
* Code the method from paper A high-order Lagrangian discontinuous Galerkin hydrodynamic method for quadratic cells using a subcell mesh stabilization on several testproblems.
* 
**********************************************************************************
*/

/**
 * @file LagrangianDG_SMS.cpp
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief 二维三阶Lagrangian DG格式
 * @version 1.01
 * @date 2022-11-28
 * 
 * @copyright Copyright (c) 2022 R.Y. Han All Rights Reserved.
 * 
 */

#include <iostream>
#include <cmath>
#include "config.h"
#include "ICBC.h"
#include "initial.h"
#include "nodal_solver.h"
#include "rhs.h"
#include "time_evo.h"
#include "plot.h"
#include <fstream>
using namespace std;

double norme();
void normux();

int main()
{
    int i, j, k, l, sk, sl;
    initial();

    double dt = 1000;
    double t = 0;
    while(t<T)
    {
        dt = choosedt(dt);
        if (t+dt >= T)
        {
            dt = T - t;
        }
        one_time_step(dt);
        t = t + dt;
        cout<<t<<endl;
        //break;
    }
    double norm;
    norm = norme();
    normux();
    plotmesh();

    ofstream f;
    const char* fn = "F:\\C++Code\\LagrangianDG_with_SMS\\output\\shockless_Noh\\norm.txt";
    f.open(fn, ios::app);
    f<<"n="<<n_element<<"\t"<<"m="<<m_element<<"\t"<<"T="<<T<<"\t"<<"Pk="<<pk<<endl;
    f<<"norm="<<norm<<endl<<endl;
    f.close();

    system("pause");
}

double norme()
{
    double temp;
    int i,j,k;
    temp = 0;
    for (i=0; i<n_element; i++)
    {
        for (j=0; j<m_element; j++)
        {
            for (k=0; k<gpn; k++)
            {
                double jt, xit, etat;
                double xt, yt;
                xit = Gausspoint_xi[k];
                etat = Gausspoint_eta[k];
                xt = o[i][j].phi_x(xit,etat);
                yt = o[i][j].phi_y(xit,etat);
                jt = o[i][j].detJacobi(xit, etat);
                double  ux, uy, e;
                ux = 0;
                uy = 0;
                e = 0;
                for (int l=0; l<pk; l++)
                {
                    ux = ux + o[i][j].ux[l] * o[i][j].Psi(l,xit,etat);
                    uy = uy + o[i][j].uy[l] * o[i][j].Psi(l,xit,etat);
                    e = e + o[i][j].tau[l] * o[i][j].Psi(l,xit,etat);
                }
                e = e - 0.5 * (ux * ux + uy * uy);
                e = max(1e-9,e);
                double e_ana;
                e_ana = ana_e(xt,yt,T);
                temp = temp + (e_ana - e) * (e_ana - e) * jt * Gaussweight[k];
            }
        }
    }
    cout<<temp<<endl;
    return temp;
}

void normux()
{
    double temp;
    int i,j,k;
    temp = 0;
    for (i=0; i<n_element; i++)
    {
        for (j=0; j<m_element; j++)
        {
            for (k=0; k<gpn; k++)
            {
                double jt, xit, etat;
                double xt, yt;
                xit = Gausspoint_xi[k];
                etat = Gausspoint_eta[k];
                xt = o[i][j].phi_x(xit,etat);
                yt = o[i][j].phi_y(xit,etat);
                jt = o[i][j].detJacobi(xit, etat);
                double  ux, uy, e;
                ux = 0;
                uy = 0;
                e = 0;
                for (int l=0; l<pk; l++)
                {
                    ux = ux + o[i][j].ux[l] * o[i][j].Psi(l,xit,etat);
                    uy = uy + o[i][j].uy[l] * o[i][j].Psi(l,xit,etat);
                }
                double ux_ana;
                double xt0, yt0;
                xt0 = o[i][j].phi_x0(xit,etat);
                yt0 = o[i][j].phi_y0(xit,etat);
                ux_ana = ini_ux(xt0,yt0);
                temp = temp + (ux_ana - ux) * (ux_ana - ux) * jt * Gaussweight[k];
            }
        }
    }
    cout<<temp<<endl;
    return ;
}
/**
 * @file plot.cpp
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief 
 * @version 0.1
 * @date 2022-12-05
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#include "plot.h"
#include "config.h"
#include <fstream>
#include <iostream>
using namespace std;

//shockless Noh
const char* fn_mesh = "F:\\C++Code\\LagrangianDG_with_SMS\\output\\shockless_Noh\\mesh.plt";
//Taylor-Green vortex
//const char* fn_mesh = "F:\\C++Code\\LagrangianDG_with_SMS\\output\\Taylor-Green_vortex\\mesh.plt";

void plotmesh()
{
    //????tecplot??????
    int i, j, r;
    remove(fn_mesh);
    ofstream f(fn_mesh);
    f<<"VARIABLES = X, Y"<<endl;
    f<<"ZONE N = "<<(n_point + 1) * (m_point + 1)<<", E = "<<n_point*m_point<<", F = FEPOINT, ET = QUADRILATERAL"<<endl;
    for (i=0; i<=n_point; i++)
    {
        for (j=0; j<=m_point; j++)
        {
            f<<"\t"<<point[i][j].x<<"\t"<<point[i][j].y<<endl;
        }
    }

    for (i=0; i<n_point; i++)
    {
        for (j=0; j<m_point; j++)
        {
            for (r=0; r<4; r++)
            {
                f<<"\t"<<subcell[i][j].vertex[r]+1;
            }
            f<<endl;
        }
    }

    f.close();
    return ;
}
#ifndef __POTENTIALGRADIENTS_TPP_
#define __POTENTIALGRADIENTS_TPP_


namespace Icarus
{
    template<typename Scalar>
    Icarus::getInnerPotentialGradients (const FullVector<Scalar>& potential, const size_t nx, const size_t ny, const size_t nz, double h, std::vector<char> discpoints,
                                     FullVector<Scalar>& xComp, FullVector<Scalar>& yComp, FullVector<Scalar>& zComp)
    {
        size_t globDim(nx*ny*nz);
        if (potential.get_dim() != globDim) LOG_ERROR("Dimension of potential vector does not equal problem dimension");
        if (xComp.get_dim() != globDim) LOG_ERROR("Dimension of vector xComp does not equal problem dimension");
        if (yComp.get_dim() != globDim) LOG_ERROR("Dimension of vector yComp does not equal problem dimension");
        if (zComp.get_dim() != globDim) LOG_ERROR("Dimension of vector zComp does not equal problem dimension");
        size_t globIdx;
        for (size_t z(1); z < nz-1; z++)
        {
            for (size_t y(1); y < ny-1 ; y++)
            {
                globIdx=z*nx*ny+y*nx+1;
                for (size_t x(1); x < nx-1; x++, globIdx++)
                {
                    if (discpoints[gloIdx]=="i")
                    {
                        zComp[globIdx]=(1.0/(2.0*h)) * (potential[globIdx + nx*ny] - potential[globIdx - nx*ny]);
                        yComp[globIdx]=(1.0/(2.0*h)) * (potential[globIdx + nx]    - potential[globIdx - nx]);
                        xComp[globIdx]=(1.0/(2.0*h)) * (potential[globIdx + 1]     - potential[globIdx - 1]);
                    }else{
                        zComp[globIdx]=0.0;
                        yComp[globIdx]=0.0;
                        xComp[globIdx]=0.0;
                    }
                }
            }
        }
    }

    template<typename Scalar>
    Icarus::getCellMidpointGradientsFEM (const FullVector<Scalar>& potential, const size_t nx, const size_t ny, const size_t nz, std::vector<char> discpoints,
                                     FullVector<Scalar>& xComp, FullVector<Scalar>& yComp, FullVector<Scalar>& zComp)
    {
        size_t globDim((nx-1)*(ny-1)*(nz-1);
        if (potential.get_dim() != globDim) LOG_ERROR("Dimension of potential vector does not equal problem dimension");
        if (xComp.get_dim() != globDim) LOG_ERROR("Dimension of vector xComp does not equal problem dimension");
        if (yComp.get_dim() != globDim) LOG_ERROR("Dimension of vector yComp does not equal problem dimension");
        if (zComp.get_dim() != globDim) LOG_ERROR("Dimension of vector zComp does not equal problem dimension");
        size_t globIdx;
        size_t vecIdx;
        std::vector<std::vector<double>> gradCoeff(8);
        std::vector<std::vector<double>> zwsp;
        for (size_t j(0); j<8, j++)
        {
            gradCoeff[j]=assembleFem::evaluate_gradient_Basis3d(0,j,0.5*h,0.5*h,0.5*h);
            zwsp=assembleFem::evaluated_gradient_Basis3d(j);
            gradCoeff[j]={zwsp[0][13],zwsp[1][13],zwsp[2][13]}
        }
        for (size_t z(0); z < nz-1; z++)
        {
            for (size_t y(0); y < ny-1 ; y++)
            {
                globIdx=z*nx*ny+y*nx;
                for (size_t x(0); x < nx-1 ; x++ , globIdx++ , vecIdx++)
                {
                    if (discpoints[gloIdx]=="i")
                    {
                        zComp[vecIdx] = gradCoeff[0][0]*potential[globIdx]          + gradCoeff[1][0]*potential[globIdx+1]
                                  + gradCoeff[2][0]*potential[globIdx+nx]       + gradCoeff[3][0]*potential[globIdx+nx+1]
                                  + gradCoeff[4][0]*potential[globIdx+nx*ny]    + gradCoeff[5][0]*potential[globIdx+nx*ny+1]
                                  + gradCoeff[6][0]*potential[globIdx+nx*ny+nx] + gradCoeff[7][0]*potential[globIdx+nx*ny+nx+1];

                        yComp[vecIdx] = gradCoeff[0][1]*potential[globIdx]          + gradCoeff[1][1]*potential[globIdx+1]
                                  + gradCoeff[2][1]*potential[globIdx+nx]       + gradCoeff[3][1]*potential[globIdx+nx+1]
                                  + gradCoeff[4][1]*potential[globIdx+nx*ny]    + gradCoeff[5][1]*potential[globIdx+nx*ny+1]
                                  + gradCoeff[6][1]*potential[globIdx+nx*ny+nx] + gradCoeff[7][1]*potential[globIdx+nx*ny+nx+1];

                        xComp[vecIdx] = gradCoeff[0][2]*potential[globIdx]          + gradCoeff[1][2]*potential[globIdx+1]
                                  + gradCoeff[2][2]*potential[globIdx+nx]       + gradCoeff[3][2]*potential[globIdx+nx+1]
                                  + gradCoeff[4][2]*potential[globIdx+nx*ny]    + gradCoeff[5][2]*potential[globIdx+nx*ny+1]
                                  + gradCoeff[6][2]*potential[globIdx+nx*ny+nx] + gradCoeff[7][2]*potential[globIdx+nx*ny+nx+1];
                    }
                }
            }
        }
    }
}


#endif // __POTENTIALGRADIENTS_TPP_

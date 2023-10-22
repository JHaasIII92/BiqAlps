#include "headers/BiqUtil.h"
#include "headers/BiqModel.h"
#include <iostream>
#include <utility>
#include <vector>
#include <array>
#include <math.h>

void FillSparseMatrix(Sparse& sMat, const double *pdData, int N)
{
    int i = 0;
    int j = 0;
    double dTmp;

    for (j = 0; j < N; ++j) 
    {
        for (i = j; i < N; ++i)
        {
            dTmp = pdData[i + j * N];
            if(fabs(dTmp) > 0.0)
            {
                sMat.push_back(BiqSparseTriple(i,j,dTmp));
            }
        }
    }
}

/// @brief 
/// @param sMat 
void PrintSparseMatrix(Sparse sMat) 
{
    std::printf("BiqUtil PrintSparseMatrix\n");
    std::printf("==========================\n");
    for(auto &it : sMat)
    {
        if(it.dVal_ != 0.0)
        {
            std::printf("M[%d, %d] = %f\n", 
                        it.i_, it.j_, it.dVal_);
        }
        else
        {
            break;
        }
    }
    std::printf("==========================\n");
}

void PrintMatrix(double * pdMat, size_t nRow, size_t nCol)
{
    size_t i = 0;
    size_t j = 0;
    std::printf("BiqUtil PrintMatrix\n");
    std::printf("===============================================================\n");
    for(size_t pos = 0; pos < nRow*nCol; ++pos)
    {
        std::printf("%8.2f ",pdMat[pos]);
        ++i;
        if(i == nCol)
        {
            std::printf("\n");
            i = 0;
        }
    }
    std::printf("===============================================================\n");
}

void print_vector(double *vec, int N) 
{
    int i;
    for (i = 0; i < N; i++) {
        printf("%24.16e\n", vec[i]);
    }
}

void print_symmetric_matrix(double *Mat, int N) 
{
    int i, j;
    double val;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            val = (i >= j) ? Mat[i + j*N] : Mat[j + i*N];
            printf("%24.16e", val);
        }
        printf("\n");
    }
}
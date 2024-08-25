#include "headers/BiqUtil.h"
#include "headers/BiqModel.h"
#include <iostream>
#include <utility>
#include <vector>
#include <array>
#include <math.h>

// use for sorting the heap on the triangle cuts
bool operator<(BiqTriInequality const &l, BiqTriInequality const &r)
{
    return l.value_ < r.value_;
}




int GreedyGains::minBestFix()
{
    // assume we are fixing to zero
    int iRet = 0;
    // if gain one is "better" then zero we fix to one
    if(gainOne_ < gainZero_) iRet = 1;

    return iRet;
}

int GreedyGains::maxBestFix()
{
    // assume we are fixing to zero
    int iRet = 0;
    // if gain one is better then zero we fix to one
    if(gainOne_ > gainZero_) iRet = 1;

    return iRet; 
}

int GreedyGains::BestFix(bool bIsMax)
{
    // wrap the other min max best fix methods
    int iRet;
    if(bIsMax)
    {
        iRet = maxBestFix();
    }
    else
    {
        iRet = minBestFix();
    }
    return iRet; 
}


int I(BiqTriTuple btt) 
{
    return std::get<0>(btt);
}
int J(BiqTriTuple btt) 
{
    return std::get<1>(btt);
}
int K(BiqTriTuple btt) 
{
    return std::get<2>(btt);
}
int TRI_TYPE(BiqTriTuple btt) 
{
    return std::get<3>(btt);
}

void FillSparseMatrix(Sparse& sMat, std::vector<double> data)
{
    int i = 0;
    int j = 0;
    double dTmp;
    int N = sqrt(data.size());

    for (j = 0; j < N; ++j) 
    {
        for (i = j; i < N; ++i)
        {
            dTmp = data.at(i + j * N);
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

void PrintMatrix(double * pdMat, int nRow, int nCol)
{
    int i = 0;
    std::printf("BiqUtil PrintMatrix\n");
    std::printf("===============================================================\n");
    for(int pos = 0; pos < nRow*nCol; ++pos)
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

void zeroOutMatrix(std::vector<double> &mat)
{
    for(auto &it: mat)
    {
        it = 0;
    }
}

void PrintMatrix(std::vector<double> mat)
{
    int N = sqrt(mat.size());

    for(int i = 0; i < N; ++ i)
    {
        for(int j = 0; j < N; ++j)
        {
            std::printf("%f ", mat.at(i + j*N));
        }
        std::printf("\n");
    }
    std::printf("\n");
}
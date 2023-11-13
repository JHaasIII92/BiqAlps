// *****************************************************************************
// *                                                                           *
// *  MLT-BiqCrunch is a shared-memory multithreaded, semidefinite-based       * 
// *  solver for binary quadratic problems. It uses a branch-and-bound         *
// *  method featuring an improved semidefinite bounding procedure, mixed      *
// *  with a polyhedral approach. MLT-BiqCrunch uses particular input files    *
// *  format to describe the combinatorial problems.                           *
// *                                                                           *
// *   Copyright (C) 2010-2017 Nathan Krislock, Jérôme Malick, Frédéric Roupin *
// *   Multi-threaded version by C.Coti, F.Butelle, E.Leclercq, F. Roupin      *
// *                                                                           *
// *                 http://www-lipn.univ-paris13.fr/BiqCrunch/                *
// *                                           *
// *****************************************************************************
//                                         *
//    This program is free software: you can redistribute it and/or modify     *
//    it under the terms of the GNU General Public License as published by     *
//    the Free Software Foundation, either version 3 of the License, or        *
//    (at your option) any later version.                                      *
//                                                                             *
//    This program is distributed in the hope that it will be useful,          *
//    but WITHOUT ANY WARRANTY; without even the implied warranty of           *
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
//    GNU General Public License for more details.                             *
//                                                                             *
//    You should have received a copy of the GNU General Public License        *
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.    *
//                                                                             *
// *****************************************************************************

#include <math.h>
#include <string.h>

#include "bb.h"
#include "biqcrunch.h"
#include "global.h"


/*
 * Evaluates each triangle inequality (i.e., cut) using the X matrix and
 * fills the Inequality array, List, with at most params.cuts inequalities that
 * are violated by at least params.gapCuts.
 *
 * Also returns the value of the cut that is violated the most by X.
 */
double getViolatedCuts(double *X, int N, Inequality *List, int *ListSize)
{
    int ii, jj, kk, type, ListCount;
    int size = 0;
    int LeastViolatedIneq = 0;
    double LeastViolatedIneqValue = -INFINITY;
    double minAllIneq = INFINITY;
    double testineqvalue;

    // Loop through all inequalities
    size = 0;
    for (type = 1; type <= 4; type++) {
        for (ii = 0; ii < N;  ii++) {
            for (jj = 0; jj < ii; jj++) {
                for (kk = 0; kk < jj; kk++) {

//                    ineqvalue     = eval_ineq(X,     N, type, ii, jj, kk);
                    testineqvalue = eval_ineq(X, N, type, ii, jj, kk);

                    // keep track of the minimum value of testineqvalue
                    minAllIneq = (testineqvalue < minAllIneq) ? testineqvalue : minAllIneq;

                    if (testineqvalue < params.gapCuts)
                    {
                        // (1) put first params.cuts violated inequalities in list and 
                        //     keep track of the least violated inequality
                        if (size < params.cuts) {

                            // add ineq to the end of List
                            List[size].type  = type;
                            List[size].i     = ii;
                            List[size].j     = jj;
                            List[size].k     = kk;
                            List[size].value = testineqvalue;

                            // update LeastViolatedIneq
                            if (size == 0 || testineqvalue > LeastViolatedIneqValue) {
                                LeastViolatedIneqValue = testineqvalue;
                                LeastViolatedIneq = size;
                            }

                            size++;
                        } 
                        else if (testineqvalue < LeastViolatedIneqValue) 
                        {
                            // (2) if you find an inequality is violated more than the 
                            //     least violated inequality, then add it to the list, 
                            //     remove the least violated inequality, and find 
                            //     the new least violated inequality in the list

                            // add ineq to the list, replacing LeastViolatedIneq
                            List[LeastViolatedIneq].type  = type;
                            List[LeastViolatedIneq].i     = ii;
                            List[LeastViolatedIneq].j     = jj;
                            List[LeastViolatedIneq].k     = kk;
                            List[LeastViolatedIneq].value = testineqvalue;

                            // update LeastViolatedIneq
                            LeastViolatedIneqValue = -INFINITY;
                            for (ListCount = 0; ListCount < size; ListCount++) 
                            {
                                // if ineq is less violated
                                if (List[ListCount].value > LeastViolatedIneqValue) 
                                {
                                    LeastViolatedIneqValue = List[ListCount].value;
                                    LeastViolatedIneq = ListCount;
                                }
                            }
                        }
                    }
                } // kk loop
            } // jj loop
        } // ii loop
    } // type loop

    *ListSize = size;

    return minAllIneq;
}


/******************** Update Inequalities ***************************/
double updateInequalities(Problem *PP, double *y, int *NumAdded, int *NumSubtracted, int count) 
{
    int ii, jj, kk, type;
    double ineqvalue; // value of the inequality:  ineqvalue >= 0 if feasible
    int ineq; // index for inequalities
    int yindex; // index for dual multipliers
    int ListCount, ListSize;
    double minAllIneq;
    int added, subtracted;
    int N = PP->n;
    int current_BobId = setget_bobId(0);

    minAllIneq = getViolatedCuts(ListNodeVars[current_BobId].X, N, ListNodeVars[current_BobId].List, &ListSize);

    // Decide which cuts that were previously added to remove
    subtracted = 0;
    yindex = PP->mB + PP->mA;
    int next_ineq = 0;
    for (ineq = 0; ineq < PP->NIneq; ineq++) 
    {
        type = ListNodeVars[current_BobId].Cuts[ineq].type;
        ii   = ListNodeVars[current_BobId].Cuts[ineq].i;
        jj   = ListNodeVars[current_BobId].Cuts[ineq].j;
        kk   = ListNodeVars[current_BobId].Cuts[ineq].k;

        ineqvalue = eval_ineq(ListNodeVars[current_BobId].X, N, type, ii, jj, kk);

        // store the dual multiplier
        ListNodeVars[current_BobId].Cuts[ineq].y = y[yindex];  yindex++;

        // remove inequality if multiplier is small 
        // and inequality is large
        if (ListNodeVars[current_BobId].Cuts[ineq].y < 1e-8 && ineqvalue > ListNodeVars[current_BobId].gradInorm / 10.0)
        {
            subtracted++;
        } 
        else // keep inequality
        {
        	ListNodeVars[current_BobId].Cuts[next_ineq].type = ListNodeVars[current_BobId].Cuts[ineq].type;
        	ListNodeVars[current_BobId].Cuts[next_ineq].i    = ListNodeVars[current_BobId].Cuts[ineq].i;
        	ListNodeVars[current_BobId].Cuts[next_ineq].j    = ListNodeVars[current_BobId].Cuts[ineq].j;
        	ListNodeVars[current_BobId].Cuts[next_ineq].k    = ListNodeVars[current_BobId].Cuts[ineq].k;
        	ListNodeVars[current_BobId].Cuts[next_ineq].y    = ListNodeVars[current_BobId].Cuts[ineq].y;
            next_ineq++;
        }
    }
    PP->NIneq -= subtracted;

    // Add List to Cuts
    added = 0;
    for (ListCount = 0; ListCount < ListSize; ListCount++) 
    { 
        // Stop if we have reached the maximum number of cuts we can add
        if (next_ineq == MaxNineqAdded) break;

        // Check if inequality is already included in Cuts
        int found_ineq = 0;
        for (ineq = 0; ineq < PP->NIneq; ineq++) 
        {
            if (ListNodeVars[current_BobId].Cuts[ineq].type == ListNodeVars[current_BobId].List[ListCount].type &&
            		ListNodeVars[current_BobId].Cuts[ineq].i == ListNodeVars[current_BobId].List[ListCount].i &&
            		ListNodeVars[current_BobId].Cuts[ineq].j == ListNodeVars[current_BobId].List[ListCount].j &&
            		ListNodeVars[current_BobId].Cuts[ineq].k == ListNodeVars[current_BobId].List[ListCount].k)
            {
                found_ineq = 1;
            }
        }
        // If inequality not already in Cuts, add it to Cuts
        if (!found_ineq)
        {
        	ListNodeVars[current_BobId].Cuts[next_ineq].type = ListNodeVars[current_BobId].List[ListCount].type;
        	ListNodeVars[current_BobId].Cuts[next_ineq].i    = ListNodeVars[current_BobId].List[ListCount].i;
        	ListNodeVars[current_BobId].Cuts[next_ineq].j    = ListNodeVars[current_BobId].List[ListCount].j;
        	ListNodeVars[current_BobId].Cuts[next_ineq].k    = ListNodeVars[current_BobId].List[ListCount].k;
        	ListNodeVars[current_BobId].Cuts[next_ineq].y    = 0.0;
            next_ineq++;
            added++;
        }
    }
    PP->NIneq += added;

    *NumAdded = added;
    *NumSubtracted = subtracted;

    return minAllIneq;
}


/************************ eval_ineq *********************************/
double eval_ineq(double *XX, int N, int type, int ii, int jj, int kk) 
{
    double ineqvalue;

    // compute the inequality
    switch (type) {
        case 1:
            ineqvalue =  XX[ii+jj*N] + XX[ii+kk*N] + XX[jj+kk*N] + 1.0;
            break;
        case 2:
            ineqvalue =  XX[ii+jj*N] - XX[ii+kk*N] - XX[jj+kk*N] + 1.0;
            break;
        case 3:
            ineqvalue = -XX[ii+jj*N] + XX[ii+kk*N] - XX[jj+kk*N] + 1.0;
            break;
        default: // case 4 :
            ineqvalue = -XX[ii+jj*N] - XX[ii+kk*N] + XX[jj+kk*N] + 1.0;
    }

    return ineqvalue;
}


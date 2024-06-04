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

#include <stdio.h>
#include <sys/stat.h>
#include <string.h>

#include "bb.h"
#include "biqcrunch.h"
#include "global.h"

FILE *final_output;
char gen_output_path[200];

/*
 * Open the output file.
 * If the output file exist it add a number a counter at the end of the filename
 * @param instance: name of the instance to use as base for the output file name
 * @return 1 if succeed in the creation of the output file, 0 otherwise
 */
int createOutputFile(char *instance) 
{
    int i, counter = 1;
    char output_path[200];
    struct stat buffer;

    sprintf(gen_output_path, "%s", instance);
    sprintf(output_path, "%s.output", instance);
    
    // Check if the file already exists, if so append _<NUMBER> to the end of the 
    // output file name
    while (stat(output_path, &buffer) == 0)
        sprintf(output_path, "%s.output_%d", instance, counter++);

    // Display final output file
    printf("Output file: %s \n", output_path);

    // Open final output file and temporary thread output files
    final_output = fopen(output_path, "w");
    if (!final_output)  return 0;

    for(i=0;i<params.nbProcs;i++) {
    	sprintf(output_path, "%s_thread%d", instance, i);
    	output[i] = fopen(output_path, "w");
    }

    for(i=0;i<params.nbProcs;i++) {
		if (!output[i])
			return 0;
    }

    // Print BiqCrunch banner
    fprintf(final_output,BANNER);

    // Report the values of the parameters used
    fprintf(final_output, "BiqCrunch Parameters:\n");
#define P(type, name, format, def_value)\
    fprintf(final_output, "%20s = "format"\n", #name, params.name);
    PARAM_FIELDS;
#undef P

    return 1;
}


/*
 * Close the output file
 */
void closeOutputFile() 
{
    fclose(final_output);
}


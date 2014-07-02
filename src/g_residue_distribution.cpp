/*
 * This file is part of g_resid_distrib
 *
 * Author: Rajendra Kumar
 * Copyright (C) 2014  Rajendra Kumar
 *
 * g_resid_distrib is a free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * g_resid_distrib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with g_resid_distrib.  If not, see <http://www.gnu.org/licenses/>.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 */

#include <iostream>
#include <vector>

#include "typedefs.h"
#include "macros.h"
#include "gstat.h"
#include "smalloc.h"
#include "copyrite.h"
#include "pdbio.h"
#include "gmx_fatal.h"
#include "xvgr.h"
#include "matio.h"
#include "index.h"
#include "tpxio.h"
#include "rmpbc.h"
#include "do_fit.h"

#include "g_resid_distrib.hpp"


void CopyRightMsg() {

    const char *copyright[] = {
            "                                                                        ",
            "               :-)  g_resid_distrib (-:                                     ",
            "                                                                        ",
            "               Author: Rajendra Kumar                                  ",
            "                                                                        ",
            "         Copyright (C) 2014  Rajendra Kumar                             ",
            "                                                                        ",
            "                                                                        ",
            "g_resid_distrib is a free software: you can redistribute it and/or modify      ",
            "it under the terms of the GNU General Public License as published by    ",
            "the Free Software Foundation, either version 3 of the License, or       ",
            "(at your option) any later version.                                     ",
            "                                                                        ",
            "g_resid_distrib is distributed in the hope that it will be useful,             ",
            "but WITHOUT ANY WARRANTY; without even the implied warranty of          ",
            "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           ",
            "GNU General Public License for more details.                            ",
            "                                                                        ",
            "You should have received a copy of the GNU General Public License       ",
            "along with g_resid_distrib.  If not, see <http://www.gnu.org/licenses/>.       ",
            "                                                                        ",
            "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS     ",
            "\"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT     ",
            "LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR   ",
            "A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT    ",
            "OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,   ",
            "SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED",
            "TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR  ",
            "PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF  ",
            "LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING    ",
            "NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      ",
            "SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.            ",
            "                                                                        ",
            "                           :-)  g_resid_distrib (-:                            ",
            "                                                                        ",
            "                                                                        "
    };
    int i = 0;
    char *str;
    for(i=0; i<35; i++) {
        str = strdup(copyright[i]);
        fprintf(stderr,"%s\n", str);
    }
}


void set_group_mass(int index_size, atom_id *index, t_topology *top, gmx_bool bTPR, atoms_group *group)	{
	int i=0;
	real *mass;
	mass = new real [index_size];

	for(i=0;i<index_size;i++)
		if(bTPR)
			mass[i] = top->atoms.atom[index[i]].m;
		else
			mass[i] = 1.0;
	group->set_mass(mass);
	delete [] mass;
}

void set_group_coord(int index_size, atom_id *index, t_topology *top, rvec *x, atoms_group *group)	{
	int i=0, j =0;;
	rvec *coord;

	coord = new rvec [index_size];

	for(i=0;i<index_size;i++)
		for(j=0;j<DIM;j++)
			coord[i][j] = x[index[i]][j];

	group->set_coord(coord);
	delete [] coord;
}


int main(int argc, char *argv[]) {
	  const char *desc[] = {
			  "This program calculates the positions of center of masses of the selected",
			  "groups along three coordinates. This tool is useful to calculate the position",
			  "of selected residues along the channel-axis of the protein channel. Further,",
			  "output files could be use to obtain the distribution of residues' position",
			  "during the MD simulations. If -fit is enabled, the molecule will be translated,",
			  "and centered at the origin."
	  };

	  gmx_bool bFit=FALSE, bM=TRUE;
	  int num_grps=4;
	  output_env_t oenv;

	  t_pargs pa[] = {
			  { "-fit", FALSE, etBOOL, {&bFit}, "Fitting frames on reference structure" },
			  { "-mwa", TRUE, etBOOL, {&bM},  "Mass weighted fitting" },
			  { "-ng", TRUE, etINT, {&num_grps},  "Number of groups of atoms or residues" }
	    };

	  t_filenm   fnm[] = {
	     { efTRX, "-f",   NULL,      ffREAD },
	     { efTPS, "-s",   NULL,      ffREAD },
	     { efNDX, "-n",   NULL,      ffOPTRD },
	     { efDAT, "-x",  "com_coord_x", ffOPTWR },
	     { efDAT, "-y",  "com_coord_y", ffOPTWR },
	     { efDAT, "-z",  "com_coord_z", ffOPTWR },
	   };


#define NFILE asize(fnm)
	  int npargs;
	  CopyRightMsg();
	  npargs = asize(pa);
	  parse_common_args(&argc,argv, PCA_CAN_TIME | PCA_TIME_UNIT | PCA_BE_NICE , NFILE,fnm,npargs,pa, asize(desc),desc,0,NULL,&oenv);

	  //topology stuffs
	  t_topology top;
	  int        ePBC;
	  t_atoms    *atoms;
	  char title[256];


	  //Trajectory stuffs
	  t_trxstatus *status;
	  matrix     box;
	  int natoms=0, nframe=0;
	  real time;

	  //Index stuffs
	  atom_id    **index;
	  int  *index_size;
	  char **grp_names;

	  //Coordinates
	  rvec *x, *xtop;

	  //PBC stuff
	  gmx_rmpbc_t  gpbc=NULL;

	  //Structure fitting stuffs
	  atom_id *fit_index;
	  int fit_index_size=0;
	  real *w_rls=NULL;
	  char *fit_name;

	  //XVG files output stuffs
	  const char *fnX, *fnY, *fnZ;
	  FILE *fX, *fY, *fZ;

	  //General variables
	  int i=0;

	  //Vraiable for analysis
	  std::vector<atoms_group> groups;
	  gmx_bool bTPR=TRUE;
	  gmx_bool bX=FALSE, bY=FALSE, bZ=FALSE;

	  //Checking Output options
	  if(opt2fn("-x",NFILE,fnm)!=NULL)
		  bX=TRUE;
	  if(opt2fn("-y",NFILE,fnm)!=NULL)
		  bY=TRUE;
	  if(opt2fn("-z",NFILE,fnm)!=NULL)
		  bZ=TRUE;
	  if(!bX && !bY && !bZ)
		  gmx_fatal(FARGS,"No output files given\n");


	  //Reading tpr file
	  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&xtop,NULL,box,FALSE);
	  atoms=&(top.atoms);

	  if (fn2bTPX(ftp2fn(efTPS,NFILE,fnm)) == FALSE)
		  bTPR=FALSE;

	  //Initialization of fitting
	  if (bFit)		{

		  if (fn2bTPX(ftp2fn(efTPS,NFILE,fnm)) == FALSE)	{
			  printf("\n\nWARNING: Not a TPR file.... Switching off mass weighted fitting...\n\n");
			  bM = FALSE;
		  }

		  fprintf(stderr, "\nChoose a group for the least squares fit\n");
		  get_index(atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&fit_index_size,&fit_index,&fit_name);
		  if (fit_index_size < 3)
			  gmx_fatal(FARGS,"Need >= 3 points to fit!\n");

		  snew(w_rls,atoms->nr);

		  if(bM)
			  for(i=0; (i<fit_index_size); i++)
				  w_rls[ fit_index[i] ]=atoms->atom[ fit_index[i] ].m;
		  else
			  for(i=0; (i<fit_index_size); i++)
				  w_rls[ fit_index[i] ]=1.0;

		  reset_x(fit_index_size,fit_index,top.atoms.nr,NULL,xtop,w_rls);
	  }

	   //Getting index
	   snew(grp_names,num_grps);
	   snew(index_size,num_grps);
	   snew(index,num_grps);
	   get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),num_grps,index_size,index,grp_names);

	   //Initialization of atom group objects and setting mass of atoms
	   for (i=0;i<num_grps;i++){
		   groups.push_back(atoms_group(index_size[i]));
		   set_group_mass(index_size[i], index[i], &top, bTPR, &groups[i]);
	   }


	   //Reading first frame
	   natoms = read_first_x(oenv,&status,ftp2fn(efTRX,NFILE,fnm),&time,&x,box);
	   if(bFit)	{
		   reset_x(fit_index_size,fit_index,top.atoms.nr,NULL,x,w_rls);
		   do_fit(natoms,w_rls,xtop,x);
	   }

	   //Opening files
	   if(bX)	{
		   fnX = opt2fn("-x",NFILE,fnm);
		   fX = ffopen(fnX,"w");
		   fprintf(fX,"#Time");
		   for(i=0;i<num_grps;i++)
			   fprintf(fX,"  %s",grp_names[i]);
		   fprintf(fX,"\n");
	   }
	   if(bY)	{
		   fnY = opt2fn("-y",NFILE,fnm);
		   fY = ffopen(fnY,"w");
		   fprintf(fY,"#Time");
		   for(i=0;i<num_grps;i++)
			   fprintf(fY,"  %s",grp_names[i]);
		   fprintf(fY,"\n");

	   }
	   if(bZ)	{
		   fnZ = opt2fn("-z",NFILE,fnm);
		   fZ = ffopen(fnZ,"w");
		   fprintf(fZ,"#Time");
		   for(i=0;i<num_grps;i++)
			   fprintf(fZ,"  %s",grp_names[i]);
		   fprintf(fZ,"\n");
	   }

	   gpbc = gmx_rmpbc_init(&top.idef,ePBC,natoms,box);

	 do{
		 gmx_rmpbc(gpbc,natoms,box,x);

		 if (bFit)	{
			 reset_x(fit_index_size,fit_index,top.atoms.nr,NULL,x,w_rls);
			 do_fit(natoms,w_rls,xtop,x);

		 }

		 if(bX)
			 fprintf(fX,"%15.5lf",time);
		 if(bY)
			 fprintf(fY,"%15.5lf",time);
		 if(bZ)
			 fprintf(fZ,"%15.5lf",time);


		 //////////////////////////////////////
		 //Functions for main calculations
		 for (i=0;i<num_grps;i++)
			 set_group_coord(index_size[i], index[i], &top, x, &groups[i]);



		 for (i=0;i<num_grps;i++)
			 groups[i].calculate_com();

		 //////////////////////////////////////

		 if(bX)	{
			for(i=0;i<num_grps;i++)
				fprintf(fX,"%15.5lf",groups[i].com_x[0]);
			fprintf(fX,"\n");
		 }
		 if(bY)	{
			for(i=0;i<num_grps;i++)
				fprintf(fY,"%15.5lf",groups[i].com_x[1]);
			fprintf(fY,"\n");
		 }

		 if(bZ)	{
			 for(i=0;i<num_grps;i++)
				fprintf(fZ,"%15.5lf",groups[i].com_x[2]);
			 fprintf(fZ,"\n");
		 }

		 nframe++;
	   } while (read_next_x(oenv,status,&time,natoms,x,box));

	 fprintf(stdout, "Thanks for using g_resid_distrib!!!\n");

	return 0;
}

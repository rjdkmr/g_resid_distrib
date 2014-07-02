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

//Included for compatibility with gromacs variables
#include "typedefs.h"
#include <iostream>
#include <vector>

#include "g_resid_distrib.hpp"

//Constructor
atoms_group::atoms_group(int num)	{
	int i =0,j;
	nr = num;
	cumul_mass=0;
	x = new real *[nr];
	for(i=0;i<nr;i++)	{
		x[i] = new real[DIM];
		for(j=0;j<DIM;j++)	{
			x[i][j] =0.0;
		}
	}
	m = new real[nr];
	for(j=0;j<nr;j++)
		m[i] =0.0;
}

int atoms_group::set_mass(real *mass)	{
	int i = 0;

	for (i=0;i<nr;i++)	{
		m[i] = mass[i];
		cumul_mass += mass[i];
	}
	return 0;
}

int atoms_group::set_coord(rvec *coord)	{
	int i = 0, j= 0;
	for (i=0;i<nr; ++i)
		for(j=0;j<DIM;++j)
			x[i][j] = coord[i][j];
	return 0;
}

void atoms_group::calculate_com(){
	int i = 0, d = 0;

	for(d=0;(d<DIM);d++)
		com_x[d] = 0;

    for(d=0;(d<DIM);d++)	{
  	  for(i=0;(i<nr);i++)
  		  com_x[d] += x[i][d] * m[i];

  	com_x[d] = com_x[d]/cumul_mass;
    }

}

int atoms_group::get_num_atoms()	{
	return nr;
}

//atoms_group::~atoms_group()	{
	//int i =0;
	//for(i=0;i<DIM;i++)
	//	delete [] x[i];
	//delete [] x;
	//delete [] m;
//}


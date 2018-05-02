/*
	Auxiliary functions to work with matrices

AUTHORS:

- Oleksandr Kazymyrov (2013-08-03): initial version

*/

/*****************************************************************************
 *       Copyright (C) 2013 Oleksandr Kazymyrov <okazymyrov@gmail.com>
 *
 *  Distributed under the terms of the GNU General Public License (GPL)
 *  as published by the Free Software Foundation; either version 2 of
 *  the License, or (at your option) any later version.
 *                  http://www.gnu.org/licenses/
 ****************************************************************************/

#include "matrix.h"

int	init_matrix(char ***m, int rows, int columns)
{
	int r=0;

	(*m)=(char**)calloc(rows,sizeof(char*));

	for(r=0;r<rows;r++)
	{
		(*m)[r]=(char*)calloc(columns,sizeof(char));
	}

	return 0;
}

int	free_matrix(char **m, int rows)
{
	int r=0;

	if(m)
	{
		for(r=0;r<rows;r++)
		{
			if(m[r])
				free(m[r]);
		}

		free(m);
	}

	return 0;
}

int	xor_vectos(char *a,char *b,char *d,int length)
{
	int i=0;

	if(!a || !b || !d)
		return -1;

	for(i=0;i<length;i++)
	{
		d[i]=a[i]^b[i];
	}

	return 0;
}

int swap_matrix_row(char **m,int r1, int r2)
{
	char *r=0;

	if(!m)
		return -1;

	r=m[r1];
	m[r1]=m[r2];
	m[r2]=r;

	return 0;
}

int	row_reduce_matrix(char **m, int rows, int columns)
{
	int c=0, r=0, found=0, delta=0;

	for(c=0;c<columns;c++)
	{
		if( c == rows )
			break;

		if( delta + c >= min(rows,columns) )
			break;

		found=0;

		if( m[c][c+delta] != 1 )
		{
			for(r=c+1;r<rows;r++)
			{
				if( (m[r][c+delta] == 1) && (m[c][c+delta] != 1) )
				{
					found=1;
					swap_matrix_row(m,r,c);
					break;
				}
			}

			if( r == rows )
			{
				delta++;
				c--;
				continue;
			}
		}
		else
			found=1;

		if( found )
		{
			for(r=0;r<rows;r++)
			{
				if( (r != c ) && ( m[r][c+delta] == 1) )
				{
					xor_vectos(m[r],m[c],m[r],columns);
				}
			}
		}
	}
	
	return 0;
}

int	rank_matrix(char **m, int rows, int columns)
{
	char **t=NULL, *zero=NULL;
	int r=0, rank=rows;

	if(!m )
		return -1;

	init_matrix(&t,rows,columns);
	zero=(char*)calloc(columns,sizeof(char));

	for(r=0;r<rows;r++)
	{
		memcpy(t[r],m[r],columns);
	}

	row_reduce_matrix(t,rows,columns);

	for(r=rows-1;r>=0;r--)
	{
		if(!memcmp(zero,t[r],columns))
			rank--;
		else
			break;
	}

	free_matrix(t,rows);
	free(zero);

	return rank;
}

int	inverse_matrix(char **m, char **inv, int rc)
{
	char **t=NULL;

	int r=0;

	if(!m || !inv)
		return -1;

	init_matrix(&t,rc,2*rc);

	for(r=0;r<rc;r++)
	{
		memcpy(t[r],m[r],rc);
		t[r][rc+r]=1;
	}

	row_reduce_matrix(t,rc,2*rc);
	
	if(t[rc-1][rc-1] == 0)
		return -2;


	for(r=0;r<rc;r++)
	{
		memcpy(inv[r],(void*)&t[r][rc],rc);
	}

	free_matrix(t,rc);

	return 0;
}

int	is_inverse(char **m, int rc)
{
	int ret=0;
	char **inv=NULL;

	init_matrix(&inv,rc,rc);

	if ( inverse_matrix(m,inv,rc) != 0)
		ret=0;
	else
		ret=1;

	free_matrix(inv,rc);

	return ret;
}

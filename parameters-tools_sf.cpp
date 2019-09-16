/* Copyright (C) 2005-2008 University of Washington
   Written by Zhirong Bao and Dan Blick
   This file is part of starrynite.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.  */

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>

#include "parameters-tools_sf.h"
#include "utility_sf.h"
#include "types_sf.h"

/* BEGIN FUNCTION PROTOTYPES */

/*
 * static void CleanUpFrame (U8 ** image_stack, FRAME_t * frame);
 */

/* END FUNCTION PROTOTYPES */


void
MakeTools (TOOL_t * const tools,  int upper,double z_factor)
{
  int i, j, k;
  int length;

  int r, slices, dist;
  double zf2, r2;
  int *radii, **xrange;

  
  tools->spheres = (int ***) malloc_exit (sizeof (int **) * upper);
  tools->s_layers = (int *) malloc_exit (sizeof (int) * upper);
  tools->s_radii = (int **) malloc_exit (sizeof (int *) * upper);
  zf2 = z_factor * z_factor;

  for (k = 2; k < upper; k += 2)
  {
    const int R = k / 2;
    const float R2 = 1.0 * R * R;

    slices = (int) (1.0 * R / z_factor);
    radii = (int *) malloc_exit (sizeof (int) * (2 * slices + 1));
    xrange = (int **) malloc_exit (sizeof (int *) * (2 * slices + 1));
    //printf("slices %d\n",slices);
    for (i = slices; i >= 0; --i)
    {
      r = (int) sqrt (R2 - zf2 * i * i);
      radii[slices - i] = r;
      xrange[slices - i] = (int *) malloc_exit (sizeof (int) * (2 * r + 1));
      r2 = 1.0 * r * r;

      for (j = r; j > 0; j--)
      {
        dist = (int) sqrt (r2 - j * j);
        xrange[slices - i][r - j] = dist;
        xrange[slices - i][r + j] = dist;
      }
      xrange[slices - i][r] = r;

      if (i != 0)
      {
        radii[slices + i] = r;
        xrange[slices + i] = (int *) malloc_exit (sizeof (int) * (2 * r + 1));
        memcpy (xrange[slices + i], xrange[slices - i],
                sizeof (int) * (2 * r + 1));
      }
    }                           /* done going through the slices */

    tools->spheres[k] = xrange;
    tools->s_layers[k] = slices;
    tools->s_radii[k] = radii;

    if (k == (k / 2) * 2 && k < upper - 1)
    {
      /* index for last element is odd, so it's taken care of by the previous k */
      tools->spheres[k + 1] = xrange;
      tools->s_layers[k + 1] = slices;
      tools->s_radii[k + 1] = radii;
    }
    else if (k > 2)
    {
      tools->spheres[k - 1] = xrange;
      tools->s_layers[k - 1] = slices;
      tools->s_radii[k - 1] = radii;
    }
  }
}



void
CleanUpTools (TOOL_t * const tools,int upper)
{
  int i;
  for (i = 2; i < upper; i += 2)
  {
    int j;
    const int slices = tools->s_layers[i];
    int **xrange = tools->spheres[i];
    for (j = slices; j >= 0; j--)
    {
      free (xrange[slices - j]);
      if (j)
        free (xrange[slices + j]);
    }
    free (xrange);
    free (tools->s_radii[i]);
  }
  free (tools->spheres);
  free (tools->s_layers);
  free (tools->s_radii);
}

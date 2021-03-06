/*  Copyright (C) 2003-2007  CAMP
 *  Please see the accompanying LICENSE file for further information. */

#include <string.h>
#include "bmgs.h"

void Z(bmgs_paste)(const T* a, const int sizea[3],
		   T* b, const int sizeb[3], const int startb[3])
{
  b += startb[2] + (startb[1] + startb[0] * sizeb[1]) * sizeb[2];
#pragma omp for schedule(static) collapse(2)
  for (int i0 = 0; i0 < sizea[0]; i0++)
    {
      for (int i1 = 0; i1 < sizea[1]; i1++)
	    {
          for (int i2 = 0; i2 < sizea[2]; i2++)
            {
              int bstart = i1  * sizeb[2] +
                          i0 * sizeb[2] * sizeb[1] + i2;
              int astart = (i0 * sizea[1]  + i1) * sizea[2] + i2;
               b[bstart] = a[astart];
            }
        }
    }
}


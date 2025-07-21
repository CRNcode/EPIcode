/* Copyright (C) 2024-5, Murad Banaji
 *
 * This file is part of EPIcode, for compartmental models in epidemiology
 *
 * EPIcode is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2, 
 * or (at your option) any later version.
 *
 * EPIcode is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with EPIcode: see the file COPYING.  If not, write to 
 * the Free Software Foundation, Inc., 59 Temple Place - Suite 330, 
 * Boston, MA 02111-1307, USA. 

 */

// Wrapper file for the "resample" routine

#include "Epi.h"

int main(int argc, char *argv[]){
  double sampling;
  if(argc<4){
    fprintf(stderr, "Insufficient parameters. These should be:\n\t(1) input file\n\t(2) output file\n\t(3) sampling rate\nEXITING.\n");exit(0);
  }
  sampling=atof(argv[3]);

  resample(argv[1], argv[2], sampling);
  
  return 0;
}

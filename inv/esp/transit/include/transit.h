/*
 * transit.h - Common headers for the Transit program.
 *
 * Copyright (C) 2003 Patricio Rojo (pato@astro.cornell.edu)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 */

#ifndef _TRANSIT_H
#define _TRANSIT_H

#include <stdarg.h>
#include <math.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <pu/sampling.h>
#include <pu/profile.h>
#include <pu/iomisc.h>
#include <pu/numerical.h>
#include <pu/xmalloc.h>
#include <strings.h>
#include <stdlib.h>
#include <stdio.h>
#include <alloca.h>
#ifdef _USE_GSL
#include <gsl/gsl_spline.h>
#endif

#define compattliversion 3

#include <flags_tr.h>
#include <constants_tr.h>


/*****   Types     *****/
#define PREC_NSAMP int		/* Type for radius and wavelength
				   indices */
#define PREC_NREC long		/* Type for record indices */
#define PREC_ZREC double	/* Type for the partition info  */
#define PREC_LNDATA double	/* Type for the line data output */
#define PREC_RES double     	/* Type for every partial result */
#define PREC_ATM double		/* Type for atmospheric data */
#define PREC_CS  double		/* Type for cross-section */

/*****   Macros   *****/
static __inline__ void
printextprogress(long wi, long wnn)
{
}

static __inline__ double
stateeqnford(_Bool mass,	/* Mass abundance? (as opposed to
				   abundance by number. */
	     double q,		/* abundance */
	     double ma,		/* Average molecular weight */
	     double mi,		/* Molecular weight of the particular
				   specie */
	     double p,		/* Pressure */
	     double t)		/* Temperature */
{
 if(mass)
   return AMU * q * ma * p / KB / t;
 return AMU * q * mi * p / KB / t;
}

#define transitassert(a,...) if(a) transiterror(TERR_CRITICAL,__VA_ARGS__)
#define transitprint(thislevel, verblevel, ...) do{                         \
  if(thislevel <= verblevel)  fprintf(stderr,__VA_ARGS__); }while(0)
#define transitacceptflag(transit,hint,flag) do{                            \
        transit|=hint&flag;hint&=~(flag);}while(0)
#define transitallocerror(nmb)                                              \
        transiterror(TERR_CRITICAL,                                         \
	             "transit:: %s: Allocation failed for %i allocation\n"  \
	             "units in line %i. Impossible to continue.\n"          \
	             ,__FILE__,nmb,__LINE__)
#define free_null(x) do{free(x);x=NULL;}while(0)

#ifdef NODEBUG_TRANSIT
#define transitDEBUG(...) ((void)0)
#define transitASSERT(...) ((void)0)
#else
#define free(x)  do{free(x);x=NULL;}while(0)
#define transitASSERT(a,...) if(a) transiterror(TERR_CRITICAL,__VA_ARGS__)
#define transitDEBUG(...) transitprint(__VA_ARGS__)
#endif

#define maxeisoname 20

extern int transit_nowarn;
extern int verblevel;              /* verbose level, greater than 10 
				      is only for debuging */
extern int maxline;
extern int version;
extern int revision;

enum isodo {unclear=0,atmfile,ignore,fixed,factor};

#include <structures_tr.h>

extern const transit_ray_solution slantpath;


/***** Prototypes *****/
#include <proto_transit.h>
#include <proto_readlineinfo.h>
#include <proto_readatm.h>
#include <proto_transitstd.h>
#include <proto_makesample.h>
#include <proto_extinction.h>
#include <proto_idxrefraction.h>
#include <proto_tau.h>
#include <proto_argum.h>
#include <proto_geometry.h>
#include <proto_observable.h>

#endif /* _TRANSIT_H */

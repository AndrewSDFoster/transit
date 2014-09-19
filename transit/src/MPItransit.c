#include <mpi.h>
#include <transit.h>

void main(int argc,      /* Number of variables                              */
          char **argv){  /* Variables                                        */

  /* Initialization of data structure's pointers. Note that s_lt is not
     assigned because the array that is going to point has not been
     initialized yet.                                                       */
  struct transit transit;
  long itr=0;
  struct timeval tv;
  double t0=0.0;
  int i;

  transitprint(1,2, "TRAN FLAG 00: This is transit\n");
  /* MPI variables:                                                         */
  int root=0, size, rank;
  double *dummy;

  double *input;  /* In*/
  /* Initialize MPI Communicator:                                          */
  MPI_Comm comm;
  MPI_Init(NULL, NULL);
  MPI_Comm_get_parent(&comm);

  /* Debug mode only:                                                       */
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  memset(&transit, 0, sizeof(struct transit));
  verblevel=2;

  /* Process the command line arguments:                                    */
  fw(processparameters, !=0, argc, argv, &transit);
  t0 = timecheck(verblevel, itr,  0, "processparameters", tv, t0);

  transitprint(1, verblevel, "verblevel: %d\n", verblevel);
  /* Silence threads with rank > 0:                                         */
  if (rank != 0)
    verblevel = 0;

  /* Accept all general hints:                                              */
  fw(acceptgenhints, !=0, &transit);
  /* Presentation:                                                          */
  printintro();

  /* No program warnings if verblevel is 0 or 1:                            */
  if(verblevel<2)
    transit_nowarn = 1;

  /* Make wavenumber binning:                                               */
  fw(makewnsample0, <0, &transit);
  t0 = timecheck(verblevel, itr,  1, "makewnsample0", tv, t0);
  if(fw_status>0)
    transitprint(7, verblevel,
                 "makewnsample() modified some of the hinted "
                 "parameters according to returned flag: 0x%lx.\n",
                 fw_status);

  /* Send the number of wavelength samples:                                 */
  int nspec = (int)transit.wns.n;
  MPI_Barrier(comm);
  MPI_Gather(&nspec, 1, MPI_INT,
             dummy,  5, MPI_INT, root, comm);
  transitprint(1, verblevel, "TRAN FLAG 56: nspec is %d\n", nspec);

  /* Get the number of parameters to receive on each iteration:             */
  MPI_Barrier(comm);
  int npars[2];
  MPI_Bcast(npars, 2, MPI_INT, root, comm);
  int ntransit = npars[0];
  int niter = npars[1];
  transitprint(1, verblevel, "TRAN FLAG 60: The number of transit params "
                             "is %d, for %d iterations\n", ntransit, niter);

  /* Send the wavenumber array:                                             */
  MPI_Barrier(comm);
  MPI_Gather(transit.wns.v, nspec, MPI_DOUBLE,
             dummy,         1,     MPI_DOUBLE, root, comm);

  /* Read Atmosphere information:                                           */
  fw(getatm, !=0, &transit);
  t0 = timecheck(verblevel, itr,  2, "getatm", tv, t0);

  /* Read line info:                                                        */
  fw(readlineinfo, !=0, &transit);
  t0 = timecheck(verblevel, itr,  3, "readlineinfo", tv, t0);

  /* Make radius binning and interpolate data to new value:               */
  fw(makeradsample, <0, &transit);
  t0 = timecheck(verblevel, itr,  4, "makeradsample", tv, t0);
  if(fw_status>0)
    transitprint(7, verblevel, "makeradsample() modified some of the hinted "
                               "parameters. Flag: 0x%lx.\n", fw_status);

  /* Hack: add an index to f_out filename to get different files
     and compare them:                                                      */
  int fout_len = (int)strlen(transit.f_out);  /* Length of tr.f_out         */
  char *dot    = strchr(transit.f_out, '.');  /* Search for '.' char        */
  int dot_pos  = fout_len - (int)strlen(dot); /* Position of dot            */
  char fout[fout_len+1];         /* Copy of tr.f_out                        */
  char str_iter[3];              /* String of iteration number              */
  strcpy(fout, transit.f_out);
  strcpy(fout+dot_pos+1, dot);
  strncpy(fout+dot_pos, "0", 1);

  /* Calculate opacity grid:                                                */
  fw(opacity, <0, &transit);
  t0 = timecheck(verblevel, itr,  5, "opacity", tv, t0);

  /* The array received from MPI:                                           */
  input = (double *)calloc(ntransit, sizeof(double));

  /* MCMC Main loop:                                                        */
  for (; niter >= 0; niter--){
    t0 = timecheck(verblevel, itr,  8, "Start loop", tv, t0);

    transitprint(1, verblevel, "TRAN FLAG 71: Prepare to receive (%d) data\n",
                               ntransit);
    /* Scatter (receive) atm. profiles from BART:                           */
    MPI_Barrier(comm);
    MPI_Scatter(dummy, 1,        MPI_DOUBLE,
                input, ntransit, MPI_DOUBLE, root, comm);
    /* Modify transit temperature and abundance profiles:                   */
    transitprint(1, verblevel, "TRAN FLAG 72: The number of layers is: %li\n"
                               "              while ntransit is:       %d\n",
                               transit.ds.at->rads.n, ntransit);
    transitprint(1, verblevel, "TRAN FLAG 72.5: INPUT=%.1f %.1f %.1f\n",
                 input[0], input[1], input[2]);
    fw(reloadatm, <0, &transit, input);

    transitprint(1, verblevel, "TRAN FLAG 73: The abundance factor is: %f\n",
                               input[100]);

    /* Compute sampling of impact parameter:                                */
    fw(makeipsample, <0, &transit);
    t0 = timecheck(verblevel, itr,  6, "makeipsample", tv, t0);
    if(fw_status>0)
      transitprint(7, verblevel, "makeipsample() modified some of the hinted "
                   "parameters. Flag: 0x%lx.\n", fw_status);
 
    /* Print sampling info:                                                 */
    if (niter == 0){
      fw(outsample, !=0, &transit);
      t0 = timecheck(verblevel, itr,  7, "outsample", tv, t0);
    }

    /* Initialize CIA:                                                      */
    fw(interpolatecia, !=0, &transit);
    t0 = timecheck(verblevel, itr,  9, "interpolatecia", tv, t0);
 
    /* Compute index of refraction:                                         */
    fw(idxrefrac, !=0, &transit);
    t0 = timecheck(verblevel, itr,  10, "idxrefrac", tv, t0);
 
    /* Calculate extinction coefficient:                                    */
    fw(extwn, !=0, &transit);
    t0 = timecheck(verblevel, itr, 11, "extwn", tv, t0);
 
    transitprint(1, verblevel, "TRAN FLAG 73.5: pre-transit calc\n");
    sleep(1);
    sprintf(str_iter, "%li", itr);
    transitprint(1, verblevel, "TRAN FLAG 73.6: pre-transit calc\n");
    sleep(1);
    strncpy(fout+dot_pos, str_iter, 1);
    transitprint(1, verblevel, "TRAN FLAG 73.7: pre-transit calc\n");
    sleep(1);
    strcpy(transit.f_out, fout);
    transitprint(1, verblevel, "TRAN FLAG 73.9: File output name: %s\n",
                 fout);

    transitprint(1, verblevel, "TRAN FLAG 74: pre-transit calc\n");
    /* Ray solutions choice:                                                */
    RaySol path=transit.ds.th->path;

    /* Calculates optical depth for eclipse                                 */
    if(path == eclipse){
      transitprint(1, verblevel, "\nCalculating eclipse:\n");

      /* Angle number                                                       */
      struct transithint *th = transit.ds.th;
      long int an = th->ann;

      /* Sets intensity grid:                                               */
      fw(intens_grid, !=0, &transit);
      for(i = 0; i < an; i++){
        /* Fills out angle index                                            */
        transit.angleIndex = i;

        fw(tau_eclipse, !=0, &transit);
        t0 = timecheck(verblevel, itr, 12, "tau eclipse", tv, t0);
  
        /* Calculates eclipse intensity:                                    */
        /* In cgs units erg/s/sr/cm                                         */
        fw(emergent_intens, !=0, &transit);
        t0 = timecheck(verblevel, itr, 13, "emergent intensity", tv, t0);
      }

      /* Calculates flux  erg/s/cm                                          */
      fw(flux, !=0, &transit);
      t0 = timecheck(verblevel, itr, 14, "flux", tv, t0);

      /* Free no longer needed memory                                       */
      freemem_intensityGrid(transit.ds.intens, &transit.pi);
    }

    /* Calculate optical depth for transit:                                 */
    else{
      transitprint(1, verblevel, "\nCalculating transit:\n");
      sleep(2);
      fw(tau, !=0, &transit);
      t0 = timecheck(verblevel, itr, 12, "tau", tv, t0); 

      /* Calculates transit modulation:                                     */
      fw(modulation, !=0, &transit);
      t0 = timecheck(verblevel, itr, 13, "modulation", tv, t0);
   }
 
    transitprint(1, verblevel, "TRAN FLAG 76: The spectrum size is: %d\n",
                                nspec);
    /* Scatter the model back to BART:                                      */
    MPI_Barrier(comm);
    MPI_Gather(transit.ds.out->o, nspec, MPI_DOUBLE,
               dummy,             1,     MPI_DOUBLE, root, comm);
    transitprint(1, verblevel, "TRAN FLAG 77: The number of layers is: %li\n",
                                transit.ds.at->rads.n);

    /* Free arrays allocated inside the loop:                               */
    free(transit.save.ext);
    freemem_cia      (transit.ds.cia, &transit.pi);
    freemem_outputray(transit.ds.out, &transit.pi);

    /* Free no longer needed memory                                         */
    //freemem_idexrefrac(transit.ds.ir,        &transit.pi);
    //freemem_extinction(transit.ds.ex,        &transit.pi);
    //freemem_tau(transit.ds.tau,              &transit.pi);

    t0 = timecheck(verblevel, itr, 14, "THE END", tv, t0);
    transitprint(1, verblevel, "----------------------------\n");
  }
  //freemem_isotopes(     transit.ds.iso, &transit.pi);
  //freemem_molecules(    transit.ds.mol, &transit.pi);
  //freemem_atmosphere(   transit.ds.at,  &transit.pi);
  //freemem_lineinfotrans(transit.ds.li,  &transit.pi);
  //freemem_transit(&transit);

  transitprint(1, verblevel, "TRAN FLAG 98: End loop.\n");
  MPI_Barrier(comm);
  MPI_Comm_disconnect(&comm);
  transitprint(1, verblevel, "TRAN FLAG 98: worker off 2\n");
  MPI_Finalize();
}

/*
  Adapted to reading NeuroSAGE filres from LabVIEW. This is the C++ version of
  earlier readneurosage.c.
  Created by Cengiz Gunay 2005/08/03 <cgunay@emory.edu>

    o Adapted from Alfonso Delgado-Reyes's readpcdx and readgenesis (07.03.2002).
    
    Any cc compiler should work ([ ] = optional)
	
	To compile:
	  In Winblows (MS Visual C++ 6.x):
		mex [-v -O] -DWIN32 -output readneurosage readneurosage.cc
	  Anywhere else:
		mex [-v -O] -output readneurosage readneurosage.cc
*/

/* For LabVIEW files are in big-endian format */
#if !defined(__BIG_ENDIAN__) && !defined(__APPLE__)
#define __BIG_ENDIAN__
#endif

// Apples can read these files without conversion
#ifdef __APPLE__
#undef __APPLE__
#endif

#define PRG_NAME "readneurosage"

// Uncomment the following if you want a lot of print outs
//#define DEBUG

#include <time.h>

#include "common_fileIO.h"
#include "common_mexutils.h"
#include "common_neurosage.cc"
#include "neurosage/Sequence.cc"
#include "neurosage/Trial.cc"
#include "neurosage/Loader.cc"

/** Checks and prints a message if the trial number chosen exceeds 
    the total number of trials in file. */
void checkTrialNumExceeds(int trial_num, int num_trials, char *filename) {
  char		tmpbuffer[BUFSIZE];  

  if (trial_num > num_trials) {
    sprintf(tmpbuffer, "\n" PRG_NAME 
	    ": Selected trial %d exceeds total number of trials %d in file %s.\n", 
	    trial_num, num_trials, filename);
    fprintf(stderr, tmpbuffer);
    assert(strlen(tmpbuffer) <= BUFSIZE);
    mexErrMsgTxt(tmpbuffer);
  }
}

void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[]) {
  FILE		*fp;
  char		*filename;
  char		tmpbuffer[BUFSIZE];  

  int		 buflen, format_num;
  mxArray	*data, *config_str;
  const mxArray *trial_array, *channel_array;
  double	*start_of_pr;
  struct header_2D_array<double> *raw = NULL;
  struct header_NS_acq_prof_v0_2 *config;
  struct header_contents *contents;
  int dims[] = { 1 };
#define STR_FIELD_NUM_TRIALS "NumTrials"
#define STR_FIELD_MAX_RATE "MaxRate"
#define STR_NUM_FIELDS 2
  const char *str_field_names[STR_NUM_FIELDS] = 
    {STR_FIELD_NUM_TRIALS, STR_FIELD_MAX_RATE};

  if (nrhs != 3) {
    mexErrMsgTxt("\nUsage: data = " PRG_NAME 
		 "('filename', trial_array, channel_array)");
  }
  else if (nlhs < 1) {
    mexErrMsgTxt("\n" PRG_NAME " has one mandatory output argument");
  }
  
  if (mxIsChar(prhs[0]) != 1) {
    mexErrMsgTxt("\n" PRG_NAME ": first argument must be a string");
  } else {
    buflen = mxGetN(prhs[0])+1;
    
    filename = (char *) mxCalloc(buflen, sizeof(char));
    mxGetString(prhs[0], filename, buflen);
  }

  int trial_num;  
  if (!mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2])) {
    mxFree(filename);
    mexErrMsgTxt("\n" PRG_NAME ": arguments 2 & 3 must be noncomplex scalar arrays.");
  } else {
    trial_array = prhs[1];
    channel_array = prhs[2];
  }

  if ((fp = fopen(filename, "rb")) == NULL) {
    sprintf(tmpbuffer, "\n" PRG_NAME ": could not open file '%s'\n", filename);
    fprintf(stderr, tmpbuffer);
    assert(strlen(tmpbuffer) <= BUFSIZE);
    mexErrMsgTxt(tmpbuffer);
  }

#ifndef NOHEAD
  format_num = guess_format(fp);
#else
  file_length = get_file_length(fp);
  format_num = FILE_FORMAT_V0_2; 
#endif

  /* TODO: do the scanning here, to be able to check for excessive
     trial_num values */

  debugPrintf(PRG_NAME ": format #%d\n", format_num);

  switch (format_num) {
  case FILE_FORMAT_RAWACQ_BUGGY:
    contents = get_contents_rawacq(fp, TRUE);
    checkTrialNumExceeds(trial_num, contents->num_trials, filename); 
    raw = get_data_rawacq(fp, trial_num, TRUE, &config);
    break;
  case FILE_FORMAT_RAWACQ:
    contents = get_contents_rawacq(fp, FALSE);
    checkTrialNumExceeds(trial_num, contents->num_trials, filename); 
    raw = get_data_rawacq(fp, trial_num, FALSE, &config);
    break;
  case FILE_FORMAT_ACQ_PROF:
    contents = get_contents_acq_prof(fp);
    checkTrialNumExceeds(trial_num, contents->num_trials, filename); 
    raw = get_data_acq_prof(fp, trial_num, &config);
    break;    
  case FILE_FORMAT_V0_2:
    plhs[0] = get_data_v0_2(fp, trial_array, channel_array);
    break;    
  case FILE_FORMAT_UNKNOWN:
  default:
    mexErrMsgTxt("File format cannot be identified. Bailing out.\n");
  }

  if (raw != NULL) {
    debugPrintf(PRG_NAME ": found %d trials \n", contents->num_trials);

    /* LabVIEW does arrays the other way! */
    data = mxCreateDoubleMatrix(raw->x, raw->y, mxREAL);

    if (data == NULL) {
      mxFree(filename);
      mexErrMsgTxt("\n" PRG_NAME ": could not create mxArray (data)");
    }
	
    start_of_pr = (double *) mxGetPr(data);
    
    memcpy(start_of_pr, raw->data, raw->x * raw->y * sizeof(double));

    plhs[0] = data;
  } 

  /* If requested, return a structure with some info */
  /*if (nlhs > 1) {
    config_str = mxCreateStructArray(1, dims, STR_NUM_FIELDS,  str_field_names);
    mxSetField(config_str, 0, STR_FIELD_NUM_TRIALS,
	       mxCreateDoubleScalar(contents->num_trials));
    mxSetField(config_str, 0, STR_FIELD_MAX_RATE,
	       mxCreateDoubleScalar(config->max_rate));
    plhs[1] = config_str;
    }*/
  
  mxFree(filename);
  fclose(fp);
  return;
}

/*
  Adapted to reading NeuroSAGE filres from LabVIEW
  by Cengiz Gunay 2005/02/04 <cgunay@emory.edu>

    Adapted for general use (a routine...)
    from Dieter Jaeger's xicanal lx_genread.c
    Alfonso Delgado-Reyes 07.03.2002

    Cengiz Gunay <cgunay@emory.edu> 03.13.2004
    Fixed memory leak of not deallocating memory for the raw data buffer.
    
    o Adapted for MATLAB 5.x and 6.x under:
        - Linux x86/PPC, September 2002
        - MS Windows, September 2002
        - Mac OS 7-9, September 2002
        - Mac OS X 10.x, September 2002
    
    Any cc compiler should work ([ ] = optional)
	
	To compile:
	  In Winblows (MS Visual C++ 6.x):
		mex [-v -O] -DWIN32 -output readneurosage readneurosage.c
	  Anywhere else:
		mex [-v -O] -output readneurosage readneurosage.c
*/

/* For LabVIEW files are in big-endian format */
#ifndef __BIG_ENDIAN__
#define __BIG_ENDIAN__
#endif

#define PRG_NAME "readneurosage"

#include "common_fileIO.h"
#include "common_mexutils.h"
#include "common_neurosage.h"

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
  double	 trial_num;
  int		 buflen, format_num;
  mxArray	*data, *config_str;
  double	*start_of_pr;
  struct header_2D_array_dbl *raw = NULL;
  struct header_NS_acq_prof *config;
  struct header_contents *contents;
  int dims[] = { 1 };
#define STR_FIELD_NUM_TRIALS "NumTrials"
#define STR_FIELD_MAX_RATE "MaxRate"
#define STR_NUM_FIELDS 2
  const char *str_field_names[STR_NUM_FIELDS] = 
    {STR_FIELD_NUM_TRIALS, STR_FIELD_MAX_RATE};

  if (nrhs != 2) {
    mexErrMsgTxt("\nUsage: [data, config] = " PRG_NAME "('filename', trialnumber)");
  }
  else if (nlhs < 1) {
    mexErrMsgTxt("\n" PRG_NAME " has one mandatory output argument");
  }
  
  if (mxIsChar(prhs[0]) != 1) {
    mexErrMsgTxt("\n" PRG_NAME ": first argument must be a string");
  }
  else {
    buflen = mxGetN(prhs[0])+1;
    
    filename = (char *) mxCalloc(buflen, sizeof(char));
    mxGetString(prhs[0], filename, buflen);
  }
  
  if (!mxIsDouble(prhs[1])) {
    mxFree(filename);
    mexErrMsgTxt("\n" PRG_NAME ": argument 2 must be a noncomplex scalar.");
  }
  else {
    trial_num = mxGetScalar(prhs[1]);
  }

  if ((fp = fopen(filename, "rb")) == NULL) {
    sprintf(tmpbuffer, "\n" PRG_NAME ": could not open file '%s'\n", filename);
    fprintf(stderr, tmpbuffer);
    assert(strlen(tmpbuffer) <= BUFSIZE);
    mexErrMsgTxt(tmpbuffer);
  }

  format_num = guess_format(fp);

  /* TODO: do the scanning here, to be able to check for excessive
     trial_num values */

  mexPrintf(PRG_NAME ": format #%d\n", format_num);

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
    header_NS_file_v0_2 file_header = get_header_v0_2(fp);
    mexPrintf(PRG_NAME ": name=%s, tag=%s\n", file_header->exp_name, file_header->exp_tag);
    break;    
  case FILE_FORMAT_UNKNOWN:
  default:
    mexErrMsgTxt("File format cannot be identified. Bailing out.\n");
  }

  mexPrintf(PRG_NAME ": found %d trials \n", contents->num_trials);

  if (raw != NULL) {
    /* LabVIEW does arrays the other way! */
    data = mxCreateDoubleMatrix(raw->x, raw->y, mxREAL);

    if (data == NULL) {
      mxFree(filename);
      mexErrMsgTxt("\n" PRG_NAME ": could not create mxArray (data)");
    }
	
    start_of_pr = (double *) mxGetPr(data);
    
    memcpy(start_of_pr, raw->data, raw->x * raw->y * sizeof(double));

    plhs[0] = data;
  } else {
    mexErrMsgTxt("\n" PRG_NAME ": error... see output above");
  }

  /* If requested, return a structure with some info */
  if (nlhs > 1) {
    config_str = mxCreateStructArray(1, dims, STR_NUM_FIELDS,  str_field_names);
    mxSetField(config_str, 0, STR_FIELD_NUM_TRIALS,
	       mxCreateDoubleScalar(contents->num_trials));
    mxSetField(config_str, 0, STR_FIELD_MAX_RATE,
	       mxCreateDoubleScalar(config->max_rate));
    plhs[1] = config_str;
  }
  
  mxFree(filename);
  return;
}

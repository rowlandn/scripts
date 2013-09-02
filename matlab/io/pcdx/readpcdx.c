/*
 o  Library routine for retrieving data from truxpac data files
    Version 1.3
    June 8, 1993
    Truxton Fulton
	
 o Edited by D. Jaeger for Xicanal April, 94

 o Adapted by Alfonso Delgado-Reyes for MATLAB 5.x and 6.x under:
	- Linux x86/PPC, November 2000
	- MS Windows, November 2000
	- Mac OS 7-9, December 2000
	- Mac OS X 10.x, August 2002
	
    Any cc compiler should work ([ ] = optional)
	
	To compile:
	  In Winblows (MS Visual C++ 6.x):
		mex [-v -O] -DWIN32 -output readpcdx readpcdx.c
	  Anywhere else:
		mex [-v -O] -output readpcdx readpcdx.c

  o Carson Roberts, Cengiz Gunay <cgunay@emory.edu> 2005/01/24
    Fixed memory leak of not deallocating memory for the raw data buffer.

*/
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#if defined(WIN32)
#include <io.h>
#else
#include <unistd.h>
#endif
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <float.h>
#include <math.h>

#include <mex.h>
#include <mat.h>

#define HEADER_SIZE		         490
#define EXTENDED_HEADER_SIZE	3629
#define STRLEN			        8192

enum program_modes
{
  pulse_mode,
  train_mode,
  unknown_mode
};

struct siu_type
{
  int		active;
  int		transition;
  char		base_program_name[61];
  char		end_program_name[61];
  enum		program_modes program_mode;
  double	delay_from_trigger;
  double	pulse_width;
  int		num_pulses;
  double	amplitude;
};
  
struct header_type
{
  double	sampling_frequency;
  double	channel_amplification[16];
  char		channel_units[16][5];
  char		trial_batch_name[80];
  int		num_samples;
  int		starting_channel;
  int		ending_channel;
  int		num_channels;
  int		trial_batch_id;
  int		trial_id;
  int		absolute_index;
  char		date_and_time[20];
  double	AD_base;
  double	AD_range;
  int		offset_from_trigger;
  char		user_1[9];
  char		user_2[9];
  char		sequence_name[61];
  struct	siu_type siu[16];
};

struct filter_type {
  double lowcut;
  double highcut;
  double notchlow;
  double notchhigh;
};

struct analog_format {
  int		type;
  int		file_type;
  int		file_idx;
  int		cross_idx;
  int		trace_no;
  char		title[240];
  int		channel;
  int		color;
  int		plot_group;
  int		select;
  int		comp_sign;
  struct	filter_type filter;
  float		gain;
  int		invert;
  float		offset;
  float		factor;
  float		overlay_pos;
  float		overlay_val;
  float		xmax;
  float		xmin;
  long		no_samples;
  float		samp_frequency;
  int		filled;
  double	*fdata;
};

struct analog_format *raw = NULL;

void civilize_header(struct header_type *header,
			         char *trial_header, 
			         char *extended_trial_header);

int	get_truxdata(char *filename, int trace, int channel);

void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[])
{
  char		*filename;
  double	trace, channel;
  int		buflen;
  mxArray	*data;
  double	*start_of_pr;
  
  if (nrhs != 3) {
    mexErrMsgTxt("\nUsage: data = readpcdx('filename', trace, channel)");
  }
   else if (nlhs != 1) {
    mexErrMsgTxt("\nreadpcdx has one output argument");
  }
  
  if (mxIsChar(prhs[0]) != 1) {
    mexErrMsgTxt("\nreadpcdx: first argument must be a string");
    return;
  }
  else {
    buflen = mxGetN(prhs[0])+5;
  
    filename = (char *) mxCalloc(buflen, sizeof(char));
    mxGetString(prhs[0], filename, buflen);
  
    if (!strstr(filename, "."))
      filename = strcat(filename, ".all");
  }
  
  if (!mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2])) {
    mexErrMsgTxt("\nreadpcdx: arguments 2 & 3 must be noncomplex scalars");
    mxFree(filename);
    return;
  }
  else {
    trace = mxGetScalar(prhs[1]);
    channel = mxGetScalar(prhs[2]);
  }
  
  if (get_truxdata(filename, (int) trace, (int) channel) >= 0) {
    data = mxCreateDoubleMatrix(raw->no_samples, 1, mxREAL);
    
    if (data == NULL)
      mexErrMsgTxt("\nreadpcdx: could not create mxArray (data)");
    
    start_of_pr = (double *) mxGetPr(data);
    
    memcpy(start_of_pr, raw->fdata, raw->no_samples*1*sizeof(double));
    plhs[0] = data;

    /* Don't forget to free memory for the raw data buffer! 
       2005/01/24 Carson Roberts, Cengiz Gunay, <cgunay@emory.edu> */
    free(raw->fdata);
  }
  else {
    mexErrMsgTxt("\nreadpcdx: error... see output above");
  }
  
  mxFree(filename);
  if (raw != NULL) free(raw);
  return;
}

int get_truxdata(char *filename, int trace, int channel)
{
  struct	header_type *header = NULL;
  
  int		ch, file_descriptor;
  char 		*ad_buffer, temps[128], 
  		    trial_header[HEADER_SIZE],
  		    extended_trial_header[EXTENDED_HEADER_SIZE];
  long		i, j, datasize, index, curr_point, fst_pt, lst_pt;
  float		prexoffs, xmax_read;
  short		v, stemp;
  
  if ((raw = (struct analog_format *)
       malloc(sizeof(struct analog_format))) == NULL) {
	fprintf(stderr, "\nreadpcdx: could not malloc data structure\n");
	return -1;
  }
  
  raw->type     = 1;
  raw->factor   = 1;
  raw->gain     = 10.0;
  raw->offset   = 0.0;
  raw->xmin     = 0.0;
  raw->xmax     = 100000.0;
  raw->fdata    = NULL;
  raw->trace_no = trace;
  raw->channel  = channel;
  raw->invert   = 0;
  
#if defined(WIN32)
	file_descriptor = open(filename, O_RDONLY | O_BINARY);
#else
	file_descriptor = open(filename, O_RDONLY);
#endif
  if (file_descriptor < 0) {
	fprintf(stderr, "\nreadpcdx: can not open file %s\n", filename);
	return -1;
  } else {
	index = 0;
	while (read(file_descriptor, 
		   trial_header, HEADER_SIZE) == HEADER_SIZE) {
			
		i = trial_header[325];
		strncpy(temps, &trial_header[326], i);
		temps[i] = 0;
			
		sscanf(temps, "%ld", &datasize);
			
		if (++index == raw->trace_no) {	
			lseek(file_descriptor, 
			      (long) (datasize-EXTENDED_HEADER_SIZE), SEEK_CUR);
			
			if (read(file_descriptor, 
				 extended_trial_header,
				 EXTENDED_HEADER_SIZE) == EXTENDED_HEADER_SIZE) {
					
				if ((header = (struct header_type *) 
				     malloc(sizeof(struct header_type))) == NULL) {
				   fprintf(stderr, "\nreadpcdx: could not malloc header\n");
				   return -1;
				}

				civilize_header(header, trial_header, extended_trial_header);
					
				lseek(file_descriptor, -datasize, SEEK_CUR);
			} else {
			    fprintf(stderr, "\nreadpcdx: could not read extended header\n");
			    break; /* EOF */
			}
				
			raw->samp_frequency = (float) header->sampling_frequency;

			prexoffs = header->offset_from_trigger/raw->samp_frequency;		
			xmax_read = prexoffs+header->num_samples/raw->samp_frequency;
			
			fprintf(stderr, 
			        "readpcdx: %s, trace %i, channel %i (%g points @ %g KHz)\n", 
                    filename, raw->trace_no, raw->channel, 
                    xmax_read, raw->samp_frequency);
			
			if (raw->xmin < prexoffs) raw->xmin = prexoffs;
			
			if (raw->xmin >= xmax_read) {
				fprintf(stderr, "readpcdx: No data available above xmin\n");
				close(file_descriptor);
				return -1;
			}
            
			if (raw->xmax > xmax_read) raw->xmax = xmax_read;
			
			fst_pt = (long) (raw->xmin*raw->samp_frequency);
			lst_pt = (long) (raw->xmax*raw->samp_frequency);
				
			raw->no_samples = (lst_pt)-(fst_pt);

			if (raw->channel < header->starting_channel || 
				raw->channel > header->ending_channel) {
				fprintf(stderr, 
				"\nreadpcdx: channel %d not available\n", raw->channel);
				close(file_descriptor);
				return -1;
			}
  						
			/* fill ad_buffer, convert to double */
			
			ad_buffer = (char *) malloc(datasize);
				
			if (ad_buffer == NULL) {
				fprintf(stderr, 
				"\nreadpcdx: cannot alloc %ld memory bytes\n", datasize);
				close(file_descriptor);
				return -1;
			}
				
			raw->fdata = (double *) malloc(raw->no_samples*sizeof(double));
			
			if (raw->fdata == NULL) {
				fprintf(stderr, "\nreadpcdx: cannot alloc data array\n");
				close(file_descriptor);
				return -1;
			}
				
			strcpy(raw->title, "Raw Trux data: ");
			strcat(raw->title, header->trial_batch_name);
			strcat(raw->title, "; Sequence: ");
			strcat(raw->title, header->sequence_name);
				
			raw->filled	= 1;
			raw->file_type	= 1;
				
			read(file_descriptor, ad_buffer, datasize);
			
			curr_point = 0;
			for (i = 0, j = 0; 
			     i < raw->no_samples && 
			     j < datasize; j += 2) {
				
				v = *((short *) (ad_buffer + j));

                #if defined(__APPLE__)      || \
                    defined(__BIG_ENDIAN__) || \
                    defined(WORDS_BIGENDIAN)
                    /* swap the bits around for big endian */
                    v = ((v & 0x00FF) << 8) + ((v & 0xFF00) >> 8);
                #endif
                                
				ch = 0x0F & v;
				
				if (ch == raw->channel) {
					curr_point++;
					if (curr_point > fst_pt) {
						stemp = 0xFFF & (v >> 4);
												
						raw->fdata[i] = ((((double) (stemp)/4096.0)* 
										header->AD_range)+ 
										header->AD_base)*1000.0;
						
						raw->fdata[i++] /= raw->gain;

					}
				}				
			}

			if (i < raw->no_samples) {
				fprintf(stderr, "\nreadpcdx: Mismatch in get_data: %d><%d\n", 
				        raw->no_samples, i);
				raw->no_samples = i;
			}

			if (raw->invert == 1)
				for (j = 0; j < raw->no_samples; j++) raw->fdata[j] *= -1;
				
			free(ad_buffer);
			free(header);
			close(file_descriptor);
			return 0;
		} else
			lseek(file_descriptor, datasize, SEEK_CUR);
	} /* while */
  }
  
  fprintf(stderr, 
          "\nreadpcdx: Trace number %i does not exist\n", raw->trace_no);
  close (file_descriptor);
  return -1;
}

void civilize_extended_header(struct header_type *header, 
				char *extended_trial_header)
{
  char	temps[STRLEN];
  int	i, j, p;

  p = 0;

  i = extended_trial_header[p];
  strncpy(temps, &extended_trial_header[p+1], i);
  temps[i] = 0;
  sprintf(header->sequence_name, "%s", temps);
  p += 61;

  for (j = 0; j < 16; j++) {
      i = extended_trial_header[p];
      strncpy(temps, &extended_trial_header[p+1], i);
      temps[i] = 0;
      header->siu[j].active = (temps[0] == 'T');
      p += 2;

      i = extended_trial_header[p];
      strncpy(temps, &extended_trial_header[p+1], i);
      temps[i] = 0;
      header->siu[j].transition = (temps[0] == 'T');
      p+=2;

      i = extended_trial_header[p];
      strncpy(temps, &extended_trial_header[p+1], i);
      temps[i] = 0;
      sprintf(header->siu[j].base_program_name, "%s", temps);
      p += 61;

      i = extended_trial_header[p];
      strncpy(temps, &extended_trial_header[p+1], i);
      temps[i] = 0;
      sprintf(header->siu[j].end_program_name, "%s", temps);
      p += 61;

      i = extended_trial_header[p];
      strncpy(temps, &extended_trial_header[p+1], i);
      temps[i] = 0;

      switch (temps[0]) {
		case 'P':
			header->siu[j].program_mode = pulse_mode;
			break;
		case 'T':
			header->siu[j].program_mode = train_mode;
			break;
		default:
			header->siu[j].program_mode = unknown_mode;
			break;
	  }
      p += 2;

      i = extended_trial_header[p];
      strncpy(temps, &extended_trial_header[p+1], i);
      temps[i] = 0;
      sscanf(temps, "%lf", &(header->siu[j].delay_from_trigger));
      p += 18;

      i = extended_trial_header[p];
      strncpy(temps, &extended_trial_header[p+1], i);
      temps[i] = 0;
      sscanf(temps, "%lf", &(header->siu[j].pulse_width));
      p += 18;

      i = extended_trial_header[p];
      strncpy(temps, &extended_trial_header[p+1], i);
      temps[i] = 0;
      sscanf(temps, "%d", &(header->siu[j].num_pulses));
      p += 9;

      i = extended_trial_header[p];
      strncpy(temps, &extended_trial_header[p+1], i);
      temps[i] = 0;
      sscanf(temps, "%lf", &(header->siu[j].amplitude));
      p += 18;
  }
}

void civilize_header(struct header_type *header,
			char *trial_header, char *extended_trial_header)
{
  char	temps[STRLEN];
  int	i, j, p;
  int	sc, ec;

  p = 0;

  i = trial_header[p];
  strncpy(temps, &trial_header[p+1], i);
  temps[i] = 0;
  sscanf(temps, "%lf", &(header->sampling_frequency));
  p += 11;

  for (j = 0; j < 16; j++) {
      i = trial_header[p];
      strncpy(temps, &trial_header[p+1], i);
      temps[i] = 0;
      sscanf(temps, "%lf", &(header->channel_amplification[j]));
      p += 9;
  }

  for (j = 0; j < 16; j++) {
      i = trial_header[p];
      strncpy(temps, &trial_header[p+1], i);
      temps[i] = 0;
      sprintf(header->channel_units[j], "%s", temps);
      p += 5;
  }

  i = trial_header[p];
  strncpy(temps, &trial_header[p+1], i);
  temps[i] = 0;
  sprintf(header->trial_batch_name, "%s", temps);
  p += 81;

  i = trial_header[p];
  strncpy(temps, &trial_header[p+1], i);
  temps[i] = 0;
  sscanf(temps, "%d", &(header->num_samples));
  p += 9;

  /* byte size of data block */

  p += 12;

  i = trial_header[p];
  strncpy(temps, &trial_header[p+1], i);
  temps[i] = 0;
  sscanf(temps, "%d", &(header->starting_channel));
  p += 3;

  i = trial_header[p];
  strncpy(temps, &trial_header[p+1], i);
  temps[i] = 0;
  sscanf(temps, "%d", &(header->ending_channel));
  p += 3;

  sc = header->starting_channel;
  ec = header->ending_channel;

  header->num_channels = 1+ec-sc;

  if (sc > ec) header->num_channels += 16;

  i = trial_header[p];
  strncpy(temps, &trial_header[p+1], i);
  temps[i] = 0;
  sscanf(temps, "%d", &(header->trial_batch_id));
  p += 9;

  i = trial_header[p];
  strncpy(temps, &trial_header[p+1], i);
  temps[i] = 0;
  sscanf(temps, "%d", &(header->trial_id));
  p += 9;

  i = trial_header[p];
  strncpy(temps, &trial_header[p+1], i);
  temps[i] = 0;
  sprintf(header->date_and_time, "%s", temps);
  p += 20;

  i = trial_header[p];
  strncpy(temps, &trial_header[p+1], i);
  temps[i] = 0;
  sscanf(temps, "%lf", &(header->AD_base));
  p += 18;

  i = trial_header[p];
  strncpy(temps, &trial_header[p+1], i);
  temps[i] = 0;
  sscanf(temps, "%lf", &(header->AD_range));
  p += 18;

  i = trial_header[p];
  strncpy(temps, &trial_header[p+1], i);
  temps[i] = 0;
  sscanf(temps, "%d", &(header->offset_from_trigger));
  p += 9;

  i = trial_header[p];
  strncpy(temps, &trial_header[p+1], i);
  temps[i] = 0;
  sprintf(header->user_1, "%s", temps);
  p += 9;

  i = trial_header[p];
  strncpy(temps, &trial_header[p+1], i);
  temps[i] = 0;
  sprintf(header->user_1, "%s", temps);
  p += 9;

  i = trial_header[p];
  strncpy(temps, &trial_header[p+1], i);
  temps[i] = 0;

  if (!strcmp(temps, "TF"))
	civilize_extended_header(header, extended_trial_header);
  else {
	sprintf(header->sequence_name, 
            "old truxpac format -- no sequence name");
	for (j = 0; j < 16; j++) header->siu[j].active = 0;
  }
  p += 3;

  /* byte size of data block --excluding extended header */

  p += 12;

  i = trial_header[p];
  strncpy(temps, &trial_header[p+1], i);
  temps[i] = 0;
  sscanf(temps, "%d", &(header->absolute_index));
  p += 9;
}

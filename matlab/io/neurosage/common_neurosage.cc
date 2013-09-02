/**
   Common NeuroSAGE/LabVIEW structures and functions header file.
   Supercedes common_neurosage.h, because this file makes use of C++ templates 
   for coping with changing data structure definitions.
   By Cengiz Gunay 2005/08/01 <cgunay@emory.edu>
*/

#ifndef PRG_NAME
#define PRG_NAME "common_neurosage_cc"
#endif

/* Because LabVIEW files are in big-endian format */
#ifndef __BIG_ENDIAN__
#define __BIG_ENDIAN__
#endif

#include "common_fileIO.h"
#include "common_mexutils.h"
#include "common_neurosage.h"

//unsigned int file_length;

template <typename scalar_type>
header_2D_array<scalar_type> *read_NS_2D_matrix(FILE *fp, char is_skip) {
  header_2D_array<scalar_type> *raw;
  char	tmpbuffer[BUFSIZE];
  int trial_samples, trial_bytes, read_length;
  double val;
  long long int llval; 
    
  raw = (header_2D_array<scalar_type> *) mex_malloc(sizeof(header_2D_array<scalar_type>));

  /* Read two int32 values for dimensions of 2D array. */
  raw->y = read_int32(fp);
  raw->x = read_int32(fp);

  /* Length of data array in bytes */
  trial_samples = raw->x * raw->y;
  trial_bytes = trial_samples * sizeof(*raw->data);

  /*debugPrintf(PRG_NAME ": reading %d x %d matrix of %d-byte scalars.\n",
    raw->y, raw->x, sizeof(*raw->data));*/

  assert(trial_samples < file_length);

  if (trial_samples == 0) {
    debugPrintf(PRG_NAME ": No data!\n");
    return raw;
  }

  if (is_skip) {		/* Skip the data */
    /*debugPrintf(PRG_NAME ": skipping %d x %d matrix of %d-byte scalars.\n",
      raw->y, raw->x, sizeof(*raw->data));*/
    /* Just fast forward to there */
    fseek(fp, trial_bytes, SEEK_CUR);
  } else {			/* Read the data */
    /*debugPrintf(PRG_NAME ": reading in %d x %d matrix of %d-byte scalars.\n",
      raw->y, raw->x, sizeof(*raw->data));*/

    /* Allocate memory for data */
    if (trial_bytes > 0)
      raw->data = (scalar_type *) mex_malloc(trial_bytes);

    /* Read the first trial (for now) */
    if ((read_length = fread(raw->data, sizeof(*raw->data), trial_samples, fp)) 
	!= trial_samples) {
      sprintf(tmpbuffer, 
	      PRG_NAME ": trial data truncated in file. Read %d samples, expected %d.\n", 
	      read_length, trial_samples);
      assert(strlen(tmpbuffer) <= BUFSIZE);
      fprintf(stderr, tmpbuffer);
      mexErrMsgTxt(tmpbuffer);
    }

    int i;
    /* Convert to little-endian for Matlab */
#if defined(__APPLE__)      || \
    defined(__BIG_ENDIAN__) || \
    defined(WORDS_BIGENDIAN)
    for (i = 0; i < trial_samples; i++) 
      bswap_ip(&raw->data[i]);
#endif
  }
  return raw;
}

template <typename scalar_type>
header_1D_array<scalar_type> *read_NS_1D_array(FILE *fp, char is_skip) {
  header_1D_array<scalar_type> *raw;
  char	tmpbuffer[BUFSIZE];
  int trial_samples, trial_bytes, read_length;
    
  raw = (header_1D_array<scalar_type> *) mex_malloc(sizeof(header_1D_array<scalar_type>));

  /* Read one int32 value for dimension of the array. */
  raw->x = read_int32(fp);

  /* Length of data array in bytes */
  trial_bytes = raw->x * sizeof(*raw->data);

  /*debugPrintf(PRG_NAME ": array of %d with %d-byte scalars.\n",
    raw->x, sizeof(*raw->data));*/

  assert(raw->x >= 0 && raw->x < file_length);

  if (is_skip) {		/* Skip the data */
    /*debugPrintf(PRG_NAME ": skipping array of %d with %d-byte scalars.\n",
      raw->x, sizeof(*raw->data));*/
    /* Just fast forward to there */
    fseek(fp, trial_bytes, SEEK_CUR);
  } else {			/* Read the data */
    /*debugPrintf(PRG_NAME ": reading in array of %d with %d-byte scalars.\n",
      raw->x, sizeof(*raw->data));*/

    /* Allocate memory for data */
    if (trial_bytes > 0)
      raw->data = (scalar_type *) mex_malloc(trial_bytes);

    /* Read the data */
    if ((read_length = fread(raw->data, sizeof(*raw->data), raw->x, fp)) != raw->x) {
      sprintf(tmpbuffer, 
	      PRG_NAME ": trial data truncated in file. Read %d samples, expected %d.\n", 
	      read_length, raw->x);
      assert(strlen(tmpbuffer) <= BUFSIZE);
      fprintf(stderr, tmpbuffer);
      mexErrMsgTxt(tmpbuffer);
    }

    /* Convert to little-endian for Matlab */
#if defined(__APPLE__)      || \
    defined(__BIG_ENDIAN__) || \
    defined(WORDS_BIGENDIAN)
    for (int i = 0; i < raw->x; i++) 
      bswap_ip(&raw->data[i]);
#endif
  }
  return raw;
}

// TODO: template makes concept of Base class useless
template <class T>
unsigned read_NS_object_array(FILE *fp, T *&objects) {
  /* Read one int32 value for dimension of the array. */
  unsigned num_items = read_int32(fp);
  assert(num_items >= 0);	// zero-length (empty) arrays possible

  /*debugPrintf(PRG_NAME ": reading in array of %d objects of size %d. ", 
    num_items, sizeof(T));*/

  //objects = (T *) mex_malloc(num_items * sizeof(T));
  objects = new T[num_items];

  /* Read the data */
  for (int i = 0; i < num_items; i++) 
    objects[i].read(fp);

  //debugPrintf("Done.\n");
  return num_items;
}

/* Dummy driver for instantiating all desired templates. */
/* TODO: remove this once they are tested and working */
void dummy() {
  assert(FALSE);
  read_NS_2D_matrix<short>(NULL, TRUE);
  read_NS_2D_matrix<int>(NULL, TRUE);
  read_NS_2D_matrix<long long int>(NULL, TRUE);
  read_NS_2D_matrix<float>(NULL, TRUE);
  read_NS_2D_matrix<double>(NULL, TRUE);

  read_NS_1D_array<short>(NULL, TRUE);
  read_NS_1D_array<int>(NULL, TRUE);
  read_NS_1D_array<long long int>(NULL, TRUE);
  read_NS_1D_array<float>(NULL, TRUE);
  read_NS_1D_array<double>(NULL, TRUE);
}

struct dict_array read_NS_dict_array(FILE *fp) {
  struct dict_array dict;

  dict.num_items = read_int32(fp);

  /* Allocate memory for channel array */
  // TODO: make classes with auto-destructors. until then, use mex_malloc to insure deallocation.
  // dict.items = new dict_item[dict.num_items];

  if (dict.num_items > 0) 
    dict.items = (struct dict_item *) 
      mex_malloc(dict.num_items * sizeof(struct dict_item));

  for (int i = 0; i < dict.num_items; i++) {
    dict.items[i].name = read_string(fp);
    dict.items[i].value = read_string(fp);
  }

  //debugPrintf("\n" PRG_NAME ": %d dictionary items found.\n", dict.num_items);

  return dict;
}

void scale_16bit2double(unsigned short *data_u16, unsigned length, 
			double range_low, double range_high, double *data_dbl) {
  int j;
  double range = range_high - range_low;

  if (range_low >= 0) {	/* If non-symmetric ranges => unsigned values */
    for (j = 0; j < length; j++) {
      data_dbl[j] = 
	((double)data_u16[j]) * range / 0xFFFFU +
	range_low;
    }
  } else {			/* If symmetric ranges => signed values  */
    for (j = 0; j < length; j++) {
      data_dbl[j] = 
	((double)((signed short)data_u16[j])) * range_high / 0x8000U;
    }
  }
}

/** Scales each row of the 16-bit matrix according to the 
    range information in each channel configuration */
struct header_2D_array<double>
*scale_16bit2double(header_NS_acq_prof_v0_2 *acq_chans, 
		    struct header_2D_array<short unsigned int> *data_u16) {
  int i, j, trial_samples, trial_bytes, row_offset, offset;
  double range_high, range_low;
  struct header_2D_array<double> *data_dbl;

  /* Length of data array in bytes */
  trial_samples = data_u16->x * data_u16->y;

  data_dbl = (struct header_2D_array<double>*) mex_malloc(sizeof(struct header_2D_array<double>));
  data_dbl->x = data_u16->x;
  data_dbl->y = data_u16->y;
  data_dbl->data = (double*) mex_malloc(trial_samples * sizeof(double));

  assert(acq_chans->num_channels == data_u16->y);

  for (i = 0; i < acq_chans->num_channels; i++) {
    range_low = (double) acq_chans->channels[i].range_low;
    range_high = (double) acq_chans->channels[i].range_high;

    /* Scale row to doubles */
    row_offset = i * data_u16->x;
    scale_16bit2double(&data_u16->data[row_offset], data_u16->x, 
		       range_low, range_high, &data_dbl->data[row_offset]);
  }

  return data_dbl;
}

struct header_contents
*get_contents_rawacq(FILE *fp, char is_buggy) {
  int trial_num;
  struct header_contents *contents;
  struct header_NS_acq_prof_v0_2 *config;
  struct header_2D_array<unsigned short> *a_chan_data;

  contents = (struct header_contents*) mex_malloc(sizeof(struct header_contents));

  /* Get the file size */
  fseek(fp, 0L, SEEK_END);
  file_length = ftell(fp);

  /* Initial rewind. */
  fseek(fp, 0L, SEEK_SET);

  debugPrintf("File length: %d\n", file_length);

  for (trial_num = 0; ftell(fp) < file_length; trial_num++) {
    config = read_NS_acq_array(fp, is_buggy);
    a_chan_data = read_NS_2D_matrix<unsigned short>(fp, TRUE);
    mxFree(config->channels);
    mxFree(config);
    mxFree(a_chan_data);
  }
  assert(trial_num > 0);

  contents->file_length = file_length;
  contents->num_trials = trial_num;

  return contents;
}

/** Read acq channel config array and then choose channel data */
struct header_2D_array<double> 
*get_data_rawacq(FILE *fp, int trial_num, char is_buggy, 
		  struct header_NS_acq_prof_v0_2 **config) {
  struct header_2D_array<double> *data;
  struct header_2D_array<short unsigned int> *a_chan_data;
  int trial_samples, trial_bytes, num_trials, read_length, skip_num;
    
  /* Get the file size */
  fseek(fp, 0L, SEEK_END);
  file_length = ftell(fp);

  /* Initial rewind. */
  fseek(fp, 0L, SEEK_SET);

  /* Skip unwanted trials */
  for (skip_num = 0; skip_num < trial_num - 1; skip_num++) {
    *config = read_NS_acq_array(fp, is_buggy);
    a_chan_data = read_NS_2D_matrix<unsigned short>(fp, TRUE);
    mxFree((*config)->channels);
    mxFree(*config);
    mxFree(a_chan_data);
  }

  *config = read_NS_acq_array(fp, is_buggy);

  /* Read the following 16-bit data matrix and convert it to a double matrix */
  data = scale_16bit2double(*config, read_NS_2D_matrix<unsigned short>(fp, FALSE));

  debugPrintf(PRG_NAME ": raw acquisition channel format, trial %d: \n" 
	    "%d channels, %d samples\n",
    trial_num, (*config)->num_channels, data->x);
    
  return data;
}

/** Read contents of file with acquisition profile format */
struct header_contents
*get_contents_acq_prof(FILE *fp) {
  int trial_num;
  struct header_contents *contents;
  struct header_NS_acq_prof_v0_2 *config;
  struct header_2D_array<unsigned short> *a_chan_data;

  contents = (struct header_contents*) mex_malloc(sizeof(struct header_contents));

  /* Get the file size */
  file_length = get_file_length(fp);

  /* Initial rewind. */
  fseek(fp, 0L, SEEK_SET);

  debugPrintf("File length: %d\n", file_length);

  for (trial_num = 0; ftell(fp) < file_length; trial_num++) {
    config = read_NS_acq_prof(fp);
    a_chan_data = read_NS_2D_matrix<unsigned short>(fp, TRUE);
    mxFree(config->channels);
    mxFree(config);
    mxFree(a_chan_data);
  }
  assert(trial_num > 0);

  contents->file_length = file_length;
  contents->num_trials = trial_num;

  return contents;
}

/** Read acq profile followed by trial data */
struct header_2D_array<double> 
*get_data_acq_prof(FILE *fp, int trial_num, 
		   struct header_NS_acq_prof_v0_2 **config) {
  struct header_2D_array<double> *data;
  struct header_2D_array<short unsigned int> *a_chan_data;
  int skip_num;
    
  /* Initial rewind. */
  fseek(fp, 0L, SEEK_SET);

  /* Skip unwanted trials */
  for (skip_num = 0; skip_num < trial_num - 1; skip_num++) {
    *config = read_NS_acq_prof(fp);
    a_chan_data = read_NS_2D_matrix<unsigned short>(fp, TRUE);
    mxFree((*config)->channels);
    mxFree(*config);
    mxFree(a_chan_data);
  }

  *config = read_NS_acq_prof(fp);

  /* Read the following 16-bit data matrix and convert it to a double matrix */
  data = scale_16bit2double(*config, read_NS_2D_matrix<unsigned short>(fp, FALSE));

  debugPrintf(PRG_NAME ": raw acquisition channel format, trial %d: \n" 
	    "%d channels, %d samples\n",
	    trial_num, (*config)->num_channels, data->x);
  /* TODO: return sample rates of each channel, and devices, etc. in a mlab structure */
    
  return data;
}

/** Guesses raw file format and returns */
enum file_formats guess_format(FILE *fp) {
  char	tmpbuffer[BUFSIZE], *profile_name;  
  struct header_2D_array<double> **data;
  struct header_2D_array<short unsigned int> *a_chan_data;
  unsigned int num_chans, trial_samples, trial_bytes, 
    num_trials, i, j, string_length;
  unsigned short rate, chan_num;
    
  /* Get the file size */
  fseek(fp, 0L, SEEK_END);
  file_length = ftell(fp);

  /* Initial rewind. */
  fseek(fp, 0L, SEEK_SET);

  /* Check for NeuroSAGE signature */
  if (is_signed(fp)) {
    return FILE_FORMAT_V0_2;
  }  

  /* Check for raw acquisition info */
  if (is_rawacq(fp)) {
    return FILE_FORMAT_RAWACQ_BUGGY;
  }

  /* Otherwise see if the first integer is followed by the above */
  /* Initial rewind. */
  fseek(fp, 0L, SEEK_SET);
  num_chans = read_int32(fp);

  if (is_rawacq(fp)) {
    sprintf(tmpbuffer, "\n" PRG_NAME ": found %d acquisition channels\n", num_chans);
    debugPrintf(tmpbuffer);
    return FILE_FORMAT_RAWACQ;
  }

  /* Finally, try to verify that the acq. prof is there. */
  fseek(fp, 0L, SEEK_SET);
  string_length = read_int32(fp);

  if (string_length < 100) {
    fseek(fp, 0L, SEEK_SET);
    profile_name = read_string(fp);
    sprintf(tmpbuffer, "\n" PRG_NAME ": found acquisition profile %s\n", profile_name);
    debugPrintf(tmpbuffer);
    return FILE_FORMAT_ACQ_PROF;    
  }

  return FILE_FORMAT_UNKNOWN;
}


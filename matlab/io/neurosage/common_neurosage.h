/**
   Common NeuroSAGE/LabVIEW structures and functions header file.
   Created by Cengiz Gunay on 2005/05/20 <cgunay@emory.edu>.
*/

unsigned int file_length;

#define SIGN_LENGTH 4
const char NeuroSAGE_signature[] = "NeSA";

/** Different file formats */
enum file_formats {
  FILE_FORMAT_UNKNOWN,
  FILE_FORMAT_RAWACQ_BUGGY,	/* Chan_num==0, rate=40000, missing number of channels */
  FILE_FORMAT_RAWACQ,		/* Same as above, but has number of channels */
  FILE_FORMAT_ACQ_PROF,		/* Contains intact acquisition profile completely. 
				   rate=max_rate, but per channel rate info known. */
  FILE_FORMAT_V0_2		/* Signed and versioned file format. Contains
				 * new sequence and trial blocks. */
};

/** NeuroSAGE units enumeration */
enum NS_units {
  UNITS_NONE
};

/** NeuroSAGE acquisition devices enumeration */
enum NS_devices {
  DEVICE_NONE
};

/** Header for reading NeuroSAGE acquisition channel clusters. */
struct header_NS_acq_chan_v0_1 {
  unsigned short	channel;
  unsigned short	rate_divisor;
  float			range_low, range_high;
  float			ext_scale;
  unsigned short	acq_device; /* see enum NS_devices */
  unsigned short	units;	/* see enum NS_units */
  double 		conversion_factor; /* This is added later */
};

struct header_NS_acq_chan_v0_2 {
  unsigned short	channel;
  unsigned short	rate_divisor;
  float			range_low, range_high;
  float			ext_scale;
  char 			*label;
  unsigned short	acq_device; /* see enum NS_devices */
  unsigned short	units;	/* see enum NS_units */
  double 		conversion_factor;
};

struct header_NS_acq_chan_v0_2 *read_NS_acq_chan_v0_2(FILE *fp) {
  struct header_NS_acq_chan_v0_2 *raw;
    
  raw = (struct header_NS_acq_chan_v0_2 *) 
    mex_malloc(sizeof(struct header_NS_acq_chan_v0_2));

  read_ip(fp, raw->channel);
  read_ip(fp, raw->rate_divisor);
  read_ip(fp, raw->range_low);
  read_ip(fp, raw->range_high);
  read_ip(fp, raw->ext_scale);
  read_ip(fp, raw->acq_device);
  read_ip(fp, raw->units);
  read_ip(fp, raw->conversion_factor);

  debugPrintf(PRG_NAME ": chan=%hd, rate_div=%hd, range=(%.2f, %.2f)"
	    ", ext_scale=%.2f, devid=%hd, units=%hd, conv_fac=%.4lf \n", 
	    raw->channel, raw->rate_divisor, raw->range_low, raw->range_high,
	    raw->ext_scale, raw->acq_device, raw->units, raw->conversion_factor);
    
  return raw;
}

/** Check for oldest format where channel numbers are always zero and
    acquisition rate is 40000 */
char is_rawacq(FILE *fp) {
  unsigned short chan_num, rate;

  chan_num = read_int16(fp);
  rate = read_int16(fp);

  if (chan_num == 0 && rate == 40000U) 
    return TRUE;
  else 
    return FALSE;
}

/** Check if file contains NeuroSAGE file header with initial signature.
 */
char is_signed(FILE *fp) {
  char magic[SIGN_LENGTH];

  /* Read first bytes */
  fread(&magic, sizeof(char), SIGN_LENGTH, fp);

  /* Compare with NeuroSAGE magic signature */
  if (!strncmp(magic, NeuroSAGE_signature, SIGN_LENGTH))
    return TRUE;
  else 
    return FALSE;
} 

/** NeuroSAGE acquisition channel array cluster */
/*struct header_NS_acq_chan_array {
  int	num_channels;
  struct header_NS_acq_chan *channels;
  };*/

/** NeuroSAGE acquisition channel array cluster */
struct header_NS_acq_prof_v0_2 {
  char *name;
  double trial_duration;
  unsigned int max_rate;
  int	num_channels;
  struct header_NS_acq_chan_v0_2 *channels;
};

struct header_NS_acq_prof_v0_2 *read_NS_acq_array(FILE *fp, char is_buggy) {
  struct header_NS_acq_prof_v0_2 *raw;
  char	tmpbuffer[BUFSIZE];
  int i, curpos, seek_forward;
    
  raw = (struct header_NS_acq_prof_v0_2 *) 
    mex_malloc(sizeof(struct header_NS_acq_prof_v0_2));

  if (is_buggy == 0)
    /* Read one int32 value for the length of 1D array. */
    raw->num_channels = read_int32(fp);
  else {			/* Otherwise we'll need to figure it out */
    curpos = ftell(fp);
    debugPrintf(PRG_NAME ": Current position at %d bytes...\n", curpos);
    /* Count structures */
    for (raw->num_channels = 0; is_rawacq(fp); raw->num_channels++) { 
      seek_forward = sizeof(struct header_NS_acq_chan_v0_2) - 2 * sizeof(short);
      /*debugPrintf(PRG_NAME ": seeking forward %d bytes...\n", seek_forward);*/
      fseek(fp, seek_forward, SEEK_CUR);
    }
    assert(raw->num_channels > 0);
    /* Rewind back to where we started from */
    fseek(fp, curpos, SEEK_SET);
  }

  /* Allocate memory for channel array */
  raw->channels = (struct header_NS_acq_chan_v0_2 *) 
    mex_malloc(raw->num_channels * sizeof(struct header_NS_acq_chan_v0_2));

  debugPrintf("\n" PRG_NAME ": %d channels found.\n", raw->num_channels);

  for (i = 0; i < raw->num_channels; i++) {
    /* Note that in the raw_acq format, the rate_divisor is the actual rate */
    raw->channels[i] = *(read_NS_acq_chan_v0_2(fp));
  }

  return raw;
}

struct header_NS_acq_prof_v0_2 *read_NS_acq_prof(FILE *fp) {
  struct header_NS_acq_prof_v0_2 *raw;
  char	tmpbuffer[BUFSIZE], *name;
  int   i, curpos, seek_forward;
  long long int foo;
  double trial_duration;
  unsigned int max_rate;
    
  raw = (struct header_NS_acq_prof_v0_2 *) 
    mex_malloc(sizeof(struct header_NS_acq_prof_v0_2));

  name = read_string(fp);
  foo = read_int64(fp);
  max_rate = read_int32(fp);

  raw = read_NS_acq_array(fp, FALSE);

  raw->name = name;
  raw->trial_duration = *((double*) &foo);
  raw->max_rate = max_rate;

  debugPrintf("\n" PRG_NAME ": Acquisition profile %s of duration %lf s,"
	    " and max rate %d.\n", raw->name, raw->trial_duration, raw->max_rate);

  return raw;
}

struct header_contents {
  int file_length;
  int num_trials;
};

/** Dictionary item */
struct dict_item {
  char *name, *value;
};

/** Dictionary array */
struct dict_array {
  int num_items; 
  struct dict_item *items;
};

/* NS file header */
struct header_NS_file_v0_2 {
  char *version, *exp_name, *exp_tag;
  time_t exp_date;
  struct dict_array dict_items;
  char *comments;
};

/** Simple templated header for reading LabVIEW 2D arrays of any scalar type. */
template <typename scalar_type>
struct header_2D_array {
  int 		x, y;
  scalar_type	*data;
};

/** Simple templated header for reading LabVIEW 1D arrays of any scalar type. */
template <typename scalar_type>
struct header_1D_array {
  int 		x;
  scalar_type	*data;
};

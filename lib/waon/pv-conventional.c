/* PV - phase vocoder : pv-conventional.c
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: pv-conventional.c,v 1.10 2007/10/14 06:23:36 kichiki Exp $
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
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>

// FFTW library
#include <fftw3.h>
// half-complex format handling routines
#include "hc.h"
#include "fft.h" // windowing()

// libsndfile
#include <sndfile.h>
#include "snd.h"

// ao device
#include <ao/ao.h>
#include "ao-wrapper.h"

// samplerate
#include <samplerate.h>

#include "memory-check.h" // CHECK_MALLOC() macro

struct default_writer_args {
  ao_device *ao;
  SNDFILE *sfout;
  SF_INFO *sfout_info;
};

struct default_reader_args {
  SNDFILE *sf;
  SF_INFO *sfinfo;
};

struct default_close_args {
  ao_device *ao;
  SNDFILE *sfout;
  SF_INFO *sfout_info;
  const char *outfile;
};

typedef int (*writer_fn)(long hop, double *l_out, double *r_out, void *user);
typedef int (*reader_fn)(long len, double *l_out, double *r_out, void *user);
typedef int (*close_fn)(long len, double *l_out, double *r_out, void *user);

int default_writer(long hop, double *l_out_src, double *r_out_src, void *user)
{
  int status;
  struct default_writer_args *args = (struct default_writer_args *)user;
  ao_device *ao = args->ao;
  SNDFILE *sfout = args->sfout;
  SF_INFO *sfout_info = args->sfout_info;

  // output
  if (sfout == NULL)
  {
    status = ao_write (ao, l_out_src, r_out_src, hop);
    status /= 4; // 2 bytes for 2 channels
  }
  else
  {
    status = sndfile_write (sfout, *sfout_info,
                            l_out_src, r_out_src, hop);
  }
  return status;
}

int default_close(long len, double *left, double *right, void* user)
{
  struct default_close_args *args = (struct default_close_args *)user;
  if (args->outfile == NULL)
  {
    ao_close (args->ao);
  }
  else
  {
    // frames left in l_out[] and r_out[]
    /*status =*/ sndfile_write (args->sfout, *args->sfout_info, left, right, len);
    sf_write_sync (args->sfout);
    sf_close (args->sfout);
  }
}

int default_reader(long len, double *left, double *right, void* user)
{
  struct default_reader_args *args = (struct default_reader_args *)user;
  return sndfile_read (args->sf, *args->sfinfo, left, right, len);
}

/** general utility routines for pv **/

/*
 * INPUT
 *  hop_res :
 *  hop_syn :
 *  l_out [hop_syn] :
 *  r_out [hop_syn] :
 *  ao, sfout, sfout_info : properties for output
 */
int
pv_play_resample (long hop_res, long hop_syn,
		  double *l_out, double *r_out,
		  writer_fn callback, void* user_args)
{
  int status = 0;

  int i;


  // samplerate conversion
  float *fl_in  = NULL;
  float *fl_out = NULL;
  double *l_out_src = NULL;
  double *r_out_src = NULL;
  SRC_DATA srdata;

  if (hop_res != hop_syn)
    {
      fl_in  = (float *)malloc (sizeof (float) * 2 * hop_syn);
      fl_out = (float *)malloc (sizeof (float) * 2 * hop_res);
      CHECK_MALLOC (fl_in,  "pv_play_resample");
      CHECK_MALLOC (fl_out, "pv_play_resample");

      srdata.input_frames  = hop_syn;
      srdata.output_frames = hop_res;
      srdata.src_ratio = (double)hop_res / (double)hop_syn;
      srdata.data_in  = fl_in;
      srdata.data_out = fl_out;

      l_out_src = (double *)malloc (sizeof (double) * hop_res);
      r_out_src = (double *)malloc (sizeof (double) * hop_res);
      CHECK_MALLOC (l_out_src, "pv_play_resample");
      CHECK_MALLOC (r_out_src, "pv_play_resample");


      // samplerate conversion (time fixed)
      for (i = 0; i < hop_syn; i ++)
	{
	  fl_in [i*2 + 0] = (float)l_out [i];
	  fl_in [i*2 + 1] = (float)r_out [i];
	}
      status = src_simple (&srdata, SRC_SINC_BEST_QUALITY, 2);
      //status = src_simple (&srdata, SRC_SINC_FASTEST, 2);
      if (status != 0)
	{
	  fprintf (stderr, "fail to samplerate conversion\n");
	  exit (1);
	}
      for (i = 0; i < hop_res; i ++)
	{
	  l_out_src [i] = (double)fl_out [i*2 + 0];
	  r_out_src [i] = (double)fl_out [i*2 + 1];
	}

      status = callback(hop_res, l_out_src, r_out_src, user_args);

      free (fl_in);
      free (fl_out);

      free (l_out_src);
      free (r_out_src);
    }
  else
    {
      // no samplerate conversion (pitch fixed)
      // output
      status = callback(hop_syn, l_out, r_out, user_args);
    }

  return (status);
}


/* estimate the superposing weight for the window with hop
 */
double
get_scale_factor_for_window (int len, long hop_syn, int flag_window)
{
  double *x = NULL;
  x = (double *)malloc (sizeof (double) * len);
  CHECK_MALLOC (x, "get_scale_factor_for_window");

  int i;
  for (i = 0; i < len; i ++)
    {
      x [i] = 1.0;
    }
  windowing (len, x, flag_window, 1.0, x);

  double acc = 0.0;
  double acc_max;
  acc_max = 0.0;
  int j;
  for (j = 0; j < hop_syn; j++)
    {
      acc = 0.0;
      for (i = 0; i < len; i += hop_syn)
	{
	  acc += x [j + i];
	}
      if (acc_max < acc) acc_max = acc;
    }

  free (x);

  // extra safety
  acc *= 1.5;

  return (acc);
}

void pv_conventional_buffer (int len, double *left, double *right, double *time, double *freq,
                             double *t_out, double *f_out, double *amp, double *ph_in,
                             double *l_ph_out, double *r_ph_out, double *omega,
                             double *l_ph_in_old, double *r_ph_in_old, fftw_plan plan,
                             fftw_plan plan_inv, double rate, double pitch_shift,
                             long hop_syn, int flag_window, int *flag_ph,
                             double *l_out, double *r_out,
                             writer_fn callback, void *user_args)
{
  long hop_ana;
  long hop_res;
  hop_res = (long)((double)hop_syn * pow (2.0, - pitch_shift / 12.0));
  hop_ana = (long)((double)hop_res * rate);

  double twopi = 2.0 * M_PI;

  int i;
  int k;

  double window_scale;
  window_scale = get_scale_factor_for_window (len, hop_syn, flag_window);


  // left channel
  apply_FFT (len, left, flag_window, plan, time, freq,
             1.0,
             amp, ph_in);
  if (flag_ph == 0)
  {
    // initialize phase
    for (k = 0; k < (len/2)+1; k ++)
    {
      l_ph_out [k] = ph_in [k] * (double)hop_syn / (double)hop_ana;
      //l_ph_out [k] = ph_in [k];

      // backup for the next step
      l_ph_in_old [k] = ph_in [k];
    }
    //*flag_ph = 1; // right channel is in the following!
  }
  else
  {
    // only for imag components who have phase
    for (k = 1; k < ((len+1)/2); k ++)
    {
      double dphi;
      dphi = ph_in [k] - l_ph_in_old [k]
          - omega [k] * (double)hop_ana;
      for (; dphi >= M_PI; dphi -= twopi);
      for (; dphi < -M_PI; dphi += twopi);

      l_ph_out [k] += dphi * (double)hop_syn / (double)hop_ana
          + omega [k] * (double)hop_syn;

      l_ph_in_old [k] = ph_in [k];
    }
  }
  polar_to_HC (len, amp, l_ph_out, 0, f_out);
  fftw_execute (plan_inv);
  // scale by len and windowing
  windowing (len, t_out, flag_window, (double)len * window_scale, t_out);
  // superimpose
  for (i = 0; i < len; i ++)
  {
    l_out [hop_syn + i] += t_out [i];
  }

  // right channel
  apply_FFT (len, right, flag_window, plan, time, freq,
             1.0,
             amp, ph_in);
  if (flag_ph == 0)
  {
    // initialize phase
    for (k = 0; k < (len/2)+1; k ++)
    {
      r_ph_out [k] = ph_in [k] * (double)hop_syn / (double)hop_ana;
      //r_ph_out [k] = ph_in [k];

      // backup for the next step
      r_ph_in_old [k] = ph_in [k];
    }
    *flag_ph = 1;
  }
  else
  {
    // only for imag components who have phase
    for (k = 1; k < ((len+1)/2); k ++)
    {
      double dphi;
      dphi = ph_in [k] - r_ph_in_old [k]
          - omega [k] * (double)hop_ana;
      for (; dphi >= M_PI; dphi -= twopi);
      for (; dphi < -M_PI; dphi += twopi);

      r_ph_out [k] += dphi * (double)hop_syn / (double)hop_ana
          + omega [k] * (double)hop_syn;

      r_ph_in_old [k] = ph_in [k];
    }
  }
  polar_to_HC (len, amp, r_ph_out, 0, f_out);
  fftw_execute (plan_inv);
  // scale by len and windowing
  //windowing (len, t_out, flag_window, (double)len, t_out);
  windowing (len, t_out, flag_window, (double)len * window_scale, t_out);
  // superimpose
  for (i = 0; i < len; i ++)
  {
    r_out [hop_syn + i] += t_out [i];
  }


  // output
  /*status =*/ pv_play_resample (hop_res, hop_syn, l_out, r_out, callback, user_args);


  /* shift acc_out by hop_syn */
  for (i = 0; i < len; i ++)
  {
    l_out [i] = l_out [i + hop_syn];
    r_out [i] = r_out [i + hop_syn];
  }
  for (i = len; i < len + hop_syn; i ++)
  {
    l_out [i] = 0.0;
    r_out [i] = 0.0;
  }


  /* for the next step */
  for (i = 0; i < (len - hop_ana); i ++)
  {
    left  [i] = left  [i + hop_ana];
    right [i] = right [i + hop_ana];
  }

}

void pv_conventional_prefilled (double *all_left, double *all_right, long all_len,
                                double rate, double pitch_shift,
                                long len, long hop_syn,
                                int flag_window, writer_fn callback, void *user_args,
                                close_fn close_callback, void *close_user_args)
{
  long hop_ana;
  long hop_res;
  hop_res = (long)((double)hop_syn * pow (2.0, - pitch_shift / 12.0));
  hop_ana = (long)((double)hop_res * rate);


  double twopi = 2.0 * M_PI;

  int i;
  int k;

  /* allocate buffers  */
  double * left = NULL;
  double * right = NULL;
  left  = (double *) malloc (sizeof (double) * len);
  right = (double *) malloc (sizeof (double) * len);
  CHECK_MALLOC (left,  "pv_conventional");
  CHECK_MALLOC (right, "pv_conventional");

  /* initialization plan for FFTW  */
  double *time = NULL;
  double *freq = NULL;
  time = (double *)fftw_malloc (len * sizeof(double));
  freq = (double *)fftw_malloc (len * sizeof(double));
  CHECK_MALLOC (time, "pv_conventional");
  CHECK_MALLOC (freq, "pv_conventional");
  fftw_plan plan;
  plan = fftw_plan_r2r_1d (len, time, freq, FFTW_R2HC, FFTW_ESTIMATE);

  double *t_out = NULL;
  double *f_out = NULL;
  f_out = (double *)fftw_malloc (len * sizeof(double));
  t_out = (double *)fftw_malloc (len * sizeof(double));
  CHECK_MALLOC (f_out, "pv_conventional");
  CHECK_MALLOC (t_out, "pv_conventional");
  fftw_plan plan_inv;
  plan_inv = fftw_plan_r2r_1d (len, f_out, t_out,
			       FFTW_HC2R, FFTW_ESTIMATE);

  double *amp = NULL;
  double *ph_in = NULL;
  amp   = (double *)malloc (((len/2)+1) * sizeof(double));
  ph_in = (double *)malloc (((len/2)+1) * sizeof(double));
  CHECK_MALLOC (amp,   "pv_conventional");
  CHECK_MALLOC (ph_in, "pv_conventional");

  double *l_ph_out = NULL;
  double *r_ph_out = NULL;
  l_ph_out    = (double *)malloc (((len/2)+1) * sizeof(double));
  r_ph_out    = (double *)malloc (((len/2)+1) * sizeof(double));
  CHECK_MALLOC (l_ph_out, "pv_conventional");
  CHECK_MALLOC (r_ph_out, "pv_conventional");


  for (i = 0; i < (len/2)+1; i ++)
  {
    ph_in [i]  = 0.0;
    l_ph_out [i] = 0.0;
    r_ph_out [i] = 0.0;
    /*l_ph_in_old [i] = 0.0;
      r_ph_in_old [i] = 0.0;*/
  }

  double *l_out = NULL;
  double *r_out = NULL;
  l_out = (double *) malloc ((hop_syn + len) * sizeof(double));
  r_out = (double *) malloc ((hop_syn + len) * sizeof(double));
  CHECK_MALLOC (l_out, "pv_conventional");
  CHECK_MALLOC (r_out, "pv_conventional");
  for (i = 0; i < (hop_syn + len); i ++)
  {
    l_out [i] = 0.0;
    r_out [i] = 0.0;
  }

  // expected frequency
  double * omega = NULL;
  omega = (double *) malloc (((len/2)+1) * sizeof(double));
  CHECK_MALLOC (omega, "pv_conventional");
  for (k = 0; k < (len/2)+1; k ++)
  {
    omega [k] = twopi * (double)k / (double)len;
  }

  double *l_ph_in_old = NULL;
  double *r_ph_in_old = NULL;
  l_ph_in_old = (double *)malloc (((len/2)+1) * sizeof(double));
  r_ph_in_old = (double *)malloc (((len/2)+1) * sizeof(double));
  CHECK_MALLOC (l_ph_in_old, "pv_conventional");
  CHECK_MALLOC (r_ph_in_old, "pv_conventional");

  for (i = 0; i < (len/2)+1; i ++)
  {
    l_ph_in_old [i] = 0.0;
    r_ph_in_old [i] = 0.0;
  }

  // first frame:
  memcpy(left, all_left, len * sizeof(double));
  memcpy(right, all_right, len * sizeof(double));

  int flag_ph = 0;
  long offset = 0;
  for (;;)
  {
    pv_conventional_buffer(len, left, right, time, freq, t_out,
                           f_out, amp, ph_in, l_ph_out, r_ph_out, omega, l_ph_in_old,
                           r_ph_in_old, plan, plan_inv, rate, pitch_shift, hop_syn,
                           flag_window, &flag_ph, l_out, r_out, callback, user_args);
   offset += len;
    if (offset >= all_len) {
      break;
    }
    memcpy(left, left + hop_ana, (len - hop_ana) * sizeof(double));
    memcpy(right, right + hop_ana, (len - hop_ana) * sizeof(double));
    memcpy(left + len - hop_ana, all_left + offset, hop_ana * sizeof(double));
    memcpy(right + len - hop_ana, all_right + offset, hop_ana * sizeof(double)); 
  }

  close_callback(len, l_out, r_out, close_user_args);

  free (l_ph_in_old);
  free (r_ph_in_old);

  free (l_out);
  free (r_out);

  free (time);
  free (freq);
  fftw_destroy_plan (plan);

  free (t_out);
  free (f_out);
  fftw_destroy_plan (plan_inv);

  free (amp);
  free (ph_in);

  free (l_ph_out);
  free (r_ph_out);

  free (omega);

  free (left);
  free (right);
}

/* standard phase vocoder
 * Ref: J.Laroche and M.Dolson (1999)
 */
void pv_conventional (const char *file, const char *outfile,
		      double rate, double pitch_shift,
		      long len, long hop_syn,
		      int flag_window)
{
  long hop_ana;
  long hop_res;
  hop_res = (long)((double)hop_syn * pow (2.0, - pitch_shift / 12.0));
  hop_ana = (long)((double)hop_res * rate);

  // open the input file
  long read_status;
  // libsndfile version
  SNDFILE *sf = NULL;
  SF_INFO sfinfo;
  memset (&sfinfo, 0, sizeof (sfinfo));
  sf = sf_open (file, SFM_READ, &sfinfo);
  if (sf == NULL)
    {
      fprintf (stderr, "fail to open %s\n", file);
      exit (1);
    }
  sndfile_print_info (&sfinfo);

  // prepare the output
  int status;
  ao_device *ao = NULL;
  SNDFILE *sfout = NULL;
  SF_INFO sfout_info;
  if (outfile == NULL)
    {
      ao = ao_init_16_stereo (sfinfo.samplerate, 1 /* verbose */);
    }
  else
    {
      sfout = sndfile_open_for_write (&sfout_info,
				      outfile,
				      sfinfo.samplerate,
				      sfinfo.channels);
      if (sfout == NULL)
	{
	  fprintf (stderr, "fail to open file %s\n", outfile);
	  exit (1);
	}
    }

  // read the first frame
  double *left = malloc(len * sizeof(double));
  double *right = malloc(len * sizeof(double));
  read_status = sndfile_read (sf, sfinfo, left, right, len);
  if (read_status != len)
    {
      exit (1);
    }

  long all_len = len;
  long offset = 0;
  for (;;) {
    /* read next segment */
    read_status = sndfile_read (sf, sfinfo,
                                left  + offset,
                                right + offset,
                                hop_ana);
    if (read_status != hop_ana)
    {
      // most likely, it is EOF.
      all_len -= hop_ana - read_status;
      break;
    }
    offset += len;
    all_len += len;
    left = realloc(left, all_len * sizeof(double));
    right = realloc(right, all_len * sizeof(double));
  }

  struct default_close_args close_args = { ao, sfout, &sfout_info, outfile };
  struct default_writer_args args = { ao, sfout, &sfout_info };
  pv_conventional_prefilled(left, right, all_len, rate, pitch_shift, len, hop_syn,
                            flag_window, default_writer, &args, default_close, &close_args);
  sf_close (sf);
}

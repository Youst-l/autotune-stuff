#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <fftw3.h>
#include <samplerate.h>
#include "hc.h"
#include "fft.h" // windowing()
#include "memory-check.h" // CHECK_MALLOC() macro

typedef int (*writer_fn)(long hop, double *l_out, double *r_out, void *user);
typedef int (*reader_fn)(long len, double *l_out, double *r_out, void *user);
typedef int (*close_fn)(long len, double *l_out, double *r_out, void *user);

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
		  double **left_out, double **right_out, long *out_len)
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

    long old_len = *out_len;
    *out_len += hop_res;
    *left_out = realloc(*left_out, *out_len * sizeof(double));
    memcpy(left_out + old_len, l_out_src, hop_res * sizeof(double));
    *right_out = realloc(*right_out, *out_len * sizeof(double));
    memcpy(right_out + old_len, r_out_src, hop_res * sizeof(double));

    free (fl_in);
    free (fl_out);

    free (l_out_src);
    free (r_out_src);
  }
  else
  {
    // no samplerate conversion (pitch fixed)
    // output
    long old_len = *out_len;
    *out_len += hop_res;
    *left_out = realloc(*left_out, *out_len * sizeof(double));
    memcpy(left_out + old_len, l_out, hop_syn * sizeof(double));
    *right_out = realloc(*right_out, *out_len * sizeof(double));
    memcpy(right_out + old_len, r_out, hop_syn * sizeof(double));
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
                             double *l_out, double *r_out, double **out_left,
                             double **out_right, long *out_len)
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
  /*status =*/ pv_play_resample (hop_res, hop_syn, l_out, r_out, out_left, out_right, out_len);


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
                                int flag_window, double **out_left, double **out_right,
                                long *out_len)
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
                           flag_window, &flag_ph, l_out, r_out, out_left, out_left, out_len);
    offset += len;
    if (offset >= all_len) {
      break;
    }
    memcpy(left, left + hop_ana, (len - hop_ana) * sizeof(double));
    memcpy(right, right + hop_ana, (len - hop_ana) * sizeof(double));
    memcpy(left + len - hop_ana, all_left + offset, hop_ana * sizeof(double));
    memcpy(right + len - hop_ana, all_right + offset, hop_ana * sizeof(double)); 
  }

  //close_callback(len, l_out, r_out, close_user_args);

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

int main() {
  double* left;
  double* right;
  long out_len;
  pv_conventional_prefilled (NULL, NULL, 0, 1.0, 0.0, 2048, 512, 3, &left, &right, &out_len);
}

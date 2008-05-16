// Define SIP wrapper for Mentor JWDB read api 

%ModuleCode
#include <ams_sim.h>
#include <malloc.h>

/* Wrapper for function that needs allocation of strings etc */
char ams_sim_get_simulation_parameter_info_wrapper   (long index, long iparam, char **name, char **val) {
  /* Allocate buffers */
  *name = malloc(LINE);
  *val = malloc(LINE);
  return ams_sim_get_simulation_parameter_info(index,iparam,*name,*val);
}  
long ams_sim_get_next_waveform_name_wrapper (long index, char ** dest, long from_where) {
  /* Allocate buffers */
  *dest = malloc(LINE);
  return ams_sim_get_next_waveform_name(index,*dest, from_where);
}   

%End

%CModule jwdb 0

  /* How to: initialize the interface */
  void ams_sim_init_file                       ();
  void ams_sim_set_connection_message          (char value);
  void ams_sim_set_data_blocksize              (long value);
  void ams_sim_set_multirun_mode               (char value);
  void ams_sim_set_jwdb_server                 (char *host, long port);

  /***********************************************************
  * Prototypes for the read API
  ***********************************************************/

  /* How to: open a file */
  long ams_sim_open_read_file                  (char *name);

  /* How to: get file informations */
  char ams_sim_is_jwdb_file                    (long index);
  char ams_sim_is_wdb_file                     (long index);
  char ams_sim_is_cou_file                     (long index);
  char ams_sim_is_compat_file                  (long index);
  
  /* How to: load a simulation */
  long ams_sim_load_simu                       (long index);
  long ams_sim_load_simu_by_index              (long index, long isimu);
  long ams_sim_load_simu_with_complex          (long index);
  long ams_sim_load_simu_with_complex_by_index (long index, long isimu);
  long ams_sim_set_simulation_xstart           (long index, double xstart);
  long ams_sim_set_simulation_xstop            (long index, double xstart);

  /* Ho to: get simulation information */
  char *ams_sim_get_x_unit_name                 (long index);
  char *ams_sim_get_creation_date               (long index);
  char *ams_sim_get_creation_time               (long index);
  long  ams_sim_get_simulation_parameter_number (long index);
char  ams_sim_get_simulation_parameter_info_wrapper   (long index, long iparam, char **name /Out,TransferBack/, char **val /Out,TransferBack/) /PyName=ams_sim_get_simulation_parameter_info/;

  /* How to: access waveforms */
  long ams_sim_get_next_waveform_index          (long index, long from_where);
  long ams_sim_get_next_waveform_name_wrapper   (long index, char ** dest /Out, TransferBack/, long from_where)  /PyName=ams_sim_get_next_waveform_name/;
  long ams_sim_get_index_to_extract_waveform    (long index, char *name);

  /* How to: get waveform information */
  int   ams_sim_is_analog_wave                  (char *WaveType /In/);
  char *ams_sim_get_wave_unit_name              (long index, long index_w);
  char  ams_sim_is_wave_complex                 (long index, long index_w);
  char  ams_sim_is_wave_analog                  (long index, long index_w);
  char  ams_sim_is_wave_digital                 (long index, long index_w);

  /* How to: set simulation boundaries for reading */
  long ams_sim_set_simulation_read_start   (long index, double xstart);
  long ams_sim_set_simulation_read_stop    (long index, double xstop);

  /* How to: read waves value */
  long ams_read_selected_waves              (long index, long *index_w, long nbw, double *x /Out/, double *y /Out/);
  long ams_read_selected_waves_complex      (long index, long *index_w, long nbw, double *x /Out/, double *y /Out/, double *yi /Out/);
  long ams_read_selected_digital_waves      (long index, long *index_w, long nbw, double *x, char **y /Out/);

  long ams_sim_read_analog_x_y_data         (long index, long index_w, double *x /Out/, double *y /Out/);
  long ams_sim_read_analog_x_y_data_complex (long index, long index_w, double *x /Out/, double *y /Out/, double *yi /Out/);
  long ams_sim_read_digital_x_y_data        (long index, long index_w, double *x /Out/, char **ydata /Out/);
  long ams_sim_read_mixed_x_y_data          (long index);

  /* How to: get waves value */
  long ams_sim_get_x_data                   (long index, double * x /Out/);
  long ams_sim_get_analog_x_y_data          (long index, long index_w, double *x /Out/, double *y /Out/, double *yi /Out/);
  long ams_sim_get_digital_x_y_data         (long index, long index_w, double *x /Out/, char **y /Out/);

  /* How to: close a file*/
  void ams_sim_close_read_file              (long index);

  /* Previous version functions */
  long ams_sim_read_analog_x_y_data_v1      (long index, long index_w, double *x /Out/, double *y /Out/);
  long ams_sim_read_digital_x_y_data_v1     (long index, long index_w, double *x /Out/, char **ydata /Out/);





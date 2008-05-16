import jwdb
#import result



#class JWDBResultSet(result.ResultSet):
#    """The JWDBResultSet class handles a JWDB waveform file used by Mentor Graphics software
#
#    """

jwdb.ams_sim_init_file()
filename = "./test/data/g3rxifvga1s_tb.wdb"
fileid = jwdb.ams_sim_open_read_file(filename)
if fileid < 0:
    raise IOError("Cannot open WDB file %s"%filename)

## Get time and date
date = jwdb.ams_sim_get_creation_date(fileid)
time = jwdb.ams_sim_get_creation_time(fileid)

## Load first simulation results
loadflag = jwdb.ams_sim_load_simu_with_complex(fileid)
while loadflag > 0:
    n = jwdb.ams_sim_get_simulation_parameter_number(fileid)

    ## Read the simulation parameters (temperature, ...)
    params = {}
    for i in range(1,n+1):
        tmp, name, val = jwdb.ams_sim_get_simulation_parameter_info(fileid, i)
        params[name] = val

    ## Get unit of X axis
    xunit = jwdb.ams_sim_get_x_unit_name(fileid)
    print "xunit", xunit

    ## Get first waveform name and identifier by putting 0 at last argument
    wid, wavename = jwdb.ams_sim_get_next_waveform_name(fileid, 0)
    while wid > 0:
        print wid, wavename

        ## Get unit of Y axis
        yunit = jwdb.ams_sim_get_wave_unit_name(fileid, wid)
        print "yunit", yunit

        ## Handle digital waveforms
        if not jwdb.ams_sim_is_analog_wave(yunit):
            pass
        ## Handle complex analog waveforms
        elif jwdb.ams_sim_is_wave_complex(fileid, wid):
            readflag, x, yr, yi = jwdb.ams_sim_read_analog_x_y_data_complex(fileid, wid)
            while readflag > 0:
                readflag, x, yr, yi = jwdb.ams_sim_read_analog_x_y_data_complex(fileid, wid)
        ## Handle real analog waveforms
        else:
            readflag, x, y = jwdb.ams_sim_read_analog_x_y_data(fileid, wid)
            while readflag > 0:
                readflag, x, y = jwdb.ams_sim_read_analog_x_y_data(fileid, wid)
                
        ## Get next waveform name and identifier 
        wid, wavename = jwdb.ams_sim_get_next_waveform_name(fileid, wid)


    ## Load next simulation results
    loadflag = jwdb.ams_sim_load_simu_with_complex(fileid)

jwdb.ams_sim_close_read_file(fileid)

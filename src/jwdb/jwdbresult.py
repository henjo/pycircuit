"""
High level JWDB waveform format reading module

By Henrik Johanson <henrik.johansson@csr.com>

"""

import jwdb
import pycircuit.result as result
import os
import numpy

## Init library
jwdb.ams_sim_init_file()

class _Simulation(object):
    def __init__(self, fileid, simid):
        self.waveforms = {}
        
        ## Read the simulation parameters (temperature, ...)
        self.params = {}
        n = jwdb.ams_sim_get_simulation_parameter_number(fileid)
        for i in range(1,n+1):
            tmp, name, val = jwdb.ams_sim_get_simulation_parameter_info(fileid, i)
            self.params[name] = val

        ## Get unit of X axis
        self.xunit = jwdb.ams_sim_get_x_unit_name(fileid)

        ## Get first waveform name and identifier by putting 0 at last argument
        wid, wavename = jwdb.ams_sim_get_next_waveform_name(fileid, 0)
        while wid > 0:
            print wid, wavename

            ## Get unit of Y axis
            self.yunit = jwdb.ams_sim_get_wave_unit_name(fileid, wid)

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

class JWDBResultSet(result.ResultSet):
    """The JWDBResultSet class handles a JWDB waveform file used by Mentor Graphics software

    Example

    >>> rs = JWDBResultSet("./test/data/g3rxifvga1s_tb.wdb")
    >>> rs
    
    """
    def __init__(self, filename):
        ## Open file
        if not os.path.isfile(filename):
            raise IOError("Cannot open WDB file %s"%filename)
        
        self.fileid = jwdb.ams_sim_open_read_file(filename)

        if self.fileid < 0:
            raise IOError("Cannot open WDB file %s"%filename)

        self.simulations = []
        self._load()

    def __del__(self):
        jwdb.ams_sim_close_read_file(self.fileid)

    def _load(self):
        ## Get time and date
        self.date = jwdb.ams_sim_get_creation_date(self.fileid)
        self.time = jwdb.ams_sim_get_creation_time(self.fileid)

        ## Load simulation results
        loadflag = jwdb.ams_sim_load_simu_with_complex(self.fileid)
        simnum = 0
        self.simulations = []
        while loadflag > 0:
            self.simulations.append(_Simulation(self.fileid))
            
            ## Load next simulation results
            loadflag = jwdb.ams_sim_load_simu_with_complex(self.fileid)


if __name__ == "__main__":
    import doctest
    doctest.testmod()

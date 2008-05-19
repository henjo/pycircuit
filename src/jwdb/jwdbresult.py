"""
High level JWDB waveform format reading module

By Henrik Johanson <henrik.johansson@csr.com>

"""

import jwdb
import pycircuit.result as result
import os
import numpy
from copy import copy
from itertools import islice, groupby
from sets import Set

analysis_types = {
    4 : 'AC',
    7 : 'DC',
    10 : 'FSST'
}

# Cartesian operator of a list
def cartesian(listList):
    if listList:
        result = []
        prod = cartesian(listList[:-1])
        for x in prod:
            for y in listList[-1]:
                result.append(x + (y,))
        return result
    return [()]

## Init library
jwdb.ams_sim_init_file()

def grouponhigh(seq, key, edge = 1):
    """grouponhigh returns sub sequences of seq limited by when they key function returns True
    
    the key function must return a boolean.

    >>> out = []
    >>> [list(subseq) for subseq in grouponhigh(range(5), lambda x: x%3 == 0)]
    [[0, 1, 2], [3, 4]]
    >>> [list(subseq) for subseq in grouponhigh([1,1,0], lambda x: x == 1, -1)]
    [[1], [1, 0]]
    >>> [list(subseq) for subseq in grouponhigh([1,1], lambda x: x == 1, -1)]
    [[1], [1]]
    
    """
    
    istart = 0
    for i, x in enumerate(seq):
        if key(x) and i > 0:
            yield islice(seq, istart, i)
            istart = i 
    yield islice(seq, istart, i+1)

    
class _Waveform(object):
    def __init__(self, fileid, simulation, id, name):
        self.fileid = fileid
        self.simulation = simulation
        self.id = id
        self.name = name

        ## Get unit of Y axis
        self.yunit = jwdb.ams_sim_get_wave_unit_name(fileid, self.id)

    def read(self):
        xlist = []
        ylist = []

        ## Load simulation
        jwdb.ams_sim_load_simu_with_complex_by_index(self.fileid, self.simulation.id)

        wname = ''
        wid = 0
        while wname != self.name and wid >= 0:
            wid, wname = jwdb.ams_sim_get_next_waveform_name(self.fileid, wid)
#        if wid < 0:
#            raise Exception('Could not find waveform %s'%self.name)

        ## Handle digital waveforms
        if not jwdb.ams_sim_is_analog_wave(self.yunit):
            pass
        ## Handle complex analog waveforms
        elif ord(jwdb.ams_sim_is_wave_complex(self.fileid, self.id)):
            readflag, x, yr, yi = jwdb.ams_sim_read_analog_x_y_data_complex(self.fileid, self.id)
            while readflag > 0:
                xlist.append(x)
                ylist.append(numpy.complex(yr, yi))
                readflag, x, yr, yi = jwdb.ams_sim_read_analog_x_y_data_complex(self.fileid, self.id)

        ## Handle real analog waveforms
        else:
            readflag, x, y = jwdb.ams_sim_read_analog_x_y_data(self.fileid, self.id)
            while readflag > 0:
                xlist.append(x)
                ylist.append(y)
                
                readflag, x, y = jwdb.ams_sim_read_analog_x_y_data(self.fileid, self.id)

        return xlist, ylist            

reserved_params = ('CKTYPE', 'NUMPTS', 'STEP', 'ANALYSIS', 'NOISEF', 'ZCHAR',
                   'Creation_Date', 'ICARLO', 'Creation_Time', 'Title')

class _Simulation(object):
    def __init__(self, fileid, simid):
        self.waveforms = {}
        self.id = simid

        ## Read the simulation parameters (temperature, ...)
        ## Try to find simulation sweep variables
        ## First comes three paramaters 
        self.params = {}
        self.variablenames = []
        self.variables = {}
        n = jwdb.ams_sim_get_simulation_parameter_number(fileid)
        for i in range(1,n+1):
            tmp, name, val = jwdb.ams_sim_get_simulation_parameter_info(fileid, i)

            if name not in reserved_params:
                self.variablenames.append(name)
                self.variables[name] = float(val)
            else:
                if name in ('Creation_Date', 'Creation_Time', 'Title'):
                    self.params[name] = val
                else:
                    self.params[name] = float(val)

        ## Get unit of X axis
        self.xunit = jwdb.ams_sim_get_x_unit_name(fileid)

        ## Get first waveform name and identifier by putting 0 at last argument
        wid, wavename = jwdb.ams_sim_get_next_waveform_name(fileid, 0)
        while wid > 0:
            waveform = _Waveform(fileid, self, wid, wavename)
            
            self.waveforms[wavename] = waveform

            wid, wavename = jwdb.ams_sim_get_next_waveform_name(fileid, wid)

    def getWaveformNames(self):
        """Get available waveform names"""
        return self.waveforms.keys()

    def getWaveform(self, name):
        """Get waveform by name"""
        return self.waveforms[name]

    def get_analysistype(self):
        """Find out analysis type"""

        if self.params.has_key('ANALYSIS'):
            analysis = int(self.params['ANALYSIS'])
            if analysis not in analysis_types:
                raise ValueError('Unknown analysis type: %d'%analysis)
            return analysis_types[analysis]
        else:
            if self.params.has_key('CKTYPE'):
                return 'Ext'
            else:
                return 'Noise'

    def __repr__(self):
        return self.__class__.__name__ + '(waveforms = ' + str(self.waveforms) + ')'

    analysis = property(get_analysistype)

class _Analysis(object):
    """An _Analysis object holds data from a simulation analysis.
    
    The analysis object contains one or more _Simulation objects and
    each one of them represents one iteration of a swept analysis,
    montecaro sample or alter statement.
    """
    

    def __init__(self, simulations):
        self.simulations = simulations

        variables = self.get_sweep_variables()

        ## Recover sweep vectors
        sweepsets = [Set([sim.variables[var] for sim in self.simulations]) \
                         for var in variables]
        self.sweepvalues = map(list, sweepsets)

        for s in self.sweepvalues: s.sort()
        

    def get_type(self):
        """Find out analysis type"""

        if len(self.simulations) < 1:
            return None
        return self.simulations[0].analysis

    def get_sweep_variables(self):
        """Return a list of sweep variables"""
        if len(self.simulations) < 1:
            return None
        return tuple(self.simulations[0].variablenames)

    def getWaveformNames(self):
        return self.simulations[0].getWaveformNames()
        
    def getWaveform(self, waveformname):
        """Get a result.Waveform object of a waveform"""

        variables = list(self.get_sweep_variables())
        sweepvalues = copy(self.sweepvalues)
        
        ## Create a dictionary between sweep values and simulations
        sweepmap = {}
        for sim in self.simulations:
            sweepmap[tuple([sim.variables[var] for var in variables])] = sim
            
        ## Load waveform data and find out if all x-vectors are the same for all sweeps
        xfirst = None
        ylist = []
        for varvalues in cartesian(sweepvalues):
            sim = sweepmap[varvalues]
            waveform = sim.getWaveform(waveformname)
            x, y = waveform.read()
            ylist.extend(y)

            if xfirst == None:
                xfirst = x
            else:
                if x != xfirst:
                    raise ValueError("All x-vectors must be the same %s != %s"%(str(xfirst), str(x)))

        variables.append(self.simulations[0].xunit)
        sweepvalues.append(xfirst)
        
        wavearray = numpy.array(ylist)
        wavearray = wavearray.reshape(map(len, sweepvalues))

        return result.Waveform(x = sweepvalues, 
                               y = wavearray,
                               xlabels = tuple(variables),
                               ylabel = waveformname,
                               xunits = (len(variables) - 1) * [''] + [self.simulations[0].xunit],
                               yunit = waveform.yunit)

    def __repr__(self):
        if len(self.simulations) > 0:
            return self.__class__.__name__ + '(' + ','.join(map(repr, self.simulations)) + ')'
        else:
            '_Analysis([])'

    type = property(get_type)


class JWDBResultSet(result.ResultSet):
    """The JWDBResultSet class handles a JWDB waveform file used by Mentor Graphics software

    Examples

    >>> dcac = JWDBResultSet("./test/data/dcac.wdb")
    >>> [a.type for a in dcac.analyses]
    ['AC', 'DC']
    >>> dcac.analyses[1].get_sweep_variables()
    []


    >>> noise = JWDBResultSet("./test/data/noise.wdb")
    >>> [a.type for a in noise.analyses]
    ['AC', 'DC', 'Noise']

    >>> sst = JWDBResultSet("./test/data/sst.wdb")
    >>> [a.type for a in sst.analyses]
    ['Ext', 'FSST']

    >>> offset = JWDBResultSet("./test/data/offset.wdb")
    >>> [a.type for a in offset.analyses]
    ['DC']
    >>> offset.analyses[0].get_sweep_variables()
    ['CARLO', 'GAINCODE']
    >>> offset.analyses[0].getWaveformNames()
    >>> offset.analyses[0].getWaveform('V(OUT<0>)')

    """
    def __init__(self, filename):
        ## Open file
        if not os.path.isfile(filename):
            raise IOError("Cannot open WDB file %s"%filename)
        
        self.fileid = jwdb.ams_sim_open_read_file(filename)

        if self.fileid < 0:
            raise IOError("Cannot open WDB file %s"%filename)

        if not ord(jwdb.ams_sim_is_jwdb_file(self.fileid)):
            raise IOError("File '%s' is not a JWDB file"%filename)

        self.__analyses = []
        self._load()

    def __del__(self):         
        jwdb.ams_sim_close_read_file(self.fileid)

    def _load(self):
        ## Get time and date
        self.date = jwdb.ams_sim_get_creation_date(self.fileid)
        self.time = jwdb.ams_sim_get_creation_time(self.fileid)

        ## Load simulation results
        simnum = 1
        simulations = []
        loadflag = jwdb.ams_sim_load_simu_with_complex(self.fileid)
        while loadflag > 0:
            simulations.append(_Simulation(self.fileid, simnum))
            
            ## Load next simulation results
            loadflag = jwdb.ams_sim_load_simu_with_complex(self.fileid)
            
            simnum += 1
            
        ## Group simulations in analyses
        simulations.sort(lambda x,y: cmp(x.analysis, y.analysis))
        for analysis, members in groupby(simulations, lambda sim: sim.analysis):
            self.__analyses.append(_Analysis(list(members)))
        
    def getResultNames(self):
        return [analysis.get_type() for analysis in self.__analyses]

    def getResult(self, name):
        for analysis in self.__analyses:
            if analysis.get_type() == name:
                return JWDBResult(analysis)


class JWDBResult(result.Result):
    def __init__(self, analysisobj):
        self._analysis = analysisobj

        result.Result.__init__(self)
        

    def getSweepDimensions(self):
        """Return the number of nested sweeps
        
        >>> dcac = JWDBResultSet("./test/data/dcac.wdb")
        >>> res = dcac.getResult('AC')
        >>> res.getSweepDimensions()
        1

        """
        return len(self._analysis.get_sweep_variables()) + 1

    def getSweepValues(self, dimension):
        """Get a numpy array of sweep values from sweep dimension dimension"""
        if dimension < len(self.__analysis.sweepvalues):
            return numpy.array(self._analysis.sweepvalues[dimension])
        else:
            ## TODO: Return x-values from waveform 
            pass

    def getSignalNames(self):
        """Returns a tuple of available signals in the result"""
        return self._analysis.getWaveformNames()

    def getSignal(self, name):
        """Returns a Waveform object if the parameter is swept and a scalar otherwise.
           If name is None it will return a dictionary of all signals keyed by
           the signal names

        >>> dcac = JWDBResultSet("./test/data/dcac.wdb")
        >>> ac = dcac.getResult('AC')
        >>> vout = result.db20(ac.getSignal('V(OUT)'))
        >>> vout.ylabel
        db20(V(OUT))
        >>> print vout.astable

        """
        return self._analysis.getWaveform(name)

if __name__ == "__main__":
    import doctest
    doctest.testmod()

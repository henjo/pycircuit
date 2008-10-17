import pycircuit.result as result
import pycircuit.waveform as waveform
import psf
import os
import re
import numpy
import operator

class PSFResultSet(result.ResultSet):
    """PSFResultSet class handles a PSF result directory

    A PSFResultSet may contain one or more runs (PSFRun objects) which in turn
    contain one or more results (PSFResult).
    """
    
    class ResultAccesser(object):
        def __init__(self, resultset):
            self._resultset = resultset
            self.__dict__.update(dict((k,None) for k in resultset.getResultNames()))
        def __getattribute__(self, attr):
            if not attr in ('__dict__', '_resultset'):
                return self._resultset.getResult(attr)
            return object.__getattribute__(self, attr)
                     
    def __init__(self, resultdir=None):
        self.resultnames = []
        self.resultdir = resultdir
        self.runObjFile = psf.PSFReader(os.path.join(resultdir,"runObjFile"), asc=True)
        self.runObjFile.open()

        runNames = self.runObjFile.getValueNames()

        self.runs = []
        self.rundict = {}
        for runName in runNames:
            run = PSFRun(runName, self.runObjFile.getValuesByName(runName), self.runObjFile.getValuePropertiesByName(runName))
            self.runs.append(run)
            self.rundict[runName] = run

        # Read log files
        for run in self.runs:
            if self.rundict.has_key(run.parentName):
                parent = self.rundict[run.parentName]
                run.setParent(parent)
                parent.addChild(run)
            run.openLogs(resultdir)

        # Add results to dictionary to allow for result introspection
        self.r = self.ResultAccesser(self)
        
    def getRuns(self):
        return self.runs

    def getRunByName(self, name):
        return self.rundict[name]

    def isfamily(self):
        """Return true if resultset is a parametric sweep

        >>> resultset=PSFResultSet('./test/resultdirs/simple')
        >>> resultset.isfamily()
        False
        >>> resultset=PSFResultSet('./test/resultdirs/parsweep/psf')
        >>> resultset.isfamily()
        True

        """
        return self.rundict.has_key("Root")

    def getResultNames(self):
        """Get name of results
        
        >>> resultset=PSFResultSet('./test/resultdirs/simple')
        >>> resultset.getResultNames()
        ('srcSweep', 'opBegin')

        >>> resultset=PSFResultSet('./test/resultdirs/parsweep/psf')
        >>> resultset.getResultNames()
        ('srcSweep', 'opBegin')

        """

        if not self.isfamily():
            return self.runs[0].getLogFile().getValueNames()
        else:
            # Find first leaf node
            for run in self.runs:
                if run.isLeaf():
                    return run.getLogFile().getValueNames()

    def getResult(self, name, run="Run1"):
        """Get result from a run by name
        
        >>> resultset=PSFResultSet('./test/resultdirs/simple')
        >>> res = resultset.getResult('srcSweep')
        >>> res.getSignalNames()
        ('VOUT', 'VIN', 'R0')

        Try attribute access

        >>> resultset=PSFResultSet('./test/resultdirs/simple')
        >>> resultset.r.srcSweep.getSignalNames()
        ('VOUT', 'VIN', 'R0')

        >>> resultset=PSFResultSet('./test/resultdirs/parsweep/psf')
        >>> res = resultset.getResult('opBegin')
        >>> res.getSignalNames()
        ('VOUT', 'VIN', 'R0')

        """
        if self.isfamily():
            return self.rundict['Root'].getResult(name, self.resultdir)
        else:
            return self.rundict[run].getResult(name, self.resultdir)

    


class PSFRun:
    def __init__(self, name, valuedict, properties):
        """
        Initiate a PSFRun

        Arguments:
        name -- Name that identifies the run
        
        """
        self.name = name
        self.logName=valuedict['logName']
        self.logDir=os.path.dirname(self.logName[0])
        self.parent=None
        self.parentName=valuedict['parent']
        self.properties = properties

        if len(valuedict['sweepVariable'])==1:
            self.sweepVariable=valuedict['sweepVariable'][0]
        elif len(valuedict['sweepVariable'])>1:
            raise ValueError("more than one sweep variable")
        else:
            self.sweepVariable=None
        
        self.sweepvalue=None
        self.logs = {}
        self.children = []

    def openLogs(self, resultdir):
        for logname in self.logName:
            if len(logname) > 0:
                self.logs[logname] = psf.PSFReader(os.path.join(resultdir, str(logname)),asc=True)
                self.logs[logname].open()                

    def addChild(self, child):
        self.children.append(child)

    def setParent(self, parent):
        if parent.sweepVariable in self.properties:
            self.sweepvalue = self.properties[parent.sweepVariable]

    def isLeaf(self):
        """Returns true if leaf node"""
        return len(self.children) == 0

    def getSweepVariables(self):
        """Get nested sweep variable names

        >>> resultset=PSFResultSet('./test/resultdirs/parsweep/psf')
        >>> run = resultset.getRunByName('Root')
        >>> run.getSweepVariables()
        ('VDC1', 'VDC2')

        """
        if self.isLeaf():
            return ()
        if not self.isLeaf():
            res = [self.sweepVariable]
            res += self.children[0].getSweepVariables()
        return tuple(res)
    
    def getLogFile(self):
        """Return a PSFReader object of log file

        >>> resultset=PSFResultSet('./test/resultdirs/simple')
        >>> run = resultset.getRunByName('Run1')
        >>> isinstance(run.getLogFile(), psf.PSFReader)
        True

        """
        for log in self.logs.keys():
            if re.search("logFile$", log):
                return self.logs[log]

    def getSweep(self):
        """Returns the sweep variable name and sweep values associated with the run
        
        >>> resultset=PSFResultSet('./test/resultdirs/parsweep/psf')
        >>> run = resultset.getRunByName('Root')
        >>> run.getSweep()
        ('VDC1', ['0.000000e+00', '1.000000e+00', '2.000000e+00'])

        """
        
        if self.isLeaf():
            raise ValueError("Run is not a parametric sweep")
        var = self.sweepVariable[0]
        return var, [child.sweepvalue for child in self.children]

    def _getNestedSweepValues(self):
        """Get all nested sweep variable values as a list of arrays with sweep values

        """
        therun = self
        sweepvalues = []

        while not therun.isLeaf():
            sweepvalues.append(numpy.array([child.sweepvalue for child in therun.children]))
            therun = therun.children[0]

        return sweepvalues

    def _findLeafs(self):
        if self.isLeaf():
            return [self]
        else:
            return reduce(operator.__add__, [child._findLeafs() for child in self.children])

    def getResultInfo(self, resultname):
        """Return information about a result such as name of data file, file format etc given a result name

        >>> resultset=PSFResultSet('./test/resultdirs/simple')
        >>> run = resultset.getRunByName('Run1')
        >>> run.getResultInfo('srcSweep')
        {'dataFile': 'srcSweep', 'description': 'Source Sweep', 'parent': '', 'format': 'PSF', 'sweepVariable': [], 'analysisType': 'dc'}

        """
        return self.getLogFile().getValuesByName(resultname)

    def getResult(self, name, dir):
        """Open result data file of result with given name located in directory dir

        >>> resultset=PSFResultSet('./test/resultdirs/simple')
        >>> run = resultset.getRunByName('Run1')
        >>> res = run.getResult('srcSweep', resultset.resultdir)
        >>> res.getSignalNames()
        ('VOUT', 'VIN', 'R0')

        >>> resultset=PSFResultSet('./test/resultdirs/parsweep/psf')
        >>> run = resultset.getRunByName('Root')
        >>> run._getNestedSweeps()
        apa
        >>> res = run.getResult('srcSweep', resultset.resultdir)
        >>> res.getSignalNames()
        ('opBegin', 'srcSweep')

        
        """
        if self.isLeaf():
            resultinfo = self.getResultInfo(name)

            if resultinfo['format'] == 'PSF':
                return PSFResult(os.path.join(dir, str(resultinfo['dataFile'])))
            else:
                raise ValueError("Cannot handle format "+resultinfo['format'])
        else:
            sweepvalues = self._getNestedSweepValues()
            leafs = self._findLeafs()
            psffiles = []
            for leaf in leafs:
                try:
                    resultinfo = leaf.getResultInfo(name)
                    psffiles.append(os.path.join(dir, leaf.logDir,
                                                 str(resultinfo['dataFile'])))
                except KeyError:
                    psffiles.append(None)

            return PSFResultFamily(self.getSweepVariables(), sweepvalues, psffiles)
        
    def __repr__(self):
        return str(self.__dict__)

    def treeString(self):
        return "\n".join([self.name]+[psf.indent(child.treeString(), 2) for child in self.children])

class PSFResult(result.Result):
    def __init__(self, psffilename=None):
        self.psfobj = psf.PSFReader(psffilename)
        self.psfobj.open()
        result.Result.__init__(self)

    def __len__(self):
        return len(self.psfobj)
        
    def getSweepDimensions(self):
        """Return the number of nested sweeps

        >>> result=PSFResult('./test/psf/srcSweep')
        >>> result.getSweepDimensions()
        1

        """
        return self.psfobj.getNSweeps()

    def getSignalNames(self):
        """Returns a tuple of available outputs in the result

        >>> result=PSFResult('./test/psf/srcSweep')
        >>> result.getSignalNames()
        ('VOUT', 'VIN', 'R0')

        >>> result=PSFResult('./test/psf/timeSweep')
        >>> result.getSignalNames()[:3]
        ('PSUP', 'INN', 'INP')

        """
        return self.psfobj.getValueNames()

    def getSweepValues(self, dimension):
        """Get a numpy array of sweep values from sweep dimension dimension

        >>> result=PSFResult('./test/psf/srcSweep')
        >>> result.getSweepValues(0)
        array([ 1.,  2.,  3.,  4.])

        >>> result=PSFResult('./test/psf/timeSweep')
        >>> result.getSweepValues(0)[:3]
        array([  0.00000000e+00,   2.00000000e-11,   5.33333333e-11])

        """
        return self.psfobj.getSweepParamValues(dimension)

    def getSignal(self, outputname=None):
        """Returns a Waveform object if the parameter is swept and a scalar otherwise. A Struct
        is given as a dictionary.
        If outputname is None it will return a dictionary of all outputs keyed by
        the output names
        
        >>> result=PSFResult('./test/psf/srcSweep')
        >>> result.getSignal()
        {'R0': Waveform([ 1.  2.  3.  4.],[-0.006 -0.004 -0.002  0.   ]), 'VIN': Waveform([ 1.  2.  3.  4.],[ 1.  2.  3.  4.]), 'VOUT': Waveform([ 1.  2.  3.  4.],[-6. -4. -2.  0.])}
        >>> result.getSignal("VOUT")
        Waveform([ 1.  2.  3.  4.],[-6. -4. -2.  0.])
        >>> result.getSignal("VIN")
        Waveform([ 1.  2.  3.  4.],[ 1.  2.  3.  4.])
        >>> result.o.VIN
        Waveform([ 1.  2.  3.  4.],[ 1.  2.  3.  4.])

        >>> result=PSFResult('./test/psf/dcOpInfo.info')
        >>> result.getSignal("R0")
        {'i': 2.5000000000000002e-06, 'res': 99999.999999999985, 'pwr': 6.2500000000000005e-07, 'v': 0.25}

        """
        if outputname == None:
            res = {}
            for output in self.getSignalNames():
                res[output] = self.getSignal(output)
            return res
        else:
            if self.psfobj.getNSweeps() > 0:
                return waveform.Waveform(self.psfobj.getSweepParamValues(0), self.psfobj.getValuesByName(outputname))
            else:
                return self.psfobj.getValuesByName(outputname)

class PSFResultFamily(result.Result):
    def __init__(self, sweepvariables, sweepvalues, psffilenames):
        self.sweepvariables = sweepvariables
        self.sweepvalues = sweepvalues
        self.psfobjects = [psf.PSFReader(psffilename) for psffilename in psffilenames]
        for po in self.psfobjects:
            po.open()
#        self.o = SignalAccesser(self)

    def getSweepDimensions(self):
        """Return the number of nested sweeps

        >>> resultset=PSFResultSet('./test/resultdirs/parsweep/psf')
        >>> result=resultset.getResult('opBegin')
        >>> result.getSweepDimensions()
        3

        """
        return len(self.sweepvariables)+self.psfobjects[0].getNSweeps()

    def getSignalNames(self):
        """Returns a tuple of available outputs in the result"""
        return self.psfobjects[0].getValueNames()

    def getSweepValues(self, dimension):
        """Get a numpy array of sweep values from the given dimension

        >>> resultset=PSFResultSet('./test/resultdirs/parsweep/psf')
        >>> result=resultset.getResult('opBegin')
        >>> result.getSweepValues(0)
        array([ 1.,  2.,  3.])

        """
        if dimension >= len(self.sweepvariables):
            return self.psfobjects[0].getSweepParamValues(dimension - len(self.sweepvariables))
        else:
            return self.sweepvalues[dimension]

    def getSignal(self, outputname=None):
        """Returns a Waveform object if the parameter is swept and a scalar otherwise. A Struct
        is given as a dictionary.
        If outputname is None it will return a dictionary of all outputs keyed by
        the output names
        
        >>> resultset=PSFResultSet('./test/resultdirs/parsweep2/psf')
        >>> result=resultset.getResult('ac-ac')
        >>> result.getSignal('net3')
        apa

        """
        if outputname == None:
            res = {}
#            for output in self.getSignalNames():
#                res[output] = self.getSignal(output)
            return res
        else:
            if self.psfobjects[0].getNSweeps() > 0:
                xvalues = self.sweepvalues + [self.psfobjects[0].getSweepParamValues(0)]
                yvalues = numpy.concatenate([psfobj.getValuesByName(outputname) for psfobj in self.psfobjects])
                yvalues = numpy.array(yvalues)
                yvalues = numpy.reshape(yvalues, map(len, xvalues))
                return waveform.Waveform(xvalues, yvalues)
            else:
                xvalues = self.sweepvalues
                yvalues = [psfobj.getValuesByName(outputname) for psfobj in self.psfobjects]
                yvalues = numpy.array(yvalues)
                yvalues = numpy.reshape(yvalues, map(len, xvalues))
                return waveform.Waveform(xvalues, yvalues)
       
if __name__ == "__main__":
    import doctest
    doctest.testmod()

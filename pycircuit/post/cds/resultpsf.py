# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import pycircuit.post.result as result
import pycircuit.post.waveform as waveform
import psf
import psfasc
import os
import re
import numpy
import operator

class PSFResultSet(result.ResultDict):
    """PSFResultSet class handles a PSF result directory

    A PSFResultSet may contain one or more runs (PSFRun objects) which in turn
    contain one or more results (PSFResult).
    """
    
    def __init__(self, resultdir=None):
        runlist = []
        self.runs = {}
        self.resultnames = []
        self.resultdir = resultdir

        ## Try to get available runs from runObjFile
        runobjfile = os.path.join(resultdir,"runObjFile")
        if os.path.exists(runobjfile):
            self.runObjFile = psf.PSFReader(runobjfile, asc=True)
            self.runObjFile.open()

            for runName in self.runObjFile.getValueNames():
                run = PSFRun(runName, self.runObjFile.getValuesByName(runName), self.runObjFile.getValuePropertiesByName(runName))
                runlist.append(run)
                self.runs[runName] = run

        ## If not present, create a single run
        else:
            run = PSFRun('Run1', {'logName': ['logFile'], 'parent': '', 'sweepVariable': [] }, {})
            self.runs['Run1'] = run
            runlist.append(run)

        # Read log files
        for run in runlist:
            if self.runs.has_key(run.parentName):
                parent = self.runs[run.parentName]
                run.setParent(parent)
                parent.addChild(run)
            run.openLogs(resultdir)

    def isfamily(self):
        """Return true if resultset is a parametric sweep

        >>> resultset=PSFResultSet('./test/resultdirs/simple')
        >>> resultset.isfamily()
        False
        >>> resultset=PSFResultSet('./test/resultdirs/parsweep/psf')
        >>> resultset.isfamily()
        True

        """
        return self.runs.has_key("Root")

    def keys(self):
        """Get name of results
        
        >>> resultset=PSFResultSet('./test/resultdirs/simple')
        >>> resultset.keys()
        ('srcSweep', 'opBegin')

        >>> resultset=PSFResultSet('./test/resultdirs/parsweep/psf')
        >>> resultset.keys()
        ('srcSweep', 'opBegin')

        """

        if not self.isfamily():
            return self.runs['Run1'].get_log().get_result_names()
        else:
            # Find first leaf node
            for run in self.runs.values():
                if run.isLeaf():
                    return run.get_log().get_result_names()

    def __getitem__(self, name):
        """Get result from a run by name
        
        >>> resultset=PSFResultSet('./test/resultdirs/simple')
        >>> res = resultset['srcSweep']
        >>> res.keys()
        ('VOUT', 'VIN', 'R0')

        Try attribute access

        >>> resultset=PSFResultSet('./test/resultdirs/simple')
        >>> resultset.r.srcSweep.keys()
        ('VOUT', 'VIN', 'R0')

        >>> resultset=PSFResultSet('./test/resultdirs/parsweep/psf')
        >>> res = resultset['opBegin']
        >>> res.keys()
        ('VOUT', 'VIN', 'R0')

        """
        if self.isfamily():
            return self.runs['Root'].getResult(name, self.resultdir)
        else:
            return self.runs['Run1'].getResult(name, self.resultdir)

class PSFRun(object):
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
                self.logs[logname] = PSFLog(os.path.join(resultdir, str(logname)))

    def addChild(self, child):
        self.children.append(child)

    def setParent(self, parent):
        if parent.sweepVariable in self.properties:
            self.sweepvalue = float(self.properties[parent.sweepVariable])

    def isLeaf(self):
        """Returns true if leaf node"""
        return len(self.children) == 0

    def getSweepVariables(self):
        """Get nested sweep variable names

        >>> resultset=PSFResultSet('./test/resultdirs/parsweep/psf')
        >>> run = resultset.runs['Root']
        >>> run.getSweepVariables()
        ('VDC1', 'VDC2')

        """
        if self.isLeaf():
            return ()
        if not self.isLeaf():
            res = [self.sweepVariable]
            res += self.children[0].getSweepVariables()
        return tuple(res)
    
    def get_log(self):
        """Return a PSFReader object of log file

        >>> resultset=PSFResultSet('./test/resultdirs/simple')
        >>> run = resultset.runs['Run1']
        >>> isinstance(run.get_log(), psf.PSFReader)
        True

        """
        for log in self.logs.keys():
            if re.search("logFile$", log):
                return self.logs[log]

    def getSweep(self):
        """Returns the sweep variable name and sweep values associated with the run
        
        >>> resultset=PSFResultSet('./test/resultdirs/parsweep/psf')
        >>> run = resultset['Root']
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
        >>> run = resultset.runs['Run1']
        >>> run.getResultInfo('srcSweep')
        {'dataFile': 'srcSweep', 'description': 'Source Sweep', 'parent': '', 'format': 'PSF', 'sweepVariable': [], 'analysisType': 'dc'}

        """
        return self.get_log().getValuesByName(resultname)

    def getResult(self, name, dir):
        """Open result data file of result with given name located in directory dir

        >>> resultset=PSFResultSet('./test/resultdirs/simple')
        >>> run = resultset.runs['Run1']
        >>> res = run.getResult('srcSweep', resultset.resultdir)
        >>> res.keys()
        ('VOUT', 'VIN', 'R0')

        >>> resultset=PSFResultSet('./test/resultdirs/parsweep/psf')
        >>> run = resultset['Root']
        >>> run._getNestedSweeps()
        apa
        >>> res = run.getResult('srcSweep', resultset.resultdir)
        >>> res.keys()
        ('opBegin', 'srcSweep')

        
        """
        if self.isLeaf():
            return self.get_log().get_result(name)
        else:
            leafs = self._findLeafs()

            firstlog = leafs[0].get_log()
            firstlogitem = firstlog[name]
            sweepvalues = self._getNestedSweepValues() + firstlogitem.get_sweep_values(firstlog)
            sweepvariables = self.getSweepVariables() + firstlogitem.get_sweep_variables(firstlog)

            psffiles = []
            for leaf in leafs:
                try:
                    logitem = leaf.get_log()[name]
                    psffiles.extend(logitem.get_datafiles(leaf.get_log()))
                except KeyError:
                    psffiles.append(None)

            print psffiles
            return PSFResultFamily(sweepvariables, sweepvalues, psffiles)
        
    def __repr__(self):
        return str(self.__dict__)

    def treeString(self):
        return "\n".join([self.name]+[psf.indent(child.treeString(), 2) for child in self.children])

class PSFLog(object):
    def __init__(self, filename):
        self.filename = filename

        self.psfobj = psf.PSFReader(filename, asc=True)
        self.psfobj.open()

        self.itemnames = []
        self.items = {}
        self.roots = []
        
        for item_name in self.psfobj.getValueNames():
            new_item = PSFLogItem(item_name, self.psfobj.getValuesByName(item_name))

            if new_item.parentname != '':
                self.items[new_item.parentname].add_child(new_item)
            else:
                self.roots.append(new_item)
                
            self.items[item_name] = new_item

    def keys(self):
        return self.itemnames

    def __getitem__(self, key):
        return self.items[key]

    def get_result_names(self):
        return [item.name for item in self.roots]

    def get_result(self, name):
        return self.items[name].get_result(self)

    def dirname(self):
        return os.path.dirname(self.filename)

class PSFLogItem(object):
    def __init__(self, name, valuedict):
        self.name = name
        self.valuedict = valuedict
        self.children = []
        self.psfobj = None

    def is_leaf(self):
        return len(self.children) == 0

    def add_child(self, item):
        self.children.append(item)

    def tree_string(self):
        return "\n".join([self.name]+[psf.indent(child.tree_string(), 2) for child in self.children])

    def get_datafile(self, log):
        return os.path.join(log.dirname(), self.valuedict['dataFile'])
    
    def get_psfobj(self, log):
        if self.psfobj is None:
            if self.valuedict['format'] == 'PSF':
                self.psfobj = create_psfreader(self.get_datafile(log))
                self.psfobj.open()
            else:
                raise ValueError("Cannot handle format " + self.valuedict['format'])

        return self.psfobj

    def get_sweep_variables(self, log):
        """Get nested sweep variable names"""
        res = []
        
        if len(self.children) > 0:
            return tuple(self.get_psfobj(log).getSweepParamNames()) + self.children[0].get_sweep_variables(log)
        else:
            return ()

    def get_sweep_values(self, log):
        """Get all nested sweep variable values as a list of arrays with sweep values

        """
        theitem = self
        sweepvalues = []

        while not theitem.is_leaf():
            sweepvalues.append(theitem.get_psfobj(log).getSweepParamValues())
            theitem = theitem.children[0]

        return sweepvalues

    def get_leafs(self):
        """Return all leafs from this log item and down in depth-first order"""
        result = []

        if self.is_leaf():
            result = [self]
        else:
            for child in self.children:
                if child.is_leaf():
                    result.append(child)
                else:
                    result.extend(child.get_leafs())
        
        return result

    def get_datafiles(self, log):
        """Return datafiles of all leafs from this log item and down in depth-first order"""
        return [item.get_datafile(log) for item in self.get_leafs()]

    def get_result(self, log):
        """Open result data file of result with given name located in directory dir"""
        if self.is_leaf():
            return PSFResult(self.get_datafile(log))
        else:
            psffiles = self.get_datafiles(log)
            return PSFResultFamily(self.get_sweep_variables(log), self.get_sweep_values(log), psffiles)

    @property
    def parentname(self):
        return self.valuedict['parent']

class PSFResult(result.ResultDict):
    def __init__(self, psffilename=None):
        self.psfobj = create_psfreader(psffilename)

        self.psfobj.open()

        super(PSFResult, self).__init__()

    def __len__(self):
        return len(self.psfobj)
        
    def getSweepDimensions(self):
        """Return the number of nested sweeps

        >>> result=PSFResult('./test/psf/srcSweep')
        >>> result.getSweepDimensions()
        1

        """
        return self.psfobj.getNSweeps()

    def keys(self):
        """Returns a tuple of available outputs in the result

        >>> result=PSFResult('./test/psf/srcSweep')
        >>> result.keys()
        ('VOUT', 'VIN', 'R0')

        >>> result=PSFResult('./test/psf/timeSweep')
        >>> result.keys()[:3]
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

    def __getitem__(self, outputname=None):
        """Returns a Waveform object if the parameter is swept and a scalar otherwise. A Struct
        is given as a dictionary.
        If outputname is None it will return a dictionary of all outputs keyed by
        the output names
        
        >>> result=PSFResult('./test/psf/srcSweep')
        >>> result[None]
        {'R0': Waveform([ 1.  2.  3.  4.],[-0.006 -0.004 -0.002  0.   ]), 'VIN': Waveform([ 1.  2.  3.  4.],[ 1.  2.  3.  4.]), 'VOUT': Waveform([ 1.  2.  3.  4.],[-6. -4. -2.  0.])}
        >>> result["VOUT"]
        Waveform([ 1.  2.  3.  4.],[-6. -4. -2.  0.])
        >>> result["VIN"]
        Waveform([ 1.  2.  3.  4.],[ 1.  2.  3.  4.])
        >>> result['VIN']
        Waveform([ 1.  2.  3.  4.],[ 1.  2.  3.  4.])

        >>> result=PSFResult('./test/psf/dcOpInfo.info')
        >>> result["R0"]
        {'i': 2.5000000000000002e-06, 'res': 99999.999999999985, 'pwr': 6.2500000000000005e-07, 'v': 0.25}

        """
        if outputname is None:
            res = {}
            for output in self.keys():
                res[output] = self[output]
            return res
        else:
            if self.psfobj.getNSweeps() > 0:
                return waveform.Waveform(self.psfobj.getSweepParamValues(0), 
                                         self.psfobj.getValuesByName(outputname),
                                         xlabels=self.psfobj.getSweepParamNames(), ylabel=outputname)
            else:
                return self.psfobj.getValuesByName(outputname)

class PSFResultFamily(result.ResultDict):
    def __init__(self, sweepvariables, sweepvalues, psffilenames):
        self.sweepvariables = sweepvariables
        self.sweepvalues = sweepvalues
        self.psfobjects = [create_psfreader(psffilename) 
                           for psffilename in psffilenames]
        for po in self.psfobjects:
            po.open()
#        self.o = SignalAccesser(self)

    def getSweepDimensions(self):
        """Return the number of nested sweeps

        >>> resultset=PSFResultSet('./test/resultdirs/parsweep/psf')
        >>> result=resultset['opBegin']
        >>> result.getSweepDimensions()
        3

        """
        return len(self.sweepvariables)+self.psfobjects[0].getNSweeps()

    def keys(self):
        """Returns a tuple of available outputs in the result"""
        return self.psfobjects[0].getValueNames()

    def getSweepValues(self, dimension):
        """Get a numpy array of sweep values from the given dimension

        >>> resultset=PSFResultSet('./test/resultdirs/parsweep/psf')
        >>> result=resultset['opBegin']
        >>> result.getSweepValues(0)
        array([ 1.,  2.,  3.])

        """
        if dimension >= len(self.sweepvariables):
            return self.psfobjects[0].getSweepParamValues(dimension - len(self.sweepvariables))
        else:
            return self.sweepvalues[dimension]

    def __getitem__(self, outputname):
        """Returns a Waveform object if the parameter is swept and a scalar otherwise. A Struct
        is given as a dictionary.
        If outputname is None it will return a dictionary of all outputs keyed by
        the output names
        
        >>> resultset=PSFResultSet('./test/resultdirs/parsweep/psf')
        >>> result=resultset['ac-ac']
        >>> result['net3']
        apa

        """
        if outputname is None:
            res = {}
#            for output in self.keys():
#                res[output] = self[output]
            return res
        else:
            if self.psfobjects[0].getNSweeps() > 0:
                xvalues = self.sweepvalues + [self.psfobjects[0].getSweepParamValues(0)]
                xlabels = self.sweepvariables + tuple(self.psfobjects[0].getSweepParamNames())
                yvalues = numpy.concatenate([psfobj.getValuesByName(outputname) for psfobj in self.psfobjects])
                yvalues = numpy.array(yvalues)
                yvalues = numpy.reshape(yvalues, map(len, xvalues))
                return waveform.Waveform(xvalues, yvalues, xlabels=xlabels, ylabel=outputname)
            else:
                xvalues = self.sweepvalues
                xlabels = self.sweepvariables
                yvalues = [psfobj.getValuesByName(outputname) for psfobj in self.psfobjects]
                yvalues = numpy.array(yvalues)
                yvalues = numpy.reshape(yvalues, map(len, xvalues))
                return waveform.Waveform(xvalues, yvalues, xlabels=xlabels, ylabel=outputname)

def create_psfreader(*args, **kwargs):
    """PSFReader factory that can use either libpsf or pure python"""
    use_libpsf = True
    
    if 'USELIBPSF' in os.environ:
        use_libpsf = bool(int(os.environ['USELIBPSF']))

    ## Try to import libpsf
    if use_libpsf:
        try:
            import psflibpsf
        except ImportError:
            use_libpsf = False
            os.environ['USELIBPSF'] = '0'

    if use_libpsf and not psfasc.is_psfasc(args[0]):
        return psflibpsf.PSFReader(*args, **kwargs)
    else:
        return psf.PSFReader(*args, **kwargs)

       
if __name__ == "__main__":
    import doctest
    doctest.testmod()

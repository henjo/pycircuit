import os
import shutil
import tempfile
import re

import cds
import skill

class InvalidSimulator(Exception):
    pass
class NetlistError(Exception):
    pass
class CellViewNotFound(Exception):
    pass

class Wrapper(object):
    def __init__(self, session, cellref):
        self.session = session
        
        wrappercellref = (cellref[0],cellref[1] + str(os.getpid()), "schematic")
        self.cellref = wrappercellref

        instid = session.dbOpenCellViewByType(*cellref)

        wrapperid = session.dbOpenCellViewByType(wrappercellref[0], 
                                                 wrappercellref[1], 
                                                 wrappercellref[2], 
                                                 "schematic", "w")
        session.dbCreateInst(wrapperid, instid, cellref[1], [0,0], "R0", 1)
        session.schCheck(wrapperid)
        session.dbSave(wrapperid)
        session.dbPurge(wrapperid)
        session.dbPurge(instid)

    def __del__(self):
        self.session.ddDeleteObj(self.session.ddGetObj(*self.cellref))

def keep_subcircuits(netlist):
    lines = []
    insubckt = 0
    for line in netlist.split('\n'):
        if line.startswith('subckt'):
            insubckt += 1

        ## Keep line if subcircuit or empty line
        if insubckt > 0 or re.match('^\s*$', line):
            lines.append(line)

        if line.startswith('ends'):
            insubckt -= 1

    return '\n'.join(lines)

def netlist_cell(session, libname, cellname, viewname, simulator='spectre', 
                 stop_view_list=('spectre', 'veriloga'),
                 switch_view_list=('simulation', 'schematic', 'spectre', 'veriloga'),
                 targetdir = None,
                 targetname = None,
                 projectdir = None,
                 write_modelincludes = False,
                 subcircuit = True):
    result = { "modelinclude_filename": None }

    if projectdir is None:
        projectdir = tempfile.mkdtemp()
        remove_projectdir_after = True
    else:
        remove_projectdir_after = False

    if targetname is None:
        targetname = cellname + '.scs'

    try:
        if targetdir is None:
            targetdir = os.curdir

        session.envSetVal('asimenv.startup', 'projectDir', skill.Symbol('string'), 
                          projectdir)
        
        if session.simulator(simulator) == skill.nil:
            raise InvalidSimulator('Invalid simulator "%s"'%simulator)

        if session.ddGetObj(libname, cellname, viewname) != skill.nil:
            session.envOption(skill.Symbol('switchViewList'), switch_view_list)
            session.envOption(skill.Symbol('stopViewList'), stop_view_list)

            if subcircuit:
                wrapper = Wrapper(session, (libname, cellname, viewname))
                libname, cellname, viewname = wrapper.cellref
                
            if session.design(libname, cellname, viewname, "r"):
                netlistfile = session.createNetlist(recreateAll=True, 
                                                    display=False)
                if netlistfile == skill.nil:
                    raise NetlistError("Could not create netlist")

                netlistdir = os.path.dirname(str(netlistfile))

                header = open(os.path.join(netlistdir, "netlistHeader")).read()
                rawnetlist = open(os.path.join(netlistdir, "netlist")).read()
                
                ## 
                if subcircuit:
                    rawnetlist = keep_subcircuits(rawnetlist)

                footer = open(os.path.join(netlistdir, "netlistFooter")).read()

                ## Write netlist 
                targetnetlist = open(os.path.join(targetdir, targetname), "w")
                targetnetlist.write(header + rawnetlist + footer)

                result["netlist_filename"] = os.path.join(targetdir, targetname)

                ## Write model definition
                if write_modelincludes:
                    modelincludes = os.path.join(targetdir, "models.scs")
                    shutil.copy(os.path.join(netlistdir, ".modelFiles"), modelincludes)
                    result["modelinclude_filename"] = modelincludes

            else:
                raise CellViewNotFound("Cellview (%s/%s/%s) not found"%(libname,cellname,viewname))

    finally:
        if remove_projectdir_after:
            shutil.rmtree(projectdir)

    return result

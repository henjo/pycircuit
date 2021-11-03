import os
import shutil
import tempfile
import re

import pycircuit.post.cds.cds as cds
import pycircuit.post.cds.skill as skill

class InvalidSimulator(Exception):
    pass
class NetlistError(Exception):
    pass
class CellViewNotFound(Exception):
    pass

class Wrapper(object):
    def __init__(self, session, cellref, 
                 stop_view_list=('spectre', 'veriloga'), 
                 switch_view_list=('simulation', 'schematic', 'spectre', 'veriloga')):
        self.session = session
        
        ## Create wrapper cellview
        wrappercellref = (cellref[0],cellref[1] + str(os.getpid()), "schematic")
        wrapperid = session.dbOpenCellViewByType(wrappercellref[0], 
                                                 wrappercellref[1], 
                                                 wrappercellref[2], 
                                                 "schematic", "w")

        self.to_be_deleted = [wrappercellref]


        # If cell to be wrapped is a config view, wrap config top cell instead
        wrapped_cell_isconfig = session.hdbIsConfig(*cellref)
        if wrapped_cell_isconfig:
            # Obtain top cell of wrapped config view
            wrapped_config = session.hdbOpen(cellref[0], cellref[1], cellref[2], "r")
            instmaster_cellref = [session.hdbGetTopLibName(wrapped_config),
                                  session.hdbGetTopCellName(wrapped_config),
                                  session.hdbGetTopViewName(wrapped_config)]
            session.hdbClose(wrapped_config)
        else:
            instmaster_cellref = cellref

        instid = session.dbOpenCellViewByType(*instmaster_cellref)

        session.dbCreateInst(wrapperid, instid, cellref[1], [0,0], "R0", 1)

        ## Save
        session.schCheck(wrapperid)
        session.dbSave(wrapperid)
        session.dbPurge(wrapperid)
        session.dbPurge(instid)

        ## Create config view
        wrapperconfig_cellref = wrappercellref[0], wrappercellref[1], "config"
        wrapperconfig = session.hdbOpen(wrapperconfig_cellref[0], wrapperconfig_cellref[1], 
                                        wrapperconfig_cellref[2], "w")
        self.to_be_deleted.append(wrapperconfig_cellref)
        
        session.hdbSetTopCellViewName(wrapperconfig, *wrappercellref)
        session.hdbSetDefaultViewListString(wrapperconfig, ' '.join(switch_view_list))
        session.hdbSetDefaultStopListString(wrapperconfig, ' '.join(stop_view_list))

        # FIXME, bind instance to wrapped cell to correct view
        session.hdbSetObjBindRule(wrapperconfig, [[instmaster_cellref[0], instmaster_cellref[1], None, None]], 
                                  [skill.Symbol('hdbcBindingRule'), [None, None, instmaster_cellref[2]]])

        session.hdbSave(wrapperconfig)
        session.hdbClose(wrapperconfig)
        
        ## Use config view for netlisting
        self.cellref = wrapperconfig_cellref

    def __del__(self):
        for cellref in self.to_be_deleted:
            self.session.ddDeleteObj(self.session.ddGetObj(*cellref))

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
                 targetdir=None,
                 targetname=None,
                 projectdir=None,
                 write_modelincludes=False,
                 write_amap=False,
                 subcircuit=True):

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
                    shutil.copytree(os.path.join(netlistdir, ".modelFiles"), modelincludes)
                    result["modelinclude_filename"] = modelincludes

                ## Copy amap and inhl directories
                if write_amap:
                    for dir in ('amap', 'ihnl'):
                        shutil.copytree(os.path.join(netlistdir, dir), 
                                        os.path.join(targetdir, dir))
            else:
                raise CellViewNotFound("Cellview (%s/%s/%s) not found"%(libname,cellname,viewname))

    finally:
        if remove_projectdir_after:
            shutil.rmtree(projectdir)

    return result

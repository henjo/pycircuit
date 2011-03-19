import os
import shutil
import tempfile

import cds
import skill

class InvalidSimulator(Exception):
    pass
class CellViewNotFound(Exception):
    pass

def netlist_cell(session, libname, cellname, viewname, simulator='spectre', 
                 stop_view_list=('spectre', 'veriloga'),
                 switch_view_list=('simulation', 'schematic', 'spectre', 'veriloga'),
                 targetdir = None,
                 targetname = None,
                 projectdir = None,
                 write_modelincludes = False):
    result = { "modelinclude_filename": None }

    if projectdir == None:
        projectdir = tempfile.mkdtemp()
        remove_projectdir_after = True
    else:
        remove_projectdir_after = False

    if targetname == None:
        targetname = cellname + '.scs'

    try:
        if targetdir == None:
            targetdir = os.curdir

        session.envSetVal('asimenv.startup', 'projectDir', skill.Symbol('string'), projectdir)
        
        if session.simulator(simulator) == skill.nil:
            raise InvalidSimulator('Invalid simulator "%s"'%simulator)

        if session.ddGetObj(libname, cellname, viewname) != skill.nil:
            session.envOption(skill.Symbol('switchViewList'), switch_view_list)
            session.envOption(skill.Symbol('stopViewList'), stop_view_list)

            if session.design(libname, cellname, viewname, "r"):
                netlistfile = session.createNetlist(recreateAll=True, display=False)

                netlistdir = os.path.dirname(netlistfile)

                header = open(os.path.join(netlistdir, "netlistHeader")).read()
                rawnetlist = open(os.path.join(netlistdir, "netlist")).read()
                footer = open(os.path.join(netlistdir, "netlistFooter")).read()

                ## Write netlist 
                targetnetlist = open(os.path.join(targetdir, targetname), "w")
                targetnetlist.write(header + rawnetlist + footer)

                result["netlist_filename"] = targetname

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
        

    

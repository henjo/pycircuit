"""
A special directive for rendering a table of parameters as defined
in pycircuit.utilities.param (Parameter)

Example
-------

    .. paramdict:: pycircuit.circuit.elements.R.instparams


"""

import re
import pycircuit.utilities.param as param
from docutils import nodes
from docutils.statemachine import ViewList
from sphinx.util.compat import Directive
import inspect

def setup(app):
    app.add_directive('paramdict', ParameterDictDirective)

class ParameterDictDirective(Directive):
    required_arguments = 1

    def run(self):
        table = nodes.table('')

        ## Create table
        group = nodes.tgroup('', cols=3)
        table.append(group)
        
        for colwidth in 10,40,5:
            group.append(nodes.colspec('', colwidth=colwidth))

        head = nodes.thead('')
        group.append(head)

        body = nodes.tbody('')
        group.append(body)

        def add_row(target, *cols):
            row = nodes.row('')
            for col in cols:
                node = nodes.paragraph(col)
                vl = ViewList()
                vl.append(col, '<apa>')
                self.state.nested_parse(vl, 0, node)
                row.append(nodes.entry('', node))

            target.append(row)

        def get_symbol(s):
            paramdict_path = s.split('.')
            
            for i in reversed(range(len(paramdict_path))):
                module = '.'.join(paramdict_path[:i])
                symbol = paramdict_path[i:]

                try:
                    m = __import__(str(module), fromlist='true')
                except ImportError:
                    continue
                else:
                    break

            parent = m
            for sym in symbol:
                print "sym", parent.__name__, sym
                parent = getattr(parent, sym)

            return parent

        add_row(head, 'Parameter', 'Description' , 'Unit')

        for param in get_symbol(self.arguments[0]):
            add_row(body, param.name, param.desc, param.unit)


        return [table]


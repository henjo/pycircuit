"""
A special directive for rendering a table of parameters as defined
in pycircuit.utilities.param (Parameter)

Example
-------

    .. parametertable:: pycircuit.circuit.elements.R.instparams


"""

import re
import pycircuit.utilities.param as param
from docutils import nodes
from docutils.statemachine import ViewList
from sphinx.util.compat import Directive
import inspect

def setup(app):
    app.add_directive('parametertable', ParameterTableDirective)

class ParameterTableDirective(Directive):
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

        def add_row(target, *column_texts):
            row = nodes.row('')
            for text in column_texts:
                if text is None:
                    text = ""
                node = nodes.paragraph('')
                vl = ViewList()
                vl.append(text, '<autosummary>')
                self.state.nested_parse(vl, 0, node)
                try:
                    if isinstance(node[0], nodes.paragraph):
                        node = node[0]
                except IndexError:
                    pass
                row.append(nodes.entry('', node))
            target.append(row)

        def get_symbol(s):
            parametertable_path = s.split('.')
            
            for i in reversed(range(len(parametertable_path))):
                module = '.'.join(parametertable_path[:i])
                symbol = parametertable_path[i:]

                try:
                    m = __import__(str(module), fromlist='true')
                except ImportError:
                    continue
                else:
                    break

            parent = m
            for sym in symbol:
                parent = getattr(parent, sym)

            return parent

        add_row(head, 'Parameter', 'Description' , 'Unit')

        for param in get_symbol(self.arguments[0]):
            add_row(body, param.name, param.desc, param.unit)

        return [table]


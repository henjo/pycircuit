# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import string

def heading1(text):
    return text + '\n' + '=' * len(text)

def heading2(text):
    return text + '\n' + '-' * len(text)

def figure(imagefile, caption=None, width="11cm"):
    s = ".. figure:: " + imagefile + "\n"
    s += "   :width: "+str(width) + "\n\n"
    if caption:
        s += "   " + caption + "\n\n"
    return s

def itemize(*items):
    """Restructured text formatted itemized list
    >>> print itemize('alpha', 'beta')
    * alpha
    * beta
    """
    
    return '\n'.join(['* ' + item for item in items])

def table(rows, header=True, headerrows = 1, vdelim=" ", padding=1, justify='right'):
    """ Outputs a list of lists as a Restructured Text Table

    - rows - list of lists
    - header - if True the first row is treated as a table header
    - vdelim - vertical delimiter betwee columns
    - padding - padding nr. of spaces are left around the longest element in the
      column
    - justify - may be left,center,right
    """

    s = ''
    
    border="=" # character for drawing the border
    justify = {'left':string.ljust,'center':string.center, 
               'right':string.rjust}[justify.lower()]

    # calculate column widhts (longest item in each col
    # plus "padding" nr of spaces on both sides)
    cols = zip(*rows)
    colWidths = [max([len(unicode(item)) + 2*padding for item in col]) 
                 for col in cols]

    # the horizontal border needed by rst
    borderline = vdelim.join([w*border for w in colWidths])

    # outputs table in rst format
    s += borderline + '\n'
    for i, row in enumerate(rows):
        if header and i < headerrows:
            justfunc = string.ljust
        else:
            justfunc = justify

        s += vdelim.join([justfunc(unicode(item),width) for (item,width) 
                          in zip(row,colWidths)]).encode('utf-8') + '\n'

        if header and headerrows == i+1: s += borderline +'\n'
    s += borderline
    return s

if __name__ == "__main__":
    import doctest
    doctest.testmod()

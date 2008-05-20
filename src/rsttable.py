import string

def toRSTtable(rows, header=True, vdelim="  ", padding=1, justify='right'):
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
    justify = {'left':string.ljust,'center':string.center, 'right':string.rjust}[justify.lower()]

    # calculate column widhts (longest item in each col
    # plus "padding" nr of spaces on both sides)
    cols = zip(*rows)
    colWidths = [max([len(unicode(item))+2*padding for item in col]) for col in cols]

    # the horizontal border needed by rst
    borderline = vdelim.join([w*border for w in colWidths])

    # outputs table in rst format
    s += borderline + '\n'
    for row in rows:
        s += vdelim.join([justify(unicode(item),width) for (item,width) in zip(row,colWidths)]).encode('utf-8') + '\n'
        if header: s += borderline +'\n'; header=False
    s += borderline
    return s

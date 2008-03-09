import skill

%%

parser skillparser:
    ignore: r'\s+'
    token FLOATNUM: r'-?[0-9]+\.[0-9e+-]*'
    token INTNUM:   r'-?[0-9]+'
    token ID:    r'[-+*/!@$%^&=.a-zA-Z0-9_]+'
    token STR:   r'"([^\\"]+|\\.)*"'
        
    ## PSFAsc grammar
    rule expr: ID                {{ return skill.Symbol(ID) }}
               | STR             {{ return eval(STR) }}
               | INTNUM          {{ return int(INTNUM) }}
               | FLOATNUM        {{ return float(FLOATNUM) }}
               | r"\("
                        {{ e = [] }}             # initialize the list
                 ( expr {{ e.append(expr) }} ) * # put each expr into the list
                 r"\)"  {{ return e }}           # return the list

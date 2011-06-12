import skill

%%

parser skillparser:
    ignore: '[ \t\n\r]+'
    token FLOATNUM: r'-?[0-9]+\.[0-9e+-]*'
    token INTNUM:   r'-?[0-9]+'
    token ID:      r'[-+*/!@$%^&=.?a-zA-Z0-9_]+'
    token OBJECT:  r'[a-z]+:0x[0-9a-f]+'
    token STR:   r'"([^\\"]+|\\.)*"'
        
    ## Skill grammar
    rule expr: ID                {{ return skill.Symbol(ID) }}
               | OBJECT          {{ return OBJECT }}
               | STR             {{ return eval(STR) }}
               | INTNUM          {{ return int(INTNUM) }}
               | FLOATNUM        {{ return float(FLOATNUM) }}
               | r"\("
                        {{ e = [] }}             # initialize the list
                 ( expr {{ e.append(expr) }} ) * # put each expr into the list
                 r"\)"  {{ return e }}           # return the list

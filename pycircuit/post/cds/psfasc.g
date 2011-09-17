import psf

def createValue(psfobj, name, typename=None, value=None):
    if psfobj.getNSweeps() == 0:
        return psf.NonSweepValue(psfobj, name=name, typeid=psfobj.types.getTypeByName(typename).id, value=value)


def updatepsf(psfobj):
        psfobj.updateHeader()
        if psfobj.getNSweeps()==0:
            psfobj.values = psf.ValuesSectionNonSweep(psfobj)
%%

parser psfasc:
    ignore: r'\s+'
    token FLOATNUM: r'-?[0-9]+\.[0-9e+-]*'
    token INTNUM:   r'-?[0-9]+'
    token STR: r'"([^\\"]+|\\.)*"'
    token EOF: r'$'
    token ARRAYLEN: r'(\*|[0-9])'

        
    ## PSFAsc grammar
    rule psfasc:                   {{ psfobj = psf.PSFReader() }}
        headersec<<psfobj>>        {{ psfobj.header = headersec }}
        [typesec<<psfobj>>
        ]                  
        [sweepsec]
        [tracesec
        ]                          {{ updatepsf(psfobj) }}
        [valuesec<<psfobj>>]
        "END"
        EOF                        {{ return psfobj }}

    ## Basic rules
    rule name: STR                 {{ return eval(STR) }}

    rule FLOAT:
        (FLOATNUM {{ return FLOATNUM }} 
        | 'Inf' {{ return 'inf' }} 
        | 'NaN' {{ return 'nan' }}
        )

    ## Header
    rule headersec<<psfobj>>:
        "HEADER"                   {{ header = psf.HeaderSection(psfobj) }}
        ( property                 {{ header.addProperty(property) }}
        ) *                        {{ return header }}

    ## Properties
    rule property: 
        name
        ( FLOAT                    {{ return psf.PropertyFloat64(name=name, value=FLOAT) }}
          | INTNUM                 {{ return psf.PropertyUInt(name=name, value=INTNUM) }}
          | STR                    {{ return psf.PropertyString(name=name, value=eval(STR)) }}
        )

    rule proplist: "PROP" r"\("    {{ props = [] }}
        ( property                 {{ props.append(property) }}
        ) *                   
        r"\)"                      {{ return props }}

    ## Type section
    rule typesec<<psfobj>>: "TYPE" {{ typesec = psf.TypeSection(psfobj) }}
        ( type<<psfobj>>           {{ psfobj.types.addType(type) }}
        ) *                        {{ return typesec }}

    rule type<<psfobj>>:                     
        name
        typedef<<psfobj>>           {{ datatypedef=psf.DataTypeDef(psfobj, name=name, datatypeid=typedef[0], structdef=typedef[1]) }}
        [ proplist                  {{ datatypedef.properties = proplist }}
        ]                           {{ return datatypedef }}

    rule typedef<<psfobj>>:
        typedefStruct<<psfobj>>     {{ return psf.TYPESTRUCT, typedefStruct }}
        | typedefFloat              {{ return psf.TYPEFLOATDOUBLE, None }}
        | typedefInt                {{ return typedefInt, None }}
        | typedefArray              {{ return psf.TYPEARRAY, typedefArray }}
        | typedefString             {{ return psf.TYPESTRING, None }}

    rule typedefFloat: "FLOAT" "DOUBLE"

    rule typedefInt: 
        "INT" 
        ( "BYTE" {{ return psf.TYPEINTBYTE }}
        | "LONG"   {{ return psf.TYPEINTLONG }}
        ) 

    rule typedefString: "STRING" "\*"

    rule typedefStruct<<psfobj>>:            
        "STRUCT\("                 {{ structdef = psf.StructDef() }}
        ( type<<psfobj>>           {{ structdef.children.append(type) }}
        ) *
        r"\)"                      {{ return structdef }}

    rule typedefArray:            
        "ARRAY" r"\(" ARRAYLEN r"\)"  
        typedef                    {{ return (typedef[0], ARRAYLEN) }}

    ## Trace section
    rule tracesec: "TRACE"
        ( tgroupdef
        ) *                        {{ return None }}
        ( tracedef
        ) *                        {{ return None }}

    rule tgroupdef:
        name "GROUP" INTNUM
        tracedef *

    rule tracedef:
        name
        STR
        [proplist]


    ## Sweep section
    rule sweepdef:
        name
        STR
        [proplist]

    rule sweepsec: "SWEEP"           
        ( sweepdef
        ) *                        {{ return None }}



    ## Value section
    rule valuesec<<psfobj>>: "VALUE" 
        ( value<<psfobj>>       {{ psfobj.values.addValue(value) }}
        ) +
        
    rule value<<psfobj>>:       {{ opttypename = None }}
        name                    
        [STR                    {{ opttypename = eval(STR) }}
        ]                       {{ value = createValue(psfobj, name, opttypename) }}
        valuedata               {{ value.setValue(valuedata) }}
        [proplist               {{ value.properties = proplist }}
        ]                       {{ return value }}

    rule nonsweptvalue:
        name
        STR
        (structarrayvalue
         | FLOAT
         | INTNUM
        )+
        [proplist]               {{ return name }}


    rule valuedata:
        structarrayvalue         {{ return structarrayvalue }}
        | FLOAT                  {{ return float(FLOAT) }}
        | INTNUM                 {{ return int(INTNUM) }}
        | STR                    {{ return eval(STR) }}

    rule structarrayvalue:
        r"\("                    {{ structval = [] }}
        (valuedata               {{ structval.append(valuedata) }}
        )*
        r"\)"                    {{ return structval }}


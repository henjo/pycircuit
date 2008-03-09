def hexdump(data):
    return " ".join(["0x%x"%ord(x) for x in data])+"\n"
def hexdumpint(data):
    return " ".join(["0x%x"%x for x in unpack(">%dI"%(len(data)/4),data)])+"\n"


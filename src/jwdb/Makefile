TARGET = jwdb.so
OFILES = sipjwdbcmodule.o
HFILES = sipAPIjwdb.h 

CC = gcc
CXX = g++
LINK = g++
CPPFLAGS = -I. -I/nfs/sw/mentor/ams//2007.2a.Linux/include -I/usr/include/python2.5
CFLAGS =  -O2 -g -fPIC -Wall -W
CXXFLAGS =  -O2 -g -fPIC -Wall -W
LFLAGS = -shared -Wl,--version-script=jwdb.exp
LIBS = -L/nfs/sw/mentor/ams//2007.2a.Linux/aol/lib -leldogwl_64
.SUFFIXES: .c .o .cpp .cc .cxx .C


.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

.C.o:
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

.c.o:
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $<

$(TARGET): $(OFILES)
	@echo '{ global: initjwdb; local: *; };' > jwdb.exp
	$(LINK) $(LFLAGS) -o $(TARGET) $(OFILES) $(LIBS)

$(OFILES): $(HFILES)

install: $(TARGET)
	@test -d $(DESTDIR)/usr/lib/python2.5/site-packages || mkdir -p $(DESTDIR)/usr/lib/python2.5/site-packages
	cp -f $(TARGET) $(DESTDIR)/usr/lib/python2.5/site-packages/$(TARGET)

clean:
	-rm -f $(TARGET)
	-rm -f sipjwdbcmodule.o
	-rm -f jwdb.exp

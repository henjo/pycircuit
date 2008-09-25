TARGET = jwdb.so
OFILES = sipjwdbcmodule.o
HFILES = sipAPIjwdb.h 

CC = gcc
CXX = g++
LINK = g++
CPPFLAGS = -I. -I/home/netcad/anacad/ams_2007.2a/include -I/home/hj03/usr.x86_64/include/python2.5
CFLAGS = -pipe -fPIC -O2 -Wall -W
CXXFLAGS = -pipe -fPIC -O2 -Wall -W
LFLAGS = -shared -Wl,-rpath,/home/netcad/anacad/ams_2007.2a/aol/lib -Wl,--version-script=jwdb.exp
LIBS = -L/home/netcad/anacad/ams_2007.2a/aol/lib -leldogwl_64
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
	@test -d $(DESTDIR)/home/hj03/usr.x86_64/lib/python2.5/site-packages || mkdir -p $(DESTDIR)/home/hj03/usr.x86_64/lib/python2.5/site-packages
	cp -f $(TARGET) $(DESTDIR)/home/hj03/usr.x86_64/lib/python2.5/site-packages/$(TARGET)
	strip $(DESTDIR)/home/hj03/usr.x86_64/lib/python2.5/site-packages/$(TARGET)

clean:
	-rm -f $(TARGET)
	-rm -f sipjwdbcmodule.o
	-rm -f jwdb.exp

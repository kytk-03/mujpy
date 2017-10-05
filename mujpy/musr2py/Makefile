prefix = /usr/local/

Linux:
	g++ -shared -Wl,-soname,MuSR2py.so -o MuSR2py.so -fPIC MuSR_td_PSI_bin.cpp MuSR2py.cpp

install: MuSR2py.so
	install -m 0755 MuSR2py.so $(prefix)/lib
#    	install -m 0644 *.png $(prefix)/share/our_program/icons
Mac:
	g++ -shared -Wl,-install_name,MuSR2py.so -o MuSR2py.so -fPIC MuSR_td_PSI_bin.cpp MuSR2py.cpp

clean:
	rm MuSR2py.so

all: pv-bare.js

pv-bare.bc: pv-bare.o fft.o hc.o
	$(CC) pv-bare.o fft.o hc.o -o pv-bare.bc

pv-bare.js: pv-bare.bc
	$(CC) pv-bare.bc /home/jdm/Downloads/libsamplerate-0.1.8/src/.libs/libsamplerate.a /home/jdm/Downloads/libsndfile-1.0.25/src/.libs/libsndfile.a /home/jdm/Downloads/fftw-3.3.3/.libs/libfftw3.a -o pv-bare.js -s INCLUDE_FULL_LIBRARY=1 -s EXPORTED_FUNCTIONS='["_pv_conventional_prefilled","_main"]'

%.o: %.c
	$(CC) -c $< -I ~/Downloads/fftw-3.3.3/api/ -I ~/Downloads/libsamplerate-0.1.8/src/

clean:
	rm -f pv-bare.bc pv-bare.o fft.o hc.o

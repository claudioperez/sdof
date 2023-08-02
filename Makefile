
src/sdof/_spectrum.so: src/_spectrum.c src/_integrate.c
	# clang -DC11THREADS -std=c11 -O3 -o thread src/spectrum.c src/_integrate.c
	gcc -std=c11 -fPIC -O3 -shared -o src/sdof/_spectrum.so src/spectrum.c src/_integrate.c -lpthread

fsdof.js:
	mkdir -p dist/
	emcc src/_integrate.c -O3 -lm -o wasm/fsdof.js \
		-s WASM=1 -s ALLOW_MEMORY_GROWTH=1 -s SINGLE_FILE=1 \
		-s EXPORTED_FUNCTIONS="['_sdof_integrate','_malloc','_free']" \
		-sINCOMING_MODULE_JS_API="['onRuntimeInitialized']" \
		-s EXPORTED_RUNTIME_METHODS="['cwrap','getValue','setValue']" 

pypa:
	python -m build

src/sdof/_integrate.%.so: src/_integrate.c Makefile
	 cc -std=c99 -pedantic -Wall -Wextra -shared -O3 src/_integrate.c -o $@  -fPIC -lm \
	    -fno-math-errno -fno-signaling-nans -fno-trapping-math \
	    -fassociative-math -ffast-math

# src/sdof/_integrate.%.so: src/_integrate.c Makefile
# 	 gcc -g -fsanitize=address -std=c99 -pedantic -Wall -Wextra -shared -Og src/_integrate.c -o $@  -fPIC -lm \
		 

thread: src/_spectrum.c src/_integrate.c
	# clang -DC11THREADS -std=c11 -O3 -o thread src/tsdof.c src/fsdof.c
	gcc -std=c11 -DHAVE_MAIN -O3 -o thread src/_spectrum.c src/_integrate.c -lpthread


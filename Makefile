
src/sdof/_tsdof.so: src/tsdof.c src/fsdof.c
	# clang -DC11THREADS -std=c11 -O3 -o thread src/tsdof.c src/fsdof.c
	gcc -std=c11 -fPIC -O3 -shared -o src/sdof/_tsdof.so src/tsdof.c src/fsdof.c -lpthread

fsdof.js:
	mkdir -p dist/
	emcc src/fsdof.c -O3 -lm -o wasm/fsdof.js \
		-s WASM=1 -s ALLOW_MEMORY_GROWTH=1 -s SINGLE_FILE=1 \
		-s EXPORTED_FUNCTIONS="['_fsdof_integrate2','_malloc','_free']" \
		-sINCOMING_MODULE_JS_API="['onRuntimeInitialized']" \
		-s EXPORTED_RUNTIME_METHODS="['cwrap','getValue','setValue']" 

pypa:
	python -m build


_fsdof.so:
	 cc -std=c99 -pedantic -Wall -Wextra -shared -O3 src/fsdof.c -o _fsdof.so  -fPIC -lm \
	    -fno-math-errno -fno-signaling-nans -fno-trapping-math \
	    -fassociative-math -ffast-math


alpha: alpha.c
	$(CC) $@ $< -lm -o alpha

thread: src/tsdof.c src/fsdof.c
	# clang -DC11THREADS -std=c11 -O3 -o thread src/tsdof.c src/fsdof.c
	gcc -std=c11 -DHAVE_MAIN -O3 -o thread src/tsdof.c src/fsdof.c -lpthread


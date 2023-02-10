

fsdof.js:
	mkdir -p dist/
	emcc src/fsdof.c -lm -o dist/fsdof.js \
		-s WASM=1 -s ALLOW_MEMORY_GROWTH=1 \
		-s EXPORTED_FUNCTIONS="['_fsdof_integrate2','_malloc']" \
		-sEXPORTED_RUNTIME_METHODS="['cwrap','getValue','setValue']"

_fsdof.so:
	 cc -std=c99 -pedantic -Wall -Wextra -shared -O3 src/fsdof.c -o _fsdof.so  -fPIC -lm \
	    -fno-math-errno -fno-signaling-nans -fno-trapping-math \
	    -fassociative-math -ffast-math


alpha: alpha.c
	$(CC) $@ $< -lm -o alpha


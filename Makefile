generalized_alpha.js:
	emcc generalized_alpha.c -lm -o dist/generalized_alpha.js -s EXPORTED_FUNCTIONS="['_generalized_alpha']"

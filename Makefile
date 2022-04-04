generalized_alpha.js:
	emcc generalized_alpha.c -lm -o dist/generalized_alpha.js -s EXPORTED_FUNCTIONS="['_generalized_alpha']"

readme: FORCE
	cat theory.md | sed 's/\\ddot U/ü/g' | sed 's/\\dot U/\\dot u/g' | pandoc -o README.md -t html

FORCE:

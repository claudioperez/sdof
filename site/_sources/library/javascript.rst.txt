

Javascript
==========


-  Install ``emscripten`` from `here <https://emscripten.org/>`__

-  run ``make``. This creates the following files:

   -  ``dist/fsdof.wasm`` - Web assembly - compiled library,
   -  ``dist/fsdof.js`` - interface to binary ``fsdof.wasm``

-  to test, you can use Python to start an HTTP server in the current
   directory as follows:

   .. code:: shell

      python -m http.server .


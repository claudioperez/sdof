// "declare" c function
const fsdof = Module.cwrap('fsdof_integrate2', 
  'number', ['number', 
             'number', 'number', 'number',
             'number', 'number', 'array', 'number', 'number']);

// Some test values from Chopra
var M = 0.2533;
var C = 0.1592;
var K = 10.0;
var dt = 0.1;
var f  =[
         0.0000,
         5.0000,
         8.6603,
        10.0000,
         8.6603,
         5.0000,
         0.0000,
         0.0000,
         0.0000,
         0.0000,
         0.0000];


function RunFSDOF(M,C,K, f,dt) {
  var n = f.length;

  var _f = new Float64Array(f);
  // change the view to Uint8, there isnt an ABI for
  // c double arrays
  var p = new Uint8Array(_f.buffer);

  // TODO: FREE!!!
  var ptr = Module._malloc(_f.byteLength*n*3);

  var u = new Float64Array(n);
  var v = new Float64Array(n);
  var a = new Float64Array(n);

  // Initialize displacement and velocity
  for (var i=0; i<3; i++)
    Module.setValue(ptr+i*_f.BYTES_PER_ELEMENT, 0.0, 'double');

  fsdof(0,M,C,K, 1.0,n,p,dt, ptr);

  for (var i=0; i<n; i++)
  {
     u[i] = Module.getValue(ptr+(3*i+0)*_f.BYTES_PER_ELEMENT, 'double');
     v[i] = Module.getValue(ptr+(3*i+1)*_f.BYTES_PER_ELEMENT, 'double');
     a[i] = Module.getValue(ptr+(3*i+2)*_f.BYTES_PER_ELEMENT, 'double');
  }
  return {u: u, v: v, a: a};
}

function TestFSDOF() {
  var out = RunFSDOF(M,C,K, f,dt);
  console.log(out.u);
  console.log(out.v);
  console.log(out.a);
}



const em_module = require('./getCellStatus.js');

console.log(em_module.ccall("getCellStatus", [0, 3])); // using ccall

console.log(em_module._getCellStatus(1, 2)); // vdirect call


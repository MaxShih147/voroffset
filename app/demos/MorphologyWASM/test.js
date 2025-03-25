// test.js

const log = (msg) => {
    console.log(msg);
    document.getElementById('log').textContent += msg + '\n';
  };
  
  MorphologyModule().then((Module) => {
    log('✔ WASM module loaded');
  
    // Dummy input: one triangle
    const verts = new Float32Array([0, 0, 0, 1, 0, 0, 0, 1, 0]);
    const inds = new Uint32Array([0, 1, 2]);
  
    const vPtr = Module._malloc(verts.length * 4);
    const iPtr = Module._malloc(inds.length * 4);
    Module.HEAPF32.set(verts, vPtr >> 2);
    Module.HEAPU32.set(inds, iPtr >> 2);
  
    const outVPtrPtr = Module._malloc(4);
    const outVSizePtr = Module._malloc(4);
    const outIPtrPtr = Module._malloc(4);
    const outISizePtr = Module._malloc(4);
  
    // Allocate MorphologyCParams manually (32 bytes assumed)
    const paramStruct = Module._malloc(48);
    const opStr = Module.allocateUTF8('dilation');
    const methodStr = Module.allocateUTF8('ours');
  
    Module.setValue(paramStruct + 0, opStr, '*');
    Module.setValue(paramStruct + 4, methodStr, '*');
    Module.setValue(paramStruct + 8, 1.0, 'double');      // dexelSize
    Module.setValue(paramStruct + 16, 8.0, 'double');      // radius
    Module.setValue(paramStruct + 24, 0, 'i8');            // radiusInMM = false
    Module.setValue(paramStruct + 25, 64, 'i32');          // numDexels
    Module.setValue(paramStruct + 32, 0.0, 'double');      // padding
  
    const handle = Module._Morphology_Create(paramStruct);
  
    const success = Module._Morphology_Run(
      handle,
      vPtr, verts.length,
      iPtr, inds.length,
      outVPtrPtr, outVSizePtr,
      outIPtrPtr, outISizePtr
    );
  
    if (success) {
      const outVPtr = Module.getValue(outVPtrPtr, '*');
      const outVSize = Module.getValue(outVSizePtr, 'i32');
      const outIPtr = Module.getValue(outIPtrPtr, '*');
      const outISize = Module.getValue(outISizePtr, 'i32');
  
      const resultVerts = new Float32Array(Module.HEAPF32.buffer, outVPtr, outVSize * 4);
      const resultInds = new Uint32Array(Module.HEAPU32.buffer, outIPtr, outISize);
  
      log(`✔ Output verts: ${resultVerts.length}, indices: ${resultInds.length}`);
    } else {
      log('✘ Morphology_Run failed');
    }
  
    Module._Morphology_Delete(handle);
  });
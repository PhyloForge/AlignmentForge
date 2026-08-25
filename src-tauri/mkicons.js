import fs from 'fs';
import zlib from 'zlib';

function makePng(width, height) {
  const signature = Buffer.from([0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A]);
  
  // IHDR
  const ihdrData = Buffer.alloc(13);
  ihdrData.writeUInt32BE(width, 0);
  ihdrData.writeUInt32BE(height, 4);
  ihdrData[8] = 8; // 8 bit depth
  ihdrData[9] = 6; // RGBA
  ihdrData[10] = 0;
  ihdrData[11] = 0;
  ihdrData[12] = 0;

  const ihdrChunk = Buffer.concat([
    Buffer.from([0x00, 0x00, 0x00, 0x0D]),
    Buffer.from('IHDR'),
    ihdrData,
    crc32(Buffer.concat([Buffer.from('IHDR'), ihdrData]))
  ]);

  // Image data
  const rawBytes = [];
  for (let y = 0; y < height; y++) {
    rawBytes.push(0); // filter
    for (let x = 0; x < width; x++) {
      rawBytes.push(30, 144, 255, 255); // RGBA
    }
  }
  const compressed = zlib.deflateSync(Buffer.from(rawBytes));
  const idatLen = Buffer.alloc(4);
  idatLen.writeUInt32BE(compressed.length, 0);
  const idatChunk = Buffer.concat([
    idatLen,
    Buffer.from('IDAT'),
    compressed,
    crc32(Buffer.concat([Buffer.from('IDAT'), compressed]))
  ]);

  const iendChunk = Buffer.from([
    0x00, 0x00, 0x00, 0x00,
    0x49, 0x45, 0x4E, 0x44,
    0xAE, 0x42, 0x60, 0x82
  ]);

  return Buffer.concat([signature, ihdrChunk, idatChunk, iendChunk]);
}

function crc32(buf) {
  let crc = 0xFFFFFFFF;
  for (let i = 0; i < buf.length; i++) {
    crc ^= buf[i];
    for (let j = 0; j < 8; j++) {
      crc = (crc >>> 1) ^ (crc & 1 ? 0xEDB88320 : 0);
    }
  }
  const res = Buffer.alloc(4);
  res.writeUInt32BE((crc ^ 0xFFFFFFFF) >>> 0, 0);
  return res;
}

if (!fs.existsSync('src-tauri/icons')) {
  fs.mkdirSync('src-tauri/icons', { recursive: true });
}

fs.writeFileSync('src-tauri/icons/32x32.png', makePng(32, 32));
fs.writeFileSync('src-tauri/icons/64x64.png', makePng(64, 64));
fs.writeFileSync('src-tauri/icons/128x128.png', makePng(128, 128));
fs.writeFileSync('src-tauri/icons/128x128@2x.png', makePng(256, 256));
fs.writeFileSync('src-tauri/icons/icon.png', makePng(512, 512));
fs.writeFileSync('src-tauri/icons/icon.ico', makePng(32, 32));
fs.writeFileSync('src-tauri/icons/icon.icns', makePng(128, 128));

console.log('Icons generated successfully with Node.js');

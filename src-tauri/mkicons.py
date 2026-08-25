#!/usr/bin/env python3
import os, struct, zlib

def make_png(width, height, color=(30, 144, 255, 255)):
    # PNG signature
    png = b'\x89PNG\r\n\x1a\n'
    # IHDR
    ihdr_data = struct.pack('>IIBBBBB', width, height, 8, 6, 0, 0, 0)
    ihdr_crc = zlib.crc32(b'IHDR' + ihdr_data) & 0xffffffff
    png += struct.pack('>I', len(ihdr_data)) + b'IHDR' + ihdr_data + struct.pack('>I', ihdr_crc)

    # Raw image scanlines
    raw_data = bytearray()
    r, g, b, a = color
    for _ in range(height):
        raw_data.append(0) # Filter byte
        for _ in range(width):
            raw_data.extend([r, g, b, a])

    compressed = zlib.compress(bytes(raw_data))
    idat_crc = zlib.crc32(b'IDAT' + compressed) & 0xffffffff
    png += struct.pack('>I', len(compressed)) + b'IDAT' + compressed + struct.pack('>I', idat_crc)

    # IEND
    iend_crc = zlib.crc32(b'IEND') & 0xffffffff
    png += struct.pack('>I', 0) + b'IEND' + struct.pack('>I', iend_crc)
    return png

os.makedirs('src-tauri/icons', exist_ok=True)

sizes = {
    'src-tauri/icons/32x32.png': (32, 32),
    'src-tauri/icons/64x64.png': (64, 64),
    'src-tauri/icons/128x128.png': (128, 128),
    'src-tauri/icons/128x128@2x.png': (256, 256),
    'src-tauri/icons/icon.png': (512, 512),
}

for path, (w, h) in sizes.items():
    with open(path, 'wb') as f:
        f.write(make_png(w, h))

# Also make dummy icon.ico and icon.icns if needed or copy icon.png
with open('src-tauri/icons/icon.ico', 'wb') as f:
    f.write(make_png(32, 32))

with open('src-tauri/icons/icon.icns', 'wb') as f:
    f.write(make_png(128, 128))

print("Icons created successfully.")

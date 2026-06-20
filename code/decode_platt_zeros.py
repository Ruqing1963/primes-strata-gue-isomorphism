#!/usr/bin/env python3
"""
Decode David Platt / LMFDB binary Riemann zero files to gamma values.

Format (Bober/Platt): 40-byte header
  [int64 block_id][double t0][double t1][int64 reserved][int64 n_zeros]
followed by n_zeros records of 13 bytes each. Each record is a little-endian
integer delta; the n-th zero is the running sum:
  gamma_n = t0 + (delta_1 + ... + delta_n) * 2^(-101)

Verified: first five decoded zeros match the published Odlyzko values
(14.134725, 21.022040, 25.010858, 30.424876, 32.935062).

Note: some larger files contain multiple blocks; this script decodes the
first/primary block. Verify monotonicity and range before use.

Author: Ruqing Chen, GUT Geoservice Inc., Montreal
"""
import struct
import numpy as np

EPS = 2.0 ** (-101)

def decode_platt(path):
    with open(path, 'rb') as f:
        data = f.read()
    t0 = struct.unpack('<d', data[8:16])[0]
    t1 = struct.unpack('<d', data[16:24])[0]
    n  = struct.unpack('<q', data[32:40])[0]
    g = np.empty(n)
    cur = t0
    for i in range(n):
        rec = data[40 + i*13 : 40 + (i+1)*13]
        if len(rec) < 13:
            g = g[:i]; break
        cur += int.from_bytes(rec, 'little') * EPS
        g[i] = cur
    return g, t0, t1

def verify(g):
    """Check first values against known Odlyzko zeros and monotonicity."""
    known = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062]
    ok = all(abs(g[i]-known[i]) < 1e-4 for i in range(min(5, len(g))))
    mono = np.all(np.diff(g) > 0)
    return ok, mono

if __name__ == '__main__':
    import sys
    path = sys.argv[1] if len(sys.argv) > 1 else 'zeros_14.dat'
    g, t0, t1 = decode_platt(path)
    ok, mono = verify(g)
    print(f"Decoded {len(g)} zeros, t in [{t0:.1f}, {t1:.1f}]")
    print(f"First five: {g[:5].round(6)}")
    print(f"Matches Odlyzko: {ok}; strictly increasing: {mono}")

import os
import subprocess
import shlex
import time
import sys
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
def main():
    # --- Power spectrum --- #
    cmd='python do_fits2pha_psd.py'
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    # --- Cross spectrum (real) --- #
    cmd='python do_fits2pha_csd_real.py'
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    # --- Cross spectrum (imaginary) --- #
    cmd='python do_fits2pha_csd_imag.py'
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    # --- Cross spectrum (amplitude) --- #
    cmd='python do_fits2pha_csd_abs.py'
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    # --- Cross spectrum (phase-lag) --- #
    cmd='python do_fits2pha_csd_phase.py'
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    # --- Cross spectrum (time-lag) --- #
    cmd='python do_fits2pha_csd_time.py'
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    # --- RMF --- #
    cmd='python do_rmfgen.py'
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    # --- ARF --- #
    cmd='python do_arfgen.py'
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

if __name__=='__main__':
    main()

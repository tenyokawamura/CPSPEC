import os
import subprocess
import shlex
import time
import sys
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
def main():
    # --- Power spectrum --- #
    cmd='python do_pyxspec_write_psds.py'
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    # --- Cross spectrum (real) --- #
    cmd='python do_pyxspec_write_csds_real.py'
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    # --- Cross spectrum (imaginary) --- #
    cmd='python do_pyxspec_write_csds_imag.py'
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    # --- Cross spectrum (amplitude) --- #
    cmd='python do_pyxspec_write_csds_abs.py'
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    # --- Cross spectrum (phase-lag) --- #
    cmd='python do_pyxspec_write_csds_phase.py'
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    # --- Cross spectrum (time-lag) --- #
    cmd='python do_pyxspec_write_csds_time.py'
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

if __name__=='__main__':
    main()

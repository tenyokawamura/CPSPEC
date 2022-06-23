import os
import subprocess
import shlex
import time
import sys
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
def main():
    # --- Power spectrum --- #
    cmd='python do_comp_psds_xspec.py'
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    # --- Cross spectrum (real) --- #
    cmd='python do_comp_csds_real_xspec.py'
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    # --- Cross spectrum (imaginary) --- #
    cmd='python do_comp_csds_imag_xspec.py'
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    # --- Cross spectrum (amplitude) --- #
    cmd='python do_comp_csds_abs_xspec.py'
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    # --- Cross spectrum (phase-lag) --- #
    cmd='python do_comp_csds_phase_xspec.py'
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

    # --- Cross spectrum (time-lag) --- #
    cmd='python do_comp_csds_time_xspec.py'
    tokens=shlex.split(cmd)
    subprocess.run(tokens)

if __name__=='__main__':
    main()

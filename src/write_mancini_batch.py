from pathlib import Path
import sys
sys.path.append('./src/')
from write_mancini_datafile_v2 import write_mancini_file

directory_path = Path("./shots/")

for subdir in directory_path.iterdir():
    cdf_file = subdir/'hyChDD.cdf'
    write_mancini_file(cdf_file)

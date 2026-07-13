import os
import glob
import shutil

class mzr_core:
    """
    A class to process observation data: creating directories, 
    moving spectral files (CaHK and H-alpha), and setting camera IDs.
    """
    
    def __init__(self, obs_dates, base_dir="."):
        if isinstance(obs_dates, str):
            self.obs_dates = [obs_dates]
        elif isinstance(obs_dates, list):
            self.obs_dates = obs_dates
        else:
            raise TypeError("Error: 'obs_dates' must be a string or a list of strings.")
            
        self.base_dir = base_dir

    def run(self):
        for date in self.obs_dates:
            print(f"[{date}] Preparing directories and files...")
            dir_c = os.path.join(self.base_dir, f"{date}C")
            dir_h = os.path.join(self.base_dir, f"{date}H")
            
            self._create_directories(dir_c, dir_h)
            self._move_data_files(dir_c, dir_h, date)
            self._create_camnum_files(dir_c, dir_h)
            print(f"[{date}] Data preparation completed.\n")

    def _create_directories(self, dir_c, dir_h):
        os.makedirs(dir_c, exist_ok=True)
        os.makedirs(dir_h, exist_ok=True)

    def _move_data_files(self, dir_c, dir_h, current_date):
        # CaHK
        c_files = glob.glob(os.path.join(self.base_dir, f"mzrc*{current_date}*"))
        if not c_files:
            c_files = glob.glob(os.path.join(self.base_dir, "mzrc*"))
        moved_c = 0
        for file_path in c_files:
            if os.path.isfile(file_path):
                shutil.move(file_path, dir_c)
                moved_c += 1
        print(f"  - Moved {moved_c} CaHK files to {dir_c}")

        # H-alpha
        h_files = glob.glob(os.path.join(self.base_dir, f"mzrh*{current_date}*"))
        if not h_files:
            h_files = glob.glob(os.path.join(self.base_dir, "mzrh*"))
        moved_h = 0
        for file_path in h_files:
            if os.path.isfile(file_path):
                shutil.move(file_path, dir_h)
                moved_h += 1
        print(f"  - Moved {moved_h} H-alpha files to {dir_h}")

    def _create_camnum_files(self, dir_c, dir_h):
        with open(os.path.join(dir_c, "camnum"), "w") as f:
            f.write("0\n")
        with open(os.path.join(dir_h, "camnum"), "w") as f:
            f.write("1\n")

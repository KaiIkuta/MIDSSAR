import os
import glob
import shutil

class mzr_analysis:
    """
    A class to process observation data: creating directories, 
    moving spectral files (CaHK and H-alpha), and setting camera IDs.
    Supports both a single date string and a list of date strings.
    """
    
    def __init__(self, obs_dates, base_dir="."):
        """
        Initialize the processor.
        
        Args:
            obs_dates (str or list): A single observation date string or a list of them 
                                     (e.g., '20260713' or ['20260713', '20260714']).
            base_dir (str): The base directory where files are located. Defaults to current directory.
        """
        # 文字列が渡された場合はリストに変換し、リストが渡された場合はそのまま保持
        if isinstance(obs_dates, str):
            self.obs_dates = [obs_dates]
        elif isinstance(obs_dates, list):
            self.obs_dates = obs_dates
        else:
            raise TypeError("Error: 'obs_dates' must be a string or a list of strings.")
            
        self.base_dir = base_dir

    def run(self):
        """Execute the entire data processing pipeline for all given dates."""
        for date in self.obs_dates:
            print(f"Starting data processing for observation date: {date}...")
            
            # Define target directory paths for the current date in loop
            dir_c = os.path.join(self.base_dir, f"{date}C")
            dir_h = os.path.join(self.base_dir, f"{date}H")
            
            self._create_directories(dir_c, dir_h)
            self._move_data_files(dir_c, dir_h, date)
            self._create_camnum_files(dir_c, dir_h)
            
            print(f"Data processing for {date} completed successfully.\n")

    def _create_directories(self, dir_c, dir_h):
        """Create directories for CaHK and H-alpha data."""
        os.makedirs(dir_c, exist_ok=True)
        os.makedirs(dir_h, exist_ok=True)
        print(f"  - Created directories: {dir_c}, {dir_h}")

    def _move_data_files(self, dir_c, dir_h, current_date):
        """Move mzrc* and mzrh* files to their respective directories."""
        
        # Process CaHK files (mzrc*)
        # ※注意: ファイル名に日付が含まれている場合（例: mzrc_20260713.fits）は、
        # "mzrc*" を f"mzrc*{current_date}*" のように変更すると日付ごとの仕分けが可能です。
        c_files = glob.glob(os.path.join(self.base_dir, "mzrc*"))
        moved_c = 0
        for file_path in c_files:
            if os.path.isfile(file_path):
                shutil.move(file_path, dir_c)
                moved_c += 1
        print(f"  - Moved {moved_c} CaHK files to {dir_c}")

        # Process H-alpha files (mzrh*)
        h_files = glob.glob(os.path.join(self.base_dir, "mzrh*"))
        moved_h = 0
        for file_path in h_files:
            if os.path.isfile(file_path):
                shutil.move(file_path, dir_h)
                moved_h += 1
        print(f"  - Moved {moved_h} H-alpha files to {dir_h}")

    def _create_camnum_files(self, dir_c, dir_h):
        """Create 'camnum' files with camera identification numbers."""
        # Write 0 for CaHK
        with open(os.path.join(dir_c, "camnum"), "w") as f:
            f.write("0\n")
        
        # Write 1 for H-alpha
        with open(os.path.join(dir_h, "camnum"), "w") as f:
            f.write("1\n")
            
        print("  - Created 'camnum' files (CaHK: 0, H-alpha: 1)")

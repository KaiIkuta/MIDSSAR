import os
import glob
import shutil
import subprocess
from mzr.core import mzr_core 

class mzr_reduction(mzr_core):
    """
    A class that wraps mzr_core and executes 
    the calibration pipeline (Bias, Flat, Sky, Anchor, and Target).
    """

    def __init__(self, obs_dates, channels=None, base_dir="."):
        self.obs_dates = [obs_dates] if isinstance(obs_dates, str) else obs_dates
        self.base_dir = base_dir
        self.channels = channels if channels else ['C', 'H']
        
        # 1. Run Data Processor
        print("Initializing pipeline: Running ObservationDataProcessor...")
        processor = mzr_core(obs_dates=self.obs_dates, base_dir=self.base_dir)
        processor.run()
        
        # 2. Setup target directories
        self.target_dirs = []
        for date in self.obs_dates:
            for ch in self.channels:
                dir_path = os.path.join(self.base_dir, f"{date}{ch}")
                if os.path.isdir(dir_path):
                    self.target_dirs.append((dir_path, ch))
                else:
                    print(f"[Warning] Directory not found and will be skipped: {dir_path}")

    def run_reduction(self):
        if not self.target_dirs:
            print("[Warning] No valid target directories found to process.")
            return

        for work_dir, ch in self.target_dirs:
            print(f"\n{'='*50}\nStarting calibration in: {work_dir} (Channel {ch})\n{'='*50}")
            
            log_files = glob.glob(os.path.join(work_dir, "mzr*.log"))
            if not log_files:
                print(f"[Error] No log file (mzr*.log) found in {work_dir}. Skipping.")
                continue
                
            log_file = log_files[0]
            print(f"Using log file: {log_file}")
            
            with open(log_file, 'r') as f:
                log_lines = f.readlines()

            # Calibration phases
            self._process_bias(work_dir, log_lines)
            self._process_high_freq_flat(work_dir, log_lines)
            self._process_sky(work_dir, log_lines)
            self._process_anchor(work_dir, log_lines)
            
            # Target (OBJECT) analysis phase
            self._process_target(work_dir, log_lines, ch)
            
            print(f"Completed processing in: {work_dir}")

    def _run_cmd(self, cmd, work_dir):
        print(f"  Running: {' '.join(cmd)}")
        subprocess.run(cmd, cwd=work_dir, check=True)

    def _remove_files(self, work_dir, pattern):
        """Helper to safely remove files matching a glob pattern."""
        files = glob.glob(os.path.join(work_dir, pattern))
        for f in files:
            try:
                os.remove(f)
            except OSError:
                pass

    def _process_bias(self, work_dir, log_lines):
        print("\n--- Processing BIAS ---")
        bias_lines = [line for line in log_lines if 'BIAS' in line]
        with open(os.path.join(work_dir, "bias"), "w") as f:
            f.writelines(bias_lines)
        self._run_cmd(["../sigclip", "bias", "bias.fits"], work_dir)

    def _process_high_freq_flat(self, work_dir, log_lines):
        print("\n--- Processing High-frequency FLAT ---")
        instflat_lines = [line for line in log_lines if 'INSTFLAT' in line]
        flat1_lines, flat2_lines = [], []
        
        for i in range(0, len(instflat_lines), 100):
            chunk = instflat_lines[i:i+100]
            if (i // 100) % 2 == 0:
                flat1_lines.extend(chunk)
            else:
                flat2_lines.extend(chunk)
                
        with open(os.path.join(work_dir, "flat1"), "w") as f: f.writelines(flat1_lines)
        with open(os.path.join(work_dir, "flat2"), "w") as f: f.writelines(flat2_lines)

        dark_cont_lines = [line for line in log_lines if 'DARK' in line and 'Cont' in line]
        with open(os.path.join(work_dir, "flatd"), "w") as f: f.writelines(dark_cont_lines)

        self._run_cmd(["../sigclip", "flatd", "flatd.fits"], work_dir)
        self._run_cmd(["../sigclip", "flat1", "flat1.fits", "flatd.fits"], work_dir)
        self._run_cmd(["../sigclip", "flat2", "flat2.fits", "flatd.fits"], work_dir)
        self._run_cmd(["../ave", "flat1.fits", "flat2.fits", "flat0.fits"], work_dir)
        self._run_cmd(["../mkflat", "flat0"], work_dir)
        
        os.rename(os.path.join(work_dir, "flat0_xy.fits"), os.path.join(work_dir, "flat.fits"))

    def _process_sky(self, work_dir, log_lines):
        print("\n--- Processing SKY ---")
        skyflat_lines = [line for line in log_lines if 'SKYFLAT' in line]
        sky_set2 = skyflat_lines[100:200] 
        with open(os.path.join(work_dir, "sky"), "w") as f: f.writelines(sky_set2)

        dark_lines = [line for line in log_lines if 'DARK' in line and 'Cont' not in line and 'FeNe' not in line]
        dark_set = dark_lines[10:20]
        with open(os.path.join(work_dir, "skyd"), "w") as f: f.writelines(dark_set)

        self._run_cmd(["../sigclip", "skyd", "skyd.fits"], work_dir)
        self._run_cmd(["../sigclip", "sky", "sky.fits", "skyd.fits", "flat.fits"], work_dir)
        self._run_cmd(["../mzwarp", "sky.fits"], work_dir)

        print("  Running 'double' 4 times...")
        for _ in range(4):
            self._run_cmd(["../double", "sky.fits", "SKY.fits", "mzwarp_dxdy.fits"], work_dir)
            shutil.copy(os.path.join(work_dir, "double_rf.fits"), os.path.join(work_dir, "flats.fits"))
            shutil.copy(os.path.join(work_dir, "double_y.dat"), os.path.join(work_dir, "flat.dat"))

    def _process_anchor(self, work_dir, log_lines):
        print("\n--- Processing ANCHOR ---")
        anchor_lines = [line for line in log_lines if 'ANCHOR' in line]
        dark_fene_lines = [line for line in log_lines if 'DARK' in line and 'FeNe' in line]

        with open(os.path.join(work_dir, "anc"), "w") as f: 
            f.writelines(anchor_lines)

        chunk_size = 10
        for i in range(0, len(anchor_lines), chunk_size):
            chunk_idx = i // chunk_size
            anc_chunk = anchor_lines[i:i+chunk_size]
            ancd_chunk = dark_fene_lines[i:i+chunk_size]
            
            if len(anc_chunk) == 0: break

            anc_name = f"anc_{chunk_idx:03d}"
            ancd_name = f"ancd_{chunk_idx:03d}"
            
            with open(os.path.join(work_dir, anc_name), "w") as f: f.writelines(anc_chunk)
            with open(os.path.join(work_dir, ancd_name), "w") as f: f.writelines(ancd_chunk)

            self._run_cmd(["../sigclip", ancd_name, f"{ancd_name}.fits"], work_dir)
            self._run_cmd(["../sigclip", anc_name, f"{anc_name}.fits", f"{ancd_name}.fits", "flat.fits"], work_dir)

    def _process_target(self, work_dir, log_lines, ch):
        """Analyze Target (OBJECT) data and generate GIFs."""
        print("\n--- Processing TARGET (OBJECT) ---")
        
        object_lines = [line for line in log_lines if 'OBJECT' in line]
        if not object_lines:
            print("[Warning] No OBJECT lines found in log. Skipping Target analysis.")
            return

        obj_file = f"OBJECT{ch}"
        with open(os.path.join(work_dir, obj_file), "w") as f:
            f.writelines(object_lines)

        obj_r = f"OBJECT_r{ch}"
        self._run_cmd(["../reduc", obj_file, obj_r, "bias.fits", "flat.fits"], work_dir)

        self._remove_files(work_dir, "spec_*.png")
        self._remove_files(work_dir, f"OBJECT_{ch}.log")

        obj_d = f"OBJECT_d{ch}"
        self._run_cmd(["../double", obj_r, obj_d, "mzwarp_dxdy.fits", "anc"], work_dir)

        spec_pngs = sorted([os.path.basename(p) for p in glob.glob(os.path.join(work_dir, "spec_*.png"))])
        if spec_pngs:
            gif_file = f"OBJECT_{ch}.gif"
            self._run_cmd(["convert", "-delay", "50", "-loop", "0"] + spec_pngs + [gif_file], work_dir)
            print(f"  - Created {gif_file}")

        self._remove_files(work_dir, "spec_*.png")


        self._run_cmd(["../merge", obj_d], work_dir)

        spec_pngs_merged = sorted([os.path.basename(p) for p in glob.glob(os.path.join(work_dir, "spec_*.png"))])
        if spec_pngs_merged:
            gif_r_file = f"OBJECT_{ch}r.gif"
            self._run_cmd(["convert", "-delay", "50", "-loop", "0"] + spec_pngs_merged + [gif_r_file], work_dir)
            print(f"  - Created {gif_r_file}")

        print(f"  - Final spectrum saved as OBJECT_{ch}.dat")

import os
import subprocess

TARGET_DIR = 'src/mzr/code' 

COMPILE_COMMANDS = [
    ["gcc", "ave.c", "-o", "ave", "-lm", "-lcfitsio", "-Wall"],
    ["gcc", "mkflat.c", "-o", "mkflat", "-lm", "-lcfitsio", "-Wall"],
    ["gcc", "mzwarp.c", "-o", "mzwarp", "-lm", "-lcfitsio", "-Wall"],
    ["gcc", "sigclip.c", "-o", "sigclip", "-lm", "-lcfitsio", "-Wall"],
    ["gcc", "reduc.c", "-o", "reduc", "-lm", "-lcfitsio", "-Wall"],
    ["gcc", "double.c", "-o", "double", "-lm", "-lcfitsio", "-Wall"],
    ["gcc", "merge.c", "-o", "merge", "-lm", "-Wall"]
]


def _compile_c_codes():
    if not os.path.exists(TARGET_DIR):
        print(f"Analysis codes are not found: {TARGET_DIR}")
        return

    print(f"Compiling codes in [{TARGET_DIR}]")

    for cmd in COMPILE_COMMANDS:
        try:
            subprocess.run(cmd, cwd=TARGET_DIR, check=True)
            print(f"Compilation Success: {' '.join(cmd)}")
        except subprocess.CalledProcessError as e:
            print(f"Error: Compilation failed ({' '.join(cmd)})\nDetails: {e}")
            raise e
        except FileNotFoundError as e:
            print("Error: 'gcc' command not found. Please check your PATH.")
            raise e

# Executed once when the module is imported
_compile_c_codes()

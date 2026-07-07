import subprocess
import os

class SigclipError(Exception):
    """Custom exception to handle errors during sigclip execution."""
    pass

def run(input_path, output_path, sigma=None, dark_image=None, flat_image=None, exec_path="../sigclip"):
    """
    A function to call the C-compiled sigclip executable from Python.

    Args:
        input_path (str): Input file name (or list file), e.g., 'bias'
        output_path (str): Output file name, e.g., 'bias.fits'
        sigma (float, optional): Clipping threshold (sigma)
        dark_image (str, optional): FITS file path for the dark image
        flat_image (str, optional): FITS file path for the flat image
        exec_path (str): Path to the compiled sigclip executable
        
    Returns:
        bool: True if execution is successful
    """
    if not os.path.isfile(exec_path):
        raise FileNotFoundError(f"Executable not found: {exec_path}\nPlease compile first using 'gcc sigclip.c -o {exec_path} -lm -lcfitsio -Wall'.")

    # Build the command
    cmd = [exec_path, str(input_path), str(output_path)]

    # Add optional arguments (following the C code's argument specifications)
    if sigma is not None:
        cmd.append(str(sigma))
    if dark_image is not None:
        cmd.append(str(dark_image))
    if flat_image is not None:
        cmd.append(str(flat_image))

    print(f"Executing: {' '.join(cmd)}")

    try:
        # Execute the C program as a subprocess
        # The C code mainly outputs logs to standard error (stderr), so we capture them
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        
        # Output upon success (Pipes C code's fprintf(stderr, ...) to Python's standard output)
        if result.stderr:
            print(result.stderr)
        if result.stdout:
            print(result.stdout)
            
        return True

    except subprocess.CalledProcessError as e:
        # Handle runtime errors (e.g., when terminated with exit(1) in C)
        print("=== Sigclip Execution Failed ===")
        print(f"Return code: {e.returncode}")
        if e.stderr:
            print(e.stderr)
        raise SigclipError("An error occurred during sigclip execution.")
"""
import sigclip

sigclip.run(
    input_path="object",
    output_path="object.fits",
    sigma=3.0,
    dark_image="dark.fits",
    flat_image="flat.fits",
    exec_path="../sigclip"
)
"""        




import subprocess
import argparse
import os
import shutil
import csv
import concurrent.futures
import glob
import re
import threading
import time

# Global variables to track progress.
progress_lock = threading.Lock()
progress_counter = 0
progress_total = 0  # Will be calculated after reading the CSV.

def call_sub(cmd, timeout=5000):
    """Runs a shell command and handles errors."""
    try:
        subprocess.run(cmd, check=True, timeout=timeout)
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
    except subprocess.TimeoutExpired as e:
        print(f"Command timed out: {e}")

def delete_event_files(base_dir, source_id, rup_id, rup_var_id, period_str):
    """
    Deletes all files in base_dir (DowntownLA) that match the event prefix for the given period string.
    """
    full_prefix = f"cs_shakemap_src_{source_id}_rup_{rup_id}_rv_{rup_var_id}_{period_str}"
    patterns = [
        f"{full_prefix}_gmpe_amps.txt",
        f"{full_prefix}_cs_amps.txt",
        f"{full_prefix}_interp.txt",
        f"{full_prefix}_interpolated.png",
        f"{full_prefix}_interpolated_marks.png"
    ]
    
    files_to_delete = []
    for pattern in patterns:
        files_to_delete.extend(glob.glob(os.path.join(base_dir, pattern)))

    for file_path in files_to_delete:
        try:
            os.remove(file_path)
            print(f"Deleted: {file_path}")
        except Exception as e:
            print(f"Error deleting {file_path}: {e}")



# def move_gmpe_files(source_id, rup_id, rup_var_id, period, gmpe):
#     """
#     Moves the generated GMPE TXT file for the given event and period into its final directory.
#     Ensures thread safety by renaming before moving.
#     """
#     base_dir = "/Users/briertonsharp/Desktop/LA_Basin/disagg_py/DowntownLA"
#     dest_base_dir = "/Users/briertonsharp/Desktop/LA_Basin/disagg_py/DowntownLA_Ruptures"

#     period_str = "pga" if period == 0.0 else f"{period}s"
#     original_filename = f"cs_shakemap_src_{source_id}_rup_{rup_id}_rv_{rup_var_id}_{period_str}_gmpe_amps.txt"
#     source_path = os.path.join(base_dir, original_filename)

#     if not os.path.exists(source_path):
#         print(f"Warning: {source_path} not found. Nothing will be moved.")
#         return

#     # Ensure thread safety by renaming the file before moving
#     unique_filename = f"cs_shakemap_src_{source_id}_rup_{rup_id}_rv_{rup_var_id}_{period_str}_{gmpe}_gmpe_amps.txt"
#     temp_path = os.path.join(base_dir, unique_filename)
#     os.rename(source_path, temp_path)  # Rename to make it unique for the GMPE

#     # Destination: /Users/.../DowntownLA_Ruptures/{source_id}/{gmpe}_{source_id}/
#     gmpe_folder = os.path.join(dest_base_dir, f"{source_id}", f"{gmpe}_{source_id}")
#     os.makedirs(gmpe_folder, exist_ok=True)  # Ensure the folder exists

#     gmpe_dest_path = os.path.join(gmpe_folder, unique_filename)

#     # Move the uniquely named file to the GMPE folder
#     shutil.move(temp_path, gmpe_dest_path)
#     print(f"Moved: {temp_path} → {gmpe_dest_path}")

def move_gmpe_files(source_id, rup_id, rup_var_id, period, gmpe):
    """
    Moves the generated GMPE TXT file for the given event and period into its structured directory,
    while renaming the file to ensure uniqueness.
    """
    base_dir = "/Users/briertonsharp/Desktop/LA_Basin/disagg_py/DowntownLA"
    dest_base_dir = "/Users/briertonsharp/Desktop/LA_Basin/disagg_py/DowntownLA_Ruptures"

    period_str = "pga" if period == 0.0 else f"{period}s"
    original_filename = f"cs_shakemap_src_{source_id}_rup_{rup_id}_rv_{rup_var_id}_{period_str}_gmpe_amps.txt"
    source_path = os.path.join(base_dir, original_filename)

    if not os.path.exists(source_path):
        print(f"Warning: {source_path} not found. Nothing will be moved.")
        return

    # Ensure thread safety by renaming the file before moving
    unique_filename = f"cs_shakemap_src_{source_id}_rup_{rup_id}_rv_{rup_var_id}_{period_str}_{gmpe}_gmpe_amps.txt"
    temp_path = os.path.join(base_dir, unique_filename)
    os.rename(source_path, temp_path)  # Rename to make it unique for the GMPE

    # Define the structured directory path
    gmpe_folder = os.path.join(dest_base_dir, f"{source_id}", f"Rup_{rup_id}", f"{gmpe}_{source_id}")
    os.makedirs(gmpe_folder, exist_ok=True)  # Ensure the full directory structure exists

    # Define the destination path
    gmpe_dest_path = os.path.join(gmpe_folder, unique_filename)

    # Move the uniquely named file to the structured directory
    shutil.move(temp_path, gmpe_dest_path)
    print(f"Moved: {temp_path} → {gmpe_dest_path}")






def plot_gmpe_maps(study, periods, source_ids, rup_ids, rup_var_ids, output_dir, gmpe_models, timeout=5000):
    for source_id in source_ids:
        for rup_id in rup_ids:
            for rup_var_id in rup_var_ids:
                for period in periods:
                    for gmpe in gmpe_models:
                        period_str = "pga" if period == 0.0 else f"{period}s"
                        filename = f"cs_shakemap_src_{source_id}_rup_{rup_id}_rv_{rup_var_id}_{period_str}_gmpe_amps.txt"
                        dest_base_dir = "/Users/briertonsharp/Desktop/LA_Basin/disagg_py/DowntownLA_Ruptures"
                        dest_folder = os.path.join(dest_base_dir, f"{source_id}", f"Rup_{rup_id}")
                        dest_path = os.path.join(dest_folder, filename)

                        if os.path.exists(dest_path):
                            print(f"Skipping {filename}, already exists at {dest_path}")
                            continue

                        print(f"Running GMPE map for Source {source_id}, Rupture {rup_id}, Variation {rup_var_id}, Period {period}, GMPE {gmpe}")
                        
                        cmd = [
                        "java", "-Xmx2G", "-cp", "opensha-cybershake-all.jar",
                        "org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator",
                        "--study", study, "--period", str(period),
                        "--source-id", str(source_id), "--rupture-id", str(rup_id),
                        "--rupture-var-id", str(rup_var_id), "--gmpe", gmpe,
                        "--spacing", "0.02",
                        "--no-plot",  # Add spacing argument here
                        "-o", output_dir
                    ]


                        call_sub(cmd, timeout=timeout)
                        move_gmpe_files(source_id, rup_id, rup_var_id, period, gmpe)
                        delete_event_files(output_dir, source_id, rup_id, rup_var_id, period_str)

import os

# def event_already_processed(source_id, rup_id, gmpe_models):
#     """
#     Checks if all GMPE files for a given source_id and rup_id already exist.
#     """
#     dest_base_dir = "/Users/briertonsharp/Desktop/LA_Basin/disagg_py/DowntownLA_Ruptures"

#     for gmpe in gmpe_models:
#         gmpe_folder = os.path.join(dest_base_dir, f"{source_id}", f"{gmpe}_{source_id}")
#         period_strs = ["pga", "0.3s", "1.0s", "3.0s"]
        
#         for period_str in period_strs:
#             filename = f"cs_shakemap_src_{source_id}_rup_{rup_id}_rv_0_{period_str}_gmpe_amps.txt"
#             file_path = os.path.join(gmpe_folder, filename)
            
#             if not os.path.exists(file_path):
#                 return False  # If any file is missing, we need to run the event
    
#     return True  # All required files exist, so we can skip processing

def event_already_processed(source_id, rup_id, gmpe_models):
    """
    Checks if all GMPE files for a given source_id and rup_id already exist.
    """
    dest_base_dir = "/Users/briertonsharp/Desktop/LA_Basin/disagg_py/DowntownLA_Ruptures"

    for gmpe in gmpe_models:
        gmpe_folder = os.path.join(dest_base_dir, f"{source_id}", f"Rup_{rup_id}", f"{gmpe}_{source_id}")
        period_strs = ["pga", "0.3s", "1.0s", "3.0s"]
        
        for period_str in period_strs:
            filename = f"cs_shakemap_src_{source_id}_rup_{rup_id}_rv_0_{period_str}_{gmpe}_gmpe_amps.txt"
            file_path = os.path.join(gmpe_folder, filename)
            
            if not os.path.exists(file_path):
                print(f"Missing file: {file_path}, event needs to be processed.")
                return False  # If any file is missing, we need to run the event
    
    print(f"Skipping Source {source_id}, Rupture {rup_id}: All GMPEs exist.")
    return True  # All required files exist, so we can skip processing



def event_runner(event, fixed_periods, gmpe_models):
    """Wrapper function to run plot_gmpe_maps() for a single event."""
    
    source_id = event["source_id"]
    rup_id = event["rupture_id"]

    if event_already_processed(source_id, rup_id, gmpe_models):
        print(f"Skipping Source {source_id}, Rupture {rup_id}: Already processed.")
        return

    print(f"Processing Source {source_id}, Rupture {rup_id}...")
    
    plot_gmpe_maps(
        study="STUDY_22_12_HF",
        periods=fixed_periods,
        source_ids=[source_id],
        rup_ids=[rup_id],
        rup_var_ids=[0],  # Force only RV = 0
        output_dir="DowntownLA",
        gmpe_models=gmpe_models
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate GMPE-based CyberShake Spectral Acceleration Maps.")
    args = parser.parse_args()
    
    fixed_periods = [0.3, 1.0, 3.0, 0.0]
    gmpe_models = ["ASK_2014", "BSSA_2014", "CB_2014", "CY_2014"]
    compiled_csv = "/Users/briertonsharp/Desktop/LA_Basin/disagg_py/DowntownLA_Ruptures/selected_events.csv"
    
    events = []
    with open(compiled_csv, "r") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            try:
                source_id = int(row["Source_ID"])
                rupture_id = int(row["Rupture_ID"])
                rup_var_ids = [int(x.strip()) for x in row["Rup_Var_ID"].split(",")]
            except ValueError:
                print(f"Skipping invalid row: {row}")
                continue
            
            events.append({
                "source_id": source_id,
                "rupture_id": rupture_id,
                "rup_var_ids": rup_var_ids
            })
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        futures = [executor.submit(event_runner, event, fixed_periods, gmpe_models) for event in events]
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print("An event failed:", e)
    
    print("All GMPE events processed.")

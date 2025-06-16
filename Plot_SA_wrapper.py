


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

# Global variables to manage staggering.
global_delay_lock = threading.Lock()
global_next_start_time = 0  # Timestamp when the next command is allowed to run.

def wait_for_stagger():
    """
    Waits until the scheduled start time for the next terminal command.
    This function ensures that each map generator command is started at least one minute apart.
    """
    global global_next_start_time
    with global_delay_lock:
        now = time.time()
        if global_next_start_time < now:
            global_next_start_time = now
        delay = global_next_start_time - now
        # Schedule next command to run 60 seconds after the current one.
        global_next_start_time += 3
    if delay > 0:
        print(f"Waiting {delay:.2f} seconds before starting next command...")
        time.sleep(delay)

def call_sub(cmd, timeout=5000):
    """Runs a shell command and handles errors."""
    try:
        subprocess.run(cmd, check=True, timeout=timeout)
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
    except subprocess.TimeoutExpired as e:
        print(f"Command timed out: {e}")

def sort_remaining_rand_mean_files(base_dir, dest_base_dir):
    """
    Scans the base_dir (DowntownLA) for any remaining rand_mean files (txt or png)
    and moves them to the appropriate destination folder based on their filename.

    Expected filename format:
      cs_shakemap_src_<source>_rup_<rupture>_rv_<rup_var>_<period>_rand_mean.<ext>
    where <period> may be something like "1.0s" or "pga".
    """
    # Regular expression to capture components: source, rupture, rupture variation, and period.
    pattern = re.compile(
        r"cs_shakemap_src_(\d+)_rup_(\d+)_rv_(\d+)_([^_]+)_rand_mean\.(txt|png)$"
    )
    
    search_pattern = os.path.join(base_dir, "*rand_mean.*")
    files = glob.glob(search_pattern)
    
    if not files:
        print("No remaining rand_mean files found to sort.")
        return
    
    for file_path in files:
        filename = os.path.basename(file_path)
        match = pattern.match(filename)
        if not match:
            print(f"Filename does not match expected pattern: {filename}")
            continue

        source_id, rup_id, rup_var_id, period_str, ext = match.groups()
        
        # Create the destination folder based on the period string.
        dest_folder = os.path.join(
            dest_base_dir, source_id, f"Rup_{rup_id}", f"RV_{rup_var_id}", period_str
        )
        os.makedirs(dest_folder, exist_ok=True)
        dest_path = os.path.join(dest_folder, filename)
        
        try:
            shutil.move(file_path, dest_path)
            print(f"Sorted and moved remaining file: {file_path} → {dest_path}")
        except Exception as e:
            print(f"Error moving {file_path} to {dest_path}: {e}")



def delete_event_files(base_dir, source_id, rup_id, rup_var_id, period_str):
    """
    Deletes all files in base_dir (DowntownLA) that match the event prefix for the given period string.
    This will match any file whose name starts with:
      cs_shakemap_src_<source_id>_rup_<rup_id>_rv_<rup_var_id>_<period_str>
    and has relevant extensions.
    """
    # Build the full prefix including the period string.
    full_prefix = f"cs_shakemap_src_{source_id}_rup_{rup_id}_rv_{rup_var_id}_{period_str}"
    
    # Patterns to delete (rand files + additional ones)
    patterns = [
        f"{full_prefix}_rand*.*",  # Existing pattern for rand files
        f"{full_prefix}_cs_amps.txt",
        f"{full_prefix}_interp.txt",
        f"{full_prefix}_interpolated.png",
        f"{full_prefix}_interpolated_marks.png"
    ]
    
    files_to_delete = []
    for pattern in patterns:
        files_to_delete.extend(glob.glob(os.path.join(base_dir, pattern)))

    if not files_to_delete:
        print(f"No associated files found for event prefix {full_prefix}.")
        return

    for file_path in files_to_delete:
        try:
            os.remove(file_path)
            print(f"Deleted: {file_path}")
        except Exception as e:
            print(f"Error deleting {file_path}: {e}")


def move_shakemap_files(source_id, rup_id, rup_var_id, period):
    """
    Moves the generated rand_mean TXT file for the given event and period from the DowntownLA folder
    into its appropriate subdirectory in DowntownLA_Ruptures. If the TXT file is not found,
    nothing is moved so that the file generation can be retried later.
    Also prints the progress percentage.
    """
    global progress_counter, progress_total

    base_dir = "/Users/briertonsharp/Desktop/LA_Basin/disagg_py/DowntownLA"
    dest_base_dir = "/Users/briertonsharp/Desktop/LA_Basin/disagg_py/DowntownLA_Ruptures"

    # For period 0.0, use "pga" instead of "0.0s"
    period_str = "pga" if period == 0.0 else f"{period}s"
    
    # Construct expected filename (only .txt is considered)
    filename_txt = f"cs_shakemap_src_{source_id}_rup_{rup_id}_rv_{rup_var_id}_{period_str}_rand_mean.txt"
    source_path = os.path.join(base_dir, filename_txt)
    
    # Check only for the TXT file.
    if not os.path.exists(source_path):
        print(f"Warning: {source_path} not found. Nothing will be moved.")
        return

    dest_folder = os.path.join(dest_base_dir, f"{source_id}", f"Rup_{rup_id}", f"RV_{rup_var_id}", period_str)
    os.makedirs(dest_folder, exist_ok=True)
    dest_path = os.path.join(dest_folder, filename_txt)

    shutil.move(source_path, dest_path)
    print(f"Moved: {source_path} → {dest_path}")

    # Update and print progress percentage.
    with progress_lock:
        progress_counter += 1
        percent = (progress_counter / progress_total) * 100
        print(f"Progress: {progress_counter}/{progress_total} ({percent:.2f}%) complete")

    # Delete any remaining associated files for this event in the base_dir.
    delete_event_files(base_dir, source_id, rup_id, rup_var_id, period_str)





def plot_SA_map(study, periods, source_ids, spatial_corr_fields, spacing, output_dir, rup_ids, rup_var_ids, no_plot=False, timeout=5000):
    for source_id in source_ids:
        for rup_id in rup_ids:
            for rup_var_id in rup_var_ids:
                for period in periods:
                    # For period 3.0, use a different study.
                    study_local = "STUDY_22_12_LF" if period == 3.0 else study

                    # Use "pga" when period is 0.0.
                    period_str = "pga" if period == 0.0 else f"{period}s"
                    
                    filename = f"cs_shakemap_src_{source_id}_rup_{rup_id}_rv_{rup_var_id}_{period_str}_rand_mean.txt"
                    dest_base_dir = "/Users/briertonsharp/Desktop/LA_Basin/disagg_py/DowntownLA_Ruptures"
                    dest_folder = os.path.join(dest_base_dir, f"{source_id}", f"Rup_{rup_id}", f"RV_{rup_var_id}", period_str)
                    dest_path = os.path.join(dest_folder, filename)

                    if os.path.exists(dest_path):
                        print(f"Skipping {filename}, already exists at {dest_path}")
                        continue

                    print(f"Running CyberShake map for Source {source_id}, Rupture {rup_id}, Variation {rup_var_id}, Period {period} using study {study_local}")
                    
                    # Wait to stagger the start of this terminal command.
                    wait_for_stagger()

                    cmd = [
                        "java", "-Xmx2G", "-cp", "opensha-cybershake-all.jar",
                        "org.opensha.sha.cybershake.maps.CyberShakeScenarioShakeMapGenerator",
                        "--study", study_local, "--period", str(period),
                        "--source-id", str(source_id), "--rupture-id", str(rup_id),
                        "--rupture-var-id", str(rup_var_id),
                        "--download-interpolated", "--spatial-corr-debug",
                        "--spatial-corr-fields", str(spatial_corr_fields),
                        "--spacing", str(spacing), "-o", output_dir
                    ]

                    if no_plot:
                        cmd.append("--no-plot")

                    call_sub(cmd, timeout=timeout)
                    move_shakemap_files(source_id, rup_id, rup_var_id, period)

    print("Processing complete.")

def event_runner(event, fixed_periods, no_plot):
    """Wrapper function to run plot_SA_map() for a single event."""
    plot_SA_map(
        study="STUDY_22_12_HF",
        periods=fixed_periods,
        source_ids=[event["source_id"]],
        spatial_corr_fields=10,
        spacing=0.02,
        output_dir="DowntownLA",
        rup_ids=[event["rupture_id"]],
        rup_var_ids=event["rup_var_ids"],
        no_plot=no_plot
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate CyberShake Spectral Acceleration Maps using CSV input and organize output concurrently."
    )
    parser.add_argument("--no-plot", action="store_true", help="Skip plotting and only generate data files")
    args = parser.parse_args()
    
    fixed_periods = [0.3, 1.0, 3.0, 0.0]
    compiled_csv = "/Users/briertonsharp/Desktop/LA_Basin/disagg_py/DowntownLA_Ruptures/selected_events_test.csv"

    # Run an extra pass to sort any remaining rand_mean files before processing events.
    downtownLA_dir = "/Users/briertonsharp/Desktop/LA_Basin/disagg_py/DowntownLA"
    dest_base_dir = "/Users/briertonsharp/Desktop/LA_Basin/disagg_py/DowntownLA_Ruptures"
    sort_remaining_rand_mean_files(downtownLA_dir, dest_base_dir)
    
    events = []
    with open(compiled_csv, "r") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            try:
                source_id = int(row["Source_ID"])
                rupture_id = int(row["Rupture_ID"])
            except ValueError:
                print(f"Skipping row with invalid Source_ID or Rupture_ID: {row}")
                continue
            
            
            try:
                rup_var_ids = [int(x.strip()) for x in row["Rup_Var_ID"].split(",")]
            except Exception as e:
                print(f"Error parsing Rup_Var_ID in row: {row}. Error: {e}")
                continue
            
            events.append({
                "source_id": source_id,
                "rupture_id": rupture_id,
                "rup_var_ids": rup_var_ids
            })
    
    print(f"Found {len(events)} events in the CSV.")
    
    # Calculate the total expected file moves:
    # For each event, each rupture variation, for each of the 4 periods.
    progress_total = sum(len(event["rup_var_ids"]) * len(fixed_periods) for event in events)
    
    max_workers = 10  # Adjust based on your computer's capacity.
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for event in events:
            print(f"Submitting event: Source {event['source_id']}, Rupture {event['rupture_id']}, Rup_Var_IDs {event['rup_var_ids']}")
            future = executor.submit(event_runner, event, fixed_periods, args.no_plot)
            futures.append(future)
        
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print("An event failed:", e)
    
    print("All events processed and remaining files sorted.")





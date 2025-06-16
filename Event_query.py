

import re
import pymysql
import glob
import os
import csv
from collections import defaultdict

# Folder containing the TXT files (each site)
folder_path = "/Users/briertonsharp/Desktop/LA_Basin/disagg_py/DowntownLA/"

# List of source IDs to loop through in the queries
source_list = [244, 242, 158, 219, 216, 245, 215, 218, 232, 184]

# Fixed list of periods to test for (in seconds)
fixed_periods = [0.2, 1.0, 3.0, 0.1]

def extract_values(file_path):
    """Extracts parameters from a text file."""
    values = {}
    with open(file_path, "r") as file:
        content = file.read()
        
        run_match = re.search(r"Run:\s*ID:\s*(\d+)", content)
        if run_match:
            values["Run_ID"] = int(run_match.group(1))
        
        erf_match = re.search(r"ERF_ID:\s*(\d+)", content)
        if erf_match:
            values["ERF_ID"] = int(erf_match.group(1))
        
        iml_match = re.search(r"for IML\s*=\s*([\d.]+)", content)
        if iml_match:
            values["IML"] = float(iml_match.group(1)) * 100  # Convert IML
        
        values["IM_Type_Component"] = "RotD50"
    
    return values

def connect_db():
    """Establishes a persistent database connection."""
    try:
        return pymysql.connect(
            host="focal.usc.edu",
            port=3306,
            user="cybershk_ro",
            password="CyberShake2007",
            db="CyberShake",
            cursorclass=pymysql.cursors.Cursor
        )
    except Exception as e:
        print("Database connection error:", e)
        return None

# Dictionary to store events found at all sites that exceed the threshold at *all* periods
final_events = {}

# Process each TXT file (each site)
txt_files = glob.glob(os.path.join(folder_path, "*.txt"))
print(f"Found {len(txt_files)} text files.")

# Open database connection once for efficiency
conn = connect_db()
if not conn:
    exit("Database connection failed.")

cursor = conn.cursor()

for file_path in txt_files:
    print(f"Processing file: {file_path}")
    params = extract_values(file_path)
    required = ["Run_ID", "ERF_ID", "IM_Type_Component", "IML"]
    if not all(k in params for k in required):
        print(f"Skipping {file_path} due to missing parameters.")
        continue

    # Temporary store for this site's events
    site_events = defaultdict(set)

    for period in fixed_periods:
        for source in source_list:
            query = f"""
            SELECT 
                P.Source_ID, 
                P.Rupture_ID, 
                P.Rup_Var_ID, 
                R.Mag
            FROM PeakAmplitudes P
            JOIN Ruptures R 
                ON P.Source_ID = R.Source_ID  
                AND P.Rupture_ID = R.Rupture_ID
            JOIN IM_Types I 
                ON P.IM_Type_ID = I.IM_Type_ID
            WHERE 
                P.Run_ID = {params["Run_ID"]}
                AND R.ERF_ID = {params["ERF_ID"]}
                AND P.Source_ID = {source}
                AND I.IM_Type_Component = '{params["IM_Type_Component"]}'
                AND I.IM_Type_Value = {period}
                AND P.IM_Value > {params["IML"]}
            ORDER BY P.IM_Value DESC;
            """

            cursor.execute(query)
            results = cursor.fetchmany(40000)  # Fetch 25000 rows at a time to avoid memory issues
            print(f"File {os.path.basename(file_path)} - Source {source}, Period {period}s returned {len(results)} rows.")

            while results:
                for row in results:
                    source_id, rupture_id, rup_var_id, mag = row
                    event_key = (source_id, rupture_id, rup_var_id, mag)
                    site_events[event_key].add(period)
                results = cursor.fetchmany(5000)  # Continue fetching in chunks

    # Filter events that appear in *all* periods
    site_events = {event: periods for event, periods in site_events.items() if len(periods) == len(fixed_periods)}

    # Update final_events to keep only events present at *all sites* for *all periods*
    if not final_events:
        final_events = site_events
    else:
        final_events = {event: periods for event, periods in final_events.items() if event in site_events}

# Close the database connection
cursor.close()
conn.close()

# Write final events to CSV
output_csv = os.path.join("/Users/briertonsharp/Desktop/LA_Basin/disagg_py/DowntownLA_Ruptures", "compiled_events.csv")
with open(output_csv, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Source_ID", "Rupture_ID", "Rup_Var_ID", "Magnitude", "Periods"])

    for event_key, periods in final_events.items():
        source_id, rupture_id, rup_var_id, mag = event_key
        period_str = ", ".join(str(p) for p in sorted(periods))
        writer.writerow([source_id, rupture_id, rup_var_id, mag, period_str])

# Count unique (Source_ID, Rupture_ID, Rup_Var_ID) combinations
unique_event_count = len(final_events)
print(f"Compiled events written to {output_csv}")
print(f"Total number of unique events: {unique_event_count}")

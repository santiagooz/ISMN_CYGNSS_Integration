import re
import datetime
import requests
import json
from requests.auth import HTTPBasicAuth
from requests.exceptions import ChunkedEncodingError
import numpy as np
import time

# Function to update input data with CYGNSS information
def update_input_data(input_data, CYGNSS_info, close_points, dist_to_station, input_vars, date):
    for var in input_vars:
        input_data[var].extend(CYGNSS_info[var][close_points])
    input_data['date'].extend([date] * len(close_points))
    input_data['dist_to_station'].extend(dist_to_station[close_points])

# Function to calculate haversine distance
def haversine(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    radius = 6371.0
    distance = radius * c
    return distance

# Function to extract values from CSV content
def extract_values(csv_content, var_list, offset):
    CYGNSS_info = {'status': False}
    start_idx = {}

    # Find starting indices for each variable in the CSV content
    for key in var_list:
        start_idx[key] = csv_content.find(key)

    # Sort indices based on their positions in the content
    sorted_start_idx = {k: v for k, v in sorted(start_idx.items(), key=lambda item: item[1])}
    start_idx_values = list(sorted_start_idx.values())
    start_idx_keys = list(sorted_start_idx.keys())

    # Check if all variables are present in the CSV content
    CYGNSS_info['status'] = not -1 in start_idx_values

    if CYGNSS_info['status']:
        # Extract values for each variable
        for j in range(len(var_list)):
            if start_idx_keys[j] != 'ddm_timestamp_utc':
                end_idx = start_idx_values[j + 1] if j < len(var_list) - 1 else len(csv_content)
                result = re.findall(r'-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?(?![^\[]*])', csv_content[start_idx_values[j]:end_idx])
                CYGNSS_info[start_idx_keys[j]] = np.array([float(val) for val in result[offset::4]])
            else:
                init_idx = start_idx_values[j] + len(start_idx_keys[j] + ', ')
                end_idx = csv_content[init_idx:].find('\n') + init_idx
                CYGNSS_info[start_idx_keys[j]] = np.array([float(val) for val in csv_content[init_idx:end_idx].split(', ')])
    return CYGNSS_info

# Function to make an OpenDAP request for CYGNSS data
def opendap_request(date, cygID, var_list, retries=3):
    time_res = 1 if (date.year >= 2019 and date.month >= 7) or date.year >= 2020 else 0.5

    D, M, YMD = str(date.day).zfill(2), str(date.month).zfill(2), str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2)
    url = f"https://opendap.earthdata.nasa.gov/collections/C2146321631-POCLOUD/granules/cyg0{cygID}.ddmi.s{YMD}-000000-e{YMD}-235959.l1.power-brcs.a31.d32.dap.csv?dap4.ce="

    for var in var_list:
        url += f"/{var};"
    url = url[:-1]

    for _ in range(retries):
        try:
            response = requests.get(url, auth=HTTPBasicAuth('santiagooz', 'Esttemilo12'))
            response.raise_for_status()  # Raise an HTTPError for bad responses
            return response
        except ChunkedEncodingError as e:
            print(f"Error: {e}")
            # Optionally, you can wait for some time before retrying
            # time.sleep(1)
        except requests.exceptions.RequestException as e:
            # Handle other request-related errors if needed
            print(f"Error: {e}")
            break  # Break out of the retry loop for non-chunked encoding errors

    # Return None or handle the case where retries are exhausted
    return None

# Load ISMN stations data from JSON file
file_path = 'H:\\Documents\\CYGNSS\\Data\\ISMN_stations.json'
with open(file_path, 'r') as json_file:
    stations_dict = json.load(json_file)

num_stations = len(stations_dict)

# Define start and end dates
start_date_str = "2022/06/30"
end_date_str = "2023/06/30"

# Convert date strings to datetime objects
start_date = datetime.datetime.strptime(start_date_str, "%Y/%m/%d")
end_date = datetime.datetime.strptime(end_date_str, "%Y/%m/%d")
number_of_days = (end_date - start_date).days

# List of variables to retrieve
var_list = ['sp_lon', 'sp_lat', 'ddm_snr', 'tx_to_sp_range', 'rx_to_sp_range', 'sp_rx_gain', 'ddm_timestamp_utc']

output_data = {'soil_moisture': []}

input_data = {var: [] for var in var_list}
input_data['date'] = []
input_data['dist_to_station'] = []

p0 = {} # stations locations
max_distance = 3  # maximum distance between ISMN sensor and SP location [km]

cntr = 0
start_time = time.time()

time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(start_time - 3*3600))
print('Process started at ' + time_str)

current_date = start_date
# Loop through each date
while current_date != end_date:
    for cygID in range(1, 9):
        req = opendap_request(current_date, cygID, var_list)

        if req is not None:
            if req.status_code == 200:
                csv_init = req.content.decode('utf-8')

                for ddmID in range(1, 5):
                    CYGNSS_info = extract_values(csv_init, var_list, ddmID - 1)
                    if CYGNSS_info['status']:
                        for i in range(num_stations):
                            p0['lat'] = float(stations_dict[str(i)]['Latitude'])
                            p0['lon'] = float(stations_dict[str(i)]['Longitude'])
                            
                            dist_to_station = haversine(CYGNSS_info['sp_lat'], CYGNSS_info['sp_lon'], p0['lat'], p0['lon'])
                            close_points = np.where(dist_to_station < max_distance)[0]
                            

                            if len(close_points) > 0:
                                print(close_points)
                                with open(stations_dict[str(i)]['file_path'], 'r') as file:
                                    lines = [line.split() for line in file]

                                dates = np.array([datetime.datetime.strptime(lines[i][0], "%Y/%m/%d") for i in range(1, len(lines))])
                                secs = np.array([int(lines[i][1][:2]) * 3600 for i in range(1, len(lines))])

                                for point in close_points:
                                    t0 = CYGNSS_info['ddm_timestamp_utc'][point]
                                    idx = np.argmax((dates == current_date) & (np.abs(secs - t0) < 1800))
                                    output_data['soil_moisture'].append(float(lines[idx + 1][2]))

                                update_input_data(input_data, CYGNSS_info, close_points, dist_to_station, var_list, current_date)

            else:
                print('404 - CYGNSS Data file not found')
        
        cntr += 1
        elapsed_time = time.time() - start_time
        remaining_time = elapsed_time / cntr * (number_of_days *8 - cntr)
        hours = int(remaining_time // 3600)
        minutes = int((remaining_time % 3600) // 60)
        seconds = int(remaining_time % 60)

        num_cases = len(output_data['soil_moisture'])
        print('Process started at ' + time_str + f' - Number of cases: {num_cases} - Time remaining: {hours} hours, {minutes} minutes, {seconds} seconds')

    input_dates = input_data['date']
    input_data['date'] = [str(date) for date in input_data['date']]
    # Save input and output data to JSON files
    with open('input_sm_data.json', 'w') as json_file:
        json.dump(input_data, json_file, indent=2)

    with open('output_sm_data.json', 'w') as json_file:
        json.dump(output_data, json_file, indent=2)
    input_data['date'] = input_dates

    current_date += datetime.timedelta(days=1)

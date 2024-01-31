import datetime
import requests
import json
from requests.auth import HTTPBasicAuth
from requests.exceptions import ChunkedEncodingError
import numpy as np
import time
import os
import netCDF4 as nc
from netCDF4 import Dataset
from scipy.interpolate import interp1d

# Function to read credentials from a file
def read_credentials(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        username = lines[0].strip()
        password = lines[1].strip()
    return username, password

# Function to calculate the width of DDM waveforms
def WF_width(DDM, K):
    peak_width = []
    for jj in range(DDM.shape[0]):
        ddm0 = DDM[jj, :, :]
        max_row_index = np.argmax(ddm0, axis=1)
        max_row_values = ddm0[range(ddm0.shape[0]), max_row_index]
        interp_function = interp1d(range(ddm0.shape[0]), max_row_values, kind='cubic')
        interp_points = np.linspace(0, ddm0.shape[0] - 1, K)
        interp_curve = interp_function(interp_points)
        max_value = np.max(interp_curve)
        threshold_value = 0.7 * max_value
        indices_below_threshold = np.where(interp_curve >= threshold_value)[0]

        # Calculate width
        if len(indices_below_threshold) >= 2:
            peak_width.append(interp_points[indices_below_threshold[-1]] - interp_points[indices_below_threshold[0]])
        else:
            peak_width.append(-1)
    return peak_width

# Function to update input data with CYGNSS information
def update_input_data(input_data, CYGNSS_info, close_points, dist_to_station, input_vars, date, ddmID):
    for var in input_vars:
        if var == 'brcs':
            ddm = CYGNSS_info[var][close_points, ddmID-1, :, :]
            input_data[var].extend([float(value) for value in WF_width(ddm, 100)])
        elif var == 'ddm_timestamp_utc':
            input_data[var].extend(CYGNSS_info[var][close_points])
        else:
            input_data[var].extend([float(value) for value in CYGNSS_info[var][close_points, ddmID-1]])

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

# Function to extract values from NetCDF content
def extract_values(response, var_list):
    CYGNSS_info = {}
    ncfile_path = 'temp.nc4'
    
    # Save the NetCDF file locally
    with open(ncfile_path, 'wb') as f:
        f.write(response.content)
    
    # Read the NetCDF file
    try:
        nc_dataset = Dataset(ncfile_path, 'r')
        for var in var_list:
            CYGNSS_info[var] = nc_dataset.variables[var][:]
        # Close the NetCDF dataset
        nc_dataset.close()
    except OSError as e:
        print(f"Error opening NetCDF file: {e}")
        CYGNSS_info = None
    
    # Delete the NetCDF file
    os.remove(ncfile_path)
    
    return CYGNSS_info

# Function to make an OpenDAP request for CYGNSS data
def opendap_request(date, cygID, var_list, retries=3):
    time_res = 1 if (date.year >= 2019 and date.month >= 7) or date.year >= 2020 else 0.5

    D, M, YMD = str(date.day).zfill(2), str(date.month).zfill(2), str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2)
    url = f"https://opendap.earthdata.nasa.gov/collections/C2146321631-POCLOUD/granules/cyg0{cygID}.ddmi.s{YMD}-000000-e{YMD}-235959.l1.power-brcs.a31.d32.dap.nc4?dap4.ce="

    for var in var_list:
        url += f"/{var};"
    url = url[:-1]
    
    username, password = read_credentials('credentials_file')
    
    for _ in range(retries):
        try:
            response = requests.get(url, auth=HTTPBasicAuth(username, password))
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
file_path = 'ISMN_stations.json'
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
var_list = ['sp_lon', 'sp_lat', 'ddm_snr', 'tx_to_sp_range', 'rx_to_sp_range', 'sp_rx_gain', 'ddm_timestamp_utc', 'brcs']

output_data = {'soil_moisture': []}

input_data = {var: [] for var in var_list}
input_data['date'] = []
input_data['dist_to_station'] = []

p0 = {} # stations locations
max_distance = 5  # maximum distance between ISMN sensor and SP location [km]

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
                CYGNSS_info = extract_values(req, var_list)
                if CYGNSS_info is not None:
                    for ddmID in range(1, 5):
                        for i in range(num_stations):
                            p0['lat'] = float(stations_dict[str(i)]['Latitude'])
                            p0['lon'] = float(stations_dict[str(i)]['Longitude'])
                            dist_to_station = haversine(CYGNSS_info['sp_lat'][:, ddmID-1], CYGNSS_info['sp_lon'][:, ddmID-1], p0['lat'], p0['lon'])
                            close_points = np.where(dist_to_station < max_distance)[0]

                            if len(close_points) > 0:
                                with open(stations_dict[str(i)]['file_path'], 'r') as file:
                                    lines = [line.split() for line in file]
                                dates = np.array([datetime.datetime.strptime(lines[i][0], "%Y/%m/%d") for i in range(1, len(lines))])
                                secs = np.array([int(lines[i][1][:2]) * 3600 for i in range(1, len(lines))])

                                for point in close_points:
                                    t0 = CYGNSS_info['ddm_timestamp_utc'][point]
                                    idx = np.argmax((dates == current_date) & (np.abs(secs - t0) < 1800))
                                    output_data['soil_moisture'].append(float(lines[idx + 1][2]))
                                update_input_data(input_data, CYGNSS_info, close_points, dist_to_station, var_list, current_date, ddmID)
            else:
                print('404 - CYGNSS Data file not found')
        
        cntr += 1
        elapsed_time = time.time() - start_time
        remaining_time = elapsed_time / cntr * (number_of_days * 8 - cntr)
        hours = int(remaining_time // 3600)
        minutes = int((remaining_time % 3600) // 60)
        seconds = int(remaining_time % 60)

        num_cases = len(output_data['soil_moisture'])
        print(f'Number of cases: {num_cases} - Time remaining: {hours} hours, {minutes} minutes, {seconds} seconds')

    input_dates = input_data['date']
    input_data['date'] = [str(date) for date in input_data['date']]
    # Save input and output data to JSON files
    with open('input_sm_data.json', 'w') as json_file:
        json.dump(input_data, json_file, indent=2, default=lambda x: float(x) if isinstance(x, np.float32) else x)

    with open('output_sm_data.json', 'w') as json_file:
        json.dump(output_data, json_file, indent=2)
    input_data['date'] = input_dates

    current_date += datetime.timedelta(days=1)
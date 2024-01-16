import os
import json

def read_stm_file(file_path):
    """Reads an STM file and returns a list of lines split into values."""
    with open(file_path, 'r') as file:
        return [line.split() for line in file]

def process_stm_line(line):
    """Processes a single line of STM file and returns a dictionary with relevant information."""
    return {
        'CSE': line[0],
        'Network': line[1],
        'Station': line[2],
        'Latitude': line[3],
        'Longitude': line[4],
        'Elevation': line[5],
        'Depth from': line[6],
        'Depth to': line[7],
    }

def process_stm_file(file_path, station_info, ctr):
    """Processes an STM file and adds station information to the dictionary."""
    lines = read_stm_file(file_path)
    
    # Skip the header line
    header = process_stm_line(lines[0])
    station_info[ctr] = {
        **header,
        'Date': [line[0] for line in lines[1:]],
        'Time': [line[1] for line in lines[1:]],
        'sm': [line[2] for line in lines[1:]],
        'file_path': file_path
    }

def process_directory(directory_path):
    """Processes all STM files in the specified directory and its subdirectories."""
    # Counter for station information entries
    ctr = 0
    # Dictionary to store station information
    station_info = {}

    # Iterate through all files in the specified directory and its subdirectories
    for root, dirs, files in os.walk(directory_path):
        for name in files:
            # Check if the file has a .stm extension
            if name.endswith(".stm"):
                # Construct the full file path
                file_path = os.path.join(root, name)
                process_stm_file(file_path, station_info, ctr)
                ctr += 1

    return station_info

def save_as_json(data, file_path):
    """Saves data as a JSON file."""
    with open(file_path, 'w') as json_file:
        json.dump(data, json_file, indent=2)

if __name__ == "__main__":
    # Specify the path to the directory containing the .stm files
    directory_path = 'H:\\Documents\\CYGNSS\\Data\\Data_separate_files_header_20220630_20230630_10792_nWPL_20240104'

    # Process the directory and get station information
    station_info = process_directory(directory_path)

    # Save the station information dictionary as a JSON file
    json_file_path = 'ISMN_stations.json'
    save_as_json(station_info, json_file_path)

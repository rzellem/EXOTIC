import csv
import subprocess

def run_tests(csv_file):
    with open(csv_file, newline='', encoding='utf-8') as file:
        reader = csv.DictReader(file)
        for row in reader:
            tic_id = row['\ufeff"TIC ID"'].strip('"')
            print(f"Running test for TIC ID: {tic_id}")
            # Call the test_exofop.py script with the TIC ID
            subprocess.run(['python', 'test_exofop_enhanced.py', tic_id], check=True)

# def run_tests(csv_file):
#     with open(csv_file, newline='', encoding='utf-8') as file:
#         reader = csv.DictReader(file)
#         # Print the fieldnames to check what headers are available
#         print(reader.fieldnames)  # Add this line to debug
#         for row in reader:
#             # Debugging line to show all available keys in the row dictionary
#             print(row.keys())  # You can remove or comment this out later
#             try:
#                 tic_id = row['TIC ID'].strip('"')  # Make sure the key name matches exactly
#                 print(f"Running test for TIC ID: {tic_id}")
#                 # Call the test_exofop.py script with the TIC ID
#                 subprocess.run(['python', 'test_exofop.py', tic_id], check=True)
#             except KeyError as e:
#                 print(f"Key error: {e} - the row keys are: {row.keys()}")
#                 continue

if __name__ == "__main__":
    csv_file = 'exofop_tess_tois.csv'
    run_tests(csv_file)

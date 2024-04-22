import sys
from exofop import ExoFOP
import json

def main(tic_id):
    # Create an instance of the ExoFOP class with the given TIC ID
    exofop_instance = ExoFOP(tic_code=tic_id)

    # Try to query ExoFOP and print the results
    try:
        data = exofop_instance.query_exofop()
        print("Raw JSON Data Retrieved:")
        print(json.dumps(data, indent=4))  # Pretty print the JSON data

        # Save raw data to a JSON file
        with open(f'tic-all-{tic_id}.json', 'w') as file:
            json.dump(data, file, indent=4)
        
        formatted_data = exofop_instance.get_formatted_data()
        print("\nFormatted Data:")
        print(json.dumps(formatted_data, indent=4))

        # Save formatted data to a JSON file
        with open(f'tic-formatted-{tic_id}.json', 'w') as file:
            json.dump(formatted_data, file, indent=4)

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python test_exofop.py <TIC ID>")
        sys.exit(1)
    main(sys.argv[1])

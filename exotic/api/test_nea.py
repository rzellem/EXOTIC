from nea import NASAExoplanetArchive
import json

def main():
    # Prompt user to enter the name of the exoplanet
    planet_name = input("Please enter the name of the exoplanet to fetch data for: ")
    
    # Create an instance of the NASAExoplanetArchive class
    nea_instance = NASAExoplanetArchive(planet=planet_name)
    
    # Try to fetch data from the NASA Exoplanet Archive and print the results
    try:
        # Attempt to resolve the name to ensure it exists in the database
        if not nea_instance.resolve_name():
            print("Exoplanet not found in the NASA Exoplanet Archive.")
            return
        
        # If the planet name is resolved, fetch the detailed information
        planet_info, is_candidate, data = nea_instance.planet_info()
        if is_candidate:
            print(f"{planet_name} is a candidate. Detailed data might not be available.")
        else:
            print("Detailed Exoplanet Data Retrieved:")
            print(json.dumps(data, indent=4))  # Assuming `data` is the final dictionary with planet info

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()

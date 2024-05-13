from nea import NASAExoplanetArchive, ExoplanetNotFoundError

def get_target():
    while True:
        target = input('Please enter the name of your exoplanet or exoplanet candidate target: ')
        if target.strip() == "":
            print("Exoplanet target may not be blank.")
        else:
            return target

def main():
    target = get_target()
    print(f"Searching for target: {target}")

    try:
        print("Trying NASAExoplanetArchive...")
        targ = NASAExoplanetArchive(planet=target)
        target = targ.planet_info()[0]
        print(f"Found target '{target}' in the NASA Exoplanet Archive")
    except ExoplanetNotFoundError as e:
        print("Caught ExoplanetNotFoundError")
        print(f"Exception message: {str(e)}")
        print("Target not found in the NASA Exoplanet Archive")
    except Exception as e:
        print("Caught unexpected exception")
        print(f"Exception message: {str(e)}")

if __name__ == "__main__":
    main()
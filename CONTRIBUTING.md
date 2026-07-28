# Contributing to EXOTIC

EXOTIC is an open source project that welcomes contributions. Please fork the repository and submit a pull request to 
the develop branch for your addition(s) to be reviewed. 

## Running the tests

Pull requests run the test suite automatically. To run it yourself, from the root of your clone:

```
python3 -m venv .venv
source .venv/bin/activate       # Windows: .venv\Scripts\activate
pip3 install -r requirements.txt -r requirements-dev.txt
pytest
```

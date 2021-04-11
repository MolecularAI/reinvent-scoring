This package contains the scoring functions for REINVENT.
To install in REINVENT's environment either install from repo or use `pip install reinvent-scoring` for the latest
official release.

For development use the enclosed environment.yml

For unit testing make sure to include a `config.json` file located in 
`\reinvent_scoring\configs\`. Use the provided example.config.json file as a template. 


Building: python setup.py sdist bdist_wheel

Upload build to test: python -m twine upload --repository testpypi dist/*

Upload build: python -m twine upload dist/*


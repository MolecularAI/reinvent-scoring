This package contains the scoring functions for REINVENT.
To install in REINVENT's environment either install from repo or use `pip install reinvent-scoring` for the latest
official release.

For development use the enclosed environment.yml


Building: python setup.py sdist bdist_wheel

Upload build to test: python -m twine upload --repository testpypi dist/*

Upload build: python -m twine upload dist/*


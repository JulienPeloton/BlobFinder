BlobFinder [![Build Status](https://travis-ci.org/JulienPeloton/BlobFinder.svg?branch=master)](https://travis-ci.org/JulienPeloton/BlobFinder)
==

Package to find point sources in sky maps.
The signal-to-noise ratio of sources is optimized with respect to the background using a matched filter,
following the work of ACT [1007.5256](https://arxiv.org/abs/1007.5256).

![ScreenShot](https://github.com/JulienPeloton/BlobFinder/blob/master/additional_files/temperature.png)

## Installation
A setup.py is provided for the installation. Just run:
```bash
python setup.py install
```
Make sure you have correct permissions (otherwise just add --user).
You can also directly use the code by updating manually your PYTHONPATH.
Just add in your bashrc:
```bash
BlobFinderPATH=/path/to/the/package
export PYTHONPATH=$PYTHONPATH:$BlobFinderPATH
```

## Example

The package contains an example using basic functionalities:
```bash
python test/test.py -setup_instrument setup_instrument.ini --plot
```

## TODO list
* Add polarisation
* Add curve-sky (input)

## License
GNU License (see the LICENSE file for details) covers all files
in the BlobFinder repository unless stated otherwise.
Some routines have been taken from this [repo](https://github.com/jeffmcm1977/CMBAnalysis_SummerSchool).

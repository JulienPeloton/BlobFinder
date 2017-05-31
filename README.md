BlobFinder
==

Package to find point sources in sky maps.

![ScreenShot](https://github.com/JulienPeloton/BlobFinder/blob/master/additional_files/temperature.png)

## Installation
We provide a setup.py for the installation. Just run:
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

We provide an example:
```bash
python test/test.py -setup_instrument setup_instrument.ini --plot
```

## TODO list
* Add polarisation
* Add curve-sky (input)

## License
GNU License (see the LICENSE file for details) covers all files
in the LaFabrique repository unless stated otherwise.

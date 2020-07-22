## Radio Source Ruler

By Tamela Maciel, 2013

**About:**

An image analysis tool to automatically measure jets and lobes of radio galaxy images.

![alt text](https://raw.githubusercontent.com/tamelamaciel/radio_source_ruler/master/thumbnail.png "radio source ruler thumbnail")


**Python libraries:**

Numpy  
Pyfits  
Shapely   
NDimage from SciPy  
Pywcs, gaussfitter, mpfit  
Matplotlib


**Key results:**

![alt text](https://raw.githubusercontent.com/tamelamaciel/radio_source_ruler/master/gallery.png "radio source gallery")

![alt text](https://raw.githubusercontent.com/tamelamaciel/radio_source_ruler/master/automated_vs_manual_size_comparison.png "accuracy of radio source ruler")

**How to run:**

python auto_size_measure.py

This script expects radio galaxy FITS files (/FITS_files/), and a database of optical ID positions ('optical_positions.dat').  
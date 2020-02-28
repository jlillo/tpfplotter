# tpfplotter
 Create paper-ready figures (1-column) overplotting the Gaia DR2 catalog to the TESS Target Pixel Files (TPF).

## Installation & Requirenments
 Clone this folder or download it to your computer. That's it!

***tpfplotter*** is written in Python2.7 (sorry!), but it will be updated to Python3 as soon as I have some time...

## Usage
Using tpfplotter is really easy, you just need to know the TIC number:

```
python tpfplotter.py 150428135 --maglim 6
```

Alternatively if you have a list of TIC values just type:

```
python tpfplotter.py list_of_tics.lis --LIST --maglim 6
```

![alt text](https://github.com/jlillo/tpfplotter/blob/master/TPF_Gaia_TIC150428135.jpg)

In case ***tpfplotter*** is not able to find the target in the Gaia catalog, you can overcome this by providing the Gaia ID and Gmag as input in the following way:

```
python tpfplotter.py 150428135 --maglim 6 --gid UCAC4 123-010026 --gmag 12.07
```

You can also save the list of Gaia sources in the TPF in a separete ascii file called Gaia_TIC*.dat with the --SAVEGAIA option. This file contains the ID, XY location on the TPF, distance to the target, Gmag and a column flag indicating if the source is inside the aperture mask (=1) or not (=0):

```
python tpfplotter.py 150428135 --maglim 6 --SAVEGAIA --sector 4
```

## Papers using tpfplotter
Several papers involving different sceince cases have already used TPF plotter. A sample of them are highlighted here:

- Aller, A., Lillo-Box, J., Jones, D., et al., 2019, arXiv e-prints, arXiv:1911.09991
- Dreizler et al., in prep

## Credits
If you use tpfplotter, please give credit to the following paper:

Aller, A., Lillo-Box, J., Jones, D., et al. (2019) "Planetary nebulae seen with TESS: Discovery of new binary central star candidates from Cycle 1," arXiv, arXiv:1911.09991 - 2019arXiv191109991A [ADS link](https://ui.adsabs.harvard.edu/abs/2019arXiv191109991A/abstract)

and add the following sentence somewhere in the paper (either footnote of acknowledgements section):
 > This work made use of \texttt{tpfplotter} by J. Lillo-Box (publicly available in www.github.com/jlillo/tpfplotter), which also made use of the python packages \texttt{astropy}, \texttt{lightkurve}, \texttt{matplotlib} and \texttt{numpy}.




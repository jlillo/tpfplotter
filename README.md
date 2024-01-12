# tpfplotter
 Create paper-ready figures (1-column) overplotting the Gaia DR3 catalog to the TESS Target Pixel Files (TPF). You can create plots for any target observed by TESS! Even if you do not have a TIC number, you can search by coordinates now (see examples below)!

![alt text](https://github.com/jlillo/tpfplotter/blob/master/logo_tpfplotter.png)


## Installation & Requirenments
 Clone this folder or download it to your computer. That's it!

 Due to the latest changes in the lightkurve package (v2 and above) you will
 need to upgrade the following packages to the corresponding versions:

```
numpy --> > 1.20.1
matplotlib --> > 3.2.1
astropy --> > 4.2
lightkurve --> > 2.0.3
```

Remember that you can upgrade your packages by just using pip as, e.g.:

```
pip install numpy --upgrade
```

***tpfplotter*** is written in both Python3.6 (tpfplotter.py) and Python2.7 (tpfplotter_py2.py). But please note that since October 2020 only the Python3 version is kept up-to-date.

Since 22-September-2022 the default Gaia catalog used is DR3. If you still want to use the DR2 catalog, please use the "--DR2" option.

## Usage
Using tpfplotter is really easy.

If you know the TIC number ():

```
python tpfplotter.py 150428135 --maglim 6
```

![alt text](https://github.com/jlillo/tpfplotter/blob/master/TPF_Gaia_TIC150428135.jpg)

Note: if the TIC is in the CTL, the mask will correspond to the pipeline mask. Otherwise, it will just show a mask obtained with tpf.create_threshold_mask(threshold=10,reference_pixel='center'):

If there is no TIC number, you can search by coordinates:

```
python tpfplotter.py TestTarget1 --COORD 166.0189667,49.15338516 --maglim 6
```

Also, if you have a list of TIC values just type:

```
python tpfplotter.py list_of_tics.lis --LIST --maglim 6
```

Your list can also specify the required options (maglim, sector, name, ra, dec) for each of the targets by including the corresponding columns. You can use the `example_list.lis` file as a template. The usage of lists in this case would be:

```
python tpfplotter.py example_list.lis --LIST
```


In case ***tpfplotter*** is not able to find the target in the Gaia catalog, you can overcome this by providing the Gaia ID and Gmag as input in the following way:

```
python tpfplotter.py 150428135 --maglim 6 --gid UCAC4 123-010026 --gmag 12.07
```

You can also save the list of Gaia sources in the TPF in a separate ascii file called Gaia_TIC*.dat with the --SAVEGAIA option. This file contains the ID, XY location on the TPF, distance to the target, Gmag and a column flag indicating if the source is inside the aperture mask (=1) or not (=0):

```
python tpfplotter.py 150428135 --maglim 6 --SAVEGAIA --sector 4
```

Note: the `--SAVEGAIA option also works with lists, saving one file per entry in the list file.`

As of June 8th, you can also plot the proper motion directions of all targets within the TPF. Simply add the "--PM" option to your command like this:

```
python tpfplotter.py 150428135 --maglim 6 --PM
```

## Papers using tpfplotter
Several papers involving different science cases have already used TPF plotter. A sample of them are highlighted here:

- Aller, A., Lillo-Box, J., Jones, D., et al., 2019, arXiv e-prints, arXiv:1911.09991  [ADS link](https://ui.adsabs.harvard.edu/abs/2020A%26A...635A.128A/abstract)
- Lillo-Box, et al., (2020), arXiv e-prints, arXiv:2010.06928. [ADS link](https://ui.adsabs.harvard.edu/abs/2020arXiv201006928L)
- Demory, et al., (2020), A&A, 642, A49. [ADS link](https://ui.adsabs.harvard.edu/abs/2020A&A...642A..49D)
- Luque, et al., (2020), arXiv e-prints, arXiv:2009.08338. [ADS link](https://ui.adsabs.harvard.edu/abs/2020arXiv200908338L)
- Bluhm, et al., (2020), A&A, 639, A132. [ADS link](https://ui.adsabs.harvard.edu/abs/2020A&A...639A.132B)
- Nowak, et al., (2020), arXiv e-prints, arXiv:2003.01140. [ADS link](https://ui.adsabs.harvard.edu/abs/2020arXiv200301140N)




## Credits
If you use tpfplotter, please give credit to the following paper:

Aller, A., Lillo-Box, J., Jones, D., et al. (2020, A&A, 635, 128) "Planetary nebulae seen with TESS: Discovery of new binary central star candidates from Cycle 1,"  [ADS link](https://ui.adsabs.harvard.edu/abs/2020A%26A...635A.128A/abstract)

and add the following sentence somewhere in the paper (either footnote or acknowledgements section):
 > This work made use of \texttt{tpfplotter} by J. Lillo-Box (publicly available in www.github.com/jlillo/tpfplotter), which also made use of the python packages \texttt{astropy}, \texttt{lightkurve}, \texttt{matplotlib} and \texttt{numpy}.

## Contributors

J. Lillo-Box, A. Aller, A. Castro (TESScut addition), D. Jones (python3 version), P. Bluhm (Giai identification), N. Espinoza (Gaia RA-DEC cross-matching), E. Jensen (fix to the orientation arrows and others).

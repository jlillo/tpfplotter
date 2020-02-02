# tpfplotter
 TESS TPF plotter with Gaia catalog

## Purpose
 Create paper-ready figures (1-column) overplotting the Gaia DR2 catalog to the TESS Target Pixel Files (TPF).

## Installation & Requirenments
 Clone this folder or download it to your computer. That's it!

tpfplotter is written in Python2.7 (sorry!), but it will be updated to Python3 as soon as I have some time...

## Usage
Using tpfplotter is really easy:

```
python tpfplotter.py 150428135 --maglim 6
```

Alternatively if you have a list of TIC values just type:

```
python tpfplotter.py list_of_tics.lis --LIST --maglim 6
```

## Papers using tpfplotter
Several papers involving different sceince cases have already used TPF plotter. A sample of them are highlighted here:

- Aller, A., Lillo-Box, J., Jones, D., et al., 2019, arXiv e-prints, arXiv:1911.09991
- Dreizler et al., in prep

## Credits
If you use tpfplotter, you have two options to give credit to it:

- If I (J. Lillo-Box) am co-author of your paper, you can use the following sentence to the acknowledgments section:
```
This work made use of \texttt{tpfplotter}, which also made use of the python packages \texttt{astropy}, \texttt{lightkurve}, \texttt{matplotlib} and \texttt{numpy}.
```

- If I am not co-author of your work that is fine but please add the following entence to the acknowledgments section:
```
This work made use of \texttt{tpfplotter} by J. Lillo-Box (publicly available in www.github.com/jlillobox/tpfplotter), which also made use of the python packages \texttt{astropy}, \texttt{lightkurve}, \texttt{matplotlib} and \texttt{numpy}.
```



import os
import sys
import time
import __future__

from lightkurve import search_targetpixelfile
from lightkurve import search_tesscut
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colorbar import Colorbar
from matplotlib import patches
import numpy as np
import matplotlib.gridspec as gridspec
from bokeh.io import export_png
from bokeh.io.export import get_screenshot_as_png
import warnings
import numpy as np
from astropy.stats import sigma_clip
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u
from astropy.visualization import SqrtStretch,LinearStretch
import astropy.visualization as stretching
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.table import Table, Column, MaskedColumn
from astropy.io import ascii
from astroquery.mast import Catalogs
import argparse

def cli():
    """command line inputs

    Get parameters from command line

    Returns
    -------
    Arguments passed by command line
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("tic", help="TIC number")
    parser.add_argument("-L", "--LIST", help="Only fit the LC", action="store_true")
    parser.add_argument("-S", "--SAVEGAIA", help="Save Gaia sources", action="store_true")
    parser.add_argument("-C", "--COORD", help="Use coordinates", default=False)
    parser.add_argument("--maglim", default=5., help="Maximum magnitude contrast respect to TIC")
    parser.add_argument("--sector", default=None, help="Select Sector if more than one")
    parser.add_argument("--gid", default=None, help="Gaia ID")
    parser.add_argument("--gmag", default=None, help="Gaia mag")
    args = parser.parse_args()
    return args

def add_gaia_figure_elements(tpf, magnitude_limit=18,targ_mag=10.):
    """Make the Gaia Figure Elements"""
    # Get the positions of the Gaia sources
    c1 = SkyCoord(tpf.ra, tpf.dec, frame='icrs', unit='deg')
    # Use pixel scale for query size
    pix_scale = 4.0  # arcseconds / pixel for Kepler, default
    if tpf.mission == 'TESS':
        pix_scale = 21.0
    # We are querying with a diameter as the radius, overfilling by 2x.
    from astroquery.vizier import Vizier
    Vizier.ROW_LIMIT = -1
    result = Vizier.query_region(c1, catalog=["I/345/gaia2"],
                                 radius=Angle(np.max(tpf.shape[1:]) * pix_scale, "arcsec"))
    no_targets_found_message = ValueError('Either no sources were found in the query region '
                                          'or Vizier is unavailable')
    too_few_found_message = ValueError('No sources found brighter than {:0.1f}'.format(magnitude_limit))
    if result is None:
        raise no_targets_found_message
    elif len(result) == 0:
        raise too_few_found_message
    result = result["I/345/gaia2"].to_pandas()
    result = result[result.Gmag < magnitude_limit]
    if len(result) == 0:
        raise no_targets_found_message
    radecs = np.vstack([result['RA_ICRS'], result['DE_ICRS']]).T
    coords = tpf.wcs.all_world2pix(radecs, 0.5) ## TODO, is origin supposed to be zero or one?
    year = ((tpf.astropy_time[0].jd - 2457206.375) * u.day).to(u.year)
    pmra = ((np.nan_to_num(np.asarray(result.pmRA)) * u.milliarcsecond/u.year) * year).to(u.arcsec).value
    pmdec = ((np.nan_to_num(np.asarray(result.pmDE)) * u.milliarcsecond/u.year) * year).to(u.arcsec).value
    result.RA_ICRS += pmra
    result.DE_ICRS += pmdec

    # Gently size the points by their Gaia magnitude
    sizes = 128.0 / 2**(result['Gmag']/targ_mag)#64.0 / 2**(result['Gmag']/5.0)
    one_over_parallax = 1.0 / (result['Plx']/1000.)
    r = (coords[:, 0]+tpf.column,coords[:, 1]+tpf.row,result['Gmag'])

    return r,result

# Plot orientation
def plot_orientation(tpf):
	"""
    Plot the orientation arrows

    Returns
    -------
    tpf read from lightkurve

	"""
	mean_tpf = np.mean(tpf.flux,axis=0)
	nx,ny = np.shape(mean_tpf)
	x0,y0 = tpf.column+int(0.9*nx),tpf.row+int(0.2*nx)
	# East
	tmp =  tpf.get_coordinates()
	ra00, dec00 = tmp[0][0][0][0], tmp[1][0][0][0]
	ra10,dec10 = tmp[0][0][0][-1], tmp[1][0][0][-1]
	theta = np.arctan((dec10-dec00)/(ra10-ra00))
	if (ra10-ra00) < 0.0: theta += np.pi
	#theta = -22.*np.pi/180.
	x1, y1 = 1.*np.cos(theta), 1.*np.sin(theta)
	plt.arrow(x0,y0,x1,y1,head_width=0.2,color='white')
	plt.text(x0+1.5*x1,y0+1.5*y1,'E',color='white')
	# North
	theta = theta +90.*np.pi/180.
	x1, y1 = 1.*np.cos(theta), 1.*np.sin(theta)
	plt.arrow(x0,y0,x1,y1,head_width=0.2,color='white')
	plt.text(x0+1.5*x1,y0+1.5*y1,'N',color='white')



def get_gaia_data(ra, dec):
    """
    Get Gaia parameters

    Returns
    -------
    RA, DEC
    """
    # Get the positions of the Gaia sources
    c1 = SkyCoord(ra, dec, frame='icrs', unit='deg')
    # We are querying with a diameter as the radius, overfilling by 2x.
    from astroquery.vizier import Vizier
    Vizier.ROW_LIMIT = -1
    result = Vizier.query_region(c1, catalog=["I/345/gaia2"],
                                 radius=Angle(10., "arcsec"))
    try:
    	result = result["I/345/gaia2"]
    except:
    	print('Not in Gaia DR2. If you know the Gaia ID and Gmag, try the options --gid and --gmag.')
    	print('Exiting without finishing...')
    	sys.exit()

    no_targets_found_message = ValueError('Either no sources were found in the query region '
                                          'or Vizier is unavailable')
    too_few_found_message = ValueError('No sources found closer than 1 arcsec to TPF coordinates')
    if result is None:
        raise no_targets_found_message
    elif len(result) == 0:
        raise too_few_found_message

    if len(result)>1:
        dist = np.sqrt((result['RA_ICRS']-ra)**2 + (result['DE_ICRS']-dec)**2)
        idx = np.where(dist == np.min(dist))[0][0]
        return result[idx]['Source'], result[idx]['Gmag']
    else:
        return result[0]['Source'], result[0]['Gmag']

def get_gaia_data_from_tic(tic):
    '''
    Get Gaia parameters

    Returns
    -----------------------
    GaiaID, Gaia_mag
    '''
    # Get the Gaia sources
    result = Catalogs.query_object('TIC'+tic, radius=.005, catalog="TIC")
    IDs = result['ID'].data.data
    k = np.where(IDs == tic)[0][0]
    GAIAs = result['GAIA'].data.data
    Gaiamags = result['GAIAmag'].data.data

    GAIA_k = GAIAs[k]
    Gaiamag_k = Gaiamags[k]

    if GAIA_k == '':
        GAIA_k = np.nan
    return GAIA_k, Gaiamag_k


def get_coord(tic):
	"""
	Get TIC corrdinates

	Returns
	-------
	TIC number
	"""
	try:
		catalog_data = Catalogs.query_object(objectname="TIC"+tic, catalog="TIC")
		ra = catalog_data[0]["ra"]
		dec = catalog_data[0]["dec"]
		return ra, dec
	except:
		print("ERROR: No gaia ID found for this TIC")


# ======================================
# 	        MAIN
# ======================================

if __name__ == "__main__":
    args = cli()
    if args.LIST:
    	if args.COORD is not False:
    		_tics = np.genfromtxt(args.tic,dtype=None)
    		tics, ras, decs = [], [], []
    		for t in _tics:
    			tics.append(str(t[0]))
    			ras.append(str(t[1]))
    			decs.append(str(t[2]))
    		ras, decs = np.array(ras), np.array(decs)
    	else:
    		_tics = np.genfromtxt(args.tic,dtype=None)
    		tics = []
    		for t in _tics: tics.append(str(t))
    else:
    	if args.COORD:
    		coords = args.COORD
    		ras, decs = np.array([coords.split(',')[0]]), np.array([coords.split(',')[1]])
    		tics = np.array([args.tic])
    	else:
    		tics = np.array([args.tic])


    for tt,tic in enumerate(tics):

        if args.COORD  is not False:
        	ra,dec = ras[tt], decs[tt]
        	print('Working on '+tic+' (ra = '+ra+', '+'dec = '+dec+') ...')
        else:
        	ra,dec = get_coord(tic)
        	print('Working on TIC'+tic+' (ra = '+str(ra)+', '+'dec = '+str(dec)+') ...')

        if args.gid != None:
        	gaia_id, mag = args.gid, np.float(args.gmag)
        else:
        	#gaia_id, mag = get_gaia_data(ra, dec)
            gaia_id, mag = get_gaia_data_from_tic(tic)



        # By coordinates -----------------------------------------------------------------
        if args.COORD  is not False:
        	                                                                             #
        	if args.sector != None:
        		tpf = search_tesscut(ra+" "+dec, sector=int(args.sector)).download(cutout_size=(12,12))     #

        	else:
        		tpf = search_tesscut(ra+" "+dec).download(cutout_size=(12,12))                             #
        	pipeline = "False"
        	print('    --> Using TESScut to get the TPF')

        # By TIC name --------------------------------------------------------------------
        else:
        	# If the target is in the CTL (short-cadance targets)...
        	try:
        		if args.sector != None:
        			tpf = search_targetpixelfile("TIC "+tic, sector=int(args.sector), mission='TESS').download()
        			a = tpf.flux        # To check it has the flux array
        			pipeline = "True"
        		else:
        			tpf = search_targetpixelfile("TIC "+tic, mission='TESS').download()
        			a = tpf.flux        # To check it has the flux array
        			pipeline = "True"

        		print("    --> Target found in the CTL!")

        	# ... otherwise if it still has a TIC number:
        	except:
        		if args.sector != None:
        			tpf = search_tesscut("TIC "+tic, sector=int(args.sector)).download(cutout_size=(12,12))
        		else:
        			tpf = search_tesscut("TIC "+tic).download(cutout_size=(12,12))
        		print("    -->  Target not in CTL. The FFI cut out was succesfully downloaded")
        		pipeline = "False"

        fig = plt.figure(figsize=(6.93, 5.5))
        gs = gridspec.GridSpec(1,3, height_ratios=[1], width_ratios=[1,0.05,0.01])
        gs.update(left=0.05, right=0.95, bottom=0.12, top=0.95, wspace=0.01, hspace=0.03)
        ax1 = plt.subplot(gs[0,0])

        # TPF plot
        mean_tpf = np.mean(tpf.flux,axis=0)
        nx,ny = np.shape(mean_tpf)
        norm = ImageNormalize(stretch=stretching.LogStretch())
        division = np.int(np.log10(np.nanmax(np.nanmean(tpf.flux,axis=0))))
        splot = plt.imshow(np.nanmean(tpf.flux,axis=0)/10**division,norm=norm, \
        				extent=[tpf.column,tpf.column+ny,tpf.row,tpf.row+nx],origin='bottom', zorder=0)

        # Pipeline aperture
        if pipeline == "True":                                           #
        	aperture_mask = tpf.pipeline_mask
        	aperture = tpf._parse_aperture_mask(aperture_mask)
        	maskcolor = 'tomato'
        	print("    --> Using pipeline aperture...")
        else:
        	aperture_mask = tpf.create_threshold_mask(threshold=10,reference_pixel='center')
        	aperture = tpf._parse_aperture_mask(aperture_mask)
        	maskcolor = 'lightgray'
        	print("    --> Using threshold aperture...")


        for i in range(aperture.shape[0]):
        	for j in range(aperture.shape[1]):
        		if aperture_mask[i, j]:
        			ax1.add_patch(patches.Rectangle((j+tpf.column, i+tpf.row),
        										   1, 1, color=maskcolor, fill=True,alpha=0.4))
        			ax1.add_patch(patches.Rectangle((j+tpf.column, i+tpf.row),
        										   1, 1, color=maskcolor, fill=False,alpha=1,lw=2))

        # Gaia sources
        r, res = add_gaia_figure_elements(tpf,magnitude_limit=mag+np.float(args.maglim),targ_mag=mag)
        x,y,gaiamags = r
        x, y, gaiamags=np.array(x)+0.5, np.array(y)+0.5, np.array(gaiamags)
        size = 128.0 / 2**((gaiamags-mag))
        plt.scatter(x,y,s=size,c='red',alpha=0.6, edgecolor=None,zorder = 10)

        # Gaia source for the target
        this = np.where(np.array(res['Source']) == int(gaia_id))[0]
        plt.scatter(x[this],y[this],marker='x',c='white',s=32,zorder = 11)

        # Legend
        add = 0
        if np.int(args.maglim) % 2 != 0:
            add = 1
        maxmag = np.int(args.maglim) + add
        legend_mags = np.linspace(-2,maxmag,np.int((maxmag+2)/2+1))
        fake_sizes = mag + legend_mags #np.array([mag-2,mag,mag+2,mag+5, mag+8])
        for f in fake_sizes:
            size = 128.0 / 2**((f-mag))
            plt.scatter(0,0,s=size,c='red',alpha=0.6, edgecolor=None,zorder = 10,label = r'$\Delta m=$ '+str(np.int(f-mag)))

        ax1.legend(fancybox=True, framealpha=0.7)

        # Source labels
        dist = np.sqrt((x-x[this])**2+(y-y[this])**2)
        dsort = np.argsort(dist)
        for d,elem in enumerate(dsort):
        	if dist[elem] < 6:
        		plt.text(x[elem]+0.1,y[elem]+0.1,str(d+1),color='white', zorder=100)

        # Orientation arrows
        plot_orientation(tpf)

        # Labels and titles
        plt.xlim(tpf.column,tpf.column+ny)
        plt.ylim(tpf.row,tpf.row+nx)
        plt.xlabel('Pixel Column Number', fontsize=16)
        plt.ylabel('Pixel Row Number', fontsize=16)
        if args.COORD is not False:                                                                                          #
        	plt.title('Coordinates '+tic+' - Sector '+str(tpf.sector), fontsize=16)# + ' - Camera '+str(tpf.camera))  #
        else:   												#
        	plt.title('TIC '+tic+' - Sector '+str(tpf.sector), fontsize=16)# + ' - Camera '+str(tpf.camera))

        # Colorbar
        cbax = plt.subplot(gs[0,1]) # Place it where it should be.
        pos1 = cbax.get_position() # get the original position
        pos2 = [pos1.x0 - 0.05, pos1.y0 ,  pos1.width, pos1.height]
        cbax.set_position(pos2) # set a new position

        cb = Colorbar(ax = cbax, mappable = splot, orientation = 'vertical', ticklocation = 'right')
        plt.xticks(fontsize=14)
        exponent = r'$\times 10^'+str(division)+'$'
        cb.set_label(r'Flux '+exponent+r' (e$^-$)', labelpad=10, fontsize=16)

        plt.savefig('TPF_Gaia_TIC'+tic+'_S'+str(tpf.sector)+'.pdf')

        # Save Gaia sources info
        if args.SAVEGAIA:
        	dist = np.sqrt((x-x[this])**2+(y-y[this])**2)
        	GaiaID = np.array(res['Source'])
        	srt = np.argsort(dist)
        	x, y, gaiamags, dist, GaiaID = x[srt], y[srt], gaiamags[srt], dist[srt], GaiaID[srt]

        	IDs = np.arange(len(x))+1
        	inside = np.zeros(len(x))

        	for i in range(aperture.shape[0]):
        		for j in range(aperture.shape[1]):
        			if aperture_mask[i, j]:
        				xtpf, ytpf = j+tpf.column, i+tpf.row
        				_inside = np.where((x > xtpf) & (x < xtpf+1) &
        						 		   (y > ytpf) & (y < ytpf+1))[0]
        				inside[_inside] = 1



        	data = Table([IDs, GaiaID, x, y, dist, dist*21., gaiamags, inside.astype('int')],
        				names=['# ID','GaiaID','x', 'y','Dist_pix','Dist_arcsec','Gmag', 'InAper'])
        	ascii.write(data, 'Gaia_TIC'+tic+'_S'+str(tpf.sector)+'.dat',overwrite=True)

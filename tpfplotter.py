from __future__ import print_function
import __future__
import os
import sys
import time
import warnings

import numpy as np
import argparse

from lightkurve import search_targetpixelfile
from lightkurve import search_tesscut
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colorbar import Colorbar
from matplotlib import patches
import matplotlib.gridspec as gridspec
from bokeh.io import export_png
from bokeh.io.export import get_screenshot_as_png

from astropy.stats import sigma_clip
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u
from astropy.visualization import SqrtStretch,LinearStretch
import astropy.visualization as stretching
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.table import Table, Column, MaskedColumn
from astropy.io import ascii
from astroquery.mast import Catalogs
from astroquery.simbad import Simbad
Simbad.add_votable_fields('pmra', 'pmdec')
from astroquery.gaia import Gaia

import warnings
warnings.filterwarnings('ignore')

# from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Tahoma'],'size':16})
# rc('text', usetex=False)


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
    parser.add_argument("-n", "--name", help="Target name to be plotted in title", default=False)
    parser.add_argument("-D2", "--DR2", help="Use Gaia DR2 catalog instead of DR3", action="store_true")
    parser.add_argument("--maglim", default=5., help="Maximum magnitude contrast respect to TIC")
    parser.add_argument("--sector", default=None, help="Select Sector if more than one")
    parser.add_argument("--gid", default=None, help="Gaia ID")
    parser.add_argument("--gmag", default=None, help="Gaia mag")
    parser.add_argument("--sradius", default=10., type=float, help="Search radius (in arcsec) for the get_gaia_data function")
    parser.add_argument("--legend", default='best', help="Legend location")
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

    if args.DR2:
        gaia_cat, catID = "I/345/gaia2", "DR2"
        print('\t --> Using Gaia DR2 as requested by user...')
    else:
        gaia_cat, catID = "I/355/gaiadr3", "DR3"

    result = Vizier.query_region(c1, catalog=[gaia_cat],
                                 radius=Angle(np.max(tpf.shape[1:]) * pix_scale, "arcsec"))
    no_targets_found_message = ValueError('Either no sources were found in the query region '
                                          'or Vizier is unavailable')
    too_few_found_message = ValueError('No sources found brighter than {:0.1f}'.format(magnitude_limit))
    if result is None:
        raise no_targets_found_message
    elif len(result) == 0:
        raise too_few_found_message
    result = result[gaia_cat].to_pandas()
    result = result[result.Gmag < magnitude_limit]
    if len(result) == 0:
        raise no_targets_found_message
    year = ((tpf.time[0].jd - 2457206.375) * u.day).to(u.year)
    pmra = ((np.nan_to_num(np.asarray(result.pmRA)) * u.milliarcsecond/u.year) * year).to(u.degree).value
    pmdec = ((np.nan_to_num(np.asarray(result.pmDE)) * u.milliarcsecond/u.year) * year).to(u.degree).value
    result.RA_ICRS += pmra
    result.DE_ICRS += pmdec
    radecs = np.vstack([result['RA_ICRS'], result['DE_ICRS']]).T
    coords = tpf.wcs.all_world2pix(radecs, 0.5) ## TODO, is origin supposed to be zero or one?

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
	x0,y0 = tpf.column+int(0.2*nx)+0.5,tpf.row+int(0.2*ny)+0.5
	# East
	tmp =  tpf.get_coordinates()
	ra00, dec00 = tmp[0][0][0][0], tmp[1][0][0][0]
	ra10,dec10 = tmp[0][0][0][-1], tmp[1][0][0][-1]
    # Each degree of RA is not a full degree on the sky if not
    # at equator; need cos(dec) factor to compensate
	cosdec = np.cos(np.deg2rad(0.5*(dec10+dec00)))
    # Reverse the order of RA arguments here relative to dec
    # args to account for handedness of RA/Dec vs. x/y coords:
	theta = np.arctan((dec10-dec00)/(cosdec*(ra00-ra10)))
	if (ra10-ra00) < 0.0: theta += np.pi
	#theta = -22.*np.pi/180.
    # If angle is small, arrows can be a bit closer to corner:
	if (abs(np.rad2deg(theta)) < 30):
		x0 -= 0.08*nx
		y0 -= 0.08*ny
	x1, y1 = 1.*np.cos(theta), 1.*np.sin(theta)
	plt.arrow(x0,y0,x1,y1,head_width=0.2,color='white')
	plt.text(x0+1.6*x1,y0+1.6*y1,'E',color='white',ha='center',va='center')
	# North
	theta = theta +90.*np.pi/180.
	x1, y1 = 1.*np.cos(theta), 1.*np.sin(theta)
	plt.arrow(x0,y0,x1,y1,head_width=0.2,color='white')
	plt.text(x0+1.6*x1,y0+1.6*y1,'N',color='white',ha='center',va='center')



def get_gaia_data(ra, dec, search_radius=10.):
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
    if args.DR2:
        gaia_cat, catID = "I/345/gaia2", "DR2"
        print('\t --> Using Gaia DR2 as requested by user...')
    else:
        gaia_cat, catID = "I/355/gaiadr3", "DR3"

    result = Vizier.query_region(c1, catalog=[gaia_cat],
                                 radius=Angle(search_radius, "arcsec"))
    try:
    	result = result[gaia_cat]
    except:
        print('Not in Gaia '+catID+'. If you know the Gaia ID and Gmag, try the options --gid and --gmag.')
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

def get_dr2_id_from_tic(tic):
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


def get_gaia_data_from_simbad(dr2ID):
    simb = Simbad.query_object('Gaia DR2 '+dr2ID)
    simbid = Simbad.query_objectids('Gaia DR2 '+dr2ID)
    if simbid == None:
        print("ERROR: TIC not found as Gaia DR2 "+str(dr2ID))
    ids = np.array(simbid['ID'].data).astype(str)
    myid = [id for id in ids if 'DR3' in id]
    if len(myid) == 0:
        myid = [id for id in ids if 'DR2' in id]
    myid = myid[0].split(' ')[2]

    query2 = "SELECT \
             TOP 1 \
             source_id, ra, dec, pmra, pmdec, parallax, phot_g_mean_mag\
             FROM gaiadr3.gaia_source\
             WHERE source_id = "+str(myid)+" \
             "
    job = Gaia.launch_job_async(query2)
    gmag = job.get_results()['phot_g_mean_mag'].data[0]

    return myid,gmag




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
        # print(catalog_data.keys())
        # print(catalog_data[0]["GAIA"])
        return ra, dec
    except:
    	print("ERROR: TIC not found in Simbad")


# ======================================
# 	        MAIN
# ======================================

if __name__ == "__main__":
    args = cli()
    # print("\n")
    print("======================")
    print("     tpfplotter       ")
    print("======================\n")
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
        	print('* Working on '+tic+' (ra = '+ra+', '+'dec = '+dec+') ...')
        else:
        	ra,dec = get_coord(tic)
        	print('* Working on TIC'+tic+' (ra = '+str(ra)+', '+'dec = '+str(dec)+') ...')

        if args.gid != None:
        	gaia_id, mag = args.gid, float(args.gmag)
        else:
            if args.COORD  is not False:
        	       gaia_id, mag = get_gaia_data(ra, dec, search_radius=args.sradius)
            else:
                dr2ID,_ = get_dr2_id_from_tic(tic)
                gaia_id, mag = get_gaia_data_from_simbad(dr2ID)
                if np.isnan(mag):
                    gaia_id, mag = get_gaia_data(ra, dec, search_radius=args.sradius)


        # By coordinates -----------------------------------------------------------------
        if args.COORD  is not False:
        	                                                                             #
        	if args.sector != None:
        		tpf = search_tesscut(ra+" "+dec, sector=int(args.sector)).download(cutout_size=(12,12))     #

        	else:
        		tpf = search_tesscut(ra+" "+dec).download(cutout_size=(12,12))                             #
        	pipeline = "False"
        	print('\t --> Using TESScut to get the TPF')

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

        		print("\t --> Target found in the CTL!")

        	# ... otherwise if it still has a TIC number:
        	except:
        		if args.sector != None:
        			tpf = search_tesscut("TIC "+tic, sector=int(args.sector)).download(cutout_size=(12,12))
        		else:
        			tpf = search_tesscut("TIC "+tic).download(cutout_size=(12,12))
        		print("\t -->  Target not in CTL. The FFI cut out was succesfully downloaded")
        		pipeline = "False"

        fig = plt.figure(figsize=(6.93, 5.5))
        gs = gridspec.GridSpec(1,3, height_ratios=[1], width_ratios=[1,0.05,0.01])
        gs.update(left=0.05, right=0.95, bottom=0.12, top=0.95, wspace=0.01, hspace=0.03)
        ax1 = plt.subplot(gs[0,0])

        # TPF plot
        mean_tpf = np.mean(tpf.flux,axis=0)
        nx,ny = np.shape(mean_tpf)
        norm = ImageNormalize(stretch=stretching.LogStretch())
        division = int(np.log10(np.nanmax(np.nanmean(tpf.flux.value ,axis=0)))) #* u.s/u.electron
        image = np.nanmean(tpf.flux,axis=0)/10**division
        splot = plt.imshow(image,norm=norm, \
        				extent=[tpf.column+0.5,tpf.column+ny+0.5,tpf.row+0.5,tpf.row+nx+0.5],origin='lower', zorder=0)

        # Pipeline aperture
        if pipeline == "True":                                           #
        	aperture_mask = tpf.pipeline_mask
        	aperture = tpf._parse_aperture_mask(aperture_mask)
        	maskcolor = 'tomato'
        	print("\t --> Using pipeline aperture...")
        else:
        	aperture_mask = tpf.create_threshold_mask(threshold=10,reference_pixel='center')
        	aperture = tpf._parse_aperture_mask(aperture_mask)
        	maskcolor = 'lightgray'
        	print("\t --> Using threshold aperture...")


        for i in range(aperture.shape[0]):
        	for j in range(aperture.shape[1]):
        		if aperture_mask[i, j]:
        			ax1.add_patch(patches.Rectangle((j+tpf.column+0.5, i+tpf.row+0.5),
        										   1, 1, color=maskcolor, fill=True,alpha=0.4))
        			ax1.add_patch(patches.Rectangle((j+tpf.column+0.5, i+tpf.row+0.5),
        										   1, 1, color=maskcolor, fill=False,alpha=1,lw=2))

        # Gaia sources
        r, res = add_gaia_figure_elements(tpf,magnitude_limit=mag+float(args.maglim),targ_mag=mag)
        # plt.figure(2)
        # plt.scatter(res.RA_ICRS,res.DE_ICRS)
        # for r,d,s in zip(res.RA_ICRS,res.DE_ICRS,res['Source']): plt.text(r,d,s)
        # plt.plot(tpf.ra,tpf.dec,'s',c='none',markeredgecolor='red')
        # plt.show()
        x,y,gaiamags = r
        x, y, gaiamags=np.array(x)+0.5, np.array(y)+0.5, np.array(gaiamags)
        size = 128.0 / 2**((gaiamags-mag))
        plt.scatter(x+0.5,y+0.5,s=size,c='red',alpha=0.6, edgecolor=None,zorder = 10)

        # Gaia source for the target
        this = np.where(np.array(res['Source']) == int(gaia_id))[0]
        plt.scatter(x[this]+0.5,y[this]+0.5,marker='x',c='white',s=32,zorder = 11)

        # Legend
        add = 0
        if int(args.maglim) % 2 != 0:
            add = 1
        maxmag = int(args.maglim) + add
        legend_mags = np.linspace(-2,maxmag,int((maxmag+2)/2+1))
        fake_sizes = mag + legend_mags #np.array([mag-2,mag,mag+2,mag+5, mag+8])
        for f in fake_sizes:
            size = 128.0 / 2**((f-mag))
            plt.scatter(0,0,s=size,c='red',alpha=0.6, edgecolor=None,
                        zorder = 10,label = r'$\Delta m=$ '+str(int(f-mag)))

        ax1.legend(fancybox=True, framealpha=0.7, loc=args.legend,fontsize=14)

        # Source labels
        dist = np.sqrt((x-x[this])**2+(y-y[this])**2)
        dsort = np.argsort(dist)
        corners = np.array([np.abs(x[this]-(tpf.column+nx)), np.abs(x[this]-tpf.column),
                                   np.abs(y[this]-(tpf.row+ny)), np.abs(y[this]-tpf.row)])
        mindist = np.min(corners)
        xmin = tpf.column + 0.05*nx
        xmax = tpf.column + 0.95*nx
        ymin = tpf.row + 0.05*ny
        ymax = tpf.row + 0.95*ny
        for d,elem in enumerate(dsort):
            if ( (x[elem]+0.5 < xmax) & (x[elem]+0.5 > xmin) & (y[elem]+0.5 < ymax) & (y[elem]+0.5 > ymin)  ):
                plt.text(x[elem]+0.1+0.5,y[elem]+0.1+0.5,str(d+1),color='white', zorder=100,fontsize=14)

        # Orientation arrows
        plot_orientation(tpf)

        # Labels and titles
        # Reverse x limits so that image plots as seen on the sky:
        plt.xlim(tpf.column+ny+0.5,tpf.column+0.5)
        plt.ylim(tpf.row+0.5,tpf.row+nx+0.5)
        plt.xlabel('Pixel Column Number', fontsize=16, zorder=200)
        plt.ylabel('Pixel Row Number', fontsize=16, zorder=200)
        if args.COORD is not False:                                                                                          #
        	plt.title('Coordinates '+tic+' - Sector '+str(tpf.sector), fontsize=16, zorder=200)# + ' - Camera '+str(tpf.camera))  #
        elif args.name is not False:
            plt.title(args.name +' - Sector '+str(tpf.sector), fontsize=16, zorder=200)
        else:   												#
        	plt.title('TIC '+tic+' - Sector '+str(tpf.sector), fontsize=16, zorder=200)# + ' - Camera '+str(tpf.camera))

        # Colorbar
        cbax = plt.subplot(gs[0,1]) # Place it where it should be.
        pos1 = cbax.get_position() # get the original position
        pos2 = [pos1.x0 - 0.05, pos1.y0 ,  pos1.width, pos1.height]
        cbax.set_position(pos2) # set a new position

        cbar_ticks = np.linspace(np.min(image), np.max(image), 8, endpoint=True)

        cb = Colorbar(ax = cbax, mappable = splot, orientation = 'vertical',
                      ticklocation = 'right')
        plt.xticks(fontsize=14)
        #cbax.set_yticklabels(["{:4.2f}".format(i) for i in cbar_ticks])
        exponent = r'$\times 10^'+str(division)+'$'
        cb.set_label(r'Flux '+exponent+r' (e$^-$/s)', labelpad=10, fontsize=16)

        plt.savefig('TPF_Gaia_TIC'+tic+'_S'+str(tpf.sector)+'.pdf')
        print('\t --> TPF plot written in file: '+'TPF_Gaia_TIC'+tic+'_S'+str(tpf.sector)+'.pdf')

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
            print('\t --> Gaia close sources saved in file: '+'Gaia_TIC'+tic+'_S'+str(tpf.sector)+'.dat')
        print("\t --> Done!\n")

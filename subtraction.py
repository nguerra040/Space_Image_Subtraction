# DOCUMENTATION AND EXPLANATION IN WRITE-UP

from wget import download

link = 'https://mast.stsci.edu/api/v0.1/Download/file?uri=https://hla.stsci.edu/cgi-bin/getdata.cgi?dataset=hst_10861_18_acs_wfc_f475w_drz.fits'
fits_name = download(link)

from wget import download
import tarfile, os, glob, gzip

link = 'http://hst.esac.esa.int/ehst-sl-server/servlet/data-action?RETRIEVAL_TYPE=PRODUCT&OBSERVATION_ID=hst_13386_t3_wfc3_ir_f160w'
tar_name = download(link)
f = tarfile.open(tar_name, mode='r|*')
working_dir_path = '/content/fits'  # CHANGE TO WHEREVER YOU WANT THE DATA TO BE EXTRACTED
f.extractall(path=working_dir_path)

os.chdir('/content/fits/HST/') # CHANGE DIR HOLDING FITS FILES 
folder_name = glob.glob('*')[0] # USE GLOB TO GET FOLDER NAME
fits_name = gzip.open(f'{folder_name}/{folder_name}_drz.fits.gz') # GUNZIP FILE

from astroquery.mast import Observations

obs_table = Observations.query_object("a2744",radius=".01 arcsec") # SEARCH
data_products_by_obs = Observations.get_product_list(obs_table) # TABLE WITH LEFT HAND SIDE BEING OBSID
print(data_products_by_obs)

satellites = []
for satellite in data_products_by_obs['obs_collection']:
  if satellites.count(satellite) == 0:
    satellites.append(satellite)
print(satellites) # SEE AVAILABLE TELESCOPES

hla_products = Observations.filter_products(data_products_by_obs,
                                         obs_collection="HLA",
                                         extension="fits") # FILTER TO ONLY HLA FITS FILES

first_hla_products = hla_products[0:5] # DOWNLOAD FIRST 5 FILES
download_hla = Observations.download_products(first_hla_products)

fits_name = download_hla['Local Path'][0] # GET PATH OF FIRST FITS FILE

import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
from wget import download
from astropy.io import fits

link = 'https://mast.stsci.edu/api/v0.1/Download/file?uri=https://hla.stsci.edu/cgi-bin/getdata.cgi?dataset=hst_10861_18_acs_wfc_f475w_drz.fits'
fits_name = download(link)

fitsfile = fits.open(fits_name)

print(fitsfile.info()) # SEE HDUS STORED
hdu = fitsfile['SCI'] # GET HDU LABELED SCI
print(hdu.header)

#PLOT DATA
mean, median, std = sigma_clipped_stats(hdu.data)
plt.imshow(hdu.data, vmin = median - 5*std, vmax = median + 5*std)
plt.colorbar()

import sep
from astropy.io import fits
from wget import download

link = 'https://mast.stsci.edu/api/v0.1/Download/file?uri=https://hla.stsci.edu/cgi-bin/getdata.cgi?dataset=hst_10861_18_acs_wfc_f475w_drz.fits'
fits_name = download(link)
hdu = fits.open(fits_name)['SCI']

data = hdu.data.byteswap().newbyteorder()
bkg = sep.Background(data) # CALCULATE BACKGROUND
data_sub_bkg = data - bkg # SUBTRACT BACKGROUND

import astroalign as aa
...
sci_aligned, footprint = aa.register(sci_data, ref_data)

from photutils.psf import create_matching_kernel
from astropy.convolution import convolve

#IDENTIFY POINT SOURCES IN REF IMAGE
bkg = sep.Background(ref_data)
sources = pd.DataFrame(sep.extract(np.array(ref_data), 500, err = bkg.globalrms)) # GET SOURCES ABOVE 500*SIGMA
final_sources = pd.DataFrame(columns = list(sources.columns))
for index, row in sources.iterrows():
  if row['b']/row['a'] > .5: #GET ROUNDER OBJECTS (1 IS PERFECT CIRCLE)
    final_sources.loc[index] = row

#CONVOLVE THE PSFS TO MATCH
for index, row in final_sources.iterrows():
  ymin, ymax, xmin, xmax = int(row['ymin']), int(row['ymax']), int(row['xmin']), int(row['xmax'])
  # CONVOLVE ONLY ACCEPTS ODD NUMBER AXES
  if (ymax - ymin)%2 == 0:
    ymax += 1
  if (xmax - xmin)%2 == 0:
    xmax += 1
  kernel = create_matching_kernel(sci_aligned[ymin:ymax,xmin:xmax],ref_data[ymin:ymax,xmin:xmax])
  sci_aligned[ymin:ymax,xmin:xmax] = convolve(sci_data[ymin:ymax,xmin:xmax], kernel)
  
  #CROP
xstart, xend, ystart, yend = 2000, 4000, 1200, 3600

w = WCS(ref.header)
dif_cropped = fits.PrimaryHDU()
dif_cropped.data = difference[ystart:yend,xstart:xend]
dif_cropped.header = ref_hdu.header
dif_cropped.header.update(w[ystart:yend,xstart:xend].to_header())

#plot
mean, median, std = sigma_clipped_stats(difference)
plt.figure(2,figsize=(11,11))
plt.subplot(1,2,1)
plt.imshow(difference, vmin = median-5*std, vmax = median+5*std)
plt.colorbar()
plt.axhline(ystart)
plt.axhline(yend)
plt.axvline(xstart)
plt.axvline(xend)
plt.title('Difference Image')

mean, median, std = sigma_clipped_stats(dif_cropped.data)
plt.subplot(1,2,2)
plt.imshow(dif_cropped.data, vmin = median-5*std, vmax = median+5*std)
plt.colorbar()
plt.title('Cropped Difference Image')
plt.show()

dbkg = sep.Background(np.array(dif_cropped.data))
objects = pd.DataFrame(sep.extract(np.array(dif_cropped.data), threshold, minarea=minarea, err = dbkg.globalrms, deblend_nthresh=deblend)) #detect

#filter out non-round objects
final_objects = pd.DataFrame(columns = list(objects.columns))
for index, row in objects.iterrows():
  if row['b']/row['a'] > ab_fraction_min:
    final_objects.loc[index] = row
final_objects.reset_index(inplace=True, drop=True)

w = WCS(dif_cropped.header) # convert from pixel to world coords
final_objects['RA (deg)']=w.pixel_to_world(final_objects['x'],final_objects['y']).ra.deg
final_objects['Dec (deg)']=w.pixel_to_world(final_objects['x'],final_objects['y']).dec.deg

fig, ax = plt.subplots(figsize=(20,20))
mean, median, std = sigma_clipped_stats(dif_cropped.data)
im = ax.imshow(dif_cropped.data, vmin = median-5*std, vmax = median+5*std)

#factor to multiply ellipses to see better
multx = int(dif_cropped.data.shape[0]*.025)
multy = int(dif_cropped.data.shape[1]*.025)
for i in range(len(final_objects)):
  e = Ellipse(xy=(final_objects['x'][i], final_objects['y'][i]),
              width=multx*final_objects['a'][i],
              height=multy*final_objects['b'][i],
              angle=final_objects['theta'][i] * 180. / np.pi)
  e.set_facecolor('none')
  e.set_edgecolor('red')
  ax.add_artist(e)
plt.show()

ra=3.574458
dec=-30.399917
c = SkyCoord(ra,dec,unit='deg')
x,y = w.world_to_pixel(c)

fits.writeto('difference.fits',dif_cropped.data,dif_cropped.header)
final_objects.to_csv('src.csv')

f = open("src.reg", "a")
f.write("# Region file format: DS9 version 4.1\n")
f.write("global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
f.write("fk5\n")
for index, row in final_objects.iterrows():
  f.write(f"circle({row['RA (deg)']},{row['Dec (deg)']},1.00\")\n")
f.close()

ra = 3.6047835350
dec = -30.4149532318
size = 1 #INITIAL CROP TO FIND CENTER OF SOURCE
radius = 2 #5x5 BOX (RADIUS=3 IS 7x7 BOX, RADIUS=4 IS 9x9 BOX)
ee = .825

#CONVERT FROM WORLD TO PIXEL
w = WCS(hdu.header)
c = SkyCoord(ra = ra, dec = dec, unit='deg')
x,y = w.world_to_pixel(c)
x,y = int(x),int(y)

time = (hdu.header['EXPSTART'] + hdu.header['EXPEND'])/2
print("MJD:", time)

#SUBTRACT BACKGROUND
bkg = sep.Background(hdu.data.byteswap().newbyteorder())
data = hdu.data.byteswap().newbyteorder() - bkg

#CROP
cropped_data = data[y-size:y+size,x-size:x+size]

#CENTER ON OBJECT
y_max, x_max = np.where(cropped_data == cropped_data.max())
x_max, y_max = int(x_max), int(y_max)
better_x, better_y = x_max+(x-size), y_max+(y-size)
object_center = pd.DataFrame({'x':[better_x],'y':[better_y]})

box = data[better_y-radius:better_y+radius+1,better_x-radius:better_x+radius+1]

cr = np.sum(box)
flux = (cr/ee)*hdu.header['PHOTFLAM']
print("Flux:",flux)

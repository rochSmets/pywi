

#------------------------------
# 
    
# -----------------------------
def main():
	import picgsfc as pr
	import numpy as np
	import matplotlib.pyplot as plt
	import matplotlib.patches as mpatches


	runby10 = pr.PICGSFCrun('by10', '/Volumes/drobo/nico/asymmetric/by10/data/')

	data   = runby10.GetN(50.0 ,specie='electrons',silent='no')

#	
#	fig = plt.figure(1,facecolor='none', linewidth=0,frameon=False, figsize=(10,10))    
#	#plt.subplots_adjust(wspace=0.5)
#	ax1=fig.add_subplot(211)
#	ax2=fig.add_subplot(212)
#	
#	
#	
#	plt.subplot(111)

	datac = data[:,:]

	plt.imshow(np.transpose(datac), # data to show \
		cmap='jet',                     # color map\
		origin='lower',                 # where is the zero\
		interpolation='nearest',        # pixels \
		extent=[0,64,0,25.6],           # range for x and y [xmin,ymin,ymax,zmax] (physical units.)\
		aspect='equal',                 # auto = aspect ratio given by 'extent'
		vmin= np.min(datac),vmax=np.max(datac))
#	
#	ax1.set_title(r'$J_z$', fontsize=25)
#	ax1.set_xlabel(r'$x/\delta_i$', fontsize=20)
#	ax1.set_ylabel(r'$y/\delta_i$', fontsize=20)
	
	plt.colorbar()
	plt.show()	
	
	
# ----------------------------






if __name__ == "__main__":
	main() 

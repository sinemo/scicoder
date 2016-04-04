import numpy as np
import pyfits
from astropy.io import fits

fname = '/home/sinemo/Desktop/fitstable/modified_last_version2.csv'
all_hi_edd = '/home/sinemo/Desktop/fitstable/all_hi_edd_catalogue.csv'

numbering_with_duplications		= np.genfromtxt(fname, dtype=int, delimiter=',', usecols=(0), skip_header=1, missing_values=None)
numbering_excluding_duplications	= np.genfromtxt(fname, dtype=int, delimiter=',', usecols=(1), skip_header=1, missing_values=None)
HIPASS_NAME 				= np.genfromtxt(fname, dtype=None, delimiter=',', usecols=(2), skip_header=1, missing_values=None)
OBS_RUN_NIGHT 				= np.genfromtxt(fname, dtype=None, delimiter=',', usecols=(3), skip_header=1, missing_values=None)
OBS_TIME 				= np.genfromtxt(fname, dtype=None, delimiter=',', usecols=(4), skip_header=1, missing_values=None)
ALT_NAME 				= np.genfromtxt(fname, dtype=None, delimiter=',', usecols=(5), skip_header=1, missing_values=None)
OBSERVATION_COMMENTS 			= np.genfromtxt(fname, dtype=None, delimiter=',', usecols=(6), skip_header=1, missing_values=None)
NOD_SHUFFLE 				= np.genfromtxt(fname, dtype=None, delimiter=',', usecols=(7), skip_header=1, missing_values=None)
NED_RA 					= np.genfromtxt(fname, dtype=None, delimiter=',', usecols=(8), skip_header=1, missing_values=None) 
NED_DEC 				= np.genfromtxt(fname, dtype=None, delimiter=',', usecols=(9), skip_header=1, missing_values=None)
USED_ORNOT				= np.genfromtxt(fname, dtype=None, delimiter=',', usecols=(10), skip_header=1, missing_values=None)
COMMENT_ON_PPXF				= np.genfromtxt(fname, dtype=None, delimiter=',', usecols=(11), skip_header=1, missing_values=None)
NED_BASIC_DATA_MORPHOLOGICAL_TYPE	= np.genfromtxt(fname, dtype=None, delimiter=',', usecols=(12), skip_header=1, missing_values=None)
MORPHOLOGY_SHORTENED			= np.genfromtxt(fname, dtype=None, delimiter=',', usecols=(13), skip_header=1, missing_values=None)
NED_REDSHIFT				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(14), skip_header=1, missing_values=None)
NED_REDSHIFT_UNCERTAINTY		= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(15), skip_header=1, missing_values=None)
NED_BASIC_DATA_MAJOR_DIAMETER		= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(16), skip_header=1, missing_values=None)
NED_BASIC_DATA_MINOR_DIAMETER		= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(17), skip_header=1, missing_values=None)
NED_BASIC_MAGNITUDE			= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(18), skip_header=1, missing_values=None)
NED_MAJOR_DIAMETER_ESO_Uppsala		= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(19), skip_header=1, missing_values=None)
NED_MINOR_MAJOR_RATIO_ESO_Uppsala	= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(20), skip_header=1, missing_values=None)
NED_POSITION_ANGLE_ESO_Uppsala		= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(21), skip_header=1, missing_values=None)
NED_MAJOR_DIAMETER_RC3			= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(22), skip_header=1, missing_values=None)
NED_MAJOR_DIAMETER_UNC_RC3		= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(23), skip_header=1, missing_values=None)
NED_MINOR_MAJOR_RATIO_RC3		= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(24), skip_header=1, missing_values=None)
NED_MINOR_MAJOR_RATIO_UNC_RC3		= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(25), skip_header=1, missing_values=None)
NED_POSITION_ANGLE_RC3			= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(26), skip_header=1, missing_values=None)
HIPASS_RV50max				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(27), skip_header=1, missing_values=None)
HIPASS_RV50min				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(28), skip_header=1, missing_values=None)
HIPASS_RV20max				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(29), skip_header=1, missing_values=None)
HIPASS_RV20min				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(30), skip_header=1, missing_values=None)
HIPASS_RVmom				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(31), skip_header=1, missing_values=None)
HIPASS_RVsp				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(32), skip_header=1, missing_values=None)
HIPASS_RVGSR				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(33), skip_header=1, missing_values=None)
HIPASS_RVLG				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(34), skip_header=1, missing_values=None)
HIPASS_RVCMB				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(35), skip_header=1, missing_values=None)
HIPASS_RV1				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(36), skip_header=1, missing_values=None)
HIPASS_RV2				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(37), skip_header=1, missing_values=None)
HIPASS_RVsp1				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(38), skip_header=1, missing_values=None)
HIPASS_RVsp2				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(39), skip_header=1, missing_values=None)
HIPASS_RVmask				= np.genfromtxt(fname, dtype=None, delimiter=',', usecols=(40), skip_header=1, missing_values=None)
HIPASS_W50max				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(41), skip_header=1, missing_values=None)
HIPASS_W50min				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(42), skip_header=1, missing_values=None)
HIPASS_W20max				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(43), skip_header=1, missing_values=None)
HIPASS_W20min				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(44), skip_header=1, missing_values=None)
HIPASS_Speak				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(45), skip_header=1, missing_values=None)
HIPASS_Sint				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(46), skip_header=1, missing_values=None)
HIPASS_RMS				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(47), skip_header=1, missing_values=None)
HIPASS_RMSclip				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(48), skip_header=1, missing_values=None)
HIPASS_RMScube				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(49), skip_header=1, missing_values=None)
HIPASS_cube				= np.genfromtxt(fname, dtype=int, delimiter=',', usecols=(50), skip_header=1, missing_values=None)
HIPASS_Sigma				= np.genfromtxt(fname, dtype=int, delimiter=',', usecols=(51), skip_header=1, missing_values=None)
HIPASS_Boxsize				= np.genfromtxt(fname, dtype=int, delimiter=',', usecols=(52), skip_header=1, missing_values=None)
HIPASS_Qual				= np.genfromtxt(fname, dtype=int, delimiter=',', usecols=(53), skip_header=1, missing_values=None)
HIPASS_nb				= np.genfromtxt(fname, dtype=int, delimiter=',', usecols=(54), skip_header=1, missing_values=None)
HIPASS_cf				= np.genfromtxt(fname, dtype=int, delimiter=',', usecols=(55), skip_header=1, missing_values=None)
HIPASS_ext				= np.genfromtxt(fname, dtype=int, delimiter=',', usecols=(56), skip_header=1, missing_values=None)
RC3_PGC					= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(57), skip_header=1, missing_values=None)
RC3_D25					= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(58), skip_header=1, missing_values=None)
RC3_e_D25				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(59), skip_header=1, missing_values=None)
RC3_R25					= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(60), skip_header=1, missing_values=None)
RC3_e_R25				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(61), skip_header=1, missing_values=None)
RC3_Ae					= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(62), skip_header=1, missing_values=None)
RC3_e_Ae				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(63), skip_header=1, missing_values=None)
RC3_A_e_in_arcsec			= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(64), skip_header=1, missing_values=None)
RC3_PA					= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(65), skip_header=1, missing_values=None)
RC3_BT					= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(66), skip_header=1, missing_values=None)
RC3_BT_code				= np.genfromtxt(fname, dtype=None, delimiter=',', usecols=(67), skip_header=1, missing_values=None)
RC3_e_BT				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(68), skip_header=1, missing_values=None)
RC3_Bmag				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(69), skip_header=1, missing_values=None)
RC3_e_Bmag				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(70), skip_header=1, missing_values=None)
RC3_BoT					= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(71), skip_header=1, missing_values=None)
RC3_Ai					= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(72), skip_header=1, missing_values=None)
RC3_Ag					= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(73), skip_header=1, missing_values=None)
RC3_W50					= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(74), skip_header=1, missing_values=None)
RC3_e_W50				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(75), skip_header=1, missing_values=None)
RC3_cz					= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(76), skip_header=1, missing_values=None)
TWOMASX_ID				= np.genfromtxt(fname, dtype=None, delimiter=',', usecols=(77), skip_header=1, missing_values=None)
TWOMASX_K_b_over_a			= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(78), skip_header=1, missing_values=None)
TWOMASX_K_PA				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(79), skip_header=1, missing_values=None)
TWOMASX_K_mag				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(80), skip_header=1, missing_values=None)
TWOMASX_K_mag_sigma			= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(81), skip_header=1, missing_values=None)
TWOMASX_K_mag_chi2			= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(82), skip_header=1, missing_values=None)
TWOMRS_E_B_V				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(83), skip_header=1, missing_values=None)
ABSOLUTE_MAGNITUDE_B			= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(84), skip_header=1, missing_values=None)
ABSOLUTE_MAGNITUDE_K			= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(85), skip_header=1, missing_values=None)
cos2_i					= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(86), skip_header=1, missing_values=None)
sin_i					= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(87), skip_header=1, missing_values=None)
INCLINATION_FROM_R25			= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(88), skip_header=1, missing_values=None)
HIPASS_W50max_z				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(89), skip_header=1, missing_values=None)
HIPASS_W50max_zs			= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(90), skip_header=1, missing_values=None)
HIPASS_W50max_zst			= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(91), skip_header=1, missing_values=None)
HIPASS_W50max_zsti			= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(92), skip_header=1, missing_values=None)
BILBI_CENTRAL_VELOCITY_DISPERSION	= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(93), skip_header=1, missing_values=None)
BILBI_ERR_CENTRAL_SIGMA			= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(94), skip_header=1, missing_values=None)
BILBI_ROTATIONAL_VELOCITY		= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(95), skip_header=1, missing_values=None)
BILBI_ERR_ROTATIONAL_VELOCITY		= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(96), skip_header=1, missing_values=None)
BILBI_POSITION_ANGLE			= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(97), skip_header=1, missing_values=None)
BILBI_STD_DEV_POS_ANG			= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(98), skip_header=1, missing_values=None)
BILBI_COMMENT_POS_ANG			= np.genfromtxt(fname, dtype=None, delimiter=',', usecols=(99), skip_header=1, missing_values=None)
BILBI_Wmax_z				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(100), skip_header=1, missing_values=None)
BILBI_Wmax_zs				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(101), skip_header=1, missing_values=None)
BILBI_Wmax_zst				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(102), skip_header=1, missing_values=None)
BILBI_Wmax_zsti				= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(103), skip_header=1, missing_values=None)
IN_THE_TABLE_OR_NOT			= np.genfromtxt(fname, dtype=float, delimiter=',', usecols=(104), skip_header=1, missing_values=None)


PGC					= np.genfromtxt(all_hi_edd, dtype=None, delimiter=',', usecols=(0), skip_header=6, missing_values=None)
Name_Profile				= np.genfromtxt(all_hi_edd, dtype=None, delimiter=',', usecols=(1), skip_header=6, missing_values=None)
Vh_av					= np.genfromtxt(all_hi_edd, dtype=float, delimiter=',', usecols=(2), skip_header=6, missing_values=None)
Wmx_av					= np.genfromtxt(all_hi_edd, dtype=float, delimiter=',', usecols=(3), skip_header=6, missing_values=None)
eW_av					= np.genfromtxt(all_hi_edd, dtype=float, delimiter=',', usecols=(4), skip_header=6, missing_values=None)
N_av					= np.genfromtxt(all_hi_edd, dtype=float, delimiter=',', usecols=(5), skip_header=6, missing_values=None)
Source1					= np.genfromtxt(all_hi_edd, dtype=None, delimiter=',', usecols=(6), skip_header=6, missing_values=None)
Tel1					= np.genfromtxt(all_hi_edd, dtype=None, delimiter=',', usecols=(7), skip_header=6, missing_values=None)
Vhel1					= np.genfromtxt(all_hi_edd, dtype=float, delimiter=',', usecols=(8), skip_header=6, missing_values=None)
Wm501					= np.genfromtxt(all_hi_edd, dtype=float, delimiter=',', usecols=(9), skip_header=6, missing_values=None)
Wcm501					= np.genfromtxt(all_hi_edd, dtype=float, delimiter=',', usecols=(10), skip_header=6, missing_values=None)
Wmx1					= np.genfromtxt(all_hi_edd, dtype=float, delimiter=',', usecols=(11), skip_header=6, missing_values=None)
e_W1					= np.genfromtxt(all_hi_edd, dtype=float, delimiter=',', usecols=(12), skip_header=6, missing_values=None)
SN1					= np.genfromtxt(all_hi_edd, dtype=float, delimiter=',', usecols=(13), skip_header=6, missing_values=None)
Flux1					= np.genfromtxt(all_hi_edd, dtype=float, delimiter=',', usecols=(14), skip_header=6, missing_values=None)
Res1					= np.genfromtxt(all_hi_edd, dtype=float, delimiter=',', usecols=(15), skip_header=6, missing_values=None)
Ns1					= np.genfromtxt(all_hi_edd, dtype=float, delimiter=',', usecols=(16), skip_header=6, missing_values=None)
Fm501					= np.genfromtxt(all_hi_edd, dtype=float, delimiter=',', usecols=(17), skip_header=6, missing_values=None)



#########################SECOND CALCULATION FOR HIPASS WIDTH CORRECTION ############################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!USE THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c=300000. #in km/s
z_HIPASS = HIPASS_RV50max/c

Wz_HIPASS = HIPASS_W50max / (1+z_HIPASS)
Wzs_HIPASS = Wz_HIPASS - 0.13*26.4

Wt_HIPASS = 5.
Wc_HIPASS = 100.

for i in range(len(HIPASS_NAME)):

	if sin_i[i] == 0:
		sin_i[i]=None
#		print HIPASS_NAME[i]

Wzst_HIPASS = np.sqrt(Wzs_HIPASS**2 + ((Wt_HIPASS**2)*(1.-2.*np.exp(-((Wzs_HIPASS/Wc_HIPASS)**2)))) - 2.*Wzs_HIPASS*Wt_HIPASS*(1. - np.exp(-((Wzs_HIPASS/Wc_HIPASS)**2))))

Wzsti_HIPASS = Wzst_HIPASS / sin_i


#for i in range(len(HIPASS_NAME)):
#	print HIPASS_NAME[i], Wz_HIPASS[i], Wzs_HIPASS[i],Wzst_HIPASS[i], np.log10(Wzsti_HIPASS[i])

#to use the new widths for the code below

HIPASS_W50max_zsti = Wzsti_HIPASS


########################### CORRECTED EDD WIDTHS #######################################


Wt_EDD = 5.	# as in Meyer's thesis. In the EDD paper they used 9 km/s for Wt, but Wc is the same
Wc_EDD = 100.
Wzs_EDD = Wcm501 #verification has been made, by choosing the values in the paper, same results came out after the correction of  Wm501
j=0
Wm501_EDD =  [None for x in range(len(RC3_PGC))]
Wzst_EDD = [None for x in range(len(RC3_PGC))]
Wzsti_EDD = [None for x in range(len(RC3_PGC))]

#print Wzst_EDD, Wzsti_EDD

#RC3_PGC = map(lambda x: (x), RC3_PGC)
#PGC = map(lambda x: (x), PGC)

for i in range(len(RC3_PGC)):
	for j, k in enumerate(PGC):
#		print k, RC3_PGC[i]
		if k == RC3_PGC[i]:
			index = j
	#		print index
			Wm501_EDD[i] = Wm501[index]
			if sin_i[i] == 0. :
				Wzst_EDD[i] = None
				Wzsti_EDD[i] = None

			else :
				Wzst_EDD[i] = np.sqrt(Wzs_EDD[index]**2 + ((Wt_EDD**2)*(1.-2.*np.exp(-((Wzs_EDD[index]/Wc_EDD)**2)))) - 2.*Wzs_EDD[index]*Wt_EDD*(np.exp(-((Wzs_EDD[index]/Wc_EDD)**2))))

				Wzsti_EDD[i] = Wzst_EDD[i] / sin_i[i]
		
#			print sin_i[i], RC3_PGC[i], PGC[index]


#print Wzsti_EDD



##################### ERROR CALCULATION ##################################

def err_calc():
	
	W_HIPASS_error = 7.5 # in km/s
	W_EDD_error = e_W1
	err_width_HIPASS = [None for x in range(len(HIPASS_NAME))]
	err_width_EDD = [None for x in range(len(PGC))]
	err_mag_K = [None for x in range(len(HIPASS_NAME))]
	err_mag_B = [None for x in range(len(HIPASS_NAME))]
	
	#print HIPASS_RV50max
	
	min_maj_rat_RC3 = [None for x in range(len(HIPASS_NAME))]
	
	for i in range(len(HIPASS_NAME)):
		min_maj_rat_RC3[i] = 1. / (10**RC3_R25[i])
	#	print min_maj_rat_RC3[i]
	#	print RC3_R25[i]
	
	min_maj_rat_RC3_sigma = ((-1/(10**(RC3_R25+RC3_e_R25)))+(1/(10**(RC3_R25-RC3_e_R25))))/2
	#print min_maj_rat_RC3_sigma
	
	z_HIPASS = HIPASS_RV50max/c
	z_EDD = Vhel1/c
	
	#### For HIPASS widths #####
	
	
	for i in range(len(HIPASS_NAME)):
	
		c1 = -3.432
		c2 = -50.
		c3 = -1./10000.
		c4 = -10.
		c5 = 0.12*0.12
		
	
		
		tt1 = (HIPASS_W50max[i]/(1.+ z_HIPASS[i])) +c1
		tt2 = 1./(1.+z_HIPASS[i])
		tt3 = np.exp(c3*(tt1**2))
		tt4 = 1./(np.sqrt(1.-((min_maj_rat_RC3[i]**2 - c5)/(1.-c5))))
	
	
		dT1_dW = 2. * tt1 * tt2
	
		dT2_dW = 2. * c2 * c3 * tt1 * tt2 * tt3
	
		dT3_dW = c4 * ((tt2 * (1.-tt3)) + (tt1 * (-2.) * c3 * tt1 * tt3 * tt2))
	
		dWf_dW = 0.5 * tt4 * (tt1**2 + 25. + c2*tt3 + c4*tt1*(1.-tt3))**(-1./2.) * (dT1_dW + dT2_dW + dT3_dW)
	
		dT4_dminmajrat = (-1./2.) * (-2.*min_maj_rat_RC3[i]) * (1.-c5)**(1./2.) * (1.-min_maj_rat_RC3[i]**2)**(-3./2.)
	
		dWf_dminmajrat = (tt1**2 + 25. + c2*tt3 + c4*tt1*(1.-tt3))**(1./2.) * dT4_dminmajrat
	
		err_width_HIPASS[i] = np.sqrt(dWf_dW**2 * W_HIPASS_error**2 + dWf_dminmajrat**2 * min_maj_rat_RC3_sigma[i]**2)
	
	print 0.434*(err_width_HIPASS/HIPASS_W50max)  #the log10 of the error!!!!!
	
	##### For EDD widths #####
	
	for i in range(len(PGC)):
	
		c1 = -3.432
		c2 = -50.
		c3 = -1./10000.
		c4 = -10.
		c5 = 0.12*0.12
		
	
		
		tt1 = (Wm501[i]/(1.+ z_EDD[i])) +c1
		tt2 = 1./(1.+z_EDD[i])
		tt3 = np.exp(c3*(tt1**2))
		tt4 = 1./(np.sqrt(1.-((min_maj_rat_RC3[i]**2 - c5)/(1.-c5))))
	
	
		dT1_dW = 2. * tt1 * tt2
	
		dT2_dW = 2. * c2 * c3 * tt1 * tt2 * tt3
	
		dT3_dW = c4 * ((tt2 * (1.-tt3)) + (tt1 * (-2.) * c3 * tt1 * tt3 * tt2))
	
		dWf_dW = 0.5 * tt4 * (tt1**2 + 25. + c2*tt3 + c4*tt1*(1.-tt3))**(-1./2.) * (dT1_dW + dT2_dW + dT3_dW)
	
		dT4_dminmajrat = (-1./2.) * (-2.*min_maj_rat_RC3[i]) * (1.-c5)**(1./2.) * (1.-min_maj_rat_RC3[i]**2)**(-3./2.)
	
		dWf_dminmajrat = (tt1**2 + 25. + c2*tt3 + c4*tt1*(1.-tt3))**(1./2.) * dT4_dminmajrat
	
		err_width_EDD[i] = np.sqrt(dWf_dW**2 * W_EDD_error[i]**2 + dWf_dminmajrat**2 * min_maj_rat_RC3_sigma[i]**2)
	
	print 0.434*(err_width_EDD/Wm501)  #the log10 of the error!!!!!
	
	###### For absolute magnitudes #######
	
	HIPASS_RVCMB_err = (10**4)*((HIPASS_Speak*1000)**(-2))+5.
	print HIPASS_RVCMB_err
	V_pec = 125 #km/s in Meyer thesis pg.105
	for i in range(len(HIPASS_NAME)): 
	
		err_mag_K[i] = np.sqrt(TWOMASX_K_mag_sigma[i]**2 + (25./HIPASS_RVCMB[i]**2)*(HIPASS_RVCMB_err[i]**2 + V_pec**2))
		err_mag_B[i] = np.sqrt(RC3_e_Bmag[i]**2 + (25./HIPASS_RVCMB[i]**2)*(HIPASS_RVCMB_err[i]**2 + V_pec**2))
	
	print err_mag[K]
	print err_mag[B]
	
######################## CULLED & UPDATED WIDTHS ####################################


outlier_list = '/home/sinemo/Desktop/fitstable/outlier_list_edd_hipass.txt'

outlier_hipass_names = np.genfromtxt(outlier_list, dtype=None, delimiter=' ', usecols=(0), missing_values=None)

#print outlier_hipass_names
#print HIPASS_NAME
culled_updated_line_widths = [None for x in range(len(HIPASS_NAME))]
culled_updated_line_widths_corrected = [None for x in range(len(HIPASS_NAME))]
culled_updated_line_width_errors = [None for x in range(len(HIPASS_NAME))]

for m in range(len(HIPASS_NAME)):

	if HIPASS_NAME[m][1:-1] in outlier_hipass_names:   #[1:-1] needed to get rid of double quotes in the HIPASS_NAME[m]
		print HIPASS_NAME[m]
		culled_updated_line_widths[m] = Wm501_EDD[m]
		culled_updated_line_widths_corrected[m] = Wzsti_EDD[m]
		culled_updated_line_width_errors[m] = err_width_EDD[m]
	elif HIPASS_NAME[m][1:-1]=='HIPASSJ2224-03' or HIPASS_NAME[m][1:-1]=='HIPASSJ0911-14':
		print HIPASS_NAME[m]
		culled_updated_line_widths[m] = None 
		culled_updated_line_widths_corrected[m] = None
		culled_updated_line_width_errors[m] = None
	else:
		culled_updated_line_widths[m] = HIPASS_W50max[m] 
		culled_updated_line_widths_corrected[m] = HIPASS_W50max_zsti[m]
		culled_updated_line_width_errors[m] = err_width_HIPASS[m]

#	print HIPASS_W50max[m], culled_updated_line_widths[m]


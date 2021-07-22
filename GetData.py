import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
from astroquery.irsa import Irsa # had to install astroquery w/ pip
import astropy.units as u
# from astroquery.simbad import Simbad

matplotlib.rcParams.update({'font.size':18})
matplotlib.rcParams.update({'font.family':'serif'})

Irsa.ROW_LIMIT = 0
Irsa.TIMEOUT = 1000


def WISE_LC(obj, alldata=False, interac=False, moreplots=False, returndata=False):
    # the WISE tables to search
    cats = ['neowiser_p1bs_psd', 'allsky_4band_p1bs_psd', 'allsky_3band_p1bs_psd', 'allsky_2band_p1bs_psd']

    table1 = Irsa.query_region(obj, catalog=cats[0], spatial='Cone', radius=3 * u.arcsec)
    table2 = Irsa.query_region(obj, catalog=cats[1], spatial='Cone', radius=3 * u.arcsec)
    table3 = Irsa.query_region(obj, catalog=cats[2], spatial='Cone', radius=3 * u.arcsec)
    table4 = Irsa.query_region(obj, catalog=cats[3], spatial='Cone', radius=3 * u.arcsec)

    table1.sort('mjd')
    table2.sort('mjd')
    table3.sort('mjd')
    table4.sort('mjd')

    df1 = table1.to_pandas()
    df2 = table2.to_pandas()
    df3 = table3.to_pandas()
    df4 = table4.to_pandas()

    # print(' found '+str(len(df1))+' + '+ str(len(df2))+' + '+ str(len(df3))+' + '+ str(len(df4))+' visits')
    totvisits = len(df1) + len(df2) + len(df3) + len(df4)

    # dump to CSV files for possible later analysis
    if not os.path.exists('data'):
        os.makedirs('data')

    # manually fix the Python3 str=>bytestr problem... boo
    df1['ph_qual'] = df1['ph_qual'].str.decode('ascii')
    df2['ph_qual'] = df2['ph_qual'].str.decode('ascii')
    df3['ph_qual'] = df3['ph_qual'].str.decode('ascii')
    df4['ph_qual'] = df4['ph_qual'].str.decode('ascii')

    df1['cc_flags'] = df1['cc_flags'].str.decode('ascii')
    df2['cc_flags'] = df2['cc_flags'].str.decode('ascii')
    df3['cc_flags'] = df3['cc_flags'].str.decode('ascii')
    df4['cc_flags'] = df4['cc_flags'].str.decode('ascii')


    df1.to_csv('data/' + obj + cats[0] + '.csv')
    df2.to_csv('data/' + obj + cats[1] + '.csv')
    df3.to_csv('data/' + obj + cats[2] + '.csv')
    df4.to_csv('data/' + obj + cats[3] + '.csv')


    #### NEED TO ADD QUALITY CUTS
    # can't add this to the latter 3 surveys...  (df1['qual_frame'] > 8)
    # removing this for now... (July 2021)
    # ok1 = (df1['ph_qual'].str[0] == 'A') & (df1['nb'] == 1) & (df1['cc_flags'].str[0:2] == '00') & (df1['w1rchi2'] < 5) & (df1['qual_frame'] > 8)
    # ok2 = (df2['ph_qual'].str[0] == 'A') & (df2['nb'] == 1) & (df2['cc_flags'].str[0:2] == '00') & (df2['w1rchi2'] < 5)
    # ok3 = (df3['ph_qual'].str[0] == 'A') & (df3['nb'] == 1) & (df3['cc_flags'].str[0:2] == '00') & (df3['w1rchi2'] < 5)
    # ok4 = (df4['ph_qual'].str[0] == 'A') & (df4['nb'] == 1) & (df4['cc_flags'].str[0:2] == '00') & (df4['w1rchi2'] < 5)


    # colors
    colors = ['#1f77b4', '#ff7f0e', '#c5b0d5', '#d62728']

    if not os.path.exists('img'):
        os.makedirs('img')

    # 1) W1 light curve
    plt.figure()
    if alldata:
        plt.scatter(df1['mjd'], df1['w1mpro'], c='k', s=8, alpha=0.5)
        plt.scatter(df2['mjd'], df2['w1mpro'], c='k', s=8, alpha=0.5)
        plt.scatter(df3['mjd'], df3['w1mpro'], c='k', s=8, alpha=0.5)
        plt.scatter(df4['mjd'], df4['w1mpro'], c='k', s=8, alpha=0.5)

    plt.errorbar(df1['mjd'], df1['w1mpro'], yerr=df1['w1sigmpro'],
                 marker='o', linestyle='none', alpha=0.5, color=colors[0])
    plt.errorbar(df2['mjd'], df2['w1mpro'], yerr=df2['w1sigmpro'],
                 marker='o', linestyle='none', alpha=0.5, color=colors[1])
    plt.errorbar(df3['mjd'], df3['w1mpro'], yerr=df3['w1sigmpro'],
                 marker='o', linestyle='none', alpha=0.5, color=colors[2])
    plt.errorbar(df4['mjd'], df4['w1mpro'], yerr=df4['w1sigmpro'],
                 marker='o', linestyle='none', alpha=0.5, color=colors[3])

    plt.ylabel('W1 (mag)')
    plt.xlabel('MJD (days)')
    plt.gca().invert_yaxis()
    plt.title(obj + ', N=' + str(totvisits))
    plt.savefig('img/'+obj + '_W1.pdf', dpi=150, bbox_inches='tight', pad_inches=0.25)
    if interac:
        plt.show()
    # else:
    #     plt.close()



    # 2) W1-W2 color light curve
    if moreplots:
        plt.figure()
        if alldata:
            plt.scatter(df1['mjd'], df1['w1mpro']-df1['w2mpro'], c='k', s=8, alpha=0.5)
            plt.scatter(df2['mjd'], df2['w1mpro']-df2['w2mpro'], c='k', s=8, alpha=0.5)
            plt.scatter(df3['mjd'], df3['w1mpro']-df3['w2mpro'], c='k', s=8, alpha=0.5)
            plt.scatter(df4['mjd'], df4['w1mpro']-df4['w2mpro'], c='k', s=8, alpha=0.5)

        plt.errorbar(df1['mjd'], df1['w1mpro'] - df1['w2mpro'],
                     yerr=np.sqrt(df1['w1sigmpro']**2 + df1['w2sigmpro']**2),
                     marker='o', linestyle='none', alpha=0.5, color=colors[0])
        plt.errorbar(df2['mjd'], df2['w1mpro'] - df2['w2mpro'],
                     yerr=np.sqrt(df2['w1sigmpro']**2 + df2['w2sigmpro']**2),
                     marker='o', linestyle='none', alpha=0.5, color=colors[1])
        plt.errorbar(df3['mjd'], df3['w1mpro'] - df3['w2mpro'],
                     yerr=np.sqrt(df3['w1sigmpro']**2 + df3['w2sigmpro']**2),
                     marker='o', linestyle='none', alpha=0.5, color=colors[2])
        plt.errorbar(df4['mjd'], df4['w1mpro'] - df4['w2mpro'],
                     yerr=np.sqrt(df4['w1sigmpro']**2 + df4['w2sigmpro']**2),
                     marker='o', linestyle='none', alpha=0.5, color=colors[3])
        plt.xlabel('MJD (days)')
        plt.ylabel('W1-W2 (mag)')
        plt.title(obj + ', N=' + str(totvisits))
        plt.savefig('img/'+obj + '_W1W2.pdf', dpi=150, bbox_inches='tight', pad_inches=0.25)
        # plt.show()
        if interac:
            plt.show()
        # else:
        #     plt.close()



        # 3) CMD
        plt.figure()
        # plt.errorbar(df1['w1mpro'] - df1['w2mpro'], df1['w1mpro'],
        #             xerr=np.sqrt(df1['w1sigmpro']**2 + df1['w2sigmpro']**2), yerr=df1['w1sigmpro'],
        #             marker='o', linestyle='none', alpha=0.5, color='#1f77b4')
        # plt.errorbar(df2['w1mpro_ep'] - df2['w2mpro_ep'], df2['w1mpro_ep'],
        #             xerr=np.sqrt(df2['w1sigmpro_ep']**2 + df2['w2sigmpro_ep']**2 ), yerr=df2['w1sigmpro_ep'],
        #             marker='o', linestyle='none', alpha=0.5, color='#ff7f0e')
        plt.errorbar(df1['w1mpro'] - df1['w2mpro'], df1['w1mpro'],
                     xerr=np.sqrt(df1['w1sigmpro']**2 + df1['w2sigmpro']**2),
                     yerr=df1['w1sigmpro'],
                     marker='o', linestyle='none', alpha=0.5, color=colors[0])
        plt.errorbar(df2['w1mpro'] - df2['w2mpro'], df2['w1mpro'],
                     xerr=np.sqrt(df2['w1sigmpro']**2 + df2['w2sigmpro']**2),
                     yerr=df2['w1sigmpro'],
                     marker='o', linestyle='none', alpha=0.5, color=colors[1])
        plt.errorbar(df3['w1mpro'] - df3['w2mpro'], df3['w1mpro'],
                     xerr=np.sqrt(df3['w1sigmpro']**2 + df3['w2sigmpro']**2),
                     yerr=df3['w1sigmpro'],
                     marker='o', linestyle='none', alpha=0.5, color=colors[2])
        plt.errorbar(df4['w1mpro'] - df4['w2mpro'], df4['w1mpro'],
                     xerr=np.sqrt(df4['w1sigmpro']**2 + df4['w2sigmpro']**2),
                     yerr=df4['w1sigmpro'],
                     marker='o', linestyle='none', alpha=0.5, color=colors[3])

        plt.ylabel('W1 (mag)')
        plt.xlabel('W1-W2 (mag)')
        plt.gca().invert_yaxis()
        plt.savefig('img/'+obj + '_cmd.pdf', dpi=150, bbox_inches='tight', pad_inches=0.25)
        if interac:
            plt.show()
        # else:
        #     plt.close()


        plt.figure()
        Per = 0.498194
        plt.errorbar((df1['mjd'] % Per)/Per, df1['w1mpro'], yerr=df1['w1sigmpro'],
                 marker='o', linestyle='none', alpha=0.5, color=colors[0])
        plt.errorbar((df2['mjd'] % Per)/Per, df2['w1mpro'], yerr=df2['w1sigmpro'],
                    marker='o', linestyle='none', alpha=0.5, color=colors[1])
        plt.errorbar((df3['mjd'] % Per)/Per, df3['w1mpro'], yerr=df3['w1sigmpro'],
                    marker='o', linestyle='none', alpha=0.5, color=colors[2])
        plt.errorbar((df4['mjd'] % Per)/Per, df4['w1mpro'], yerr=df4['w1sigmpro'],
                    marker='o', linestyle='none', alpha=0.5, color=colors[3])

        plt.ylabel('W1 (mag)')
        plt.xlabel('Phase (P=0.498194 days)')
        plt.gca().invert_yaxis()
        plt.title(obj + ', N=' + str(totvisits))
        plt.savefig('img/'+obj + '_phase.pdf', dpi=150, bbox_inches='tight', pad_inches=0.25)
        if interac:
            plt.show()
        # else:
        #     plt.close()



        # RA,Dec to make sure not a blend, etc
        plt.figure()
        plt.scatter(df1['ra'], df1['dec'],
                    marker='o', alpha=0.5, c=color[0])
        plt.scatter(df2['ra'], df2['dec'],
                    marker='o', alpha=0.5, c=colors[1])
        plt.scatter(df3['ra'], df3['dec'],
                    marker='o', alpha=0.5, c=colors[2])
        plt.scatter(df4['ra'], df4['dec'],
                    marker='o', alpha=0.5, c=colors[3])
        plt.xlabel('RA (deg)')
        plt.ylabel('Dec (deg)')
        plt.savefig('img/' + obj + '_radec.pdf', dpi=150, bbox_inches='tight', pad_inches=0.25)
        # plt.close()

    if returndata:
        return df1, df2, df3, df4
    else:
        return


if __name__ == "__main__":
    import sys
    WISE_LC(str(sys.argv[1:]), interac=True, alldata=True, moreplots=True)
